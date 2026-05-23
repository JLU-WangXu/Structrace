# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = Path(__file__).resolve().parent
MESSAGE = "npj SB"
START_FREQ = 5
RADIUS_ANGSTROM = 2.0
SIGMAS = [0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0]
N_REPEATS = 50
RANDOM_SEED = 20260518
CASES = ["6MRR", "8HFE", "8VC8"]


def message_to_bits(message: str) -> str:
    return "".join(f"{byte:08b}" for byte in list(message.encode("utf-8", "ignore")) + [0])


def atom_key(line: str) -> tuple[str, str, str, str, str]:
    return (
        line[12:16].strip(),
        line[17:20].strip(),
        line[21:22].strip(),
        line[22:26].strip(),
        line[26:27].strip(),
    )


def read_atom_records(pdb_path: Path) -> list[dict[str, object]]:
    records = []
    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            records.append(
                {
                    "key": atom_key(line),
                    "atom": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21:22].strip(),
                    "resseq": line[22:26].strip(),
                    "coord": np.array(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                        dtype=float,
                    ),
                    "b": float(line[60:66]),
                }
            )
        except (ValueError, IndexError):
            continue
    return records


def ca_records(records: list[dict[str, object]]) -> list[dict[str, object]]:
    return [record for record in records if record["atom"] == "CA"]


def matched_ca_arrays(
    original_records: list[dict[str, object]], query_records: list[dict[str, object]]
) -> tuple[list[tuple[str, str, str, str, str]], np.ndarray, np.ndarray, np.ndarray]:
    query_by_key = {record["key"]: record for record in ca_records(query_records)}
    keys = []
    original_coords = []
    query_coords = []
    b_factors = []
    for original_record in ca_records(original_records):
        query_record = query_by_key.get(original_record["key"])
        if query_record is None:
            continue
        keys.append(original_record["key"])
        original_coords.append(original_record["coord"])
        query_coords.append(query_record["coord"])
        b_factors.append(float(original_record["b"]))
    if len(keys) < 3:
        raise ValueError("Too few matched CA atoms")
    return keys, np.vstack(original_coords), np.vstack(query_coords), np.asarray(b_factors, dtype=float)


def selected_indices(b_factors: np.ndarray, bit_length: int) -> np.ndarray:
    top_n = max(35, 2 * (math.ceil(bit_length / 3) + START_FREQ))
    ranked = sorted(range(len(b_factors)), key=lambda index: float(b_factors[index]), reverse=True)
    return np.asarray(ranked[:top_n], dtype=int)


def decode_from_ca_coords(
    original_ca_coords: np.ndarray,
    query_ca_coords: np.ndarray,
    selected: np.ndarray,
    bit_length: int,
) -> str:
    coords0 = original_ca_coords[selected]
    coords1 = query_ca_coords[selected]
    f0 = [np.fft.fft(coords0[:, 0]), np.fft.fft(coords0[:, 1]), np.fft.fft(coords0[:, 2])]
    f1 = [np.fft.fft(coords1[:, 0]), np.fft.fft(coords1[:, 1]), np.fft.fft(coords1[:, 2])]
    bits = []
    for i in range(bit_length):
        freq = START_FREQ + i // 3
        delta = np.abs(f1[i % 3][freq]) - np.abs(f0[i % 3][freq])
        bits.append("1" if delta > 0 else "0")
    return "".join(bits)


def bits_to_text_until_null(bits: str) -> str:
    out = bytearray()
    for i in range(0, len(bits), 8):
        chunk = bits[i : i + 8]
        if len(chunk) < 8:
            break
        byte = int(chunk, 2)
        if byte == 0:
            break
        out.append(byte)
    return out.decode("utf-8", errors="replace")


def bitacc(decoded: str, expected: str) -> float:
    return sum(a == b for a, b in zip(decoded, expected)) / len(expected) * 100.0


def rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((coords_a - coords_b) ** 2, axis=1))))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expected_bits = message_to_bits(MESSAGE)
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []

    for case_id in CASES:
        original_pdb = ROOT / case_id / f"{case_id}_original.pdb"
        watermarked_pdb = ROOT / case_id / f"{case_id}_watermarked_npj_SB.pdb"
        original_atoms = read_atom_records(original_pdb)
        watermarked_atoms = read_atom_records(watermarked_pdb)
        ca_keys, original_ca, watermarked_ca, b_factors = matched_ca_arrays(original_atoms, watermarked_atoms)
        selected = selected_indices(b_factors, len(expected_bits))

        all_atom_coords = np.vstack([record["coord"] for record in watermarked_atoms])
        ca_key_to_index = {key: index for index, key in enumerate(ca_keys)}
        atom_to_ca_index = np.full(len(watermarked_atoms), -1, dtype=int)
        for atom_index, record in enumerate(watermarked_atoms):
            if record["atom"] == "CA" and record["key"] in ca_key_to_index:
                atom_to_ca_index[atom_index] = ca_key_to_index[record["key"]]

        center_indices = rng.integers(0, len(ca_keys), size=N_REPEATS)

        for repeat, center_ca_index in enumerate(center_indices):
            center_coord = watermarked_ca[int(center_ca_index)]
            distances = np.linalg.norm(all_atom_coords - center_coord, axis=1)
            perturbed_atom_mask = distances <= RADIUS_ANGSTROM
            perturbed_atom_indices = np.flatnonzero(perturbed_atom_mask)
            perturbed_ca_indices = np.unique(atom_to_ca_index[perturbed_atom_indices])
            perturbed_ca_indices = perturbed_ca_indices[perturbed_ca_indices >= 0]
            center_key = ca_keys[int(center_ca_index)]

            for sigma in SIGMAS:
                disturbed_ca = watermarked_ca.copy()
                atom_noise = rng.normal(0.0, sigma, size=(len(perturbed_atom_indices), 3))
                for local_row, atom_index in enumerate(perturbed_atom_indices):
                    ca_index = atom_to_ca_index[atom_index]
                    if ca_index >= 0:
                        disturbed_ca[ca_index] = disturbed_ca[ca_index] + atom_noise[local_row]

                decoded = decode_from_ca_coords(original_ca, disturbed_ca, selected, len(expected_bits))
                global_ca_rmsd = rmsd(original_ca, disturbed_ca)
                local_ca_rmsd = (
                    rmsd(watermarked_ca[perturbed_ca_indices], disturbed_ca[perturbed_ca_indices])
                    if len(perturbed_ca_indices) > 0
                    else 0.0
                )
                max_ca_displacement = (
                    float(
                        np.max(
                            np.linalg.norm(
                                disturbed_ca[perturbed_ca_indices] - watermarked_ca[perturbed_ca_indices],
                                axis=1,
                            )
                        )
                    )
                    if len(perturbed_ca_indices) > 0
                    else 0.0
                )

                rows.append(
                    {
                        "case_id": case_id,
                        "attack": "local_spherical_gaussian_distortion",
                        "radius_angstrom": RADIUS_ANGSTROM,
                        "sigma_angstrom": sigma,
                        "repeat": repeat + 1,
                        "random_seed": RANDOM_SEED,
                        "center_atom": "CA",
                        "center_chain": center_key[2],
                        "center_resseq": center_key[3],
                        "center_resname": center_key[1],
                        "perturbed_atom_count": int(len(perturbed_atom_indices)),
                        "perturbed_ca_count": int(len(perturbed_ca_indices)),
                        "perturbed_atom_fraction": float(len(perturbed_atom_indices) / len(watermarked_atoms)),
                        "perturbed_ca_fraction": float(len(perturbed_ca_indices) / len(ca_keys)),
                        "decoded_bits": decoded,
                        "decoded_message": bits_to_text_until_null(decoded),
                        "bitacc_percent": bitacc(decoded, expected_bits),
                        "exact_recovery": decoded == expected_bits,
                        "global_ca_rmsd_angstrom": global_ca_rmsd,
                        "local_ca_rmsd_angstrom": local_ca_rmsd,
                        "max_local_ca_displacement_angstrom": max_ca_displacement,
                    }
                )

    with (OUT_DIR / "local_distortion_sigma_results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (OUT_DIR / "local_distortion_sigma_results.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    with (OUT_DIR / "local_distortion_sigma_summary.txt").open("w", encoding="utf-8") as handle:
        for case_id in CASES:
            for sigma in SIGMAS:
                group = [row for row in rows if row["case_id"] == case_id and row["sigma_angstrom"] == sigma]
                values = [row["bitacc_percent"] for row in group]
                exact = sum(bool(row["exact_recovery"]) for row in group)
                handle.write(
                    f"{case_id} sigma={sigma:g} A: n={len(group)} "
                    f"BitAcc_mean={np.mean(values):.3f}% BitAcc_min={np.min(values):.3f}% "
                    f"BitAcc_max={np.max(values):.3f}% exact={exact}/{len(group)}\n"
                )
    print(f"Wrote {len(rows)} rows to {OUT_DIR / 'local_distortion_sigma_results.csv'}")


if __name__ == "__main__":
    main()
