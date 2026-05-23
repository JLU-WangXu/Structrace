# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent
MESSAGE = "npj SB"
START_FREQ = 5
SIGMAS = [0.0010, 0.0012, 0.0014, 0.0016, 0.0018, 0.0020]
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
    f0 = [np.fft.fft(coords0[:, axis]) for axis in range(3)]
    f1 = [np.fft.fft(coords1[:, axis]) for axis in range(3)]
    bits = []
    for index in range(bit_length):
        freq = START_FREQ + index // 3
        delta = np.abs(f1[index % 3][freq]) - np.abs(f0[index % 3][freq])
        bits.append("1" if delta > 0 else "0")
    return "".join(bits)


def bits_to_text_until_null(bits: str) -> str:
    out = bytearray()
    for index in range(0, len(bits), 8):
        chunk = bits[index : index + 8]
        if len(chunk) < 8:
            break
        byte = int(chunk, 2)
        if byte == 0:
            break
        out.append(byte)
    return out.decode("utf-8", errors="replace")


def bitacc(decoded: str, expected: str) -> float:
    return sum(observed == target for observed, target in zip(decoded, expected)) / len(expected) * 100.0


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
        ca_key_to_index = {key: index for index, key in enumerate(ca_keys)}
        ca_indices_in_atom_table = [
            atom_index
            for atom_index, record in enumerate(watermarked_atoms)
            if record["atom"] == "CA" and record["key"] in ca_key_to_index
        ]
        ca_row_for_atom_index = {
            atom_index: ca_key_to_index[watermarked_atoms[atom_index]["key"]]
            for atom_index in ca_indices_in_atom_table
        }

        for sigma in SIGMAS:
            for repeat in range(1, N_REPEATS + 1):
                atom_noise = rng.normal(0.0, sigma, size=(len(watermarked_atoms), 3))
                disturbed_ca = watermarked_ca.copy()
                for atom_index, ca_row in ca_row_for_atom_index.items():
                    disturbed_ca[ca_row] = disturbed_ca[ca_row] + atom_noise[atom_index]

                decoded = decode_from_ca_coords(original_ca, disturbed_ca, selected, len(expected_bits))
                clean_to_noisy_ca_rmsd = rmsd(watermarked_ca, disturbed_ca)
                original_to_noisy_ca_rmsd = rmsd(original_ca, disturbed_ca)
                selected_ca_rmsd = rmsd(watermarked_ca[selected], disturbed_ca[selected])
                max_ca_displacement = float(np.max(np.linalg.norm(disturbed_ca - watermarked_ca, axis=1)))

                rows.append(
                    {
                        "case_id": case_id,
                        "attack": "gaussian_coordinate_noise_all_atoms",
                        "sigma_angstrom": sigma,
                        "repeat": repeat,
                        "random_seed": RANDOM_SEED,
                        "target": "all_protein_atoms",
                        "protein_atom_count": len(watermarked_atoms),
                        "ca_count": len(ca_keys),
                        "selected_ca_count": len(selected),
                        "decoded_bits": decoded,
                        "decoded_message": bits_to_text_until_null(decoded),
                        "bitacc_percent": bitacc(decoded, expected_bits),
                        "ber_percent": 100.0 - bitacc(decoded, expected_bits),
                        "exact_recovery": decoded == expected_bits,
                        "clean_to_noisy_ca_rmsd_angstrom": clean_to_noisy_ca_rmsd,
                        "original_to_noisy_ca_rmsd_angstrom": original_to_noisy_ca_rmsd,
                        "selected_ca_rmsd_angstrom": selected_ca_rmsd,
                        "max_ca_displacement_angstrom": max_ca_displacement,
                    }
                )

    with (OUT_DIR / "gaussian_noise_all_atoms_sigma_0.001_to_0.002_results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (OUT_DIR / "gaussian_noise_all_atoms_sigma_0.001_to_0.002_results.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    with (OUT_DIR / "gaussian_noise_all_atoms_sigma_0.001_to_0.002_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write(
            "Interpolated low-amplitude global Gaussian coordinate noise on all protein ATOM records. "
            "Noise is sampled independently per coordinate from N(0, sigma^2); "
            "sigma values follow the local coordinate distortion experiment.\n\n"
        )
        for case_id in CASES:
            for sigma in SIGMAS:
                group = [row for row in rows if row["case_id"] == case_id and row["sigma_angstrom"] == sigma]
                values = np.array([row["bitacc_percent"] for row in group], dtype=float)
                exact = sum(bool(row["exact_recovery"]) for row in group)
                rmsds = np.array([row["clean_to_noisy_ca_rmsd_angstrom"] for row in group], dtype=float)
                handle.write(
                    f"{case_id} sigma={sigma:g} A: n={len(group)} "
                    f"BitAcc_mean={values.mean():.3f}% BitAcc_min={values.min():.3f}% "
                    f"BitAcc_max={values.max():.3f}% exact={exact}/{len(group)} "
                    f"clean_to_noisy_CA_RMSD_mean={rmsds.mean():.4f} A\n"
                )
    print(f"Wrote {len(rows)} rows to {OUT_DIR / 'gaussian_noise_all_atoms_sigma_0.001_to_0.002_results.csv'}")


if __name__ == "__main__":
    main()
