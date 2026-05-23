# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


ATTACK_ROOT = Path(__file__).resolve().parents[1]
CASES_ROOT = ATTACK_ROOT / "3cases"
OUT_DIR = Path(__file__).resolve().parent

MESSAGE = "npj SB"
START_FREQ = 5
CASES = ["6MRR", "8HFE", "8VC8"]

SIGMA_ANGSTROM = 0.05
N_TARGET_CA = 12
N_REPEATS = 100
RANDOM_SEED = 20260518


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
        original_pdb = CASES_ROOT / case_id / f"{case_id}_original.pdb"
        watermarked_pdb = CASES_ROOT / case_id / f"{case_id}_watermarked_npj_SB.pdb"
        original_atoms = read_atom_records(original_pdb)
        watermarked_atoms = read_atom_records(watermarked_pdb)
        ca_keys, original_ca, watermarked_ca, b_factors = matched_ca_arrays(original_atoms, watermarked_atoms)

        selected = selected_indices(b_factors, len(expected_bits))
        selected_set = set(int(index) for index in selected)
        nonselected = np.asarray([index for index in range(len(ca_keys)) if index not in selected_set], dtype=int)
        target_count = min(N_TARGET_CA, len(selected), len(nonselected))
        if target_count < 1:
            raise ValueError(f"{case_id} does not have enough non-selected CA atoms")

        for repeat in range(1, N_REPEATS + 1):
            for target_group, target_pool in [
                ("selected_CA", selected),
                ("nonselected_CA", nonselected),
            ]:
                target_indices = rng.choice(target_pool, size=target_count, replace=False)
                attacked_ca = watermarked_ca.copy()
                attacked_ca[target_indices] += rng.normal(0.0, SIGMA_ANGSTROM, size=(target_count, 3))

                decoded = decode_from_ca_coords(original_ca, attacked_ca, selected, len(expected_bits))
                target_displacements = np.linalg.norm(attacked_ca[target_indices] - watermarked_ca[target_indices], axis=1)

                rows.append(
                    {
                        "case_id": case_id,
                        "attack": "selected_vs_nonselected_ca_gaussian_perturbation",
                        "target_group": target_group,
                        "sigma_angstrom": SIGMA_ANGSTROM,
                        "target_ca_count": int(target_count),
                        "repeat": repeat,
                        "random_seed": RANDOM_SEED,
                        "ca_count": len(ca_keys),
                        "selected_ca_count": len(selected),
                        "nonselected_ca_count": len(nonselected),
                        "target_ca_indices": ";".join(str(int(index)) for index in target_indices),
                        "decoded_bits": decoded,
                        "decoded_message": bits_to_text_until_null(decoded),
                        "bitacc_percent": bitacc(decoded, expected_bits),
                        "ber_percent": 100.0 - bitacc(decoded, expected_bits),
                        "exact_recovery": decoded == expected_bits,
                        "clean_to_attacked_ca_rmsd_angstrom": rmsd(watermarked_ca, attacked_ca),
                        "target_ca_rmsd_angstrom": rmsd(watermarked_ca[target_indices], attacked_ca[target_indices]),
                        "max_target_ca_displacement_angstrom": float(np.max(target_displacements)),
                    }
                )

    with (OUT_DIR / "selected_vs_nonselected_ca_perturbation_results.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    (OUT_DIR / "selected_vs_nonselected_ca_perturbation_results.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    with (OUT_DIR / "selected_vs_nonselected_ca_perturbation_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write(
            "Random Gaussian perturbation of watermark-bearing selected C-alpha atoms versus non-selected "
            "C-alpha atoms. Each repeat perturbs the same number of C-alpha atoms with identical sigma.\n\n"
        )
        handle.write(
            f"sigma={SIGMA_ANGSTROM:g} A; target_ca_count={N_TARGET_CA}; repeats={N_REPEATS}; "
            f"message={MESSAGE!r}; start_freq={START_FREQ}\n\n"
        )
        for case_id in CASES:
            for target_group in ["nonselected_CA", "selected_CA"]:
                group = [row for row in rows if row["case_id"] == case_id and row["target_group"] == target_group]
                values = np.asarray([float(row["bitacc_percent"]) for row in group], dtype=float)
                exact = sum(bool(row["exact_recovery"]) for row in group)
                handle.write(
                    f"{case_id} {target_group}: n={len(group)} "
                    f"BitAcc_median={np.median(values):.3f}% BitAcc_mean={np.mean(values):.3f}% "
                    f"IQR={np.percentile(values, 25):.3f}-{np.percentile(values, 75):.3f}% "
                    f"min={np.min(values):.3f}% max={np.max(values):.3f}% exact={exact}/{len(group)}\n"
                )

    print(f"Wrote {len(rows)} rows to {OUT_DIR / 'selected_vs_nonselected_ca_perturbation_results.csv'}")


if __name__ == "__main__":
    main()
