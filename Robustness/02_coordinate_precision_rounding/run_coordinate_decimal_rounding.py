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
DECIMALS = [3, 2, 1, 0]
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


def rounded_coord(value: float, decimals: int) -> float:
    return float(np.round(value, decimals=decimals))


def round_pdb_coordinates(input_pdb: Path, output_pdb: Path, decimals: int) -> None:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for line in input_pdb.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("ATOM"):
            try:
                x = rounded_coord(float(line[30:38]), decimals)
                y = rounded_coord(float(line[38:46]), decimals)
                z = rounded_coord(float(line[46:54]), decimals)
                line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            except (ValueError, IndexError):
                pass
        lines.append(line)
    output_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expected_bits = message_to_bits(MESSAGE)
    rows = []

    for case_id in CASES:
        original_pdb = ROOT / case_id / f"{case_id}_original.pdb"
        watermarked_pdb = ROOT / case_id / f"{case_id}_watermarked_npj_SB.pdb"
        original_atoms = read_atom_records(original_pdb)
        watermarked_atoms = read_atom_records(watermarked_pdb)
        ca_keys, original_ca, clean_watermarked_ca, b_factors = matched_ca_arrays(original_atoms, watermarked_atoms)
        selected = selected_indices(b_factors, len(expected_bits))

        for decimals in DECIMALS:
            rounded_pdb = OUT_DIR / "rounded_pdb" / case_id / f"{case_id}_watermarked_round_{decimals}dec.pdb"
            round_pdb_coordinates(watermarked_pdb, rounded_pdb, decimals)
            rounded_atoms = read_atom_records(rounded_pdb)
            _, _, rounded_ca, _ = matched_ca_arrays(original_atoms, rounded_atoms)

            decoded = decode_from_ca_coords(original_ca, rounded_ca, selected, len(expected_bits))
            bit_accuracy = bitacc(decoded, expected_bits)
            rows.append(
                {
                    "case_id": case_id,
                    "attack": "coordinate_decimal_rounding",
                    "decimals": decimals,
                    "coordinate_step_angstrom": 10.0 ** (-decimals),
                    "rounded_pdb": str(rounded_pdb),
                    "ca_count": len(ca_keys),
                    "selected_ca_count": len(selected),
                    "decoded_bits": decoded,
                    "decoded_message": bits_to_text_until_null(decoded),
                    "bitacc_percent": bit_accuracy,
                    "ber_percent": 100.0 - bit_accuracy,
                    "exact_recovery": decoded == expected_bits,
                    "clean_to_rounded_ca_rmsd_angstrom": rmsd(clean_watermarked_ca, rounded_ca),
                    "original_to_rounded_ca_rmsd_angstrom": rmsd(original_ca, rounded_ca),
                    "max_ca_rounding_displacement_angstrom": float(
                        np.max(np.linalg.norm(rounded_ca - clean_watermarked_ca, axis=1))
                    ),
                }
            )

    with (OUT_DIR / "coordinate_decimal_rounding_results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (OUT_DIR / "coordinate_decimal_rounding_results.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    with (OUT_DIR / "coordinate_decimal_rounding_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write(
            "Deterministic decimal rounding of all protein ATOM coordinates in the watermarked PDB files. "
            "The original files are standard PDB-style 3-decimal coordinate files; 3 decimals therefore serve "
            "as the baseline precision.\n\n"
        )
        for case_id in CASES:
            for decimals in DECIMALS:
                row = next(row for row in rows if row["case_id"] == case_id and row["decimals"] == decimals)
                handle.write(
                    f"{case_id} decimals={decimals}: BitAcc={row['bitacc_percent']:.3f}% "
                    f"BER={row['ber_percent']:.3f}% exact={row['exact_recovery']} "
                    f"clean_to_rounded_CA_RMSD={row['clean_to_rounded_ca_rmsd_angstrom']:.5f} A\n"
                )

    print(f"Wrote {len(rows)} rows to {OUT_DIR / 'coordinate_decimal_rounding_results.csv'}")


if __name__ == "__main__":
    main()
