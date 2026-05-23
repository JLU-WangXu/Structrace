# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "Global_translation"
MESSAGE = "npj SB"
START_FREQ = 5
TRANSLATIONS = [
    ("tx_0.1A", (0.1, 0.0, 0.0)),
    ("tx_0.5A", (0.5, 0.0, 0.0)),
    ("tx_1A", (1.0, 0.0, 0.0)),
    ("tx_5A", (5.0, 0.0, 0.0)),
    ("tx_10A", (10.0, 0.0, 0.0)),
    ("tx_50A", (50.0, 0.0, 0.0)),
    ("tx_100A", (100.0, 0.0, 0.0)),
]
CASES = ["6MRR", "8HFE", "8VC8"]


def message_to_bits(message: str) -> str:
    return "".join(f"{byte:08b}" for byte in list(message.encode("utf-8", "ignore")) + [0])


def extract_ca_from_lines(lines: list[str]) -> list[dict[str, object]]:
    data = []
    for line in lines:
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                data.append(
                    {
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                        "b": float(line[60:66]),
                        "line": line,
                    }
                )
            except (ValueError, IndexError):
                continue
    return data


def select_top_bfactor_ca(records: list[dict[str, object]], bit_length: int) -> list[dict[str, object]]:
    top_n = max(35, 2 * (math.ceil(bit_length / 3) + START_FREQ))
    return sorted(records, key=lambda record: float(record["b"]), reverse=True)[:top_n]


def decode_watermark(original_pdb: Path, query_pdb: Path, bit_length: int) -> str:
    lines0 = original_pdb.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    lines1 = query_pdb.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    cas0 = select_top_bfactor_ca(extract_ca_from_lines(lines0), bit_length)
    cas1_dict = {str(record["line"]): record for record in extract_ca_from_lines(lines1)}
    cas1 = []
    for record0 in cas0:
        line0 = str(record0["line"])
        if line0 in cas1_dict:
            cas1.append(cas1_dict[line0])
            continue
        found = False
        for line1, record1 in cas1_dict.items():
            if line1[17:26] == line0[17:26]:
                cas1.append(record1)
                found = True
                break
        if not found:
            raise ValueError(f"Cannot match CA atom {line0[17:26]!r}")

    coords0 = np.array([[record["x"], record["y"], record["z"]] for record in cas0])
    coords1 = np.array([[record["x"], record["y"], record["z"]] for record in cas1])
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


def translate_pdb(input_pdb: Path, output_pdb: Path, vector: tuple[float, float, float]) -> None:
    dx, dy, dz = vector
    output_lines = []
    for line in input_pdb.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                x = float(line[30:38]) + dx
                y = float(line[38:46]) + dy
                z = float(line[46:54]) + dz
                line = line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
            except ValueError:
                pass
        output_lines.append(line)
    output_pdb.write_text("".join(output_lines), encoding="utf-8")


def bitacc(decoded: str, expected: str) -> float:
    return sum(a == b for a, b in zip(decoded, expected)) / len(expected) * 100.0


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expected_bits = message_to_bits(MESSAGE)
    rows = []
    for case_id in CASES:
        case_dir = OUT_DIR / case_id
        case_dir.mkdir(parents=True, exist_ok=True)
        original = ROOT / case_id / f"{case_id}_original.pdb"
        watermarked = ROOT / case_id / f"{case_id}_watermarked_npj_SB.pdb"
        for label, vector in TRANSLATIONS:
            out_pdb = case_dir / f"{case_id}_watermarked_npj_SB_global_translation_{label}.pdb"
            translate_pdb(watermarked, out_pdb, vector)
            decoded = decode_watermark(original, out_pdb, len(expected_bits))
            rows.append(
                {
                    "case_id": case_id,
                    "attack": "global_translation",
                    "strength_label": label,
                    "translation_x_angstrom": vector[0],
                    "translation_y_angstrom": vector[1],
                    "translation_z_angstrom": vector[2],
                    "query_pdb": str(out_pdb),
                    "decoded_bits": decoded,
                    "decoded_message": bits_to_text_until_null(decoded),
                    "bitacc_percent": bitacc(decoded, expected_bits),
                    "exact_recovery": decoded == expected_bits,
                }
            )

    with (OUT_DIR / "global_translation_results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (OUT_DIR / "global_translation_results.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    with (OUT_DIR / "global_translation_summary.txt").open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                f"{row['case_id']} {row['strength_label']}: decoded={row['decoded_message']!r} "
                f"BitAcc={row['bitacc_percent']:.3f}% exact={row['exact_recovery']}\n"
            )
    for row in rows:
        print(
            f"{row['case_id']} {row['strength_label']}: decoded={row['decoded_message']!r}, "
            f"BitAcc={row['bitacc_percent']:.3f}%, exact={row['exact_recovery']}"
        )


if __name__ == "__main__":
    main()
