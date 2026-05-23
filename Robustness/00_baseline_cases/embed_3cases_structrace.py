# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
import json
import math
import shutil
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent

MESSAGE = "npj SB"
START_FREQ = 5
STRENGTH = 0.02

CASES = [
    ("6MRR", ROOT / "PDB" / "6mrr_fixed.pdb"),
    ("8HFE", ROOT / "PDB" / "8HFE_fixed.pdb"),
    ("8VC8", ROOT / "PDB" / "8vc8_fixed.pdb"),
]


def message_to_bits(message: str) -> str:
    # Match Structrace high-density text mode: UTF-8 bytes followed by NULL terminator.
    msg_bytes = list(message.encode("utf-8", "ignore")) + [0]
    return "".join(f"{byte:08b}" for byte in msg_bytes)


def extract_ca_from_lines(lines: list[str]) -> list[dict[str, object]]:
    """Strictly mirror Structrace Watermark.py: use ATOM CA only, excluding HETATM solvent/ligands."""
    data: list[dict[str, object]] = []
    for line in lines:
        if not (line.startswith("ATOM") and line[12:16].strip() == "CA"):
            continue
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


def required_top_n(bits_len: int) -> int:
    needed_freqs = math.ceil(bits_len / 3)
    return max(35, 2 * (needed_freqs + START_FREQ))


def select_top_bfactor_ca(ca_records: list[dict[str, object]], bits_len: int) -> list[dict[str, object]]:
    top_n = required_top_n(bits_len)
    if len(ca_records) < top_n:
        raise ValueError(f"Need at least {top_n} CA atoms, found {len(ca_records)}")
    return sorted(ca_records, key=lambda record: float(record["b"]), reverse=True)[:top_n]


def embed_watermark(pdb_path: Path, watermark_bits: str, output_path: Path) -> dict[str, object]:
    lines = pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    all_ca = extract_ca_from_lines(lines)
    selected = select_top_bfactor_ca(all_ca, len(watermark_bits))

    x = np.array([record["x"] for record in selected])
    y = np.array([record["y"] for record in selected])
    z = np.array([record["z"] for record in selected])
    fx, fy, fz = np.fft.fft(x), np.fft.fft(y), np.fft.fft(z)
    fx1, fy1, fz1 = fx.copy(), fy.copy(), fz.copy()
    n = len(fx1)

    for i, bit_char in enumerate(watermark_bits):
        freq = START_FREQ + i // 3
        arr = [fx1, fy1, fz1][i % 3]
        mag, phi = np.abs(arr[freq]), np.angle(arr[freq])
        mag_new = mag + (STRENGTH if bit_char == "1" else -STRENGTH)
        arr[freq] = mag_new * np.exp(1j * phi)
        arr[n - freq] = mag_new * np.exp(-1j * phi)

    x1 = np.fft.ifft(fx1).real
    y1 = np.fft.ifft(fy1).real
    z1 = np.fft.ifft(fz1).real

    new_lines = lines.copy()
    line_to_indices: dict[str, list[int]] = {}
    for idx, line in enumerate(lines):
        line_to_indices.setdefault(line, []).append(idx)

    used_indices: set[int] = set()
    selected_rows = []
    for i, record in enumerate(selected):
        line = str(record["line"])
        idx_in_lines = None
        for candidate in line_to_indices[line]:
            if candidate not in used_indices:
                idx_in_lines = candidate
                used_indices.add(candidate)
                break
        if idx_in_lines is None:
            continue

        xs, ys, zs = f"{x1[i]:8.3f}", f"{y1[i]:8.3f}", f"{z1[i]:8.3f}"
        new_lines[idx_in_lines] = line[:30] + xs + ys + zs + line[54:]
        selected_rows.append(
            {
                "rank": i + 1,
                "line_index": idx_in_lines,
                "chain": line[21].strip(),
                "resseq": line[22:26].strip(),
                "icode": line[26].strip(),
                "atom": line[12:16].strip(),
                "bfactor": float(record["b"]),
                "x_original": float(record["x"]),
                "y_original": float(record["y"]),
                "z_original": float(record["z"]),
                "x_watermarked": float(x1[i]),
                "y_watermarked": float(y1[i]),
                "z_watermarked": float(z1[i]),
            }
        )

    output_path.write_text("".join(new_lines), encoding="utf-8")

    encoded_ca = extract_ca_from_lines(new_lines)
    encoded_ca = extract_ca_from_lines(new_lines)
    original_coords = np.array([[record["x"], record["y"], record["z"]] for record in all_ca])
    encoded_coords = np.array([[record["x"], record["y"], record["z"]] for record in encoded_ca])
    global_ca_rmsd = float(np.sqrt(np.mean(np.sum((encoded_coords - original_coords) ** 2, axis=1))))
    selected_disp = np.sqrt((x1 - x) ** 2 + (y1 - y) ** 2 + (z1 - z) ** 2)
    selected_ca_rmsd = float(np.sqrt(np.mean(selected_disp**2)))

    return {
        "selected_rows": selected_rows,
        "global_ca_rmsd": global_ca_rmsd,
        "selected_ca_rmsd": selected_ca_rmsd,
        "max_selected_ca_displacement": float(np.max(selected_disp)),
        "top_n": len(selected),
    }


def decode_watermark(original_pdb_path: Path, encoded_pdb_path: Path, bit_length: int) -> str:
    """Strictly mirror Structrace Watermark.py decode_watermark for CA records."""
    lines0 = original_pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    lines1 = encoded_pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    all_ca0 = extract_ca_from_lines(lines0)
    cas0 = select_top_bfactor_ca(all_ca0, bit_length)
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
            raise ValueError("Cannot find corresponding CA atom in encoded PDB")

    x0 = np.array([record["x"] for record in cas0])
    y0 = np.array([record["y"] for record in cas0])
    z0 = np.array([record["z"] for record in cas0])
    x1 = np.array([record["x"] for record in cas1])
    y1 = np.array([record["y"] for record in cas1])
    z1 = np.array([record["z"] for record in cas1])
    f0 = [np.fft.fft(x0), np.fft.fft(y0), np.fft.fft(z0)]
    f1 = [np.fft.fft(x1), np.fft.fft(y1), np.fft.fft(z1)]

    decoded = []
    for i in range(bit_length):
        freq = START_FREQ + i // 3
        delta = np.abs(f1[i % 3][freq]) - np.abs(f0[i % 3][freq])
        decoded.append("1" if delta > 0 else "0")
    return "".join(decoded)


def bits_to_text_until_null(bits: str) -> str:
    output = bytearray()
    for i in range(0, len(bits), 8):
        byte_bits = bits[i : i + 8]
        if len(byte_bits) < 8:
            break
        byte = int(byte_bits, 2)
        if byte == 0:
            break
        output.append(byte)
    return output.decode("utf-8", errors="replace")


def write_selected_atoms(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    watermark_bits = message_to_bits(MESSAGE)
    results = []

    for case_id, pdb_path in CASES:
        case_dir = OUT_DIR / case_id
        case_dir.mkdir(parents=True, exist_ok=True)
        original_copy = case_dir / f"{case_id}_original.pdb"
        output_path = case_dir / f"{case_id}_watermarked_npj_SB.pdb"
        selected_csv = case_dir / f"{case_id}_selected_ca_atoms.csv"

        shutil.copy2(pdb_path, original_copy)
        metrics = embed_watermark(original_copy, watermark_bits, output_path)
        decoded_bits = decode_watermark(original_copy, output_path, len(watermark_bits))
        decoded_message = bits_to_text_until_null(decoded_bits)
        exact_recovery = decoded_bits == watermark_bits
        bitacc = sum(a == b for a, b in zip(decoded_bits, watermark_bits)) / len(watermark_bits) * 100.0
        write_selected_atoms(selected_csv, metrics["selected_rows"])

        results.append(
            {
                "case_id": case_id,
                "case_dir": str(case_dir),
                "source_original_pdb": str(pdb_path),
                "original_pdb": str(original_copy),
                "watermarked_pdb": str(output_path),
                "selected_ca_csv": str(selected_csv),
                "message": MESSAGE,
                "watermark_bits": watermark_bits,
                "payload_bits_with_null": len(watermark_bits),
                "start_freq": START_FREQ,
                "strength": STRENGTH,
                "top_n_selected_ca": metrics["top_n"],
                "selected_atom_policy": "ATOM records only; CA atoms only; ranked by original PDB B-factor; HETATM/water/solvent excluded",
                "global_ca_rmsd_angstrom": metrics["global_ca_rmsd"],
                "selected_ca_rmsd_angstrom": metrics["selected_ca_rmsd"],
                "max_selected_ca_displacement_angstrom": metrics["max_selected_ca_displacement"],
                "decoded_bits": decoded_bits,
                "decoded_message": decoded_message,
                "bitacc_percent": bitacc,
                "exact_recovery": exact_recovery,
            }
        )

    (OUT_DIR / "embed_results.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    with (OUT_DIR / "embed_results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)
    with (OUT_DIR / "embed_summary.txt").open("w", encoding="utf-8") as handle:
        for row in results:
            handle.write(
                f"{row['case_id']}: output={row['watermarked_pdb']} decoded_message={row['decoded_message']!r} "
                f"BitAcc={row['bitacc_percent']:.3f}% exact={row['exact_recovery']} "
                f"global_CA_RMSD={row['global_ca_rmsd_angstrom']:.12f} Å "
                f"selected_CA_RMSD={row['selected_ca_rmsd_angstrom']:.12f} Å\n"
            )

    for row in results:
        print(
            f"{row['case_id']}: decoded={row['decoded_message']!r}, "
            f"BitAcc={row['bitacc_percent']:.3f}%, exact={row['exact_recovery']}"
        )


if __name__ == "__main__":
    main()
