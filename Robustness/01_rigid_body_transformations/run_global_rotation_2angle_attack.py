# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from run_global_rotation_attack import (
    CASES,
    MESSAGE,
    OUT_DIR,
    align_pdb_to_reference,
    bitacc,
    bits_to_text_until_null,
    decode_watermark,
    escaped_text,
    message_to_bits,
    read_coordinate_centroid,
    rotation_matrix,
    safe_console_text,
)


ANGLES_DEG = list(range(0, 361, 30))
ROOT = Path(__file__).resolve().parent


def rotate_pdb_two_angles(input_pdb: Path, output_pdb: Path, angle_x_deg: float, angle_z_deg: float) -> None:
    center = read_coordinate_centroid(input_pdb)
    matrix = rotation_matrix("z", angle_z_deg) @ rotation_matrix("x", angle_x_deg)
    output_lines = []
    for line in input_pdb.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                rotated = center + matrix @ (coord - center)
                line = line[:30] + f"{rotated[0]:8.3f}{rotated[1]:8.3f}{rotated[2]:8.3f}" + line[54:]
            except ValueError:
                pass
        output_lines.append(line)
    output_pdb.write_text("".join(output_lines), encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expected_bits = message_to_bits(MESSAGE)
    rows = []
    for case_id in CASES:
        case_dir = OUT_DIR / case_id / "two_angle_grid"
        case_dir.mkdir(parents=True, exist_ok=True)
        original = ROOT / case_id / f"{case_id}_original.pdb"
        watermarked = ROOT / case_id / f"{case_id}_watermarked_npj_SB.pdb"
        for angle_x in ANGLES_DEG:
            for angle_z in ANGLES_DEG:
                label = f"rx_{angle_x}deg_rz_{angle_z}deg"
                attacked_original = case_dir / f"{case_id}_original_global_rotation_{label}.pdb"
                attacked_watermarked = case_dir / f"{case_id}_watermarked_npj_SB_global_rotation_{label}.pdb"
                aligned_original = case_dir / f"{case_id}_original_global_rotation_{label}_aligned.pdb"
                aligned_watermarked = case_dir / f"{case_id}_watermarked_npj_SB_global_rotation_{label}_aligned.pdb"

                rotate_pdb_two_angles(original, attacked_original, angle_x, angle_z)
                rotate_pdb_two_angles(watermarked, attacked_watermarked, angle_x, angle_z)
                original_alignment_rmsd = align_pdb_to_reference(original, attacked_original, aligned_original)
                watermarked_alignment_rmsd = align_pdb_to_reference(original, attacked_watermarked, aligned_watermarked)
                decoded = decode_watermark(aligned_original, aligned_watermarked, len(expected_bits))

                rows.append(
                    {
                        "case_id": case_id,
                        "attack": "global_rotation_2angle_aligned",
                        "attack_mode": "combined_xz_rotation_then_kabsch_alignment",
                        "angle_x_deg": angle_x,
                        "angle_z_deg": angle_z,
                        "strength_label": label,
                        "attacked_original_pdb": str(attacked_original),
                        "attacked_query_pdb": str(attacked_watermarked),
                        "aligned_original_pdb": str(aligned_original),
                        "aligned_query_pdb": str(aligned_watermarked),
                        "original_alignment_ca_rmsd": original_alignment_rmsd,
                        "watermarked_alignment_ca_rmsd": watermarked_alignment_rmsd,
                        "decoded_bits": decoded,
                        "decoded_message": bits_to_text_until_null(decoded),
                        "decoded_message_escaped": escaped_text(bits_to_text_until_null(decoded)),
                        "bitacc_percent": bitacc(decoded, expected_bits),
                        "exact_recovery": decoded == expected_bits,
                    }
                )

    with (OUT_DIR / "global_rotation_2angle_results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (OUT_DIR / "global_rotation_2angle_results.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    with (OUT_DIR / "global_rotation_2angle_summary.txt").open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                f"{row['case_id']} rx {row['angle_x_deg']} deg rz {row['angle_z_deg']} deg: "
                f"decoded={row['decoded_message_escaped']!r} BitAcc={row['bitacc_percent']:.3f}% "
                f"exact={row['exact_recovery']} align_RMSD={row['watermarked_alignment_ca_rmsd']:.6f}\n"
            )
    for row in rows:
        print(
            f"{row['case_id']} rx {row['angle_x_deg']} deg rz {row['angle_z_deg']} deg: "
            f"decoded={safe_console_text(str(row['decoded_message']))!r}, "
            f"BitAcc={row['bitacc_percent']:.3f}%, exact={row['exact_recovery']}"
        )


if __name__ == "__main__":
    main()
