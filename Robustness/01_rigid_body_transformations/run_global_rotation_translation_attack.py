# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

from run_global_rotation_attack import (
    CASES,
    MESSAGE,
    START_FREQ,
    bitacc,
    bits_to_text_until_null,
    escaped_text,
    message_to_bits,
    rotation_matrix,
)


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "Global_rotation_translation"
ANGLES_DEG = list(range(0, 361, 30))
TRANSLATION_DISTANCES = [0.0, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0]


def atom_key(line: str) -> tuple[str, str, str, str, str]:
    return (
        line[12:16].strip(),
        line[17:20].strip(),
        line[21:22].strip(),
        line[22:26].strip(),
        line[26:27].strip(),
    )


def read_ca_records(pdb_path: Path) -> list[dict[str, object]]:
    records = []
    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                records.append(
                    {
                        "key": atom_key(line),
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


def matched_ca_arrays(original_pdb: Path, watermarked_pdb: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    original_records = read_ca_records(original_pdb)
    watermarked_by_key = {record["key"]: record for record in read_ca_records(watermarked_pdb)}
    original_coords = []
    watermarked_coords = []
    b_factors = []
    for original_record in original_records:
        watermarked_record = watermarked_by_key.get(original_record["key"])
        if watermarked_record is None:
            continue
        original_coords.append(original_record["coord"])
        watermarked_coords.append(watermarked_record["coord"])
        b_factors.append(float(original_record["b"]))
    if len(original_coords) < 3:
        raise ValueError(f"Too few matched CA atoms for {original_pdb}")
    return np.vstack(original_coords), np.vstack(watermarked_coords), np.asarray(b_factors, dtype=float)


def selected_indices(b_factors: np.ndarray, bit_length: int) -> np.ndarray:
    top_n = max(35, 2 * (math.ceil(bit_length / 3) + START_FREQ))
    ranked = sorted(range(len(b_factors)), key=lambda index: float(b_factors[index]), reverse=True)
    return np.asarray(ranked[:top_n], dtype=int)


def apply_rigid_attack(coords: np.ndarray, angle_x_deg: float, angle_z_deg: float, translation_x: float) -> np.ndarray:
    center = coords.mean(axis=0)
    matrix = rotation_matrix("z", angle_z_deg) @ rotation_matrix("x", angle_x_deg)
    translation = np.array([translation_x, 0.0, 0.0], dtype=float)
    return (coords - center) @ matrix.T + center + translation


def kabsch_transform(mobile_coords: np.ndarray, reference_coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mobile_centroid = mobile_coords.mean(axis=0)
    reference_centroid = reference_coords.mean(axis=0)
    mobile_centered = mobile_coords - mobile_centroid
    reference_centered = reference_coords - reference_centroid
    covariance = mobile_centered.T @ reference_centered
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[2, 2] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    translation = reference_centroid - mobile_centroid @ rotation
    return rotation, translation


def apply_transform(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return coords @ rotation + translation


def coordinate_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((coords_a - coords_b) ** 2, axis=1))))


def kabsch_align(mobile_coords: np.ndarray, reference_coords: np.ndarray) -> tuple[np.ndarray, float]:
    rotation, translation = kabsch_transform(mobile_coords, reference_coords)
    aligned = apply_transform(mobile_coords, rotation, translation)
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - reference_coords) ** 2, axis=1))))
    return aligned, rmsd


def decode_from_coords(original_coords: np.ndarray, query_coords: np.ndarray, selected: np.ndarray, bit_length: int) -> str:
    coords0 = original_coords[selected]
    coords1 = query_coords[selected]
    f0 = [np.fft.fft(coords0[:, 0]), np.fft.fft(coords0[:, 1]), np.fft.fft(coords0[:, 2])]
    f1 = [np.fft.fft(coords1[:, 0]), np.fft.fft(coords1[:, 1]), np.fft.fft(coords1[:, 2])]
    bits = []
    for i in range(bit_length):
        freq = START_FREQ + i // 3
        delta = np.abs(f1[i % 3][freq]) - np.abs(f0[i % 3][freq])
        bits.append("1" if delta > 0 else "0")
    return "".join(bits)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expected_bits = message_to_bits(MESSAGE)
    rows = []

    for case_id in CASES:
        original_pdb = ROOT / case_id / f"{case_id}_original.pdb"
        watermarked_pdb = ROOT / case_id / f"{case_id}_watermarked_npj_SB.pdb"
        original_coords, watermarked_coords, b_factors = matched_ca_arrays(original_pdb, watermarked_pdb)
        selected = selected_indices(b_factors, len(expected_bits))

        for translation_distance in TRANSLATION_DISTANCES:
            for angle_x in ANGLES_DEG:
                for angle_z in ANGLES_DEG:
                    attacked_original = apply_rigid_attack(
                        original_coords,
                        angle_x,
                        angle_z,
                        translation_distance,
                    )
                    attacked_watermarked = apply_rigid_attack(
                        watermarked_coords,
                        angle_x,
                        angle_z,
                        translation_distance,
                    )
                    rotation, translation = kabsch_transform(attacked_original, original_coords)
                    aligned_original = apply_transform(attacked_original, rotation, translation)
                    aligned_watermarked = apply_transform(attacked_watermarked, rotation, translation)
                    original_alignment_rmsd = coordinate_rmsd(aligned_original, original_coords)
                    watermarked_alignment_rmsd = coordinate_rmsd(aligned_watermarked, original_coords)
                    decoded = decode_from_coords(
                        aligned_original,
                        aligned_watermarked,
                        selected,
                        len(expected_bits),
                    )
                    rows.append(
                        {
                            "case_id": case_id,
                            "attack": "global_rotation_translation",
                            "attack_mode": "combined_xz_rotation_x_translation_then_kabsch_alignment",
                            "angle_x_deg": angle_x,
                            "angle_z_deg": angle_z,
                            "translation_x_angstrom": translation_distance,
                            "original_alignment_ca_rmsd": original_alignment_rmsd,
                            "watermarked_alignment_ca_rmsd": watermarked_alignment_rmsd,
                            "decoded_bits": decoded,
                            "decoded_message": bits_to_text_until_null(decoded),
                            "decoded_message_escaped": escaped_text(bits_to_text_until_null(decoded)),
                            "bitacc_percent": bitacc(decoded, expected_bits),
                            "exact_recovery": decoded == expected_bits,
                        }
                    )

    with (OUT_DIR / "global_rotation_translation_results.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (OUT_DIR / "global_rotation_translation_results.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    with (OUT_DIR / "global_rotation_translation_summary.txt").open("w", encoding="utf-8") as handle:
        for case_id in CASES:
            case_rows = [row for row in rows if row["case_id"] == case_id]
            bitacc_values = [float(row["bitacc_percent"]) for row in case_rows]
            exact_count = sum(bool(row["exact_recovery"]) for row in case_rows)
            handle.write(
                f"{case_id}: n={len(case_rows)} BitAcc_min={min(bitacc_values):.3f}% "
                f"BitAcc_max={max(bitacc_values):.3f}% exact={exact_count}/{len(case_rows)}\n"
            )

    print(f"Wrote {len(rows)} rows to {OUT_DIR / 'global_rotation_translation_results.csv'}")


if __name__ == "__main__":
    main()
