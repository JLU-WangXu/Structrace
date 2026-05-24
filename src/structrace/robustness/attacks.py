from __future__ import annotations

from pathlib import Path
from typing import Tuple, Union

import numpy as np


def _rewrite_atom_coord(line: str, coord: np.ndarray) -> str:
    x, y, z = (f"{float(value):8.3f}" for value in coord)
    return line[:30] + x + y + z + line[54:]


def translate_pdb(input_pdb: Union[str, Path], output_pdb: Union[str, Path], vector: Tuple[float, float, float]) -> None:
    vector_arr = np.asarray(vector, dtype=float)
    out_lines = []
    for line in Path(input_pdb).read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            try:
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                out_lines.append(_rewrite_atom_coord(line, coord + vector_arr))
                continue
            except ValueError:
                pass
        out_lines.append(line)
    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    Path(output_pdb).write_text("".join(out_lines), encoding="utf-8")


def round_pdb_coordinates(input_pdb: Union[str, Path], output_pdb: Union[str, Path], decimals: int = 3) -> None:
    out_lines = []
    for line in Path(input_pdb).read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            try:
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                coord = np.round(coord, decimals=decimals)
                out_lines.append(_rewrite_atom_coord(line, coord))
                continue
            except ValueError:
                pass
        out_lines.append(line)
    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    Path(output_pdb).write_text("".join(out_lines), encoding="utf-8")


def add_gaussian_noise_to_atoms(
    input_pdb: Union[str, Path],
    output_pdb: Union[str, Path],
    sigma: float,
    *,
    seed: int = 1,
    atom_records: Tuple[str, ...] = ("ATOM",),
) -> None:
    rng = np.random.default_rng(seed)
    out_lines = []
    for line in Path(input_pdb).read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True):
        if line.startswith(atom_records):
            try:
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                out_lines.append(_rewrite_atom_coord(line, coord + rng.normal(0.0, sigma, size=3)))
                continue
            except ValueError:
                pass
        out_lines.append(line)
    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    Path(output_pdb).write_text("".join(out_lines), encoding="utf-8")
