from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np


@dataclass(frozen=True)
class CAAtom:
    """Parsed C-alpha atom record used for deterministic watermarking."""

    line_index: int
    line: str
    chain: str
    resseq: str
    icode: str
    bfactor: float
    coord: np.ndarray

    @property
    def residue_key(self) -> Tuple[str, str, str]:
        return (self.chain, self.resseq, self.icode)


def read_pdb_lines(path: Union[str, Path]) -> List[str]:
    return Path(path).read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)


def write_pdb_lines(path: Union[str, Path], lines: List[str]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("".join(lines), encoding="utf-8")


def extract_ca_atoms(lines: List[str]) -> List[CAAtom]:
    atoms: List[CAAtom] = []
    for idx, line in enumerate(lines):
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            coord = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
            bfactor = float(line[60:66])
        except (ValueError, IndexError):
            continue
        atoms.append(
            CAAtom(
                line_index=idx,
                line=line,
                chain=line[21:22].strip(),
                resseq=line[22:26].strip(),
                icode=line[26:27].strip(),
                bfactor=bfactor,
                coord=coord,
            )
        )
    return atoms


def selected_top_bfactor_atoms(ca_atoms: List[CAAtom], top_n: int) -> List[CAAtom]:
    if len(ca_atoms) < top_n:
        raise ValueError(f"Need at least {top_n} C-alpha atoms, found {len(ca_atoms)}.")
    return sorted(ca_atoms, key=lambda atom: atom.bfactor, reverse=True)[:top_n]


def required_top_n(bit_length: int, start_freq: int = 5, minimum: int = 35) -> int:
    if bit_length <= 0:
        raise ValueError("bit_length must be positive.")
    if start_freq < 1:
        raise ValueError("start_freq must be >= 1.")
    needed_freqs = int(np.ceil(bit_length / 3))
    return max(minimum, 2 * (needed_freqs + start_freq))


def replace_atom_coordinates(lines: List[str], atoms: List[CAAtom], coords: np.ndarray) -> List[str]:
    if len(atoms) != len(coords):
        raise ValueError("atoms and coords must have the same length.")
    out = list(lines)
    for atom, coord in zip(atoms, coords):
        x, y, z = (f"{float(value):8.3f}" for value in coord)
        line = out[atom.line_index]
        out[atom.line_index] = line[:30] + x + y + z + line[54:]
    return out


def ca_coordinate_array(atoms: List[CAAtom]) -> np.ndarray:
    return np.array([atom.coord for atom in atoms], dtype=float)
