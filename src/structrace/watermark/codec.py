from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import numpy as np

from .payload import bits_to_text, text_to_bits, validate_bits
from .pdb import (
    ca_coordinate_array,
    extract_ca_atoms,
    read_pdb_lines,
    replace_atom_coordinates,
    required_top_n,
    selected_top_bfactor_atoms,
    write_pdb_lines,
)


@dataclass(frozen=True)
class EmbedResult:
    input_pdb: str
    output_pdb: str
    bit_length: int
    top_n_selected_ca: int
    global_ca_rmsd: float
    selected_ca_rmsd: float
    max_selected_ca_displacement: float


@dataclass(frozen=True)
class DecodeResult:
    master_pdb: str
    query_pdb: str
    bit_length: int
    decoded_bits: str
    decoded_text: Optional[str]
    bit_accuracy: Optional[float]
    exact_recovery: Optional[bool]


def _embed_selected_coordinates(coords: np.ndarray, bits: str, start_freq: int, strength: float) -> np.ndarray:
    axes = [np.fft.fft(coords[:, axis]) for axis in range(3)]
    n = len(coords)
    for i, char in enumerate(bits):
        freq = start_freq + i // 3
        axis = i % 3
        arr = axes[axis]
        mag = np.abs(arr[freq])
        phase = np.angle(arr[freq])
        mag_new = mag + (strength if char == "1" else -strength)
        arr[freq] = mag_new * np.exp(1j * phase)
        arr[n - freq] = mag_new * np.exp(-1j * phase)
    return np.column_stack([np.fft.ifft(axis).real for axis in axes])


def embed_bits(
    input_pdb: Union[str, Path],
    bits: str,
    output_pdb: Union[str, Path],
    *,
    start_freq: int = 5,
    strength: float = 0.02,
    minimum_top_n: int = 35,
) -> EmbedResult:
    bits = validate_bits(bits)
    lines = read_pdb_lines(input_pdb)
    ca_atoms = extract_ca_atoms(lines)
    top_n = required_top_n(len(bits), start_freq=start_freq, minimum=minimum_top_n)
    selected = selected_top_bfactor_atoms(ca_atoms, top_n)
    selected_coords = ca_coordinate_array(selected)
    embedded_coords = _embed_selected_coordinates(selected_coords, bits, start_freq, strength)
    out_lines = replace_atom_coordinates(lines, selected, embedded_coords)
    write_pdb_lines(output_pdb, out_lines)

    out_ca_atoms = extract_ca_atoms(out_lines)
    global_rmsd = float(np.sqrt(np.mean(np.sum((ca_coordinate_array(out_ca_atoms) - ca_coordinate_array(ca_atoms)) ** 2, axis=1))))
    selected_displacements = np.linalg.norm(embedded_coords - selected_coords, axis=1)
    selected_rmsd = float(np.sqrt(np.mean(np.sum((embedded_coords - selected_coords) ** 2, axis=1))))
    return EmbedResult(
        input_pdb=str(input_pdb),
        output_pdb=str(output_pdb),
        bit_length=len(bits),
        top_n_selected_ca=top_n,
        global_ca_rmsd=global_rmsd,
        selected_ca_rmsd=selected_rmsd,
        max_selected_ca_displacement=float(np.max(selected_displacements)),
    )


def embed_text(
    input_pdb: Union[str, Path],
    text: str,
    output_pdb: Union[str, Path],
    *,
    add_null: bool = True,
    start_freq: int = 5,
    strength: float = 0.02,
    minimum_top_n: int = 35,
) -> EmbedResult:
    return embed_bits(
        input_pdb,
        text_to_bits(text, add_null=add_null),
        output_pdb,
        start_freq=start_freq,
        strength=strength,
        minimum_top_n=minimum_top_n,
    )


def decode_bits(
    master_pdb: Union[str, Path],
    query_pdb: Union[str, Path],
    bit_length: int,
    *,
    expected_bits: Optional[str] = None,
    start_freq: int = 5,
    minimum_top_n: int = 35,
) -> DecodeResult:
    master_lines = read_pdb_lines(master_pdb)
    query_lines = read_pdb_lines(query_pdb)
    master_ca = extract_ca_atoms(master_lines)
    query_ca = extract_ca_atoms(query_lines)
    query_by_residue = {atom.residue_key: atom for atom in query_ca}
    top_n = required_top_n(bit_length, start_freq=start_freq, minimum=minimum_top_n)
    selected_master = selected_top_bfactor_atoms(master_ca, top_n)
    selected_query = []
    for atom in selected_master:
        try:
            selected_query.append(query_by_residue[atom.residue_key])
        except KeyError as exc:
            raise ValueError(f"Missing matching C-alpha atom for residue {atom.residue_key}.") from exc

    master_coords = ca_coordinate_array(selected_master)
    query_coords = ca_coordinate_array(selected_query)
    master_fft = [np.fft.fft(master_coords[:, axis]) for axis in range(3)]
    query_fft = [np.fft.fft(query_coords[:, axis]) for axis in range(3)]
    decoded = []
    for i in range(bit_length):
        freq = start_freq + i // 3
        axis = i % 3
        delta = np.abs(query_fft[axis][freq]) - np.abs(master_fft[axis][freq])
        decoded.append("1" if delta > 0 else "0")
    decoded_bits = "".join(decoded)
    accuracy = None
    exact = None
    if expected_bits is not None:
        expected_bits = validate_bits(expected_bits)
        compare_len = min(len(expected_bits), len(decoded_bits))
        accuracy = sum(a == b for a, b in zip(expected_bits[:compare_len], decoded_bits[:compare_len])) / compare_len
        exact = decoded_bits == expected_bits
    decoded_text = bits_to_text(decoded_bits) if bit_length >= 8 else None
    return DecodeResult(
        master_pdb=str(master_pdb),
        query_pdb=str(query_pdb),
        bit_length=bit_length,
        decoded_bits=decoded_bits,
        decoded_text=decoded_text,
        bit_accuracy=accuracy,
        exact_recovery=exact,
    )


def decode_text(
    master_pdb: Union[str, Path],
    query_pdb: Union[str, Path],
    bit_length: int,
    *,
    expected_text: Optional[str] = None,
    start_freq: int = 5,
    minimum_top_n: int = 35,
) -> DecodeResult:
    expected_bits = text_to_bits(expected_text) if expected_text is not None else None
    return decode_bits(
        master_pdb,
        query_pdb,
        bit_length,
        expected_bits=expected_bits,
        start_freq=start_freq,
        minimum_top_n=minimum_top_n,
    )
