from __future__ import annotations

import numpy as np


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"RMSD arrays must have the same shape, got {a.shape} and {b.shape}.")
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=-1))))


def bit_accuracy(expected_bits: str, decoded_bits: str) -> float:
    if not expected_bits:
        raise ValueError("expected_bits must not be empty.")
    n = min(len(expected_bits), len(decoded_bits))
    if n == 0:
        return 0.0
    return sum(a == b for a, b in zip(expected_bits[:n], decoded_bits[:n])) / n


def bit_error_rate(expected_bits: str, decoded_bits: str) -> float:
    return 1.0 - bit_accuracy(expected_bits, decoded_bits)
