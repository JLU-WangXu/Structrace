"""Fourier-domain watermark embedding and decoding."""

from .codec import DecodeResult, EmbedResult, decode_bits, decode_text, embed_bits, embed_text
from .pdb import CAAtom, extract_ca_atoms, read_pdb_lines

__all__ = [
    "CAAtom",
    "DecodeResult",
    "EmbedResult",
    "decode_bits",
    "decode_text",
    "embed_bits",
    "embed_text",
    "extract_ca_atoms",
    "read_pdb_lines",
]
