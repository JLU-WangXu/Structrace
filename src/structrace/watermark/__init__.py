"""Fourier-domain watermark embedding and decoding."""

from .codec import DecodeResult, EmbedResult, decode_text, embed_text

__all__ = [
    "DecodeResult",
    "EmbedResult",
    "decode_text",
    "embed_text",
]
