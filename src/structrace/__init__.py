"""StrucTrace package API."""

from .watermark import (
    DecodeResult,
    EmbedResult,
    decode_text,
    embed_text,
)

__all__ = [
    "DecodeResult",
    "EmbedResult",
    "decode_text",
    "embed_text",
]

__version__ = "0.1.0"
