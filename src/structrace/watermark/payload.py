from __future__ import annotations


def text_to_bits(text: str, add_null: bool = True) -> str:
    data = text.encode("utf-8")
    bits = "".join(f"{byte:08b}" for byte in data)
    return bits + ("00000000" if add_null else "")


def bits_to_text(bits: str, stop_at_null: bool = True) -> str:
    if len(bits) % 8:
        bits = bits[: len(bits) - (len(bits) % 8)]
    data = bytearray(int(bits[i : i + 8], 2) for i in range(0, len(bits), 8))
    if stop_at_null and 0 in data:
        data = data[: data.index(0)]
    return data.decode("utf-8", errors="replace")


def validate_bits(bits: str) -> str:
    bits = bits.strip()
    if not bits or any(char not in "01" for char in bits):
        raise ValueError("bits must be a non-empty string containing only 0 and 1.")
    return bits
