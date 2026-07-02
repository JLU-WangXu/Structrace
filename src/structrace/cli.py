from __future__ import annotations

import argparse
import json
from typing import List, Optional

from .watermark import decode_text, embed_text


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="structrace", description="Embed and decode StrucTrace text watermarks.")
    sub = parser.add_subparsers(dest="command", required=True)

    embed = sub.add_parser("embed", help="Embed a text watermark into a PDB file.")
    embed.add_argument("input_pdb", help="Original PDB structure.")
    embed.add_argument("--text", required=True, help="Text watermark to embed.")
    embed.add_argument("-o", "--output", required=True, help="Output watermarked PDB path.")
    embed.add_argument("--start-freq", type=int, default=5, help="Starting Fourier frequency index.")
    embed.add_argument("--strength", type=float, default=0.02, help="Watermark modulation strength.")

    decode = sub.add_parser("decode", help="Decode a text watermark from a watermarked PDB file.")
    decode.add_argument("master_pdb", help="Original PDB structure.")
    decode.add_argument("query_pdb", help="Watermarked PDB structure.")
    decode.add_argument("--bits", type=int, default=4, help="Payload bit length to decode. Defaults to 4 bits.")
    decode.add_argument("--expected-text", help="Optional expected text for accuracy checking.")
    decode.add_argument("--start-freq", type=int, default=5, help="Starting Fourier frequency index.")

    return parser


def main(argv: Optional[List[str]] = None) -> None:
    args = build_parser().parse_args(argv)

    if args.command == "embed":
        result = embed_text(args.input_pdb, args.text, args.output, start_freq=args.start_freq, strength=args.strength)
        print(json.dumps(result.__dict__, indent=2))
        return

    if args.command == "decode":
        result = decode_text(
            args.master_pdb,
            args.query_pdb,
            args.bits,
            expected_text=args.expected_text,
            start_freq=args.start_freq,
        )
        print(json.dumps(result.__dict__, indent=2))
        return


if __name__ == "__main__":
    main()
