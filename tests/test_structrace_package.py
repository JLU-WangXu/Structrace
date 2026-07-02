from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from structrace.watermark import decode_text, embed_text


def test_watermark_text_roundtrip(tmp_path: Path) -> None:
    master = Path("Robustness/00_baseline_cases/6MRR/6MRR_original.pdb")
    output = tmp_path / "6MRR_watermarked.pdb"
    embed = embed_text(master, "npj SB", output)
    decoded = decode_text(master, output, expected_text="npj SB")
    assert embed.global_ca_rmsd > 0
    assert decoded.decoded_text == "npj SB"
    assert decoded.bit_accuracy == 1.0
    assert decoded.exact_recovery is True


def test_watermark_text_auto_decode_until_null(tmp_path: Path) -> None:
    master = Path("Robustness/00_baseline_cases/6MRR/6MRR_original.pdb")
    output = tmp_path / "6MRR_watermarked.pdb"
    embed_text(master, "npj SB", output)
    decoded = decode_text(master, output)
    assert decoded.decoded_text == "npj SB"


def test_cli_decode_bits_default_and_explicit(tmp_path: Path) -> None:
    master = Path("Robustness/00_baseline_cases/6MRR/6MRR_original.pdb")
    output = tmp_path / "6MRR_watermarked.pdb"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "structrace",
            "embed",
            str(master),
            "--text",
            "npj SB",
            "-o",
            str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    default_decode = subprocess.run(
        [sys.executable, "-m", "structrace", "decode", str(master), str(output)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert json.loads(default_decode.stdout)["bit_length"] == 4

    explicit_decode = subprocess.run(
        [sys.executable, "-m", "structrace", "decode", str(master), str(output), "--bits", "56"],
        check=True,
        capture_output=True,
        text=True,
    )
    explicit_payload = json.loads(explicit_decode.stdout)
    assert explicit_payload["bit_length"] == 56
    assert explicit_payload["decoded_text"] == "npj SB"
