from __future__ import annotations

from pathlib import Path

from structrace.security import LedgerEvent, build_ledger, decrypt_file, encrypt_file, verify_ledger
from structrace.watermark import decode_text, embed_text
from structrace.watermark.payload import text_to_bits


def test_watermark_text_roundtrip(tmp_path: Path) -> None:
    master = Path("Robustness/00_baseline_cases/6MRR/6MRR_original.pdb")
    output = tmp_path / "6MRR_watermarked.pdb"
    embed = embed_text(master, "npj SB", output)
    decoded = decode_text(master, output, len(text_to_bits("npj SB")), expected_text="npj SB")
    assert embed.global_ca_rmsd > 0
    assert decoded.decoded_text == "npj SB"
    assert decoded.bit_accuracy == 1.0
    assert decoded.exact_recovery is True


def test_tier2_encrypt_decrypt_roundtrip(tmp_path: Path) -> None:
    source = Path("Robustness/00_baseline_cases/6MRR/6MRR_original.pdb")
    encrypted = tmp_path / "asset.pdbenc"
    decrypted = tmp_path / "asset.pdb"
    encrypt_file(source, encrypted, machine_id="TEST_MACHINE")
    decrypt_file(encrypted, decrypted, machine_id="TEST_MACHINE")
    assert source.read_bytes() == decrypted.read_bytes()


def test_tier3_ledger_verification(tmp_path: Path) -> None:
    ledger = tmp_path / "ledger.json"
    build_ledger(
        [
            LedgerEvent("asset-001", "mint_asset_record", "owner_lab", "registry", "provenance_registration"),
            LedgerEvent("asset-001", "grant_license", "owner_lab", "institute_A", "view_and_verify"),
        ],
        output_json=ledger,
    )
    assert verify_ledger(ledger) is True
