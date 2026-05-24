from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CASE = ROOT / "case_1BIL"
OUT = ROOT / "hardware_bound_access_audit.csv"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    encrypted = CASE / "1BIL_encrypted.pdbenc"
    encrypted_txt = CASE / "1BIL_encrypted.txt"
    decrypted = CASE / "1BIL_decrypted.pdb"
    decrypted_txt = CASE / "1BIL_decrypted_re.txt"

    rows = [
        {
            "tier": "Tier 2",
            "case_id": "1BIL",
            "artifact": path.name,
            "role": role,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path, role in (
            (encrypted, "encrypted_binary_asset"),
            (encrypted_txt, "encrypted_text_export"),
            (decrypted, "authorized_decrypted_pdb"),
            (decrypted_txt, "authorized_decrypted_text_export"),
        )
    ]

    rows.append(
        {
            "tier": "Tier 2",
            "case_id": "1BIL",
            "artifact": "integrity_check",
            "role": "decrypted_pdb_matches_text_export",
            "bytes": decrypted.stat().st_size,
            "sha256": str(sha256_file(decrypted) == sha256_file(decrypted_txt)),
        }
    )
    rows.append(
        {
            "tier": "Tier 2",
            "case_id": "1BIL",
            "artifact": "integrity_check",
            "role": "encrypted_pdbenc_matches_text_export",
            "bytes": encrypted.stat().st_size,
            "sha256": str(sha256_file(encrypted) == sha256_file(encrypted_txt)),
        }
    )

    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
