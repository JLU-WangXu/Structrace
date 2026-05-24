from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
LEDGER_JSON = ROOT / "drm_ledger.json"
LEDGER_CSV = ROOT / "drm_ledger.csv"


def hash_record(record: dict[str, str], previous_hash: str) -> str:
    payload = json.dumps(
        {"previous_hash": previous_hash, **record},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def main() -> None:
    events = [
        {
            "timestamp_utc": "2026-05-23T00:00:00Z",
            "asset_id": "StrucTrace-6IYC-abstract-watermark",
            "event": "mint_asset_record",
            "actor": "owner_lab",
            "counterparty": "registry",
            "permission": "provenance_registration",
        },
        {
            "timestamp_utc": "2026-05-23T00:05:00Z",
            "asset_id": "StrucTrace-6IYC-abstract-watermark",
            "event": "grant_license",
            "actor": "owner_lab",
            "counterparty": "authorized_institute",
            "permission": "view_and_verify",
        },
        {
            "timestamp_utc": "2026-05-23T00:10:00Z",
            "asset_id": "StrucTrace-6IYC-abstract-watermark",
            "event": "revoke_license",
            "actor": "owner_lab",
            "counterparty": "expired_client",
            "permission": "access_revoked",
        },
    ]

    previous = "GENESIS"
    ledger = []
    for index, event in enumerate(events, start=1):
        record = {"index": str(index), **event, "recorded_at_utc": datetime.now(timezone.utc).isoformat()}
        current = hash_record(record, previous)
        ledger.append({**record, "previous_hash": previous, "record_hash": current})
        previous = current

    LEDGER_JSON.write_text(json.dumps(ledger, indent=2), encoding="utf-8")
    with LEDGER_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ledger[0].keys()))
        writer.writeheader()
        writer.writerows(ledger)


if __name__ == "__main__":
    main()
