from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Union


@dataclass(frozen=True)
class LedgerEvent:
    asset_id: str
    event: str
    actor: str
    counterparty: str
    permission: str
    timestamp_utc: Optional[str] = None


def _hash_record(record: Dict[str, str], previous_hash: str) -> str:
    payload = json.dumps({"previous_hash": previous_hash, **record}, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def build_ledger(events: List[LedgerEvent], output_json: Optional[Union[str, Path]] = None, output_csv: Optional[Union[str, Path]] = None) -> List[Dict[str, str]]:
    previous = "GENESIS"
    rows = []
    for index, event in enumerate(events, start=1):
        record = asdict(event)
        record["timestamp_utc"] = record["timestamp_utc"] or datetime.now(timezone.utc).isoformat()
        record["index"] = str(index)
        current = _hash_record(record, previous)
        rows.append({**record, "previous_hash": previous, "record_hash": current})
        previous = current
    if output_json is not None:
        Path(output_json).parent.mkdir(parents=True, exist_ok=True)
        Path(output_json).write_text(json.dumps(rows, indent=2), encoding="utf-8")
    if output_csv is not None:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        with Path(output_csv).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    return rows


def verify_ledger(ledger: Union[List[Dict[str, str]], str, Path]) -> bool:
    if not isinstance(ledger, list):
        ledger = json.loads(Path(ledger).read_text(encoding="utf-8"))
    previous = "GENESIS"
    for row in ledger:
        record = {key: value for key, value in row.items() if key not in {"previous_hash", "record_hash"}}
        if row["previous_hash"] != previous:
            return False
        if _hash_record(record, previous) != row["record_hash"]:
            return False
        previous = row["record_hash"]
    return True
