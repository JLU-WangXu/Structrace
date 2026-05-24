from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Union

from .common import sha256_file


@dataclass(frozen=True)
class ProvenanceRecord:
    case_id: str
    master_pdb: str
    query_pdb: str
    master_sha256: str
    query_sha256: str
    payload: str
    verification_mode: str = "reference-guided"


def build_public_registry(records: List[Tuple[str, Union[str, Path], Union[str, Path], str]], output_csv: Union[str, Path]) -> List[ProvenanceRecord]:
    registry = [
        ProvenanceRecord(
            case_id=case_id,
            master_pdb=str(master_pdb),
            query_pdb=str(query_pdb),
            master_sha256=sha256_file(master_pdb),
            query_sha256=sha256_file(query_pdb),
            payload=payload,
        )
        for case_id, master_pdb, query_pdb, payload in records
    ]
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    with Path(output_csv).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(registry[0]).keys()))
        writer.writeheader()
        writer.writerows(asdict(record) for record in registry)
    return registry


def verify_public_registry(registry_csv: Union[str, Path]) -> List[Dict[str, Union[str, bool]]]:
    results = []
    with Path(registry_csv).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            query_hash = sha256_file(row["query_pdb"])
            results.append(
                {
                    "case_id": row["case_id"],
                    "query_sha256": query_hash,
                    "registry_sha256": row["query_sha256"],
                    "integrity_match": query_hash == row["query_sha256"],
                    "payload": row["payload"],
                    "decision": "accepted_for_reference_guided_decoding" if query_hash == row["query_sha256"] else "rejected_integrity_mismatch",
                }
            )
    return results
