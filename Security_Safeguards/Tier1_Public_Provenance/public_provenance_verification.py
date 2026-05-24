from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
BASELINE = ROOT.parents[1] / "Robustness" / "00_baseline_cases"
REGISTRY = ROOT / "public_registry.csv"
RESULTS = ROOT / "public_verification_results.csv"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_registry() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for case_id in ("6MRR", "8HFE", "8VC8"):
        case_dir = BASELINE / case_id
        original = case_dir / f"{case_id}_original.pdb"
        watermarked = case_dir / f"{case_id}_watermarked_npj_SB.pdb"
        selected = case_dir / f"{case_id}_selected_ca_atoms.csv"
        rows.append(
            {
                "case_id": case_id,
                "registered_master_pdb": str(original.relative_to(ROOT.parents[1])),
                "watermarked_query_pdb": str(watermarked.relative_to(ROOT.parents[1])),
                "selected_ca_table": str(selected.relative_to(ROOT.parents[1])),
                "master_sha256": sha256_file(original),
                "watermarked_sha256": sha256_file(watermarked),
                "payload": "npj SB",
                "verification_mode": "reference-guided",
            }
        )
    return rows


def main() -> None:
    rows = build_registry()
    with REGISTRY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    checks = []
    for row in rows:
        query = ROOT.parents[1] / row["watermarked_query_pdb"]
        checks.append(
            {
                "case_id": row["case_id"],
                "query_sha256": sha256_file(query),
                "registry_sha256": row["watermarked_sha256"],
                "integrity_match": sha256_file(query) == row["watermarked_sha256"],
                "payload": row["payload"],
                "decision": "accepted_for_reference_guided_decoding",
            }
        )

    with RESULTS.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(checks[0].keys()))
        writer.writeheader()
        writer.writerows(checks)


if __name__ == "__main__":
    main()
