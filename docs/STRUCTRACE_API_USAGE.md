# StrucTrace API and CLI usage

## Install

From the repository root:

```bash
pip install -e .
```

## Python API

### Embed text into a PDB

```python
from structrace.watermark import embed_text

result = embed_text(
    "Robustness/00_baseline_cases/6MRR/6MRR_original.pdb",
    "npj SB",
    "tmp/6MRR_watermarked.pdb",
)
print(result.global_ca_rmsd)
```

### Decode text from a watermarked PDB

```python
from structrace.watermark import decode_text
from structrace.watermark.payload import text_to_bits

bits = len(text_to_bits("npj SB"))
result = decode_text(
    "Robustness/00_baseline_cases/6MRR/6MRR_original.pdb",
    "tmp/6MRR_watermarked.pdb",
    bits,
    expected_text="npj SB",
)
print(result.decoded_text, result.bit_accuracy, result.exact_recovery)
```

### Apply robustness perturbations

```python
from structrace.robustness import round_pdb_coordinates, translate_pdb, add_gaussian_noise_to_atoms

round_pdb_coordinates("tmp/6MRR_watermarked.pdb", "tmp/6MRR_round3.pdb", decimals=3)
translate_pdb("tmp/6MRR_watermarked.pdb", "tmp/6MRR_tx10.pdb", (10.0, 0.0, 0.0))
add_gaussian_noise_to_atoms("tmp/6MRR_watermarked.pdb", "tmp/6MRR_noise.pdb", sigma=0.001, seed=1)
```

### Tier 1: public provenance

```python
from structrace.security import build_public_registry, verify_public_registry

build_public_registry(
    [
        (
            "6MRR",
            "Robustness/00_baseline_cases/6MRR/6MRR_original.pdb",
            "Robustness/00_baseline_cases/6MRR/6MRR_watermarked_npj_SB.pdb",
            "npj SB",
        )
    ],
    "tmp/public_registry.csv",
)
print(verify_public_registry("tmp/public_registry.csv"))
```

### Tier 2: hardware-bound access

```python
from structrace.security import encrypt_file, decrypt_file, machine_fingerprint

machine_id = machine_fingerprint()
encrypt_file("input.pdb", "input.pdbenc", machine_id=machine_id)
decrypt_file("input.pdbenc", "input.decrypted.pdb", machine_id=machine_id)
```

### Tier 3: DRM ledger

```python
from structrace.security import LedgerEvent, build_ledger, verify_ledger

events = [
    LedgerEvent("asset-001", "mint_asset_record", "owner_lab", "registry", "provenance_registration"),
    LedgerEvent("asset-001", "grant_license", "owner_lab", "institute_A", "view_and_verify"),
]
build_ledger(events, output_json="tmp/ledger.json", output_csv="tmp/ledger.csv")
print(verify_ledger("tmp/ledger.json"))
```

## CLI examples

### Embed and decode

```bash
structrace embed Robustness/00_baseline_cases/6MRR/6MRR_original.pdb --text "npj SB" -o tmp/6MRR_cli_watermarked.pdb
structrace decode Robustness/00_baseline_cases/6MRR/6MRR_original.pdb tmp/6MRR_cli_watermarked.pdb --bits 56 --expected-text "npj SB"
```

### Robustness perturbations

```bash
structrace attack round tmp/6MRR_cli_watermarked.pdb --decimals 3 -o tmp/6MRR_round3.pdb
structrace attack translate tmp/6MRR_cli_watermarked.pdb --vector 10 0 0 -o tmp/6MRR_tx10.pdb
structrace attack noise tmp/6MRR_cli_watermarked.pdb --sigma 0.001 --seed 1 -o tmp/6MRR_noise.pdb
```

### Tier 1

```bash
structrace tier1 build-registry \
  --record 6MRR Robustness/00_baseline_cases/6MRR/6MRR_original.pdb Robustness/00_baseline_cases/6MRR/6MRR_watermarked_npj_SB.pdb "npj SB" \
  -o tmp/public_registry.csv

structrace tier1 verify tmp/public_registry.csv -o tmp/public_verification.csv
```

### Tier 2

```bash
structrace tier2 fingerprint
structrace tier2 encrypt input.pdb -o input.pdbenc --machine-id DEMO_MACHINE
structrace tier2 decrypt input.pdbenc -o input.decrypted.pdb --machine-id DEMO_MACHINE
structrace tier2 audit input.pdb input.decrypted.pdb -o tmp/audit.csv
```

### Tier 3

```bash
structrace tier3 build-ledger \
  --event asset-001 mint_asset_record owner_lab registry provenance_registration \
  --event asset-001 grant_license owner_lab institute_A view_and_verify \
  --json tmp/ledger.json --csv tmp/ledger.csv

structrace tier3 verify tmp/ledger.json
```
