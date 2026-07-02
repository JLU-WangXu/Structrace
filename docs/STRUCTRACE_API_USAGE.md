# StrucTrace API and CLI usage

StrucTrace exposes two user-facing watermark functions:

- `embed_text(input_pdb, text, output_pdb)`: add a text watermark to a PDB structure.
- `decode_text(master_pdb, query_pdb)`: decode the text watermark from the original and watermarked structures.

## Install

From the repository root:

```bash
pip install -e .
```

After publication to PyPI, users can install with:

```bash
pip install structrace
```

## Python API

### Embed a watermark

```python
from structrace.watermark import embed_text

result = embed_text(
    "input.pdb",
    "07022026CHIWANGTEST",
    "input_watermarked.pdb",
)

print(result.output_pdb)
print(result.global_ca_rmsd)
```

### Decode a watermark

```python
from structrace.watermark import decode_text

result = decode_text(
    "input.pdb",
    "input_watermarked.pdb",
)

print(result.decoded_text)
```

## CLI

### Embed a watermark

```bash
python -m structrace embed input.pdb \
  --text "07022026CHIWANGTEST" \
  -o input_watermarked.pdb
```

### Decode a watermark

```bash
python -m structrace decode input.pdb input_watermarked.pdb --bits 160
```

The `--bits` option can be set explicitly. If it is omitted, the CLI decodes 4 bits by default.
