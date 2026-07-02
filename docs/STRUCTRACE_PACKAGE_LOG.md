# StrucTrace package implementation log

## Summary

The installable Python package `structrace` exposes the core StrucTrace watermark workflow as two user-facing functions:

- `embed_text(input_pdb, text, output_pdb)`
- `decode_text(master_pdb, query_pdb)`

The command-line interface mirrors these two operations with `structrace embed` and `structrace decode`.

## Implemented package structure

```text
src/structrace/
  watermark/
    pdb.py        # PDB parsing, C-alpha extraction and coordinate replacement
    payload.py    # text <-> bit conversion
    codec.py      # FFT embedding and reference-guided decoding
  cli.py          # structrace embed/decode command-line entry point
```

Additional manuscript-validation modules may remain in the repository, but they are not part of the public quick-start API.

## Packaging

`pyproject.toml` defines a standard Python package:

```bash
pip install -e .
```

After installation, the CLI is available as:

```bash
python -m structrace --help
```

## Design decisions

- `embed_text()` embeds a UTF-8 text payload and appends a null terminator.
- `decode_text()` is reference-guided: the watermarked query PDB is decoded relative to the original master PDB.
- The CLI `decode --bits` option can be set explicitly; if omitted, it defaults to 4 bits.

## Validation performed

The package was smoke-tested with repository PDB examples:

- `embed_text()` and `decode_text()` round-trip on the 6MRR baseline case.
- CLI `decode` supports an explicit `--bits` value and a 4-bit default.
- CLI `python -m structrace embed` and `python -m structrace decode` round-trip.
