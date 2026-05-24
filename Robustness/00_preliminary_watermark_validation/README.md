# Preliminary watermark validation

This directory preserves the early, small-scale StrucTrace watermarking tests that were previously stored under `Watermark/`. The contents are kept as historical algorithm-validation examples and are separate from the main robustness operating-regime experiments.

## Contents

- `scripts/legacy_watermark_fft.py`: legacy FFT-based watermark embedding and decoding script.
- `natural_proteins/originals/`: original natural-protein PDB examples.
- `natural_proteins/watermarked/`: corresponding early watermarked natural-protein examples when available.
- `designed_proteins/rosetta/`: original and watermarked Rosetta-designed protein examples.
- `designed_proteins/rfdiffusion/`: original and watermarked RFdiffusion-related protein examples.

## Scope

Only core method and result files were retained:

- PDB inputs and watermarked PDB outputs.
- Text recovery/result records.
- The legacy watermarking script.

Intermediate screenshots, binding `.npz` files and `.DS_Store` metadata were intentionally excluded because they do not affect the watermarking method or the reported reliability conclusions.
