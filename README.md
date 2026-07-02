# StrucTrace: A Universal Fourier Watermark for Traceable Biomolecular Structures

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![GROMACS](https://img.shields.io/badge/GROMACS-2024.6-green.svg)](https://www.gromacs.org/)

**Authors:** <i>Xu Wang</i><sup>&dagger;</sup>, <i>Chi Wang</i><sup>&dagger;</sup>, <i>Tin-Yeh Huang</i>, <i>Yiquan Wang</i>, <i>Siyuan Jiang</i>, <i>Yafei Yuan</i><sup>&#42;</sup>  
<sup>&dagger;</sup>These authors contributed equally to this work.  
<sup>&#42;</sup>Corresponding authors.

StrucTrace is a reference-guided Fourier-domain watermarking framework for traceable biomolecular structures. It embeds provenance payloads into flexible C-alpha coordinate regions as a post-processing step, preserving atomic-scale structural fidelity while enabling public provenance verification, hardware-bound access control and auditable digital-rights management.

## Key Results

- Universal post-processing watermarking for deposited PDB files, cryo-EM structures and AI-designed proteins.
- 100% bit recovery in fidelity benchmarks, with side-chain RMSD no greater than 0.0015 Angstrom.
- High-capacity embedding of a 6,584-bit abstract into the human gamma-secretase cryo-EM structure.
- Defined recovery boundaries under rigid-body transformations, coordinate rounding, all-atom noise and local C-alpha perturbation.
- Three-tier safeguard model: public provenance, hardware-bound access and digital-rights management.

## Figures

![Figure 1: Universal Fourier-domain watermarking workflow](Figs/Figs/figure1.png)

**Figure 1. Universal Fourier-domain watermarking workflow.** StrucTrace selects flexible C-alpha atoms, embeds payload bits by mid-frequency FFT amplitude modulation and reconstructs watermarked coordinates with minimal geometric perturbation.

![Figure 2: Universality and structural integrity](Figs/Figs/figure2.png)

**Figure 2. Structural fidelity and universality.** Large-scale RMSD and Rosetta-energy analyses show negligible perturbation, and 50 ns MD simulations preserve dynamic behaviour across natural, Rosetta-designed and RFdiffusion-generated proteins.

![Figure 3: High-density information embedding](Figs/Figs/figure3.png)

**Figure 3. High-capacity embedding.** A 6,584-bit abstract was embedded into the human gamma-secretase tetramer while retaining coordinate-level agreement with the original cryo-EM structure.

![Figure 4: Watermark recovery under perturbation](Figs/Figs/figure4.png)

**Figure 4. Robustness and tamper boundary.** Reference-guided decoding is invariant to rigid-body transformations and standard PDB precision, while all-atom noise, local distortion and selected C-alpha perturbation define the recovery boundary.

![Figure 5: Multi-tiered security ecosystem](Figs/Figs/figure5.png)

**Figure 5. Three-tier safeguards.** Tier 1 provides public provenance verification, Tier 2 provides hardware-bound encrypted access, and Tier 3 provides auditable digital-rights management.

## Install

After PyPI publication, install StrucTrace with:

```bash
pip install structrace
```

For local development, clone and install the package in editable mode:

```bash
git clone https://github.com/JLU-WangXu/Structrace.git
cd Structrace
pip install -e .
```

The command-line interface is then available as:

```bash
python -m structrace --help
```

If your Python scripts directory is on `PATH`, you can also use:

```bash
structrace --help
```

External software used for manuscript validation includes DSSP, GROMACS 2024.6, Rosetta and Foldseek.

## Publishing

Build and validate the package locally:

```bash
python -m pip install --upgrade build twine
python -m build
python -m twine check dist/*
```

Upload to TestPyPI first:

```bash
python -m twine upload --repository testpypi dist/*
```

If the TestPyPI package installs and imports correctly, upload the same version to PyPI:

```bash
python -m twine upload dist/*
```

## Quick Start

### CLI

Embed a text watermark into a structure:

```bash
python -m structrace embed Robustness/00_baseline_cases/6MRR/6MRR_original.pdb \
  --text "npj SB" \
  -o tmp/6MRR_watermarked.pdb
```

Decode the watermark from the original and watermarked structures:

```bash
python -m structrace decode Robustness/00_baseline_cases/6MRR/6MRR_original.pdb \
  tmp/6MRR_watermarked.pdb \
  --bits 56
```

The `--bits` option can be set explicitly. If it is omitted, the CLI decodes 4 bits by default.

### Python API

```python
from structrace.watermark import decode_text, embed_text

master = "Robustness/00_baseline_cases/6MRR/6MRR_original.pdb"
query = "tmp/6MRR_watermarked.pdb"

embed_result = embed_text(master, "npj SB", query)
decode_result = decode_text(master, query)

print(embed_result.global_ca_rmsd)
print(decode_result.decoded_text)
```

## Package Functions

### Watermarking

- `embed_text(input_pdb, text, output_pdb)`: encode UTF-8 text as bits and embed it.
- `decode_text(master_pdb, query_pdb)`: recover the embedded UTF-8 text watermark from the original and watermarked structures.

## Repository Map

- `src/structrace/`: installable Python package.
- `Robustness/`: robustness and preliminary watermark-validation examples.
- `Security_Safeguards/`: Tier 1, Tier 2 and Tier 3 safeguard demonstrations.
- `High-density_embedding_stress_test/`: gamma-secretase high-capacity embedding case.
- `Molecular_Dynamics_Validation/`: MD inputs, structures and analysis outputs.
- `Bit_RMSD/` and `Bit_Rosetta_energy/`: large-scale fidelity benchmark tables.
- `Figs/Figs/`: revised manuscript figures.
- `docs/`: package implementation log and extended API examples.

## Citation

```bibtex
@article{Wang2026StrucTrace,
  title  = {StrucTrace: A Universal Fourier Watermark for Traceable Biomolecular Structures},
  author = {Wang, Xu and Wang, Chi and Huang, Tin-Yeh and Wang, Yiquan and Jiang, Siyuan and Yuan, Yafei},
  journal = {npj structural biology},
  year   = {2026},
  note   = {Accepted}
}
```
