# StrucTrace: A Universal Fourier Watermark for Traceable Biomolecular Structures

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![PyPI](https://img.shields.io/pypi/v/structrace.svg)](https://pypi.org/project/structrace/)
[![GROMACS](https://img.shields.io/badge/GROMACS-2024.6-green.svg)](https://www.gromacs.org/)

**Authors:** <i>Xu Wang</i><sup>&dagger;</sup>, <i>Chi Wang</i><sup>&dagger;</sup>, <i>Tin-Yeh Huang</i>, <i>Yiquan Wang</i>, <i>Siyuan Jiang</i>, <i>Yafei Yuan</i><sup>&#42;</sup>  
<sup>&dagger;</sup>These authors contributed equally to this work.  
<sup>&#42;</sup>Corresponding authors.

StrucTrace is a reference-guided Fourier-domain watermarking framework for traceable biomolecular structures. It embeds provenance payloads into flexible C-alpha coordinate regions as a post-processing step, preserving atomic-scale structural fidelity while enabling public provenance verification, hardware-bound access control and auditable digital-rights management.

## Install

Install the released package from PyPI:

```bash
pip install structrace
```

For local development, clone the repository and install in editable mode:

```bash
git clone https://github.com/JLU-WangXu/Structrace.git
cd Structrace
pip install -e .
```

Check the command-line interface:

```bash
python -m structrace --help
```

## Quick Start

### Add a watermark

```bash
python -m structrace embed input.pdb \
  --text "07022026CHIWANGTEST" \
  -o input_watermarked.pdb
```

This creates a new PDB file and leaves the original structure unchanged. The embedded text can contain English, Chinese, numbers and common symbols because StrucTrace encodes text as UTF-8 bytes.

### Decode a watermark

Decoding is reference-guided: provide the original structure and the watermarked structure.

```bash
python -m structrace decode input.pdb input_watermarked.pdb --bits 160
```

For the example text `07022026CHIWANGTEST`, use `160` bits because the payload has 19 ASCII characters plus the default null terminator:

```text
(19 + 1) * 8 = 160 bits
```

If `--bits` is omitted, the CLI decodes 4 bits by default. This is useful for quick bit-level inspection, but full text recovery should use the correct payload bit length.

### Python API

```python
from structrace.watermark import decode_text, embed_text

input_pdb = "input.pdb"
watermarked_pdb = "input_watermarked.pdb"
payload = "07022026CHIWANGTEST"

embed_result = embed_text(input_pdb, payload, watermarked_pdb)
decode_result = decode_text(input_pdb, watermarked_pdb)

print(embed_result.global_ca_rmsd)
print(decode_result.decoded_text)
```

The Python API writes StrucTrace payload-length metadata into the watermarked PDB and can usually recover the full text with `decode_text(input_pdb, watermarked_pdb)`. For legacy files or manually controlled decoding, pass an explicit bit length:

```python
decode_result = decode_text(input_pdb, watermarked_pdb, bit_length=160)
```

## Usage Guidance

### Choosing watermark text

- Use short identifiers for routine provenance, such as dates, lab IDs, manuscript IDs or accession tags.
- Prefer stable payloads such as `07022026CHIWANGTEST`, `LAB_ASSET_0001` or `project:sample:version`.
- Chinese text is supported, but it uses more bits because UTF-8 Chinese characters usually require 3 bytes, or 24 bits, per character.
- Avoid very long text for small proteins. Longer payloads require more selected C-alpha atoms.

### Calculating `--bits`

The CLI `decode` command needs the payload bit length for full recovery:

```text
bits = (len(payload.encode("utf-8")) + 1) * 8
```

The extra `+ 1` is the null terminator added by `embed_text()`.

Examples:

```text
"npj SB"                  -> (6 + 1) * 8  = 56 bits
"07022026CHIWANGTEST"    -> (19 + 1) * 8 = 160 bits
"zhongwen"               -> (8 + 1) * 8  = 72 bits
```

You can calculate it in Python:

```bash
python -c "s='07022026CHIWANGTEST'; print((len(s.encode('utf-8')) + 1) * 8)"
```

### Recommended structure inputs

- Use PDB files containing protein `ATOM` records with C-alpha (`CA`) atoms.
- Keep residue numbering, chain IDs and C-alpha records consistent between the original and watermarked files.
- Decode against the exact original structure used for embedding.
- Multi-chain protein structures are supported if the PDB keeps standard chain and residue fields.
- Structures with more C-alpha atoms can carry longer payloads more comfortably.

### Output interpretation

`embed` prints metrics such as:

- `bit_length`: number of embedded payload bits.
- `top_n_selected_ca`: number of selected C-alpha atoms used by the Fourier embedding.
- `global_ca_rmsd`: global C-alpha RMSD between original and watermarked structures.
- `selected_ca_rmsd`: RMSD on selected C-alpha atoms.
- `max_selected_ca_displacement`: largest selected-atom displacement.

For a successful text decode, check:

```json
"decoded_text": "your payload"
```

If `--expected-text` is supplied, the output also reports `bit_accuracy` and `exact_recovery`.

```bash
python -m structrace decode input.pdb input_watermarked.pdb \
  --bits 160 \
  --expected-text "07022026CHIWANGTEST"
```

### Practical recommendations

- Keep an archived copy of the original PDB. Decoding is reference-guided and needs the original structure.
- Use explicit `--bits` for CLI text recovery.
- For public examples, use short ASCII payloads so readers can reproduce the result easily.
- For non-English payloads, compute `--bits` using UTF-8 byte length, not character count.
- Do not edit residue IDs, chain IDs or C-alpha coordinates by unrelated tools before decoding.
- Standard PDB coordinate precision is supported by the workflow, but aggressive coordinate noise, atom deletion or residue renumbering can reduce recovery accuracy.

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

External software used for manuscript validation includes DSSP, GROMACS 2024.6, Rosetta and Foldseek.

## Publishing

See [docs/PYPI_RELEASE.md](docs/PYPI_RELEASE.md) for the package release checklist.

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
