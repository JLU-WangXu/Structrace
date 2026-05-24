# StrucTrace: A universal Fourier watermark for traceable biomolecular structures

StrucTrace is a post-processing framework for embedding provenance information into three-dimensional biomolecular coordinate files. It was developed for traceable structural assets in generative protein design, cryo-EM and other high-value structural biology workflows, where authorship, controlled distribution and auditability must be added without compromising atomic-scale utility.

Unlike generation-coupled watermarking methods, StrucTrace does not require retraining or modifying a protein design model. It selects thermodynamically flexible C-alpha regions, modulates mid-frequency Fourier coefficients with a deterministic key-guided procedure, and reconstructs a watermarked coordinate file with minimal structural deviation. Decoding is reference-guided: a query structure is aligned to a registered master structure before differential Fourier-domain recovery of the embedded payload.

## Highlights

- Universal post-processing watermarking for existing PDB files, cryo-EM models and AI-designed structures.
- Deterministic Fourier-domain embedding in flexible C-alpha coordinate regions identified by DSSP and B-factor analysis.
- 100% bit recovery in fidelity benchmarks, with side-chain RMSD no greater than 0.0015 Angstrom.
- High-capacity embedding demonstrated by encoding a 6,584-bit research abstract into the human gamma-secretase cryo-EM structure.
- Defined operating regime under rigid-body transformations, coordinate rounding, all-atom noise and local C-alpha perturbations.
- Layered deployment model for public provenance, hardware-bound access control and digital rights management.

## Scientific Rationale

Generative protein design and high-resolution structure determination are making biomolecular structures increasingly valuable as digital intellectual property. Existing watermarking approaches are often tied to generation-time control and can introduce coordinate deviations on the Angstrom scale, which is problematic for docking, molecular dynamics and other precision-sensitive analyses.

StrucTrace addresses this gap by decoupling provenance encoding from structure generation. The method embeds information after a structure has been produced or experimentally determined, allowing deposited archives, unpublished models and proprietary structures to be marked retrospectively. The intended use is high-fidelity, reference-guided provenance verification rather than unrestricted blind watermark detection.

## Method Overview

![Figure 1: Universal Fourier-domain watermarking workflow](Figs/Figs/figure1.png)

**Figure 1. Universal Fourier-domain watermarking workflow of StrucTrace.** Cartesian coordinates are extracted from an input PDB structure, provenance metadata are encoded as a binary payload, and selected coordinate vectors are transformed with FFT. Thermodynamically permissive C-alpha atoms are selected from flexible regions using secondary-structure assignment and B-factor analysis. Payload bits are embedded by modulating mid-frequency amplitudes, preserving low-frequency global fold information and high-frequency local structural detail.

## Structural Fidelity and Universality

Across large-scale benchmarks, StrucTrace introduced only negligible coordinate perturbations. Median C-alpha RMSD values remained in the 0.0005 to 0.001 Angstrom range across payload densities, and folding free energy changes remained centered near zero. In comparisons with generation-coupled watermarking methods, StrucTrace achieved 100% bit accuracy while maintaining side-chain RMSD values at or below 0.0015 Angstrom.

The method was further evaluated on three representative structural asset classes: an experimentally resolved natural protein, a physics-based Rosetta design and an RFdiffusion-generated heme-binding protein. Fifty-nanosecond molecular dynamics simulations showed that watermarked structures retained backbone RMSD and residue-level RMSF profiles closely matching their original counterparts.

![Figure 2: Universality and structural integrity](Figs/Figs/figure2.png)

**Figure 2. Universality and consistent structural integrity across biomolecular asset classes.** Large-scale RMSD and folding-energy analyses indicate negligible geometric and thermodynamic perturbation. Molecular dynamics trajectories and representative structural views show close agreement between original and watermarked structures for 8HFE, 6MRR and 8VC8.

## High-Capacity Embedding

To test payload capacity under experimentally realistic constraints, StrucTrace was applied to the 2.6 Angstrom cryo-EM structure of the human gamma-secretase complex (PDB: 6IYC). A complete research abstract containing 6,584 bits was embedded across automatically identified high-B-factor C-alpha atoms. The resulting structure retained coordinate-level agreement with the original, with a global heavy-atom RMSD of 0.017 Angstrom and 100% bit recovery.

![Figure 3: High-density information embedding](Figs/Figs/figure3.png)

**Figure 3. High-density information embedding in a complex cryo-EM target.** The human gamma-secretase tetramer was encoded with a full-text abstract distributed across 4,400 high-B-factor C-alpha atoms. The watermarked and original structures remain visually indistinguishable, and the measured deviation is far below the experimental resolution.

## Decoding Robustness and Operating Boundary

StrucTrace is designed for reference-guided provenance verification. Under this workflow, rigid-body translations and rotations are normalized by alignment to the registered master structure and therefore preserve exact payload recovery. Standard three-decimal PDB coordinate precision also preserves exact recovery, whereas coarser rounding progressively reduces decoding accuracy.

All-atom coordinate noise and local perturbation experiments define a more nuanced recovery boundary. Bit accuracy remains high as coordinate noise approaches 0.002 Angstrom, but exact payload recovery requires sufficient coordinate precision. Local distortions can be tolerated in many cases, whereas direct perturbation of watermark-bearing C-alpha atoms substantially reduces recovery. StrucTrace should therefore be interpreted as a high-fidelity provenance watermark with a defined white-box tamper boundary.

![Figure 4: Watermark recovery under perturbation](Figs/Figs/figure4.png)

**Figure 4. Watermark recovery under coordinate perturbations.** Rigid-body transformations preserve recovery after reference-guided alignment. Coordinate rounding, all-atom noise, local distortions and targeted C-alpha perturbations define the practical decoding boundary of the watermark.

## Security Architecture

The coordinate-level watermark is designed to operate within a broader security architecture. Tier 1 supports public provenance verification by matching a query structure to a master database and recovering the embedded signature. Tier 2 adds hardware-bound access control through encrypted structural files and device-specific authorization. Tier 3 extends the framework to digital rights management, where structural assets can be linked to auditable licensing and transaction records.

![Figure 5: Multi-tiered security ecosystem](Figs/Figs/figure5.png)

**Figure 5. Technical workflows for the multi-tiered security ecosystem.** The architecture combines public verification, AES-256 hardware-bound access control and tokenized digital rights management for biomolecular structural assets.

## Repository Contents

- `Robustness/`: robustness analyses and preliminary watermark-validation examples.
- `Security_Safeguards/`: three-tier safeguard examples for public provenance, hardware-bound access and digital rights management.
- `High-density_embedding_stress_test/`: scripts and structures for the gamma-secretase high-capacity embedding experiment.
- `Molecular_Dynamics_Validation/`: molecular dynamics inputs, structures, trajectories and analysis outputs for 8HFE, 6MRR and 8VC8.
- `Rosetta_energy/`: Rosetta energy outputs used for thermodynamic stability analysis.
- `Figs/Figs/`: revised figures corresponding to the latest manuscript version.

## Installation

Clone the repository:

```bash
git clone https://github.com/JLU-WangXu/Structrace.git
cd Structrace
```

Create a Python environment and install the core Python dependencies used by the analysis scripts:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install numpy pandas scipy biopython requests tqdm
```

Some validation workflows require external scientific software:

- DSSP for secondary-structure assignment.
- GROMACS 2024.6 for molecular dynamics validation.
- Rosetta for folding-energy calculations.
- Foldseek for fast reference retrieval in public provenance workflows.

## Reproducibility

The repository provides representative data and scripts for the major analyses described above:

- High-capacity embedding: `High-density_embedding_stress_test/run_watermark_6iyc.py`
- MD protocol: `Molecular_Dynamics_Validation/configs/protocol.md`
- GROMACS parameter files: `Molecular_Dynamics_Validation/configs/`
- Rosetta energy summaries: `Rosetta_energy/`
- Preliminary watermark examples: `Robustness/00_preliminary_watermark_validation/`
- Three-tier security examples: `Security_Safeguards/`

## Citation

If you use StrucTrace in academic work, please cite:

```bibtex
@article{Wang2025StrucTrace,
  title   = {StrucTrace: A universal Fourier watermark for traceable biomolecular structures},
  author  = {Wang, Xu and Wang, Chi and Huang, Tin-Yeh and Wang, Yiquan and Jiang, Siyuan and Yuan, Yafei},
  year    = {2026},
  note    = {Manuscript under review}
}
```

## Code Availability

The StrucTrace package and analysis scripts are provided for academic and noncommercial research use through this repository. For licensing, commercial use or controlled-access deployment, please contact the corresponding authors listed in the manuscript.
