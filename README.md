# Structrace
# StrucTrace: A Universal Fourier Watermark for Traceable Biomolecular Structures

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/release/python-380/)
[![GROMACS](https://img.shields.io/badge/GROMACS-2024.6-green.svg)](https://www.gromacs.org/)

**Authors:** Xu Wang¹, Chi Wang¹, Tin-Yeh Huang¹, Yiquan Wang¹, Yafei Yuan*¹
**Affiliation:** ¹School of Life Sciences, Tsinghua University

---

## 📖 Overview

[cite_start]**StrucTrace** is a universal, post-processing watermarking framework designed to transform biomolecular structures (Cryo-EM, AlphaFold, RFdiffusion, etc.) into traceable, secure, and auditable digital assets [cite: 76-84].

As AI-generated proteins become indistinguishable from natural ones, protecting digital biomolecular intellectual property (IP) is critical. [cite_start]Unlike generation-coupled methods that introduce notable structural deviations ($RMSD \approx 2.0 \AA$), **StrucTrace** leverages **Fourier-domain geometric modulation** to embed data [cite: 85-86]. [cite_start]It achieves **100% bit recovery** while maintaining **functionally imperceptible structural perturbations** ($scRMSD \le 0.0015 \AA$), ensuring atomic-level precision for downstream applications like molecular dynamics and docking[cite: 87, 107].

### Framework Workflow
The algorithm utilizes a deterministic workflow that decouples protection from generation. [cite_start]It targets thermodynamically flexible regions (identified via DSSP and B-factor analysis) and modulates mid-frequency coefficients to distribute information globally without affecting the core fold [cite: 115-119].

![Figure 1: Framework Overview](Figs/Figs/Fig1.png)
> **Figure 1: The universal Fourier-domain watermarking framework.** (a) Schematic of the complete pipeline. (b) Deterministic atom selection integrates DSSP and B-factor analysis to isolate flexible $C_{\alpha}$ atoms (red). (c) [cite_start]FFT-based encoding modulates mid-frequency bands to embed data while preserving global structural integrity [cite: 210-216].

---

## 🚀 Key Features

* [cite_start]**Universal Compatibility**: Validated on diverse asset classes: Natural proteins (e.g., 8HFE), Physics-based designs (Rosetta), and AI-generated structures (RFdiffusion)[cite: 126].
* [cite_start]**Zero-Distortion**: Median $C_\alpha$-RMSD $< 0.004 \AA$ and side-chain RMSD $\le 0.0015 \AA$, effectively invisible to standard structural analysis tools[cite: 119, 123].
* [cite_start]**Thermodynamic Stability**: Watermarked structures exhibit thermodynamic neutrality ($\Delta\Delta G \approx 0$)[cite: 125].
* [cite_start]**High Capacity**: Capable of embedding extensive metadata (e.g., >6,000 bits) into complex structures without compromising experimental resolution[cite: 130].
* [cite_start]**Multi-Tiered Security**: Integrates public provenance tracking, hardware-bound access control, and blockchain-based DRM[cite: 135].

---

## 📊 Universality & Stability Validation

[cite_start]To validate the universality of StrucTrace across distinct eras of structural biology, we performed **50 ns Molecular Dynamics (MD) simulations** on three representative test cases [cite: 384-387]:
1.  **Natural Protein:** Human Norepinephrine Transporter (PDB: **8HFE**).
2.  **Physics-based Design:** *De novo* protein designed via Rosetta (PDB: **6MRR**).
3.  **AI-Generated:** Heme-binder generated via RFdiffusion (PDB: **8VC8**).

![Figure 2: MD Validation](Figs/Figs/Fig2.jpg)
> **Figure 2: Universality and consistent structural integrity.** (a) Global RMSD distributions showing negligible perturbation ($<0.0015\AA$). (b) Thermodynamic neutrality ($\Delta\Delta G \approx 0$). (c-e) [cite_start]Comparative 50 ns MD trajectories confirming that watermarked structures (blue) retain dynamic stability statistically indistinguishable from original structures (green) [cite: 295-303].

---

## 🧬 High-Capacity Embedding (Stress Test)

StrucTrace can transform PDB files into self-contained data capsules. [cite_start]We demonstrated this by embedding a **complete research abstract (6,584 bits)** into the 2.6 Å Cryo-EM structure of the human $\gamma$-secretase tetramer (PDB: 6IYC) [cite: 389-391].

![Figure 3: High Capacity](Figs/Figs/Fig3.jpg)
> **Figure 3: High-density information embedding in a complex Cryo-EM target.** The payload (6,584 bits) was distributed across 4,400 high B-factor atoms. [cite_start]The resulting structure maintained a global RMSD of **0.017 Å**, two orders of magnitude smaller than the experimental resolution, with 100% bit recovery [cite: 331-335].

---

## 🛡️ Multi-Tiered Security Ecosystem

[cite_start]StrucTrace establishes a scalable infrastructure for bio-asset management through a three-tiered architecture [cite: 135-138]:

* [cite_start]**Tier 1 (Public Provenance):** Lightweight signature verification using high-speed alignment tools (Foldseek)[cite: 395].
* [cite_start]**Tier 2 (Hardware-Bound Access):** AES-256 encryption bound to physical device identifiers ($ID_{machine}$) to prevent unauthorized exfiltration [cite: 397-402].
* [cite_start]**Tier 3 (Digital Rights Management):** Blockchain-based tokenization for commercial licensing and immutable audit trails [cite: 405-407].

![Figure 4: Security Architecture](Figs/Figs/Fig4.png)
> **Figure 4: Technical workflows for the multi-tiered security ecosystem.** Top: Public verification workflow. Middle: Hardware-bound mechanism (Tier 2) validating unique device IDs. [cite_start]Bottom: Tokenized transaction flow on the blockchain (Tier 3) [cite: 376-381].

---

## 🔬 Reproducibility: GROMACS MD Protocol

To ensure reproducibility of our stability results, we provide the exact molecular dynamics protocol used in the paper.

[cite_start]**Simulation Parameters [cite: 311-318]:**
* **Software:** GROMACS v2024.6
* **Force Field:** AMBER99SB-ILDN
* **Water Model:** TIP3P (Cubic box)
* **Ions:** Neutralized with 0.15 M NaCl ($Na^+/Cl^-$)
* **Ensembles:** NVT (V-rescale) & NPT (Parrinello-Rahman)
* **Duration:** 50 ns production run per system

### Command Line Workflow
```bash
# 1. Topology Generation
gmx pdb2gmx -f 8HFE_watermarked.pdb -o protein.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh

# 2. Box Definition & Solvation
gmx editconf -f protein.gro -o protein_box.gro -c -d 1.0 -bt cubic
gmx solvate -cp protein_box.gro -cs spc216.gro -o protein_SOL.gro -p topol.top

# 3. Ionization (Physiological conditions: 0.15 M)
gmx grompp -f em.mdp -c protein_SOL.gro -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -p topol.top -o system.gro -pname NA -nname CL -neutral -conc 0.15

# 4. Energy Minimization (Steepest Descent)
gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -ntmpi 1 -nt 50

# 5. Equilibration (NVT & NPT)
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -ntmpi 1 -nt 50
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt -ntmpi 1 -nt 50

# 6. Production Run (50 ns)
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
gmx mdrun -v -deffnm md -ntmpi 1 -nt 50
