# Supplementary Information: Molecular Dynamics Simulation Protocol

**Study:** StrucTrace: A universal Fourier watermark for traceable biomolecular structures  
**Software:** GROMACS v2024.6 [1]  
**Force Field:** AMBER99SB-ILDN [2]  
**Water Model:** TIP3P [3]

---

## 1. System Setup and Topology Generation

Molecular dynamics (MD) simulations were performed to validate the dynamic stability and structural integrity of the watermarked biomolecules. [cite_start]The simulation systems were prepared using the **GROMACS v2024.6** software package.

### 1.1. Topology Preparation
[cite_start]The atomic coordinates for the target structures—including the natural membrane transporter (**8HFE**), the de novo designed protein (**6MRR**), and the AI-generated binder (**8VC8**)—were processed using the `pdb2gmx` module.
* [cite_start]**Force Field:** The **AMBER99SB-ILDN** force field was employed to describe protein parameters.
* [cite_start]**Water Model:** The **TIP3P** water model was selected for solvation.

### 1.2. Solvation and Ionization
* [cite_start]**Simulation Box:** Structures were centered in a **cubic** box with periodic boundary conditions applied in all directions.
* [cite_start]**Neutralization:** The system charge was neutralized by the addition of sodium ($Na^+$) and chloride ($Cl^-$) ions to simulate physiological conditions.

## 2. Energy Minimization

[cite_start]To remove steric clashes and relax the high-energy contacts inherent in the initial coordinates, the system underwent energy minimization using the **Steepest Descent** algorithm. Minimization was considered converged when the maximum force ($F_{max}$) dropped below the defined threshold (default: $1000 \, kJ \cdot mol^{-1} \cdot nm^{-1}$).

## 3. Equilibration Protocol

[cite_start]Following minimization, the systems were equilibrated in two phases to stabilize temperature and pressure.

### 3.1. NVT Equilibration (Canonical Ensemble)
The system was heated to the target temperature (300 K) under a constant number of particles, volume, and temperature (NVT).
* **Thermostat:** A modified Berendsen (V-rescale) thermostat was typically used to maintain temperature.
* **Restraints:** Position restraints were applied to the protein backbone to allow solvent relaxation around the solute.

### 3.2. NPT Equilibration (Isothermal-Isobaric Ensemble)
The system was subsequently equilibrated under constant number of particles, pressure, and temperature (NPT) to stabilize the system density.
* **Barostat:** Pressure coupling was applied to maintain a reference pressure of 1 bar.
* **Duration:** Sufficient equilibration time was allocated to ensure density convergence.

## 4. Production Molecular Dynamics

Unrestrained production MD runs were conducted for both the original (wild-type) and watermarked structures to allow for direct comparative analysis.

* [cite_start]**Duration:** **50 ns** per trajectory for each system.
* **Integration:** The equations of motion were integrated using the leap-frog algorithm with a time step of 2 fs.
* **Constraints:** All bonds involving hydrogen atoms were constrained using the LINCS algorithm.
* **Electrostatics:** Long-range electrostatic interactions were treated using the Particle Mesh Ewald (PME) method.

## 5. Trajectory Analysis

[cite_start]The resulting trajectories were analyzed to quantify structural deviations and flexibility changes.

* **Global Stability (RMSD):** The Backbone Root Mean Square Deviation (RMSD) was calculated relative to the starting structure to assess global structural stability over the 50 ns simulation.
* **Local Flexibility (RMSF):** Per-residue Root Mean Square Fluctuation (RMSF) was computed to localize potential changes in flexibility induced by the watermarking process.

---

## Appendix: Reproducible Command Line Workflow

The following commands outline the exact workflow used for the simulations (example using 8HFE).

```bash
# 1. Topology Generation
# Force Field: AMBER99SB-ILDN, Water: TIP3P
gmx pdb2gmx -f 8HFE_watermarked.pdb -o protein.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh

# 2. Box Definition
gmx editconf -f protein.gro -o protein_box.gro -c -d 1.0 -bt cubic

# 3. Solvation
gmx solvate -cp protein_box.gro -cs spc216.gro -o protein_SOL.gro -p topol.top

# 4. Ionization (Physiological conditions)
gmx grompp -f em.mdp -c protein_SOL.gro -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -p topol.top -o system.gro -pname NA -nname CL -neutral -conc 0.15

# 5. Energy Minimization
gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -ntmpi 1 -nt 50

# 6. Equilibration (NVT & NPT)
# NVT
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -ntmpi 1 -nt 50
# NPT
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt -ntmpi 1 -nt 50

# 7. Production Run (50 ns)
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
gmx mdrun -v -deffnm md -ntmpi 1 -nt 50

# 8. Analysis
# Backbone RMSD
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns
# Per-residue RMSF
gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res