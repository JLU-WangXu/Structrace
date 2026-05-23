# Robustness

This directory contains the core robustness experiments supporting the StrucTrace decoding operating regime. The original workspace at `E:\Watermarker\Attack` included many generated attacked PDB files, intermediate plotting variants and planning records. Those files were intentionally excluded here to keep the repository focused and auditable. The retained files are the baseline structures, selected C-alpha tables, experiment scripts, final result tables, summaries, statistics and figure outputs needed to understand and reproduce the conclusions.

## Layout

- `00_baseline_cases/`: original and watermarked PDB files for 6MRR, 8HFE and 8VC8, plus selected C-alpha atom tables and embedding summary.
- `01_rigid_body_transformations/`: translation, rotation and combined rotation-translation tests.
- `02_coordinate_precision_rounding/`: coordinate decimal precision tests.
- `03_low_amplitude_all_atom_noise/`: fine all-atom Gaussian noise sweep near the recovery boundary.
- `04_local_coordinate_distortion/`: local C-alpha-centered coordinate distortion tests.
- `05_selected_vs_nonselected_ca_perturbation/`: selected versus non-selected C-alpha perturbation tests.
- `06_integrated_figures/`: assembled robustness figures and captions.

## Core Experiments

### 1. Rigid-body transformations

Tests whether reference-guided decoding is invariant to global coordinate-frame changes.

Retained outputs include `global_translation_results.csv`, `global_rotation_2angle_results.csv` and `global_rotation_translation_results.csv`. Across the tested transformations, payload recovery remained exact with 100% bit accuracy after alignment to the registered master structure.

### 2. Coordinate precision rounding

Tests whether the watermark survives reduced coordinate precision. Standard three-decimal PDB precision preserved exact recovery, whereas coarser rounding progressively reduced bit accuracy.

Key file: `02_coordinate_precision_rounding/results/coordinate_decimal_rounding_results.csv`.

### 3. Low-amplitude all-atom noise

Defines the fine recovery boundary under global stochastic coordinate noise. Exact recovery remained complete at 0.001 Angstrom noise and declined as noise approached 0.002 Angstrom, while bit accuracy remained high.

Key file: `03_low_amplitude_all_atom_noise/results/gaussian_noise_all_atoms_sigma_0.001_to_0.002_results.csv`.

### 4. Local coordinate distortion

Tests local C-alpha-centered perturbation regions with Gaussian coordinate noise. These experiments quantify structure-dependent bit-accuracy distributions under local distortions.

Key file: `04_local_coordinate_distortion/results/local_distortion_sigma_results.csv`.

### 5. Selected versus non-selected C-alpha perturbation

Tests the white-box tamper boundary by perturbing watermark-bearing selected C-alpha atoms versus the same number of non-selected C-alpha atoms. Perturbing non-selected atoms preserved recovery, whereas perturbing selected atoms reduced bit accuracy to approximately 55-60%.

Key files: `05_selected_vs_nonselected_ca_perturbation/results/selected_vs_nonselected_ca_perturbation_results.csv` and `selected_vs_nonselected_ca_perturbation_statistics.csv`.

## Notes

The experiments decode the 56-bit UTF-8 payload `npj SB` plus a null terminator. Generated attacked PDB files were omitted because they are deterministic outputs of the retained scripts and account for most of the original workspace size without changing the reported conclusions.
