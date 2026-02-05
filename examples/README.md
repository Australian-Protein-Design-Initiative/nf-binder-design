# Example Configurations

This directory contains example configurations for the different workflows available in this repository.
These examples are setup for quick test runs on a local GPU workstation. For real production runs, or running on HPC (SLURM etc), you will need to adjust parameters (e.g., number of designs/trajectories) and the `nextflow.config`.

## RFdiffusion Workflows

Run with `--method rfd` or `--method rfd_partial`.

- **`pdl1-rfd`**: Standard binder design against PD-L1.
  - Workflow: RFdiffusion → ProteinMPNN → AlphaFold2 (initial guess) → Boltz-2 (optional refolding).
- **`pdl1-rfd-partial`**: Partial diffusion refinement.
  - Workflow: Refines existing designs using RFdiffusion partial diffusion.
- **`egfr-rfd-hypermpnn`**: EGFR binder design.
  - Uses **HyperMPNN** weights for the sequence generation step.

## BindCraft Workflows

Run with `--method bindcraft`.

- **`pdl1-bindcraft`**: PD-L1 binder design.
  - Standard protein binder design using BindCraft.
- **`egfr-bindcraft`**: EGFR binder design.

## BoltzGen Workflows

Run with `--method boltzgen`.

- **`boltzgen-protein`**: Protein binder design.
  - Protocol: `protein-anything`
  - Target: 1G13
- **`boltzgen-cyclic-peptide`**: Cyclic peptide binder design.
  - Protocol: `peptide-anything`
  - Target: 2VSM
- **`boltzgen-small-molecule-binder`**: Small molecule binder design.
  - Protocol: `protein-small_molecule`
  - Target: PFOA (Perfluorooctanoic acid) binding protein
