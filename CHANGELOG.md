# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.3] - 2025-07-11

### Added
- A plugin system for filtering designs based on calculated metrics.
- `--rfd_filters` parameter to apply filters to RFDiffusion backbones (e.g., `--rfd_filters "rg<25"`).
- Initial filter plugin for radius of gyration (`rg`).
- Integration of BindCraft-derived scoring of designs, as `extra_scores.tsv`. This adds metrics including:
  - Interface score, shape complementarity, dG, and dSASA.
  - Secondary structure percentages (helix, sheet, loop) for the binder and interface.
  - Unrelaxed and relaxed clash scores.
  - Hotspot and target RMSD.
  - Sequence-based metrics like extinction coefficient.
- `--pmpnn_omit_aas` flag to exclude certain amino acids from designs (default 'CX').
- `--rfd_compress_trajectories` flag to gzip RFDiffusion traj/*.pdb files (default `true`).

### Changed
- Filtering results are now saved to a subdirectory named after the pipeline step (in `filtering/rfdiffusion/`).
- The score merging script (`merge_scores.py`) was improved to be more robust and handle an arbitrary number of score files.
- Use more lightweight container for `get_contigs.nf` and `renumber_residues.nf` modules.

### Fixed
- The `gpu_device` parameter is now correctly passed to the partial diffusion process, allowing proper GPU selection.

## [0.1.2] - 2025-06-19

### Added
- Added the `boltz_pulldown.nf` protocol.
- Allow gpu device to be selected.
- Support relaxation in ProteinMPNN (`pmpnn_relax_cycles` can now be non-zero)

### Changed
- Change `dl_binder_design` output filenaming (include "_mpnn{n}" suffix)
- Change RFDiffusion config handling when unspecified
- Docs: Add APPTAINER_TMPDIR and NXF_APPTAINER_CACHEDIR advice to M3-specific docs

### Fixed
- Fix `get_contigs` for multi-chain targets.
- Fix failure that occurred when using `--rfd_batch_size` > 1.

## [0.1.1] - 2025-04-04

### Changed
- Update to use containers from Github package registry

## [0.1] - 2025-03-24

### Added
- Initial version with RFDiffusion->ProteinMPNN->af2_initial_guess binder
  design and partial diffusion pipelines.