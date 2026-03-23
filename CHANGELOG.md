# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- New `nci_gadi.config` configuration profile for NCI Gadi HPC cluster.
- New `--method rfd3` workflow for RFDiffusion3-based binder design using `RosettaCommons/foundry`.

### Fixed
- Boltz: `BOLTZ`, `BOLTZ_COMPARE_COMPLEX`, and `BOLTZ_COMPARE_BINDER_MONOMER` tee `boltz predict` to `.boltz_predict_console.log` and exit 1 if the log contains `ran out of memory, skipping batch` (Boltz may otherwise exit 0 and leave outputs missing).

### Changed
- `merge_scores.py`: drop from the right any column that exists in the current left before each merge so the result has no `_x`/`_y` suffixes; treat `.cif` as path-like (use basename for merge key) in addition to `.pdb`.
- `trim_to_contigs.py`: `parse_contigs()` now supports RFD3 v3 contig format (comma-separated, `/0` as separate element, e.g. `A18-132,/0,65-120`) in addition to v1 style.
- Major project restructure shifting individual workflows into `workflows/`, each launched via a single `main.nf` entry point with the `--method` flag.
  - `--method rfd` for RFdiffusion binder design (previously `main.nf`)
  - `--method rfd_partial` for partial diffusion (previously `partial.nf`)
  - `--method bindcraft` for BindCraft (previously `bindcraft.nf`)
  - `--method boltzgen` for BoltzGen (previously `boltzgen.nf`)
  - `--method boltz_pulldown` for Boltz Pulldown (previously `boltz_pulldown.nf`)
- Modules reorganised into `modules/local/` with workflow-specific subdirectories (`rfd/`, `bindcraft/`, `boltzgen/`, `common/`).
- Extracted common Boltz-2 refolding and scoring logic into `BOLTZ_REFOLD_SCORING` subworkflow (`subworkflows/local/boltz_refold_scoring.nf`).
- `filter_designs.py` + `filters.d/rg.py` now support Rg calculation for mmCIF inputs (including `.cif.gz`) via a `gemmi` fallback path.

### Removed
- `bin/rfd3/normalise_combined_tsv.py` (replaced by csvtk + `merge_scores.py` in the RFD3 combine step).
- `rmsd4all.py`: remove `--max-ca-for-tm-score` and `--pair-timeout` options.

## [0.1.5] - 2026-01-28

### Added
- New BoltzGen pipeline.
- Refold binder designs with Boltz, with post-AF2ig filtering and 
  RMSD analysis of predicted complex and binder monomer.
- Added DOI (Zenodo) badge to `README.md`, added `CITATION.cff`.
- Some parameter validation for `bindcraft.nf`.
- Write `params.json` to output directory.

### Changed
- Move 'filtering' result folder to 'rfdiffison/filtered'
- Don't output redundant .tsv files to the results directory.
- Changed default `--pmpnn_relax_cycles` from 0 to 3
- Made default queue size 1, for single local GPU mode.
- Added `m3-bdi.config`, site specific for M3/MASSIVE HPC cluster.
- Use the `nf-binder-design-utils` container instead of `mdanalysis`.
- Update config and docs to use -profile for site-specifc configurations

### Fixed
- Fixed Quarto rendering permissions issues (copy Qmd to work folder).

## [0.1.4] - 2025-08-12

### Added
- BindCraft end-to-end workflow with basic HTML report
- GPU allocation heuristics for local multi-GPU workstations.
- HyperMPNN weights download script (`models/download_hypermpnn_weights.sh`).
- New runnable examples in `examples/`.

### Changed
- RFDiffusion `--hotspot_res` no longer requires brackets in `main.nf`.
- Minor tweaks to site-specific configs for `m3` and `mlerp`.

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