# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- Trimmed README.md, testing section moved to `docs/docs/extra/development.md`, general docs cleaunp and corrections.
- Docs: note Nextflow version compatibility. On Nextflow `26.04+` (new strict parser default), set `NXF_SYNTAX_PARSER=v1` to use the legacy parser.
- `manifest.nextflowVersion` now bounds the supported range to `!>=23.04.0, <26.10` (hard failure outside this range).

### Added
- Spartan HPC platform configs `spartan-a100.config` (gpu-a100-short) and `spartan-l40s.config` (gpu-l40s) for University of Melbourne Spartan.
- Agent skill at `.agents/skills/nf-binder-design/` for AI-assisted pipeline setup and execution.
- RFD workflow docs: table of built-in and bind-mounted HyperMPNN `--pmpnn_weights` checkpoints in the `proteinmpnn_dl_binder_design` container.
- New `--method rfd3` workflow for RFDiffusion3-based binder design using `RosettaCommons/foundry` (RF3 batching, Boltz full-refold scoring, optional FoldSeek on refolded designs).
- Germinal antibody/nanobody design workflow (`--method germinal`).
- FoldSeek structural search (`--do_foldseek`) for the `rfd`, `bindcraft`, `boltzgen`, and `rfd3` workflows. Searches designed binder chains against structural databases (default: CATH50) to identify known folds and annotate results with CATH hierarchy descriptions. Supports local or remote search, gzip output, and optional HTML reports.
- `bin/complex_sasa.py`: per-residue delta SASA for target chains when a binder is removed from a complex, with optional site sums, batch PDB input, and `--min-change-percent` column pruning.
- nf-test `tests/pipeline/compilation.nf.test`: launches `rfd`, `rfd_partial`, and `rfd3` in `-preview` mode to verify every workflow/module compiles; run per Nextflow version with `NXF_VER=<version> nf-test ...` to guard against version-specific DSL parser regressions.

### Fixed
- Nextflow 24.04.3 compatibility: the `rfd3` workflow and `boltz_refold_core` subworkflow no longer trigger the "Variable already defined in the process scope" DSL parser error on Nextflow 24.04.3 (a bug fixed in the 24.10 parser rewrite). Channel/path local variables use plain assignments instead of `def` in these workflow bodies.
- `rfd`: configurable `RFDIFFUSION` via `rfd_command` and `rfd_model_directory_path` (Pawsey config overrides in `pawsey_setonix.config`; GPU behaviour uses existing `require_gpu` and `gpu_devices`).

## [0.2.0] - 2026-05-06

### Added
- New `nci_gadi.config` configuration profile for NCI Gadi HPC cluster.
- Initial [nf-test](https://www.nf-test.com/) test scaffold (`nf-test.config`, `tests/`) with a process test for `UNIQUE_ID` and `RFDIFFUSION`; documented in `docs/docs/testing.md`.
- BoltzGen: support for list-valued `entities[].file.path` and multiple entities so all referenced files are staged as Nextflow `path()` inputs. Referenced config YAMLs (`.yaml`/`.yml` in `entities[].file.path`, eg for nanobody scaffolds) and the PDB/CIF files they reference internally are collected and staged so BoltzGen’s per-generation random selection over those configs is preserved.
- Versioned documentation using [mike](https://github.com/jimporter/mike); docs are now deployed for `main` (alias: `latest`), `develop` (alias: `dev`), and version tags.

### Fixed
- Boltz refold RMSD and ipSAE: `BOLTZ_COMPARE_COMPLEX` and `BOLTZ_COMPARE_BINDER_MONOMER` now use Boltz complex chain **A** = binder and **B** = target (from `create_boltz_yaml.py` IDs), while the input design keeps `${binder_chain}` / `${target_chain}`. Previously the same chain letters were used on both structures, which mis-superposed RFD3-style complexes (target A, binder B) and inflated `rmsd_target_aligned_binder` / ruined aligned PDBs; `ipsae.py` now receives `--binder-chain A --target-chain B` for Boltz outputs.
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

### Fixed
- BindCraft (`workflows/bindcraft.nf`): omitting `--hotspot_res` no longer fails in `validateHotspotRes` with `Unknown method invocation 'trim' on Boolean type`.
- BindCraft (`workflows/bindcraft.nf`): allow explicit `--hotspot_res=""` to pass through as a no-hotspot BindCraft setting, and normalise empty hotspot list entries before writing settings.
- Declare `bindcraft_batch_size` default in `nextflow.config` so Nextflow does not warn when parsing `modules/local/bindcraft/bindcraft.nf` (unrelated to `-profile m3`).
- `rmsd4all.py`: cap worker processes to the number of pairs so single-pair comparisons (e.g. RFD3_RMSD with one design vs one refold) no longer spawn a large Pool and appear to hang; sequential path is used for one pair with progress logged.
- `rmsd4all.py`: add `--max-structural-iterations` (default 100). Biotite's refinement loop uses `max_iterations=inf` by default and only stops when anchors stabilize; with 3di the anchor set can fail to converge (oscillate) so the loop never exits. Capping iterations fixes the hang; 0 = no limit.
- `rmsd4all.py`: fix use of `array_length` (method) as if it were an attribute; use `len()` for atom counts.

### Removed
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