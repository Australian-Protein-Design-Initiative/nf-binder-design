# Germinal Workflow

Parallel [Germinal](https://github.com/SantiagoMille/germinal) execution for antibody and nanobody design across multiple GPUs.

## Overview

The `--method germinal` workflow runs Germinal trajectories in parallel across multiple GPUs — ideal for HPC clusters or multi-GPU workstations. Configuration is supplied via a Germinal Hydra YAML file (combined or partial config with Hydra default groups).

## Key Differences

Unlike a single long Germinal run that loops until stopping criteria are met, this pipeline:

- Runs a **fixed number of trajectories** (`--germinal_n_traj`)
- Splits work into **parallel batches** (`--germinal_batch_size`)
- Merges per-batch CSVs and structure folders into a single output tree

If you want a specific number of accepted designs, run a small pilot (`--germinal_n_traj 100` or more) to estimate acceptance rate, then scale up. The Germinal documentation and paper detail the specific parameter sweeps you may want to try (via different config files) to improve the acceptance rate.

## Command-line Options

See available options with `--method germinal` and no config:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method germinal
```

## Example Usage

```bash
#!/bin/bash

DATESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_DIR=/path/to/runs/germinal/test-protenix2

nextflow run /path/to/nf-binder-design-germinal/main.nf \
  --method germinal \
  --germinal_config "${RUN_DIR}/configs/pdl1_vhh_protenix.yaml" \
  --germinal_pdb_dir "${RUN_DIR}/pdbs" \
  --germinal_experiment_name pdl1_vhh \
  --germinal_n_traj 2 \
  --germinal_batch_size 1 \
  --outdir "${RUN_DIR}/results/nf-germinal" \
  -profile local \
  -resume \
  -with-report "results/logs/report_${DATESTAMP}.html" \
  -with-trace "results/logs/trace_${DATESTAMP}.txt"
```

For SLURM on M3 BDI, use `-profile slurm,m3_bdi` with `--slurm_account=yt41`.

Partial configs (e.g. `configs/config.yaml` with `defaults: [run: vhh, target: pdl1, ...]`) are supported: the workflow copies built-in Hydra config groups from the container and places your config alongside them.

## Key Parameters

| Flag | Description |
|------|-------------|
| `--germinal_config` | Germinal Hydra config YAML (required) |
| `--germinal_pdb_dir` | Directory with target PDB and optional `nb.pdb` scaffold (default: `../pdbs` relative to config) |
| `--germinal_experiment_name` | Output subdirectory name under `results/` |
| `--germinal_n_traj` | Total trajectory attempts across all batches |
| `--germinal_batch_size` | Trajectories per parallel batch |
| `--germinal_max_passing_designs` | Max accepted designs per batch (default: high) |
| `--germinal_max_hallucinated_trajectories` | Max hallucinated trajectories per batch (default: high) |
| `--gpu_devices` | GPU devices for local multi-GPU runs, e.g. `--gpu_devices=0,1` |

## Output Structure

```
results/germinal/
├── all_trajectories.csv      # merged across all batches
├── accepted_designs.csv      # merged from per-batch accepted/designs.csv
├── failure_counts.csv
├── config/
│   └── final_config.yaml     # resolved Hydra config (batch 0 only)
├── accepted/
│   └── structures/           # flattened accepted PDBs from all batches
├── trajectories/
│   ├── designs.csv
│   └── structures/
├── redesign_candidates/
│   ├── designs.csv
│   └── structures/
└── batches/
    └── 0/
        └── pdl1_vhh/         # per-batch Germinal output (includes final_config.yaml)
```

## Hotspot numbering

`target.target_hotspots` in the Germinal config (e.g. `"A37,A39,A41"`) uses **1-indexed positions relative to the start of each target chain** — the Nth residue in that chain as loaded, not arbitrary PDB/mmCIF auth residue numbers.

Germinal's hotspot proximity filter (`find_nearby_residues_from_pdb`) maps each hotspot to pose residue `chain_start + index - 1`. Residue indices must therefore be contiguous from 1 with **no gaps**. If your input PDB is numbered from a non-1 start, has numbering gaps (e.g. missing residues left as holes in the numbering), or otherwise does not match 1…N sequential order, hotspot values will not refer to the residues you intend.

**Recommendation:** renumber the input target PDB so each chain is numbered sequentially from 1 with no gaps, then choose hotspots against that renumbered structure (e.g. in ChimeraX / PyMOL / Mol*). Use `bin/renumber_chains.py`:

```bash
uv run bin/renumber_chains.py input/target.pdb -o pdbs/target.pdb
```

After renumbering, set `target_hotspots` (and `hotspot_residue` if used) to match the new 1-based sequential numbers.

## Notes

- The nanobody scaffold (`nb.pdb`) is copied from the container into `pdb_dir` if missing.
- `target.target_pdb_path` in the config should be relative (e.g. `pdbs/pdl1.pdb`) so paths resolve from the process working directory.
- Stopping criteria per batch: whichever of `max_trajectories`, `max_hallucinated_trajectories`, or `max_passing_designs` is reached first.
