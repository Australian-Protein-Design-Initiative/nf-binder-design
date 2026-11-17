# Germinal Workflow — Full Reference

## Overview

The `--method germinal` workflow runs [Germinal](https://github.com/SantiagoMille/germinal) trajectories in parallel across multiple GPUs — ideal for HPC clusters or multi-GPU workstations. Configuration is supplied via a Germinal Hydra YAML file (combined or partial config with Hydra default groups).

Unlike a single long Germinal run that loops until stopping criteria are met, this pipeline:

- Runs a **fixed number of trajectories** (`--germinal_n_traj`)
- Splits work into **parallel batches** (`--germinal_batch_size`)
- Merges per-batch CSVs and structure folders into a single output tree

If you want a specific number of accepted designs, run a small pilot (`--germinal_n_traj 10`) to estimate acceptance rate, then scale up.

## Parameters

Use `--method germinal --help` for the definitive list.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--germinal_config` | Germinal Hydra config YAML (required) | — |
| `--germinal_pdb_dir` | Directory with target PDB and optional `nb.pdb` scaffold | `../pdbs` relative to config |
| `--germinal_experiment_name` | Output subdirectory name under `results/` | Required |
| `--germinal_n_traj` | Total trajectory attempts across all batches | Required |
| `--germinal_batch_size` | Trajectories per parallel batch | 1 |
| `--germinal_max_passing_designs` | Max accepted designs per batch | High default |
| `--germinal_max_hallucinated_trajectories` | Max hallucinated trajectories per batch | High default |
| `--gpu_devices` | GPU devices for local multi-GPU runs (e.g. `"0,1"`) | Auto |

Partial configs (e.g. `configs/config.yaml` with `defaults: [run: vhh, target: pdl1, ...]`) are supported: the workflow copies built-in Hydra config groups from the container and places your config alongside them.

## Example

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method germinal \
  --germinal_config configs/pdl1_vhh.yaml \
  --germinal_pdb_dir pdbs \
  --germinal_experiment_name pdl1_vhh \
  --germinal_n_traj 4 \
  --germinal_batch_size 1 \
  --outdir results \
  -profile local -resume
```

For local dual-GPU:

```bash
nextflow run main.nf \
  -c nextflow.dual-gpu.config \
  --method germinal \
  --germinal_config configs/pdl1_vhh.yaml \
  --germinal_pdb_dir pdbs \
  --germinal_experiment_name pdl1_vhh \
  --germinal_n_traj 4 \
  --germinal_batch_size 1 \
  --gpu_devices=0,1 \
  --outdir results \
  -profile local -resume
```

For SLURM on M3 BDI: `-profile slurm,m3_bdi` with `--slurm_account=yt41`.

See `examples/pdl1-germinal/` for a complete example.

## Output Structure

```
results/germinal/
├── all_trajectories.csv
├── accepted_designs.csv
├── failure_counts.csv
├── config/final_config.yaml
├── accepted/structures/
├── trajectories/
├── redesign_candidates/
└── batches/0/{experiment_name}/
```

## Refold Scoring (Optional)

`--refold_with` folds the accepted designs again with one or more independent predictors, after Germinal's own AF2/AF3 pass and its accept/reject gate. It is a held-out second opinion, not a replacement, and Germinal's own scores are what the selection reads.

```bash
--refold_with boltz,protenix
--refold_sort_by='-<metric>'
--refold_max 50
--use_msa_server
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--refold_with` | Comma-separated engines: `af2`, `af2_mono`, `boltz`, `rf3`, `protenix`, `af3`, `openfold3`, `chai1`, `esmfold2`, `esmfold2_fast` | None (no refolding) |
| `--refold_filter` | Threshold gate on the `accepted_designs.csv` columns, semicolon-separated. Per-design, so it streams. Check the column names in that file before writing one: a metric it does not have stops the run | None |
| `--refold_sort_by` | Rank by one `accepted_designs.csv` metric, `-` prefix for descending. Does not stream. **Write it with an `=`** | None |
| `--refold_max` | Cap. With `--refold_sort_by` the best N, without it the first N to be scored | No limit |
| `--refold_binder_chain` | Binder chain id in the accepted complexes | `B` |
| `--refold_target_chain` | Target chain id in the accepted complexes | `A` |

Germinal is one of the two workflows with no chain derivation of its own, which is why `--refold_binder_chain` / `--refold_target_chain` matter here.

Outputs go to `results/germinal/refold/`: `refold_scores.tsv` (engine-prefixed columns `<tool>_iptm_max`, `<tool>_ipsae_max`, `<tool>_target_aligned_binder_rmsd_all`, `<tool>_interface_dG`, ... plus `consensus_z`), `refold_sample_scores.tsv` and `refold_selection.tsv`.

## Hotspot numbering

`target.target_hotspots` uses **1-indexed positions relative to the start of each target chain** (Nth residue in that chain as loaded), not arbitrary PDB auth residue numbers. Germinal's hotspot proximity filter maps hotspots via sequential indices (`chain_start + index - 1`), so numbering must be contiguous from 1 with **no gaps**.

Renumber the input target PDB to 1…N per chain (e.g. `bin/renumber_chains.py`), then set `target_hotspots` / `hotspot_residue` against that renumbered structure.

## Notes

- The nanobody scaffold (`nb.pdb`) is copied from the container into `pdb_dir` if missing.
- `target.target_pdb_path` in the config should be relative (e.g. `pdbs/pdl1.pdb`) so paths resolve from the process working directory.
- Stopping criteria per batch: whichever of `max_trajectories`, `max_hallucinated_trajectories`, or `max_passing_designs` is reached first.
