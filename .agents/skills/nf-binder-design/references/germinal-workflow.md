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

## Notes

- The nanobody scaffold (`nb.pdb`) is copied from the container into `pdb_dir` if missing.
- `target.target_pdb_path` in the config should be relative (e.g. `pdbs/pdl1.pdb`) so paths resolve from the process working directory.
- Stopping criteria per batch: whichever of `max_trajectories`, `max_hallucinated_trajectories`, or `max_passing_designs` is reached first.
