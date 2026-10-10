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
  -resume
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

## Structure prediction with AlphaFold3 (`af3`)

Germinal's filter stage re-predicts each design with an external structure
predictor. The default is `protenix`, but AlphaFold3 (`af3`) is also an option. 

Select AlphaFold3 in the Germinal config:

```yaml
structure_model: "af3"   # default is "protenix"
```

This requires the `-af3` germinal image, which the workflow already pins
(`germinal:20260611-104bbdd7-af3`). It bundles AlphaFold3 3.0.4, but contains **neither the model parameters nor the public sequence and template structure databases**. Both must be bind-mounted by you via the `nextflow.config`.

### Config settings

| Setting | Value | Provided by |
|---------|-------|-------------|
| `af3_repo_path` | `"/app/alphafold"` | **In the container — leave as-is** |
| `af3_sif_path` | `"/usr/local/bin/singularity"` | **Dummy wrapper, in the container — leave as-is** |
| `af3_model_dir` | `"/root/models"` | Mount point — **you bind the weights here** |
| `af3_db_dir` | `"/root/public_databases"` | Mount point — **you bind the databases here** |
| `msa_db_dir` | `null` | Only used by `msa_mode: "local"` |
| `use_metagenomic_db` | `false` | — |
| `cache_binder_msa` | `false` | Requires `msa_mode: "colabfold"` |

All five paths are already set to these values in the container's own
`configs/run/*.yaml`, so a config derived from those needs no changes.

`af3_sif_path` looks odd because upstream Germinal shells out to
`singularity exec … <image> python /root/alphafold3/run_alphafold.py`. There is
no nested Apptainer inside the container, so the image provides a shim at
`/usr/local/bin/singularity` that discards the image argument and runs AF3 from
its bundled virtualenv. Both it and `af3_repo_path` should be left as shipped.

!!! warning "Do not point `af3_model_dir` / `af3_db_dir` at host paths"
    These two must stay as the in-container mount points above. The shim
    implements Germinal's `--bind src:dest` as `rm -rf dest && ln -sfn src dest`,
    and `/root` is read-only under Apptainer, so a host path fails. Leave the
    config pointing at `/root/models` and `/root/public_databases`, and supply
    the real directories as bind mounts instead.

### Bind mounts

Add both mounts to the `GERMINAL` process, e.g. in a `nextflow.config` in your
run directory:

```groovy
process {
    withName: /^(GERMINAL|GERMINAL_PROCESS)$/ {
        containerOptions = '-B /path/to/af3_weights:/root/models -B /path/to/af3_databases:/root/public_databases'
    }
}
```

**Weights** — a directory holding exactly one `af3.bin.zst` (or `af3.bin`).
Obtain them under the AlphaFold3 Model Parameters Terms of Use; see
[AlphaFold3 weights](fold.md#alphafold3-weights) for `models/download_af3_weights.sh`
and the licensing conditions.

**Databases** — the full AlphaFold3 public database set (~630 GB), as produced by
AlphaFold3's own `fetch_databases.sh`. All nine entries are required *even though
Germinal supplies the MSAs itself*, because `run_alphafold.py` validates every
default database path against `--db_dir` before it reads the input. An empty or
missing directory fails with:

```
FileNotFoundError: ${DB_DIR}/bfd-first_non_consensus_sequences.fasta
  with ${DB_DIR} not found in any of ['/root/public_databases']
```

The directory must contain `bfd-first_non_consensus_sequences.fasta`,
`mgy_clusters_2022_05.fa`, `mmcif_files/`,
`nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta`,
`pdb_seqres_2022_09_28.fasta`, `rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta`,
`rnacentral_active_seq_id_90_cov_80_linclust.fasta`, `uniprot_all_2021_04.fa`
and `uniref90_2022_05.fa`.

### What AF3 actually runs

`msa_mode` controls the MSA, which **Germinal** generates before calling AF3 —
`"target"` (the default) builds one for the target chain only and leaves the
binder single-sequence. AF3's own MSA search is therefore unused, but its data
pipeline still runs and performs a **template search for both chains** against
`pdb_seqres` and `mmcif_files`.

On M3 the databases are already installed cluster-wide — see
[AlphaFold3 databases on M3](../extra/m3-hpc-examples.md#alphafold3-databases).

A complete worked example is in `examples/pdl1-germinal/`
(`configs/pdl1_vhh_af3.yaml` and `run-af3.sh`).

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

## External folding validation (optional) {#external-folding-validation}

Germinal already re-predicts every design with one external predictor and gates
on its scores — that is the filter stage described above, and it is what decides
which designs reach `accepted/`. This section is about an **additional,
independent re-fold afterwards**, with predictors Germinal did not use.

It is optional and sits outside the workflow today: you run it as a separate
`--method fold_pulldown` job over the designs a Germinal run produced. It may
become an option on the Germinal workflow itself in future.

### Running it

Build a FASTA of the target and a FASTA of the binder sequences (from
`accepted/designs.csv`, or the `trajectory_sequence` column of
`all_trajectories.csv`), then fold every pair with whichever predictors you want
to compare:

```bash
nextflow run main.nf --method fold_pulldown \
    --targets target.fasta --binders binders.fasta \
    --methods af3,protenix,boltz,esmfold2 \
    --create_target_msa true --create_binder_msa false \
    --msa_method mmseqs2_colabfold --use_remote_server true \
    --outdir results_refold
```

`--create_target_msa true --create_binder_msa false` is Germinal's
`msa_mode: "target"` — an MSA for the target chain only, binder left
single-sequence — so the MSA conditioning matches what the Germinal run used.

### Matching Germinal's AF3 settings

If AlphaFold3 is one of the predictors and you intend to compare its numbers
against the ones in `all_trajectories.csv`, the default `fold_pulldown` settings
are **not** a like-for-like comparison. Germinal's AF3 call supplies no
cross-chain pairing and omits the `templates` key, so AF3 searches the PDB for
templates on both chains (see [What AF3 actually runs](#what-af3-actually-runs));
`fold_pulldown` defaults to inference-only with a paired MSA and no templates.

`-profile af3_germinal_parity` switches all three, and bundles them with the MSA
and sampling settings Germinal uses:

```bash
nextflow run main.nf --method fold_pulldown \
    --targets target.fasta --binders binders.fasta --methods af3 \
    -profile slurm,m3,af3_germinal_parity \
    --af3_model_dir /path/to/af3_weights \
    --af3_db_dir /mnt/datasets/alphafold3/3.0.0 \
    --outdir results_refold
```

It composes with any platform profile, in either order — Nextflow merges
profiles in the order they are declared in `nextflow.config`, not the order you
list them, and the feature profiles are declared last. That matters because this
profile raises the AF3 walltime: with the data pipeline on, AF3 runs a template
search before inference, which the platform formulas (sized for inference only)
do not allow for.

See [Germinal AF3 parity](fold.md#germinal-af3-parity) for the full mapping.

### Comparing the numbers

- **Structure selection.** Germinal keeps the **worst** of its sampled
  structures by `ranking_score` (`af3_structure_select_mode: "worst"`) and
  reports every metric from that one, while `fold_pulldown` keeps and scores all
  of them. A per-sample score from one is not comparable to a Germinal number
  from the other until you reduce it the same way.

No rescaling is needed beyond that: `fold_pulldown` already computes `plddt` as
the mean over all atom pLDDTs and `pae` as the mean of the whole PAE matrix,
which are Germinal's definitions, and runs ipSAE at the same 10/10 cutoffs.

## Notes

- The nanobody scaffold (`nb.pdb`) is copied from the container into `pdb_dir` if missing.
- `target.target_pdb_path` in the config should be relative (e.g. `pdbs/pdl1.pdb`) so paths resolve from the process working directory.
- Stopping criteria per batch: whichever of `max_trajectories`, `max_hallucinated_trajectories`, or `max_passing_designs` is reached first.
