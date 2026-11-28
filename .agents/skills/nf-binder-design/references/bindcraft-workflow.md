# BindCraft Workflow — Full Reference

## Overview

The `--method bindcraft` workflow runs BindCraft trajectories in parallel across multiple GPUs. Unlike vanilla BindCraft (which runs indefinitely until finding N accepted designs), this pipeline runs a **fixed number of trajectories** and stops — giving predictable runtime and efficient GPU utilisation on HPC clusters.

For testing new targets and parameters, first run a small number of trajectories (`--bindcraft_n_traj 10`, `--bindcraft_n_traj 100`, or `--bindcraft_n_traj 300`) to assess the acceptance rate, then knowing the ratio of accepted designs to total trajectories, do a larger run to generate (approximately) the desired number of accepted designs. The acceptance rate is the number of accepted designs (in `results/bindcraft/accepted/`) divided by the total number of trajectories (eg `--bindcraft_n_traj` or `tail -n +2 results/bindcraft/trajectory_stats.csv | wc -l`). Acceptance rates lower than 1% will require substantial compute to generate enough binders for typical screening by assay — iterate on the target and hotspots to increase the acceptance rate in this case.

## GPU & VRAM Requirements

GPU memory is the primary bottleneck for BindCraft:

| GPU VRAM | Max complex size (target + binder) | Example GPUs              |
|----------|------------------------------------|---------------------------|
| 12 GB    | ~200–250 residues                  | RTX 3060 (12 GB), RTX 4070        |
| 16 GB    | ~300–350 residues                  | RTX 4080, T4                  |
| 24 GB    | ~400–450 residues                  | RTX 3090, RTX 4090, A5000, L4 |
| 32 GB    | ~550 residues                      | A100 (40 GB), V100 (32 GB)        |
| 48 GB    | ~750 residues                      | A6000, RTX 6000 Ada, L40S      |
| 80 GB    | ~950 residues                      | A100 (80 GB), H100        |

- **Speed**: An H100 is roughly 4× faster than an A100. As a rough guide on an H100: a 900-residue trajectory takes 2–3 hours; a 250-residue trajectory takes ~5 minutes.
- **CPU/RAM**: A single CPU core is sufficient. At least 40 GB RAM is recommended to avoid OOM from model compilation or PyRosetta.
- **Parallelisation**: BindCraft cannot split a single trajectory across multiple GPUs, but you can run multiple `--bindcraft_batch_size` batches in parallel (each on its own GPU). On a local workstation add `-profile local,multi_gpu`, which assigns each `BINDCRAFT` task its own card and caps concurrency at the card count. Under SLURM leave `--gpu_devices` unset and let the scheduler assign GPUs.
- **Run time**: "Easy" targets may only need ~100 trajectories; difficult targets can require 1,000–10,000. Some targets yield no successful designs.


## Parameters

Use `--method bindcraft --help` for the definitive list.

### Input / Output

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input_pdb` | Target PDB file, directory of PDBs, or quoted glob (e.g. `'input/*.pdb'`). `--bindcraft_n_traj` is run per input PDB | Required |
| `--outdir` | Output directory | `results` |

### Design Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--target_chains` | Target chain IDs, comma-separated (e.g. `"A"` or `"A,B"`) | Required |
| `--hotspot_res` | Hotspot residues (e.g. `"A56,A125"`) | Optional |
| `--hotspot_subsample` | Random proportion of hotspots per trajectory (0.0–1.0) | 1.0 (all) |
| `--binder_length_range` | Binder length range (e.g. `"55-120"`) | Required |
| `--contigs` | Alternative to `--target_chains`: RFdiffusion-style contigs (e.g. `"[A18-132/0]"`) | None |

You can use either `--target_chains` or `--contigs` to define the target binding surface. `--target_chains` is simpler (uses entire chain), while `--contigs` allows specifying specific residue ranges.

### Execution Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--bindcraft_n_traj` | Number of trajectories per input PDB | Required |
| `--bindcraft_batch_size` | Trajectories per batch (each batch runs on one GPU) | 1 |
| `--gpu_devices` | Comma-separated GPU IDs for multi-GPU (local profile only, e.g. `"0,1"`) | Auto |

### BindCraft Presets

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--bindcraft_advanced_settings_preset` | Settings preset from BindCraft `settings_advanced/` | `"default_4stage_multimer"` |
| `--bindcraft_filters_preset` | Filter preset from BindCraft `settings_filters/` | `"default_filters"` |

Preset names are without the `.json` extension. The experimental success rate of alternative presets may not have been as rigorously validated as `default_4stage_multimer`.

### SLURM

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--slurm_account` | SLURM account/project ID | None |

## Target Preparation and BindCraft Tips

To maximize success with BindCraft, carefully prepare your target:

- **Input Structure**: Even small variations of the "same" target (experimental vs. predicted vs. alternative trimming) can change *in silico* success rates significantly. Using an experimental structure as a template for a new AlphaFold2 or Boltz-2 prediction can help fill in missing loops.
- **Trimming**: Trim the target to only essential domains to save memory and time. Trim at realistic points like domain boundaries or hinge residues (e.g., Gly, Pro). "Unrealistic" trims (e.g., splitting GPCR helices) are supported *only if* they preserve the local structural context and don't expose core hydrophobic residues.
- **Hotspots**: Targeting a radial patch of surface residues on secondary structures is recommended, especially sites containing hydrophobic residues (F, Y, W, I, L, M). Single residues also work in well-defined binding sites.
- **Off-Target Binding**: If your hotspots are ignored because there is a significantly better binding site nearby, try:
  - Mutating the unwanted binding site residues to lysines
  - Trimming away the off-target region entirely
  - Pre-blocking the off-target site with another binder (generate a binder for the off-target site first, and use the complex as input)

## Estimating Trajectory Count

BindCraft acceptance rates vary by target. A suggested approach:

1. Start with a small run (`--bindcraft_n_traj 100`) to assess the acceptance rate
2. Calculate the ratio of accepted designs to total trajectories (the report headline collapses MPNN duplicates so multiple accepted sequences from one trajectory count once)
3. Scale up for a larger run to generate approximately the desired number of accepted designs

Multi-PDB runs tag each stats row with a `Target` column (input structure filename including extension, e.g. `PDL1.pdb`; the HTML report strips the extension for display) and the HTML report includes per-target accept rates.

## Output Structure

```
results/
└── bindcraft/
    ├── accepted/results/Accepted/     # Accepted design PDBs
    ├── batches/                       # Per-batch results
    │   ├── 0/results/
    │   └── 1/results/
    ├── bindcraft_report.html          # Summary report
    ├── failure_csv.csv                # Combined failure data
    ├── final_design_stats.csv         # Combined design statistics
    ├── mpnn_design_stats.csv          # MPNN statistics
    └── trajectory_stats.csv           # Trajectory statistics
```

## Refold Scoring (Optional)

`--refold_with` folds the accepted designs again with one or more independent predictors, after BindCraft's own AF2 passes and its accept/reject gate. It is a held-out second opinion, not a replacement, and BindCraft's own scores are what the selection reads.

```bash
--refold_with boltz,protenix
--refold_filter 'i_pTM>=0.7'
--refold_sort_by='-i_pTM'
--refold_max 50
--use_msa_server
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--refold_with` | Comma-separated engines: `af2`, `af2_mono`, `boltz`, `rf3`, `protenix`, `af3`, `openfold3`, `chai1`, `esmfold2`, `esmfold2_fast` | None (no refolding) |
| `--refold_filter` | Threshold gate on the `final_design_stats.csv` columns, semicolon-separated. Per-design, so it streams. A metric that table does not have stops the run | None |
| `--refold_sort_by` | Rank by one metric, `-` prefix for descending. Does not stream. **Write it with an `=`** | None |
| `--refold_max` | Cap. With `--refold_sort_by` the best N, without it the first N to be scored | No limit |
| `--refold_binder_chain` | Binder chain id in the accepted complexes | `B` |
| `--refold_target_chain` | Target chain id in the accepted complexes | `A` |

BindCraft is one of the two workflows with no chain derivation of its own, which is why `--refold_binder_chain` / `--refold_target_chain` matter here.

Outputs go to `results/bindcraft/refold/`: `refold_scores.tsv` (engine-prefixed columns `<tool>_iptm_max`, `<tool>_ipsae_max`, `<tool>_target_aligned_binder_rmsd_all`, `<tool>_interface_dG`, ... plus `consensus_z`), `refold_sample_scores.tsv` and `refold_selection.tsv`.

## Complete Example

### Local Workstation

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method bindcraft \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56" \
  --binder_length_range "55-120" \
  --bindcraft_n_traj 4 \
  --bindcraft_batch_size 1 \
  --bindcraft_advanced_settings_preset "default_4stage_multimer" \
  -profile local \
  -resume
```

### HPC Cluster (SLURM)

Use `-profile slurm` instead of `-profile local`, add `--slurm_account`, and optionally add a site-specific profile or `-c` config:

```bash
-profile slurm \
--slurm_account=ab12
```
