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
- **Parallelisation**: BindCraft cannot split a single trajectory across multiple GPUs, but you can run multiple `--bindcraft_batch_size` batches in parallel (each on its own GPU). Increase `--bindcraft_batch_size` and/or `--gpu_devices` to use multiple GPUs. **On local workstations**, multi-GPU execution is unreliable — use a single GPU unless on a SLURM cluster where the scheduler handles GPU assignment. For local multi-GPU, a special Nextflow config is required (see `examples/` in the pipeline repo) and `--gpu_devices` must be set. Race conditions in GPU allocation are common.
- **Run time**: "Easy" targets may only need ~100 trajectories; difficult targets can require 1,000–10,000. Some targets yield no successful designs.


## Parameters

Use `--method bindcraft --help` for the definitive list.

### Input / Output

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input_pdb` | Target PDB file (single file, not a glob) | Required |
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
| `--bindcraft_n_traj` | Total number of trajectories | Required |
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
2. Calculate the ratio of accepted designs to total trajectories
3. Scale up for a larger run to generate approximately the desired number of accepted designs

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

## Complete Example

### Local Workstation

```bash
DATESTAMP=$(date +%Y%m%d_%H%M%S)

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
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
```

### HPC Cluster (SLURM)

Use `-profile slurm` instead of `-profile local`, add `--slurm_account`, and optionally add a site-specific profile or `-c` config:

```bash
-profile slurm \
--slurm_account=ab12
```
