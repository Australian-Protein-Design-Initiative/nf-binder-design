# BoltzGen Workflow — Full Reference

## Overview

The `--method boltzgen` workflow uses the BoltzGen generative model for binder design. It supports multiple protocols and parallelises design batches across GPUs.

**Pipeline steps**: Design → Inverse Folding → Folding → Design Folding → (Merge batches) → Analysis & Filtering

## Design Scale: Testing vs Production

BoltzGen runs a **fixed number of designs** (`--num_designs`) then filters to a smaller `--budget`. Choose `--num_designs` based on whether you are validating the setup or running a full campaign.

| Run type | `--num_designs` | Purpose |
|----------|-----------------|---------|
| **Smoke test** | 2–10 | Verify YAML, paths, protocol, and pipeline wiring |
| **Pilot / parameter tuning** | 100–1,000 | Assess hit rate, filtering thresholds, and hotspot/YAML choices |
| **Production** | **40,000–60,000** | Full design campaign for experimental screening |

For a new target or YAML config, **always start small** (smoke test, then pilot) before scaling to production. Check how many designs pass filtering in `boltzgen/filtered/final_ranked_designs/` and inspect merged metrics under `boltzgen/merged/` to decide whether hotspots, `--alpha`, `--budget`, or the YAML need adjustment.

**Production notes:**

- Use `--batch_size` and `--devices` (or SLURM multi-GPU) to parallelise batches across GPUs; increase `--batch_size` only if VRAM allows.
- Set `--budget` to the number of diverse, high-quality designs you want for wet-lab ordering (often tens to low hundreds), not the raw `--num_designs`.
- A production run at 40k–60k designs is compute-heavy — confirm pilot results and storage (`work/`, containers, merged outputs) before launching on HPC.

Example production-scale command (adjust `--budget` and `--batch_size` for your cluster):

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method boltzgen \
  --config_yaml /abs/path/to/design.yaml \
  --outdir results \
  --protocol protein-anything \
  --num_designs 50000 \
  --batch_size 10 \
  --budget 100 \
  -profile slurm \
  --slurm_account=ab12 \
  -resume
```

See `examples/pdl1-boltzgen/` and `examples/boltzgen-nanobody/` for small local test runs (`--num_designs 4` or similar).

## YAML Configuration

BoltzGen requires a YAML configuration file (`--config_yaml`) that defines the input structure, constraints, and design objectives. Refer to the [BoltzGen documentation](https://github.com/HannesStark/boltzgen?tab=readme-ov-file#how-to-make-a-design-specification-yaml) for preparing configuration files.

**Important**: BoltzGen uses 1-based residue index numbering, always starting at 1, regardless of the numbering in your PDB/mmCIF file. Use `bin/renumber_chains.py` to renumber your input PDB to sequential numbering from 1 to simplify hotspot residue selection in molecular viewers.

**Important**: Use **absolute paths** in BoltzGen YAML config files — BoltzGen runs inside a container with a different working directory:

```yaml
# Wrong
protein:
  filepath: input/target.pdb

# Correct
protein:
  filepath: /home/user/project/input/target.pdb
```

## Parameters

Use `--method boltzgen --help` for the definitive list.

### Core Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--config_yaml` | **(Required)** Path to BoltzGen YAML configuration file | Required |
| `--outdir` | Output directory | `results` |
| `--design_name` | Name prefix for outputs | Basename of config file |
| `--protocol` | Protocol type (see below) | `protein-anything` |

### Protocols

| Protocol | Description |
|----------|-------------|
| `protein-anything` | Standard protein binder design |
| `peptide-anything` | Peptide binder (including cyclic peptides) |
| `protein-small_molecule` | Small molecule binding protein design |
| `nanobody-anything` | Nanobody design |

### Design Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--num_designs` | Total number of designs to generate (see **Design Scale** above; default is not production scale) | 100 |
| `--inverse_fold_num_sequences` | Sequences per backbone during inverse folding | 1 |
| `--batch_size` | Designs per batch | 10 |
| `--devices` | Number of GPU devices | 1 |
| `--num_workers` | DataLoader workers | Auto |

### Filtering Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--budget` | Final designs to keep after diversity optimisation and filtering | 10 |
| `--alpha` | Trade-off: 0.0=quality-only, 1.0=diversity-only | Auto |
| `--filter_biased` | Remove amino-acid composition outliers | `true` |
| `--metrics_override` | Per-metric inverse-importance weights (e.g. `'plip_hbonds_refolded=4'`) | None |
| `--additional_filters` | Extra hard filters (e.g. `'design_ALA>0.3'`) | None |
| `--size_buckets` | Max designs per size range (e.g. `'10-20:5' '20-30:10'`) | None |
| `--refolding_rmsd_threshold` | RMSD-based filter threshold (lower = stricter) | Auto |

## Output Structure

```
results/
├── params.json                           # Run parameters
└── boltzgen/
    ├── batches/                           # Independent design batches
    │   └── <batch>/
    │       ├── design/
    │       ├── inverse_folding/
    │       ├── folding/
    │       └── design_folding/
    ├── merged/                            # All batches merged (≈ boltzgen run output)
    └── filtered/
        └── final_ranked_designs/          # Final ranked designs
```

## Re-running Filtering

To re-filter existing results with different parameters without re-running the full pipeline:

```bash
nextflow run boltzgen_filter.nf \
  --run results/boltzgen/merged \
  --budget 20 \
  --alpha 0.05
```

`--config_yaml` and `--protocol` are auto-detected from `params.json` if not specified.

## Complete Examples

### Protein Binder (protein-anything)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method boltzgen \
  --config_yaml my_design.yaml \
  --outdir results \
  --design_name my_protein_binder \
  --protocol protein-anything \
  --num_designs 4 \
  --batch_size 1 \
  --budget 2 \
  -profile local \
  -resume \
  -with-report results/logs/report_$(date +%Y%m%d_%H%M%S).html \
  -with-trace results/logs/trace_$(date +%Y%m%d_%H%M%S).txt
```

### Nanobody Design (nanobody-anything)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method boltzgen \
  --config_yaml nanobody.yaml \
  --outdir results \
  --design_name pdl1-nano \
  --protocol nanobody-anything \
  --num_designs 4 \
  --batch_size 2 \
  --budget 2 \
  --devices 2 \
  --num_workers 2 \
  --alpha 0.2 \
  -profile local \
  -resume \
  -with-report results/logs/report_$(date +%Y%m%d_%H%M%S).html \
  -with-trace results/logs/trace_$(date +%Y%m%d_%H%M%S).txt
```

### HPC Cluster

Add `-profile slurm` and `--slurm_account`:

```bash
-profile slurm \
--slurm_account=ab12
```
