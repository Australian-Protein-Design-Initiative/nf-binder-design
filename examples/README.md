# Example configurations

Runnable examples for each `--method`, using tracked `run*.sh` scripts in this directory.
Each example is set up for short test runs on a local GPU workstation; adjust design counts, batch sizes, and `nextflow.config` for production or HPC (SLURM).

Run scripts from the example directory, e.g. `cd examples/pdl1-rfd && ./run-local-simple.sh`.

Workflow documentation: https://australian-protein-design-initiative.github.io/nf-binder-design/

## RFdiffusion (`--method rfd`)

### `pdl1-rfd`

PD-L1 binder design: RFdiffusion → ProteinMPNN → AlphaFold2 initial guess → optional Boltz-2 refolding.

| Script | Description |
|--------|-------------|
| `run-local-simple.sh` | Minimal local run |
| `run-local-full.sh` | Additional filters, refolding, and reporting options |
| `run-dual-gpu.sh` | Local run with dual-GPU `nextflow.config` |
| `run-local-full-dual-gpu.sh` | Full options with dual-GPU config |
| `run-m3-full.sh` | SLURM/M3-style run with refold, filters, and FoldSeek |

### `pdl1-rfd-partial`

Partial diffusion refinement of existing binder designs (`--method rfd_partial`).

| Script | Description |
|--------|-------------|
| `run-local.sh` | Local partial diffusion with refolding options |
| `run-dual-gpu.sh` | Dual-GPU local run |

### `egfr-rfd-hypermpnn`

EGFR binder design using HyperMPNN weights for inverse folding.

| Script | Description |
|--------|-------------|
| `run.sh` | Local run |
| `run-dual-gpu.sh` | Dual-GPU local run |

## RFdiffusion3 (`--method rfd3`)

### `pdl1-rfd3`

PD-L1 binder design with RFDiffusion3 → MPNN → RosettaFold3.

| Script | Description |
|--------|-------------|
| `run-local-simple.sh` | Params mode (`--contigs`, `--hotspot_res`) |
| `run-local-config.sh` | Config-file mode (`--rfd3_config`) |
| `run-local-full.sh` | Extended options (refold, filters, etc.) |
| `run-dual-gpu.sh` | Dual-GPU local run |

## BindCraft (`--method bindcraft`)

### `pdl1-bindcraft`

PD-L1 protein binder design.

| Script | Description |
|--------|-------------|
| `run.sh` | Local run |
| `run-gpu0.sh` | Pin to GPU 0 |
| `run-dual-gpu.sh` | Dual-GPU local run |
| `run-m3.sh` | SLURM/M3 example |

### `egfr-bindcraft`

EGFR binder design.

| Script | Description |
|--------|-------------|
| `run.sh` | Local run |
| `run-dual-gpu.sh` | Dual-GPU local run |

### `mdm2-bindcraft-peptide`

Peptide binder design for MDM2 (based on [Filus et al., 2025](https://doi.org/10.1101/2025.07.23.666285)).

| Script | Description |
|--------|-------------|
| `run.sh` | Local run |
| `run-dual-gpu.sh` | Dual-GPU local run |

## BoltzGen (`--method boltzgen`)

### `pdl1-boltzgen`

PD-L1 protein binder (`protein-anything`) with optional FoldSeek.

| Script | Description |
|--------|-------------|
| `run-local.sh` | Local run |

### `boltzgen-protein`

Protein binder against 1G13 (`protein-anything`).

| Script | Description |
|--------|-------------|
| `run-local.sh` | Local run |
| `run-local-dual-gpu.sh` | Dual-GPU local run |
| `run-filter.sh` | Run with custom filtering |

### `boltzgen-cyclic-peptide`

Cyclic peptide binder against 2VSM (`peptide-anything`).

| Script | Description |
|--------|-------------|
| `run-local.sh` | Local run |
| `run-local-dual-gpu.sh` | Dual-GPU local run |

### `boltzgen-small-molecule-binder`

Small-molecule binder for PFOA (`protein-small_molecule`).

| Script | Description |
|--------|-------------|
| `run-local.sh` | Local run |
| `run-dual-gpu.sh` | Dual-GPU local run |

### `boltzgen-nanobody`

Nanobody binder design (`nanobody-anything`).

| Script | Description |
|--------|-------------|
| `run-local.sh` | Local run |

## Boltz Pulldown (`--method boltz_pulldown`)

No tracked example directory yet. See the [Boltz Pulldown workflow documentation](https://australian-protein-design-initiative.github.io/nf-binder-design/workflows/boltz-pulldown/).
