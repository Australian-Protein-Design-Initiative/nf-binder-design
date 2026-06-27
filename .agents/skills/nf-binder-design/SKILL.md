---
name: nf-binder-design
description: >-
  Runs the nf-binder-design Nextflow pipeline for de novo protein binder design.
  Covers RFdiffusion (--method rfd), partial diffusion (--method rfd_partial),
  RFdiffusion3 (--method rfd3), BindCraft (--method bindcraft), Germinal
  (--method germinal), BoltzGen (--method boltzgen), Boltz Pulldown
  (--method boltz_pulldown), and FoldSeek (--method foldseek or --do_foldseek).
  Use when the user wants to design protein binders, nanobodies, or peptides,
  set up or run nf-binder-design, configure HPC/SLURM, or troubleshoot pipeline
  errors.
priority: 5
---

# nf-binder-design

You are an expert at running the `nf-binder-design` Nextflow pipeline for de novo protein binder design.

Your job is to help the user choose the right workflow, construct correct `nextflow run` commands, set up HPC execution, and interpret outputs.

**GitHub**: https://github.com/Australian-Protein-Design-Initiative/nf-binder-design
**Docs**: https://australian-protein-design-initiative.github.io/nf-binder-design/
**Examples**: `examples/pdl1-rfd/`, `examples/pdl1-rfd3/`, `examples/pdl1-bindcraft/`, `examples/pdl1-boltzgen/`, `examples/pdl1-germinal/`

## When to Use This Skill

Activate when the user mentions nf-binder-design, protein binder design, RFdiffusion, RFDiffusion3, ProteinMPNN, AF2 initial guess, Boltz-2, BindCraft, Germinal, BoltzGen, or FoldSeek.

Examples:

- "Design me a binder targeting PDL1"
- "Run nf-binder-design on HPC"
- "Create nanobodies with BoltzGen"
- "Design de novo binders with RFDiffusion3"
- "Troubleshoot RFdiffusion pipeline"

## Quick Start

Create a separate working directory per run:

```bash
mkdir -p {project}/runs/{method}/{target}_hotspots_A56
cd {project}/runs/{method}/{target}_hotspots_A56
```

Invoke via a single entry point:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method <method> [options]
```

From a local git clone: `nextflow run /path/to/nf-binder-design/main.nf --method <method> [options]`

Method-specific help: add `--method <method> --help`

### First-time setup

If the user does not have a working Nextflow installation, follow the ordered checklist in `references/setup-and-hpc.md`:

1. `java -version` (need Java 17+)
2. `nextflow info` (install Nextflow if missing — [Seqera install guide](https://docs.seqera.io/nextflow/install))
3. Check/pull pipeline (`~/.nextflow/assets/...` or `nextflow pull`)
4. Smoke test: `--method rfd --help`
5. Choose profile: `-profile local` (workstation GPU) vs `-profile slurm` / site profile (HPC)

→ **Full install, profile, and HPC details**: `references/setup-and-hpc.md`

## Utility Scripts

Helper scripts live in the pipeline `bin/` directory. Run with [uv](https://docs.astral.sh/uv/):

**Git clone** — from the repo root:

```bash
cd /path/to/nf-binder-design
uv run bin/get_contigs.py --help
```

**After `nextflow pull`** — scripts are in the Nextflow assets cache:

```bash
uv run ~/.nextflow/assets/Australian-Protein-Design-Initiative/nf-binder-design/bin/get_contigs.py --help
```

Use the `bin/` from the same pipeline version you intend to run. See `references/setup-and-hpc.md` for the full utility script table.

| Script | Use for |
|--------|---------|
| `get_contigs.py` | Extract contigs from a cropped PDB (`rfd`, `rfd3`) |
| `trim_to_contigs.py` | Trim a structure to contigs (`bindcraft`, `boltzgen`) |
| `renumber_chains.py` | Renumber residues from 1 (`boltzgen`) |

→ **Target preparation checklist**: `references/target-preparation.md`

## Available Methods

| Method | Description | Key Use Case |
|--------|-------------|-------------|
| `rfd` | RFdiffusion → ProteinMPNN → AF2 initial guess → Boltz-2 refolding | Standard de novo binder design |
| `rfd_partial` | Partial diffusion refinement | Optimise/diversify existing binder designs |
| `rfd3` | RFdiffusion3 → MPNN → RosettaFold3 → optional Boltz refold | De novo design with RFDiffusion3 |
| `bindcraft` | BindCraft in parallel | End-to-end design with built-in validation |
| `germinal` | Germinal in parallel | Antibody and nanobody design via Hydra YAML |
| `boltzgen` | BoltzGen generative model | Protein, peptide, nanobody, or small-molecule binders |
| `boltz_pulldown` | Boltz-2 multimer predictions | Validate designed binders (AlphaPulldown-like) |
| `foldseek` | FoldSeek structural search | Annotate designs against CATH/PDB databases |

Add `--do_foldseek` to `rfd`, `rfd3`, `bindcraft`, or `boltzgen` to run FoldSeek on outputs inline.

## Choosing a Method

1. **New binder from scratch?** → `rfd` (established) or `bindcraft` (end-to-end) or `boltzgen` (protein/peptide/nanobody/small-molecule)
2. **Antibody or nanobody?** → `germinal` (Hydra YAML) or `boltzgen` with `nanobody-anything`
3. **Refine existing designs?** → `rfd_partial`
4. **Validate binder sequences?** → `boltz_pulldown`
5. **Annotate structural similarity?** → `--do_foldseek` or `--method foldseek`

## Common Flags (All Methods)

```bash
--outdir results
-profile local          # Workstation with local GPU
-profile slurm          # SLURM HPC (generic)
-profile slurm,m3       # SLURM + site-specific (see setup-and-hpc.md)
-profile nci_gadi       # NCI Gadi (PBS Pro)
-resume
-with-report results/logs/report_$(date +%Y%m%d_%H%M%S).html
-with-trace results/logs/trace_$(date +%Y%m%d_%H%M%S).txt
-params-file params.json
```

**Profiles:** See `references/setup-and-hpc.md` for the full decision tree (local vs SLURM vs PBS vs HyperQueue vs custom `-c` config).

## Method Examples

One minimal example per method. See the matching `references/*-workflow.md` for advanced options, filters, refolding, and HPC examples.

### rfd (RFdiffusion)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --input_pdb 'input/target.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs 10 \
  -profile local -resume
```

→ `references/rfd-workflow.md`

### rfd_partial (Partial Diffusion)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd_partial \
  --input_pdb 'existing_designs/*.pdb' \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_partial_per_binder 10 \
  --rfd_partial_T "5,15" \
  -profile local -resume
```

Binder is chain A in `rfd` output; residue numbering is 1–N. → `references/rfd-workflow.md`

### rfd3 (RFdiffusion3)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd3 \
  --input_pdb input/target.pdb \
  --outdir results \
  --contigs "A17-131,/0,50-120" \
  --hotspot_res "A56,A115" \
  --rfd3_n_designs 10 \
  -profile local -resume
```

Config mode (atom-level hotspots): `--rfd3_config input/rfd3_config.json`. → `references/rfd3-workflow.md`

### bindcraft (BindCraft)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method bindcraft \
  --input_pdb 'input/target.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56" \
  --binder_length_range "55-120" \
  --bindcraft_n_traj 100 \
  --bindcraft_batch_size 1 \
  -profile local -resume
```

→ `references/bindcraft-workflow.md`

### germinal (Germinal)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method germinal \
  --germinal_config configs/target_vhh.yaml \
  --germinal_pdb_dir pdbs \
  --germinal_experiment_name target_vhh \
  --germinal_n_traj 4 \
  --germinal_batch_size 1 \
  --outdir results \
  -profile local -resume
```

→ `references/germinal-workflow.md`

### boltzgen (BoltzGen)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method boltzgen \
  --config_yaml design.yaml \
  --outdir results \
  --protocol protein-anything \
  --num_designs 100 \
  --batch_size 10 \
  --budget 20 \
  -profile local -resume
```

Protocols: `protein-anything`, `peptide-anything`, `protein-small_molecule`, `nanobody-anything`. → `references/boltzgen-workflow.md`

### boltz_pulldown (Boltz Pulldown)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method boltz_pulldown \
  --input_fasta 'binders/*.fasta' \
  --target_fasta input/target.fasta \
  --outdir results \
  -profile local -resume
```

→ `references/boltz-pulldown-workflow.md`

### foldseek (FoldSeek)

Standalone search on existing structures:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method foldseek \
  --input_pdbs 'results/designs/*.pdb' \
  --outdir results/foldseek \
  -profile local -resume
```

Or add `--do_foldseek` to a design workflow. → `references/foldseek-workflow.md`

## Key Outputs

| Method | Primary outputs |
|--------|-----------------|
| `rfd` / `rfd_partial` | `combined_scores.tsv`, `binders.fasta`, `af2_initial_guess/pdbs/`, `boltz_refold/` (if enabled) |
| `rfd3` | `rfd3/combined_scores.tsv`, `rfdiffusion3/output/`, `rosettafold3/output/` |
| `bindcraft` | `bindcraft/accepted/`, `bindcraft/bindcraft_report.html`, `bindcraft/final_design_stats.csv` |
| `germinal` | `germinal/accepted_designs.csv`, `germinal/accepted/structures/` |
| `boltzgen` | `boltzgen/filtered/final_ranked_designs/`, `boltzgen/merged/` |
| `boltz_pulldown` | `boltz_pulldown/boltz_pulldown.tsv`, `boltz_pulldown_report.html` |
| `foldseek` | `foldseek_results.tsv`, `foldseek_results_annotated.tsv` (CATH databases) |

## Critical Gotchas

- **Quote glob patterns** — `--input_pdb 'input/*.pdb'` (not unquoted globs)
- **Contig syntax** — `rfd`: `"[A18-132/0 65-120]"` (v1); `rfd3`: `"A17-131,/0,50-120"` (v3, comma-separated)
- **~60 GB+ storage** for containers; **NVIDIA GPU** required
- **Rosetta licence** — RFdiffusion and BindCraft use PyRosetta (non-commercial only)
- **One pipeline per directory** — exclusive lock on `.nextflow/cache/`; use separate run dirs
- **Glycans near hotspots** — check N-X-S/T sequons, PDB LINK records, UniProt PTM annotations
- **VRAM** — 12 GB → ~250 residues total; 24 GB → ~400; 32 GB → ~550 (target + max binder)
- **BoltzGen YAML paths** — use absolute paths in YAML configs (container working dir differs)

## Setup and Troubleshooting

→ `references/setup-and-hpc.md` — installation, HPC/SLURM, environment variables
→ `references/nextflow-troubleshooting.md` — monitoring, failure diagnosis, common errors

## Reference Files

Read on demand — do not load all upfront:

| File | Contents |
|------|----------|
| `references/rfd-workflow.md` | RFdiffusion, partial diffusion, filters, refolding, RMSD |
| `references/rfd3-workflow.md` | RFD3 params/config modes, MPNN, RF3, Boltz full-refold |
| `references/bindcraft-workflow.md` | BindCraft parameters, presets, VRAM |
| `references/germinal-workflow.md` | Germinal Hydra config, parallelisation |
| `references/boltzgen-workflow.md` | BoltzGen YAML, protocols, filtering |
| `references/boltz-pulldown-workflow.md` | Boltz Pulldown MSA options |
| `references/foldseek-workflow.md` | `--do_foldseek` flag and `--method foldseek` |
| `references/target-preparation.md` | Target viability and structure preparation |
| `references/setup-and-hpc.md` | Installation and HPC configuration |
| `references/nextflow-troubleshooting.md` | Monitoring and error diagnosis |
