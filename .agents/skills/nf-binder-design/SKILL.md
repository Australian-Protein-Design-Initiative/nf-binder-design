---
name: nf-binder-design
description: >-
  Runs the nf-binder-design Nextflow pipeline for de novo protein binder design.
  Covers RFdiffusion (--method rfd), partial diffusion (--method rfd_partial),
  RFdiffusion3 (--method rfd3), BindCraft (--method bindcraft), Germinal
  (--method germinal), BoltzGen (--method boltzgen),
  Fold Pulldown (--method fold_pulldown),
  Fold (--method fold), and FoldSeek (--method foldseek or --do_foldseek).
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
| `rfd` | RFdiffusion → ProteinMPNN → AF2 initial guess | Standard de novo binder design |
| `rfd_partial` | Partial diffusion refinement | Optimise/diversify existing binder designs |
| `rfd3` | RFdiffusion3 → MPNN → RosettaFold3 | De novo design with RFDiffusion3 |
| `bindcraft` | BindCraft in parallel | End-to-end design with built-in validation |
| `germinal` | Germinal in parallel | Antibody and nanobody design via Hydra YAML |
| `boltzgen` | BoltzGen generative model | Protein, peptide, nanobody, or small-molecule binders |
| `fold_pulldown` | Multi-model target × binder co-fold | Validate designed binders (AlphaPulldown-like) with AF2/Boltz/RF3/Protenix/AF3/OpenFold3/Chai-1/ESMFold2 |
| `fold` | Multi-method structure prediction | Fold FASTA complexes with AF2/Boltz/RF3/Protenix/AF3/OpenFold3/Chai-1/ESMFold2 |
| `foldseek` | FoldSeek structural search | Annotate designs against CATH/PDB databases |

Add `--do_foldseek` to `rfd`, `rfd3`, `bindcraft`, or `boltzgen` to run FoldSeek on outputs inline.

### Refold scoring (`--refold_with`)

Every design workflow - `rfd`, `rfd_partial`, `rfd3`, `rfd4_proteina`, `bindcraft`, `germinal` - takes `--refold_with`, which folds selected designs again with one or more independent predictors. It is a held-out second opinion **in addition to** each workflow's own cofold (AF2 initial guess in `rfd`/`rfd_partial`, RF3 in `rfd3`, BindCraft's and Germinal's own AF2 passes), not a replacement for it. That first cofold's scores are what the refold selects on. Off by default.

| Flag | Meaning |
|------|---------|
| `--refold_with` | Comma-separated engines: `af2`, `af2_mono`, `boltz`, `rf3`, `protenix`, `af3`, `openfold3`, `chai1`, `esmfold2`, `esmfold2_fast`. Unset means no refolding |
| `--refold_filter` | Threshold gate on the first cofold's scores, same syntax as `--rfd_filters`, eg `'ipsae>=0.5;pae_interaction<=10'`. Per-design, so it streams |
| `--refold_sort_by` | Rank by one metric, `-` prefix for descending. Does **not** stream. **Must be written with an `=`**: `--refold_sort_by='-ipsae'`, or Nextflow parses the leading `-` as its own option |
| `--refold_max` | Cap. With `--refold_sort_by`, the best N; without it, the first N to be scored |
| `--refold_binder_chain` / `--refold_target_chain` | Only used by `bindcraft` and `germinal` (defaults binder `B`, target `A`); other workflows derive the chains |
| `--refold_target_fasta` | Refold against a different target than the designs were built against - the full-length protein when the design used a truncated construct. FASTA of **exactly one record**. Unset means the target chain read back out of the design structures |
| `--templates` / `--target_msa` / `--create_target_msa` | The refold target's template structures, an A3M built elsewhere, or an MSA search. These carry no `--refold_` prefix: in a design workflow the fold stack only ever runs for the refold |

With a target longer than the design's, `<tool>_complex_rmsd_all` and `<tool>_complex_tm_score_all` come back as `nan` (different residue counts); the `_pruned` variants and `<tool>_target_aligned_binder_rmsd_*`, the pose metric, still work.

Changing `--refold_filter` / `--refold_sort_by` / `--refold_max` on a `-resume` needs a fresh `--outdir`: the refold tables are gathered with `collectFile(storeDir:)`, which reuses the file already there, so the scores and the selection report end up describing different design sets with no error.

Results go to `<outdir>/<method>/refold/`: `refold_scores.tsv` (merged into `combined_scores.tsv` as engine-prefixed columns `<tool>_iptm_max`, `<tool>_ipsae_max`, `<tool>_target_aligned_binder_rmsd_all`, `<tool>_interface_dG`, ...), `refold_sample_scores.tsv` and `refold_selection.tsv` (every design's selection verdict, rank and sort value).

**Renamed/removed.** The whole `--full_refold_*` family is now `--refold_*` (`--full_refold_filter_sort` → `--refold_sort_by`). Deprecated but working, with a warning: `--refold_af2ig_filters` (→ `--refold_filter`), `--refold_oracles` (→ `--refold_with`), `--refold_oracles_max` (→ `--refold_max`). Removed, and setting one stops the run: `--refold_create_target_msa` (→ `--create_target_msa`), `--refold_use_msa_server` (→ `--use_msa_server`), `--refold_target_templates` (→ `--templates`), `--refold_alignment` (→ `--target_msa`), `--boltz_refold_batch_size`. `--refold_target_fasta` is **not** removed - it kept its name, see the table above.

**`rfd3` behaviour change.** `--full_refold_filter_sort` defaulted to `pair_pae_min`, so `--refold_max` implicitly took the best designs. `--refold_sort_by` defaults to unset, so a cap alone now takes the first N to finish. Add `--refold_sort_by='pair_pae_min'` to keep the old behaviour.

## Choosing a Method

1. **New binder from scratch?** → `rfd` (established) or `bindcraft` (end-to-end) or `boltzgen` (protein/peptide/nanobody/small-molecule)
2. **Antibody or nanobody?** → `germinal` (Hydra YAML) or `boltzgen` with `nanobody-anything`
3. **Refine existing designs?** → `rfd_partial`
4. **Validate binder sequences?** → `fold_pulldown` (use `--methods boltz` for a Boltz-only run)
5. **Fold arbitrary FASTAs?** → `fold`
6. **Annotate structural similarity?** → `--do_foldseek` or `--method foldseek`

## Common Flags (All Methods)

```bash
--outdir results
-profile local          # Workstation with local GPU
-profile slurm          # SLURM HPC (generic)
-profile slurm,m3       # SLURM + site-specific (see setup-and-hpc.md)
-profile nci_gadi       # NCI Gadi (PBS Pro)
-resume
-params-file params.json
```

Reports are written automatically under `${outdir}/logs/` (`report_*.html`, `trace_*.txt`, `timeline_*.html`, `dag_*.html`) — no `-with-report` / `-with-trace` needed.

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

### fold_pulldown (multi-model pulldown)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method fold_pulldown \
  --targets targets.fasta \
  --binders binders.fasta \
  --methods boltz,rf3,protenix \
  --create_target_msa true \
  --msa_method jackhmmer_af2 \
  --outdir results \
  -profile slurm,m3 -resume
```

### fold (multi-method structure prediction)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method fold \
  --input 'input/*.fasta' \
  --methods af2,boltz,rf3,protenix \
  --msa_method jackhmmer_af2 \
  --outdir results \
  -profile slurm,m3 -resume
```

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
| `rfd` / `rfd_partial` | `combined_scores.tsv`, `binders.fasta`, `af2_initial_guess/pdbs/`, `rfd/refold/` (if `--refold_with`) |
| `rfd3` | `rfd3/combined_scores.tsv`, `rfdiffusion3/output/`, `rosettafold3/output/`, `rfd3/refold/` (if `--refold_with`) |
| `bindcraft` | `bindcraft/accepted/`, `bindcraft/bindcraft_report.html`, `bindcraft/final_design_stats.csv`, `bindcraft/refold/` (if `--refold_with`) |
| `germinal` | `germinal/accepted_designs.csv`, `germinal/accepted/structures/`, `germinal/refold/` (if `--refold_with`) |
| `boltzgen` | `boltzgen/filtered/final_ranked_designs/`, `boltzgen/merged/` |
| `fold_pulldown` | `fold_pulldown/fold_pulldown_scores.tsv`, `fold_pulldown_summary.tsv`, `fold_pulldown_report.html` |
| `fold` | `fold/fold_scores.tsv`, `fold/predictions/`, optional `engens/` |
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
| `references/rfd-workflow.md` | RFdiffusion, partial diffusion, filters, refold scoring, RMSD |
| `references/rfd3-workflow.md` | RFD3 params/config modes, MPNN, RF3, refold scoring |
| `references/bindcraft-workflow.md` | BindCraft parameters, presets, VRAM |
| `references/germinal-workflow.md` | Germinal Hydra config, parallelisation |
| `references/boltzgen-workflow.md` | BoltzGen YAML, protocols, filtering |
| `references/foldseek-workflow.md` | `--do_foldseek` flag and `--method foldseek` |
| `references/target-preparation.md` | Target viability and structure preparation |
| `references/setup-and-hpc.md` | Installation and HPC configuration |
| `references/nextflow-troubleshooting.md` | Monitoring and error diagnosis |
