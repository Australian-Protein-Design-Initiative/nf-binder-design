---
name: nf-binder-design
description: >-
  Runs the nf-binder-design Nextflow pipeline for de novo protein binder design.
  Covers RFdiffusion (--method rfd), partial diffusion (--method rfd_partial),
  RFdiffusion3 (--method rfd3), RFD4-Proteina (--method rfd4_proteina),
  BindCraft (--method bindcraft), BindCraft2 (--method bc2), Germinal
  (--method germinal), BoltzGen (--method boltzgen),
  Fold Pulldown (--method fold_pulldown),
  Fold (--method fold), and FoldSeek (--method foldseek or --do_foldseek_search).
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
**Examples**: `examples/pdl1-rfd/`, `examples/pdl1-rfd3/`, `examples/pdl1-bindcraft/`, `examples/pdl1-bc2/`, `examples/pdl1-bc2-mini/`, `examples/bindcraft2-vhh/`, `examples/pdl1-boltzgen/`, `examples/pdl1-germinal/`

## When to Use This Skill

Activate when the user mentions nf-binder-design, protein binder design, RFdiffusion, RFDiffusion3, ProteinMPNN, AF2 initial guess, Boltz-2, BindCraft, BindCraft2/BC2, Germinal, BoltzGen, or FoldSeek.

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
| `rfd4_proteina` | RFD4-Proteina flow matching, sequence and structure in one pass | De novo design with no separate inverse-folding step; the only arm64/GB10 generator |
| `bindcraft` | BindCraft in parallel | End-to-end design with built-in validation |
| `bc2` | BindCraft2 in parallel (alias `bindcraft2`) | End-to-end design with binder modalities (VHH, scFv, Fab, cyclic peptide) and composable design properties |
| `germinal` | Germinal in parallel | Antibody and nanobody design via Hydra YAML |
| `boltzgen` | BoltzGen generative model | Protein, peptide, nanobody, or small-molecule binders |
| `fold_pulldown` | Multi-model target × binder co-fold | Validate designed binders (AlphaPulldown-like) with AF2/Boltz/RF3/Protenix/AF3/OpenFold3/Chai-1/ESMFold2 |
| `fold` | Multi-method structure prediction | Fold FASTA complexes with AF2/Boltz/RF3/Protenix/AF3/OpenFold3/Chai-1/ESMFold2 |
| `foldseek` | FoldSeek structural search | Annotate designs against CATH/PDB databases |

Add `--do_foldseek_search` to `rfd`, `rfd3`, `rfd4_proteina`, `bindcraft`, `bc2`, or `boltzgen` to run FoldSeek on outputs inline.

FoldSeek is used for two different jobs, with separate flags and output directories:
`--do_foldseek_search` searches a **database** for homologs ("has this fold been seen before?"),
while `--do_foldseek_cluster` clusters the run's designs **against each other** ("how many
distinct binders did this run make?"), adding `foldseek_cluster_*` columns to
`combined_scores.tsv`. Clustering is currently wired into `rfd4_proteina`.

### Refold scoring (`--refold_with`)

Every design workflow - `rfd`, `rfd_partial`, `rfd3`, `rfd4_proteina`, `bindcraft`, `germinal` - takes `--refold_with`, which folds selected designs again with one or more independent predictors. It is a held-out second opinion **in addition to** each workflow's own cofold (AF2 initial guess in `rfd`/`rfd_partial`, RF3 in `rfd3`, BindCraft's and Germinal's own AF2 passes), not a replacement for it. That first cofold's scores are what the refold selects on. Off by default.

`rfd4_proteina` is the exception: it co-designs sequence and structure in one pass and has no cofold, so the refold is the only validation and the selection can only read the structural metrics from `--rfd4_filters`.

| Flag | Meaning |
|------|---------|
| `--refold_with` | Comma-separated engines: `af2`, `af2_mono`, `boltz`, `rf3`, `protenix`, `af3`, `openfold3`, `chai1`, `esmfold2`, `esmfold2_fast`. Unset means no refolding |
| `--refold_filter` | Threshold gate on the first cofold's scores, same syntax as `--rfd_filters`, eg `'ipsae>=0.5;pae_interaction<=10'`. Per-design, so it streams |
| `--refold_sort_by` | Rank by one metric, `-` prefix for descending. Does **not** stream. **Must be written with an `=`**: `--refold_sort_by='-ipsae'`, or Nextflow parses the leading `-` as its own option |
| `--refold_max` | Cap. With `--refold_sort_by`, the best N; without it, the first N to be scored |
| `--refold_binder_chain` / `--refold_target_chain` | Only used by `bindcraft` and `germinal` (defaults binder `B`, target `A`); other workflows derive the chains |
| `--refold_target_fasta` | Refold against a different target than the designs were built against - the full-length protein when the design used a truncated construct. FASTA of **exactly one record**. Unset means the target chain read back out of the design structures |
| `--templates` / `--target_msa` / `--create_target_msa` | The refold target's template structures, an A3M built elsewhere, or an MSA search. These carry no `--refold_` prefix: in a design workflow the fold stack only ever runs for the refold |

Sampling and aggregation, shared with `fold` / `fold_pulldown`:

| Flag | Default | Meaning |
|------|---------|---------|
| `--n_predictions` | 5 | Structures per design per engine. Run as N separate single-sample jobs with distinct seeds (42, 43, ...), not one job of N samples: samples from one seed share the trunk and differ only in diffusion noise, separate seeds are independent trajectories |
| `--refold_<engine>_batch_size` | 1 | Samples per fold job for that engine (`--refold_boltz_batch_size`, `--refold_rf3_batch_size`, ...; also `--protenix_batch_size`). `--n_predictions 6 --refold_boltz_batch_size 3` is 2 seeds of 3 samples |
| `--fold_seed_list` | `true` | Hand the whole seed list to engines that take several seeds per invocation (`esmfold2`, `esmfold2_fast`, `openfold3`, `protenix`, `af3` via `--af3_seeds`), saving repeated model loads. Boltz and RF3 run one job per seed regardless. `false` disables grouping everywhere |
| `--refold_structure_stat` | `selected` | How structural and interface columns reduce across samples. `selected` takes them all from the one sample ranked best by `--consensus_metric`, so they describe a single coherent structure; `min`/`mean`/`max`/`first` apply that operator per column instead. Confidence columns always get direct aggregates (`_max`, `_mean`, `n`) |
| `--consensus_metric` | `ipsae` | The metric that ranks samples for `selected`, and the basis of `consensus_z` |
| `--refold_metrics` | iptm/ipsae/plddt/ptm/pae | Which summary metrics reach `combined_scores.tsv` |
| `--refold_binder_monomer` | `true` | Also predict the binder alone, giving `<tool>_monomer_vs_design_rmsd_all` and `<tool>_monomer_vs_complex_rmsd_all` |
| `--refold_monomer_single_sequence` | `true` | Fold those monomers without an MSA |
| `--refold_interface_scoring` | `true` | BindCraft interface metrics on the refolded complexes (`<tool>_interface_dG`, `<tool>_interface_dSASA`, ...) |
| `--refold_dockq` | `true` | DockQ and its components, interface lDDT and a clash count for the target-aligned comparison (`<tool>_target_aligned_binder_dockq`, `_dockq_fnat`, `_dockq_irmsd`, `_interface_lddt`, ...). DockQ's receptor is the target and its ligand the binder, which differs from the reference implementation when the binder is the larger chain |

`refold_sample_scores.tsv` keeps the unaggregated per-sample values; `<tool>_selected_batch` and `<tool>_selected_model` in `combined_scores.tsv` name the sample the structural metrics came from.

**Engine parameter names are prefixed for the refold where a prefixed copy exists**, and that prefixed name is then the only one read: `--refold_boltz_batch_size`, `--refold_boltz_seed`, the whole `--refold_rf3_*` set, `--refold_esmfold2_*`, `--refold_chai1_*`, `--refold_openfold3_*`, `--refold_af3_batch_size`. Plain `--boltz_*` / `--rf3_*` do nothing to the refold - `rfd3`'s `--rf3_batch_size` is MPNN designs per RF3 task, a different meaning in the same flat namespace. Where no prefixed copy is declared the plain name is read, which is the case for `--protenix_*`, `--af2_*` and `--af3_seeds`. Full table in `references/fold-stack.md`.

With a target longer than the design's, `<tool>_complex_rmsd_all` and `<tool>_complex_tm_score_all` come back as `nan` (different residue counts); the `_pruned` variants and `<tool>_target_aligned_binder_rmsd_*`, the pose metric, still work.

Changing `--refold_filter` / `--refold_sort_by` / `--refold_max` on a `-resume` needs a fresh `--outdir`: the refold tables are gathered with `collectFile(storeDir:)`, which reuses the file already there, so the scores and the selection report end up describing different design sets with no error.

Results go to `<outdir>/<method>/refold/`: `refold_scores.tsv` (merged into `combined_scores.tsv` as engine-prefixed columns `<tool>_iptm_max`, `<tool>_ipsae_max`, `<tool>_target_aligned_binder_rmsd_all`, `<tool>_target_aligned_binder_dockq`, `<tool>_interface_dG`, ...), `refold_sample_scores.tsv` and `refold_selection.tsv` (every design's selection verdict, rank and sort value).

**Renamed/removed.** Deprecated but working, with a warning: `--refold_af2ig_filters` (→ `--refold_filter`), `--refold_oracles` (→ `--refold_with`), `--refold_oracles_max` (→ `--refold_max`). Removed, and setting one stops the run with its replacement named: `--refold_create_target_msa` (→ `--create_target_msa`), `--refold_use_msa_server` (→ `--use_msa_server`), `--refold_target_templates` (→ `--templates`), `--refold_alignment` (→ `--target_msa`), `--refold_oracle_metrics` (→ `--refold_metrics`), `--refold_oracle_target_id` (→ `--refold_target_id`), `--boltz_refold_batch_size` (→ `--refold_boltz_batch_size`, now samples per fold job). `--refold_target_fasta` is **not** removed - it kept its name, see the table above.

The old `--full_refold_*` family (`--full_refold_with`, `--full_refold_filter_sort`, ...) is now `--refold_*` (`--full_refold_filter_sort` → `--refold_sort_by`), but unlike the names above these are not declared anywhere, so setting one draws only Nextflow's generic unknown-parameter warning and the run carries on with the default. Check for them by hand when porting an old command line.

**`rfd3` behaviour change.** `--full_refold_filter_sort` defaulted to `pair_pae_min`, so `--refold_max` implicitly took the best designs. `--refold_sort_by` defaults to unset, so a cap alone now takes the first N to finish. Add `--refold_sort_by='pair_pae_min'` to keep the old behaviour.

## Choosing a Method

1. **New binder from scratch?** → `rfd` (established) or `bindcraft` / `bc2` (end-to-end) or `boltzgen` (protein/peptide/nanobody/small-molecule) or `rfd4_proteina` (sequence and structure co-designed in one pass)
2. **Antibody or nanobody?** → `germinal` (Hydra YAML), `bc2` with `--bc2_modality VHH|scFv|Fab`, or `boltzgen` with `nanobody-anything`
3. **Refine existing designs?** → `rfd_partial`
4. **Validate binder sequences?** → `fold_pulldown` (use `--methods boltz` for a Boltz-only run)
5. **Fold arbitrary FASTAs?** → `fold`
6. **Annotate structural similarity?** → `--do_foldseek_search` or `--method foldseek`
7. **Measure design diversity?** → `--do_foldseek_cluster`

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

### rfd4_proteina (RFD4-Proteina)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd4_proteina \
  --input_pdb input/target.cif \
  --contigs "50-120, /0, A17-131" \
  --hotspot_res "A56,A115" \
  --rfd4_n_designs 20 \
  --outdir results \
  -profile local -resume
```

Condition-spec mode, which is the full conditioning vocabulary: `--rfd4_config pdl1_binder.yaml`, with no `--input_pdb`. Residue numbers in contigs and conditions are mmCIF `label_seq_id`, not author numbering.

In parameter mode the generated spec names the designed segment `BINDER` and excludes cysteines from it (`C_SEQ_PROB_CYS`). `--contigs` may name it instead (`"BINDER: 50-120, /0, A17-131"`). To allow cysteines, copy the generated spec, drop that condition, and rerun with `--rfd4_config`. → `references/rfd4-proteina-workflow.md`

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

### bc2 (BindCraft2)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method bc2 \
  --input_pdb 'input/target.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56" \
  --binder_length_range "55-120" \
  --bc2_n_traj 100 \
  --bc2_batch_size 1 \
  -profile local -resume
```

Config mode: `--bc2_config campaign.json`. Any BindCraft2 setting: `--bc2_set 'KEY=VALUE;KEY=VALUE'`.
→ `references/bc2-workflow.md`

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

Binders can instead come from existing complexes, which then become the reference each
prediction is compared against (DockQ, interface lDDT, target-aligned binder RMSD):

```bash
  --complex_structures 'complexes/*.pdb' --binder_chains B
```

`--targets` is optional with it - without one the targets are read back out of the
complexes. `--target_chains` picks the superposition reference when a complex has more
than two chains; `--target_chains auto` matches it against `--targets` by sequence.
Results land in `fold_pulldown/complex_rmsd.tsv`.

Samples per complex, job splitting and seeds: `--n_predictions` defaults to 5 here.
→ `references/fold-stack.md`

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

→ `references/fold-stack.md` for `--n_predictions`, `--*_batch_size` and seeds.

### foldseek (FoldSeek)

Standalone search on existing structures:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method foldseek \
  --input_pdbs 'results/designs/*.pdb' \
  --outdir results/foldseek \
  -profile local -resume
```

Or add `--do_foldseek_search` to a design workflow. → `references/foldseek-workflow.md`

## Key Outputs

| Method | Primary outputs |
|--------|-----------------|
| `rfd` / `rfd_partial` | `combined_scores.tsv`, `binders.fasta`, `af2_initial_guess/pdbs/`, `rfd/refold/` (if `--refold_with`) |
| `rfd3` | `rfd3/combined_scores.tsv`, `rfdiffusion3/output/`, `rosettafold3/output/`, `rfd3/refold/` (if `--refold_with`) |
| `rfd4_proteina` | `rfd4_proteina/combined_scores.tsv`, `rfd4_proteina/designs/` (mmCIF), `rfd4_proteina/batches/batch<n>/` (raw rfd4-generate run dir, inline metrics and refolds), `rfd4_proteina/binders.fasta`, `rfd4_proteina/filter/` (if `--rfd4_filters`), `rfd4_proteina/refold/` (if `--refold_with`) |
| `bindcraft` | `bindcraft/accepted/`, `bindcraft/bindcraft_report.html`, `bindcraft/final_design_stats.csv`, `bindcraft/refold/` (if `--refold_with`) |
| `bc2` | `bc2/ranked/` (mmCIF), `bc2/bc2_report.html`, `bc2/summary.csv`, `bc2/ranked.csv` |
| `germinal` | `germinal/accepted_designs.csv`, `germinal/accepted/structures/`, `germinal/refold/` (if `--refold_with`) |
| `boltzgen` | `boltzgen/filtered/final_ranked_designs/`, `boltzgen/merged/` |
| `fold_pulldown` | `fold_pulldown/fold_pulldown_scores.tsv`, `fold_pulldown_summary.tsv`, `fold_pulldown_report.html`, `complex_rmsd.tsv` (if `--complex_structures`) |
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
| `references/rfd4-proteina-workflow.md` | RFD4-Proteina condition specs, contigs, filters, arm64/GB10 notes |
| `references/bindcraft-workflow.md` | BindCraft parameters, presets, VRAM |
| `references/bc2-workflow.md` | BindCraft2 config/params modes, `--bc2_set`, modalities, properties |
| `references/germinal-workflow.md` | Germinal Hydra config, parallelisation |
| `references/boltzgen-workflow.md` | BoltzGen YAML, protocols, filtering |
| `references/fold-stack.md` | Seeds, samples and batch sizes shared by `fold`, `fold_pulldown` and the refold |
| `references/foldseek-workflow.md` | `--do_foldseek_search` flag and `--method foldseek` |
| `references/target-preparation.md` | Target viability and structure preparation |
| `references/setup-and-hpc.md` | Installation and HPC configuration |
| `references/nextflow-troubleshooting.md` | Monitoring and error diagnosis |
