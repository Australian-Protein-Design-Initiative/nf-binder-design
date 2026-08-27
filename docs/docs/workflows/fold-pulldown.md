# Fold Pulldown

Multi-model target × binder pulldown, in the spirit of
[Boltz Pulldown](boltz-pulldown.md) but using the shared fold engines
(AlphaFold2, Boltz-2, RosettaFold3, Protenix).

## Overview

For every target–binder pair the pipeline co-folds the complex with each
selected predictor (optionally with multiple seeds / samples per model), then
writes:

| File | Contents |
|------|----------|
| `fold_pulldown/fold_pulldown_scores.tsv` | One row per predicted structure (canonical fold scores + `target`, `binder`) |
| `fold_pulldown/fold_pulldown_summary.tsv` | One row per `(target, binder, tool)` with mean/median/max/sd for `iptm` and `ipsae`, within-tool z-scores, and cross-tool `consensus_z` |
| `fold_pulldown/fold_pulldown_report.html` | Simple Quarto overview (boxplots, heatmaps, top hits, cross-tool Spearman) |

MSA cost is **O(N_targets + N_binders)**, not O(N × M): each sequence is searched
once. Binders are treated as having no useful homologs (query-only MSA unless
`--create_binder_msa` is set). There is **no cross-chain MSA pairing** — this is
intentional for designed binders.

## Command-line options

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method fold_pulldown --help
```

### Required

| Flag | Description |
|------|-------------|
| `--targets` | FASTA of target sequences |
| `--binders` | FASTA of binder sequences |

### Common options

| Flag | Default | Description |
|------|---------|-------------|
| `--methods` | `boltz` | Comma-separated: `af2`, `boltz`, `rf3`, `protenix` |
| `--msa_method` | `jackhmmer_af2` | `jackhmmer_af2` or `mmseqs2_colabfold` |
| `--create_target_msa` | `false` | Build MSA for each target |
| `--create_binder_msa` | `false` | Build MSA for each binder (usually leave off for de novo binders) |
| `--n_predictions` | unset | Samples per complex per method (engine defaults if unset) |
| `--skip_engens` | `true` | EnGens is off by default (would emit N × M reports) |

AF2 needs the 2021 DB snapshot with `uniprot/` (default `--af2_db_path` points at
`alphafold_20211129`). Target MSAs for AF2 are built once (jackhmmer dir, or
ColabFold/mmseqs2 a3m materialised into AF2 per-chain files) and assembled into a
multimer tree with a query-only binder chain plus a `features.pkl` (the predict
stage loads that pickle; it does not rebuild features from the raw MSA files).

## Example

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method fold_pulldown \
  --targets targets.fasta \
  --binders binders.fasta \
  --methods boltz,rf3,protenix \
  --create_target_msa true \
  --msa_method jackhmmer_af2 \
  --n_predictions 5 \
  --outdir results \
  -profile slurm,m3
```

## Interpreting scores

Different predictors have different absolute score scales. The summary table
provides **within-tool z-scores** (`iptm_z`, `ipsae_z`) so complexes can be
compared on a common footing inside each model, and a **`consensus_z`**
(mean of per-tool z for that complex) for a simple cross-model ranking.

For custom statistics (Mann–Whitney, mixed models, score calibration), use
`fold_pulldown_scores.tsv` / `fold_pulldown_summary.tsv` directly — the HTML
report stays deliberately simple.

## Related

- [Boltz Pulldown](boltz-pulldown.md) — Boltz-2-only predecessor
- [Fold](fold.md) — fold individual FASTA complexes (no target × binder fan-out)
