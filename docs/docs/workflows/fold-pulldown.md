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
| `fold_pulldown/fold_pulldown_summary.tsv` | One row per `(target, binder, tool)` with mean/median/max/sd for `iptm` and `ipsae`, per-target within-tool z-scores, cross-tool `consensus_z`, and the `z_basis` / `n_pool` / `z_pool_small` provenance of that z |
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
| `--methods` | `boltz` | Comma-separated: `af2`, `af2_mono`, `boltz`, `rf3`, `protenix` |
| `--msa_method` | `jackhmmer_af2` | `jackhmmer_af2` or `mmseqs2_colabfold` |
| `--create_target_msa` | `false` | Build MSA for each target |
| `--create_binder_msa` | `false` | Build MSA for each binder (usually leave off for de novo binders) |
| `--n_predictions` | unset | Samples per complex per method (engine defaults if unset) |
| `--skip_engens` | `true` | EnGens is off by default (would emit N × M reports) |
| `--consensus_metric` | `ipsae` | `ipsae` or `iptm`; which metric's per-tool z-scores are averaged into `consensus_z` |
| `--z_stat` | `max` | `max` or `mean`; the per-complex statistic over samples that gets standardised |
| `--z_scope` | `target` | `target` standardises within `(target, tool)`; `global` pools all targets into one distribution |
| `--min_pool` | `10` | Pools holding fewer complexes than this are flagged `z_pool_small` in the summary and warned about on stderr |

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
compared on a common footing inside each model, and a **`consensus_z`** (the mean
over tools of one metric's per-tool z for that complex) for a cross-model ranking.
Both z-columns are always written, whichever metric `consensus_z` is built from.

Three defaults decide that ranking, and each can be reverted:

- **`--consensus_metric ipsae`.** ipSAE ([Dunbrack
  2025](https://doi.org/10.1101/2025.02.10.637595)) normalises by interface size and
  isolates the interface from whole-complex confidence, which `iptm` and `ptm` mix
  together. Pass `--consensus_metric iptm` to rank on ipTM instead.
- **`--z_stat max`.** One diffusion sample per complex is noisy, so the statistic
  worth standardising is the best of the samples rather than their average. This is
  also what makes `--n_predictions > 1` pay for itself. Pass `--z_stat mean` for the
  average.
- **`--z_scope target`.** Raw co-folding scores are not comparable across targets,
  and target difficulty generally varies more than design quality does within one
  target, so pooling several targets into one distribution makes a complex rank
  partly on which target it was paired with. Pass `--z_scope global` to pool them.

Each summary row records `z_basis` (metric, statistic and scope, e.g.
`ipsae_max/target`), `n_pool`, the number of complexes its z-score was computed
over, and `z_pool_small`, `True` where that pool was smaller than `--min_pool`.
Watch those two: with *k* complexes in a pool the largest possible absolute z is
(*k*−1)/√*k*, so a two-complex pool can only ever report ±0.707 and the z-score
carries the ordering and nothing else. A saturated small-pool z looks like a
mediocre one. The flag is also printed as an stderr warning, but stderr from a
Nextflow task lands in the work directory, so the column is what a downstream
consumer should filter on.

For custom statistics (Mann–Whitney, mixed models, score calibration), use
`fold_pulldown_scores.tsv` / `fold_pulldown_summary.tsv` directly — the HTML
report stays deliberately simple.

## Related

- [Boltz Pulldown](boltz-pulldown.md) — Boltz-2-only predecessor
- [Fold](fold.md) — fold individual FASTA complexes (no target × binder fan-out)
