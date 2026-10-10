# Fold Pulldown

[![Fold Pulldown workflow](../images/fold_pulldown_metro_map.svg)](../images/fold_pulldown_metro_map.svg){target="_blank" rel="noopener" title="Open the full-size diagram in a new tab"}

An AlphaPulldown-like target × binder pulldown across the shared fold engines
(AlphaFold2, Boltz-2, RosettaFold3, Protenix, AlphaFold3, OpenFold3, ESMFold2,
ESMFold2-Fast).

!!! info "Previously `--method boltz_pulldown`"

    This workflow began as Boltz Pulldown, which could only co-fold with
    Boltz-2. It has since been expanded to run any of the fold engines above,
    and gained per-chain template matching, multiple samples per pair and a
    cross-engine score summary. `--method boltz_pulldown` has been removed; use
    `--method fold_pulldown --methods boltz` for the equivalent Boltz-only run.

## Overview

For every target–binder pair the pipeline co-folds the complex with each
selected predictor (optionally with multiple seeds / samples per model), then
writes:

| File | Contents |
|------|----------|
| `fold_pulldown/fold_pulldown_scores.tsv` | One row per predicted structure (canonical fold scores + `target`, `binder`) |
| `fold_pulldown/fold_pulldown_summary.tsv` | One row per `(target, binder, tool)` with mean/median/max/sd for `iptm` and `ipsae`, per-target within-tool z-scores, cross-tool `consensus_z`, and the `z_basis` / `n_pool` / `z_pool_small` provenance of that z |
| `fold_pulldown/fold_pulldown_report.html` | Simple Quarto overview (boxplots, heatmaps, top hits, cross-tool Spearman) |

## Setup

Setup is the same as for [Fold](fold.md#setup). Engine containers and their
weights are pulled automatically, except the AlphaFold3 weights
([AlphaFold3 weights](fold.md#alphafold3-weights)). Local sequence databases are
needed only for `--msa_method jackhmmer_af2`, or for `--msa_method mmseqs2_colabfold`
without `--use_remote_server true`. See [Fold databases](../extra/fold-databases.md).

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

## Multiple sequence alignments

MSA cost is **O(N_targets + N_binders)**, not O(N × M): each sequence is searched
once. Binders are treated as having no useful homologs (query-only MSA unless
`--create_binder_msa` is set). The pipeline never builds a joint paired
alignment; each chain is rendered on its own, and an engine pairs rows only in
the cases below.

FASTA header ids are used as per-chain filenames and are sanitised (non
`[a-zA-Z0-9_.-]` characters become `_`) before use. The pipeline will fail early with a warning is there are duplicate ids within `--targets` or `--binders`, or an id is present in both `--targets` and `--binders` (after sanitisation). Since output complexes are named with `<target>_and_<binder>`, ids that contain `_and_` may be rejected if they would be ambiguous.

### Taxonomic pairing

Pairing is done by the fold engine, from per-chain files, and only when **both**
chains have homologs that share a taxonomy id. For _de novo_ binder 'pulldowns' where we typically don't use a binder MSA (`--create_binder_msa false`), the taxonomy is not used. The pairing process uses the same heuristics as in
[Fold](fold.md#paired-msas-how-each-engine-differs): RF3 `TaxID=` a3m,
Protenix species-mnemonic a3m, Boltz `key,sequence` CSV, and the AF3 /
OpenFold3 equivalents.

| Condition | Pairing |
|-----------|---------|
| `--create_binder_msa false` (the default) | None. Chain B is the binder sequence alone. Boltz is also given `msa: empty` for chain B, including when `--use_msa_server` is set, so the server does not fetch a binder MSA either. |
| `--create_target_msa false` | None for the MSAs this pipeline builds: chain A is a single sequence too. The exception is Boltz with both `--use_msa_server` and `--create_binder_msa true`, which fetches and pairs both chains itself. |
| Both `--create_target_msa` and `--create_binder_msa`, `--msa_method jackhmmer_af2`, and `boltz`, `rf3`, `protenix`, `af3`, `openfold3` or `esmfold2` | The engine pairs shared taxids. Jackhmmer headers carry `TaxID=` / species mnemonics; each chain is still written independently and the engine matches them. |
| Same flags, but `af2` or `af2_mono` | None. Those engines always write the binder chain query-only, whatever `--create_binder_msa` is. If they are the only selected `--methods`, the binder search is skipped. |
| Both MSAs requested, `--msa_method mmseqs2_colabfold` | None. ColabFold headers carry no taxonomy, so the renders stay unpaired. Boltz with `--use_msa_server` and `--create_binder_msa true` still fetches and pairs its own MSAs. |
| `esmfold2_fast`, or `esmfold2` with `--esmfold2_single_sequence` | None. Those runs fold from sequence alone. |

## Template structures

When a target's structure is known, pass it with `--templates` (a directory or
glob of `.pdb` / `.cif` files, e.g. the target PDBs themselves). Templates are
matched to target chains automatically by sequence; see
[Templates](fold.md#templates). In a pulldown:

- Only the **target** chain is templated by default. `--binder_templates true`
  also matches binder chains against the same files, e.g. when a complex
  structure includes a known binder.
- Boltz holds templated chains close to their template by default
  (`--boltz_template_force`, 1.0 Å), because the chain structures are known and
  the binding pose is what is being predicted. Pass `--boltz_template_force false`
  to let Boltz depart from them.
- Matching cutoffs are the same as in `--method fold`
  ([Match cutoffs](fold.md#match-cutoffs)): identity and coverage of at least
  0.3, at least 40 aligned residues, and the best 4 templates per chain.
- Every engine except ESMFold2 uses templates, most of them up to 4 per chain
  ([how each engine uses them](fold.md#how-each-engine-uses-them)).

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
| `--methods` | `boltz` | Comma-separated: `af2`, `af2_mono`, `boltz`, `rf3`, `protenix`, `af3`, `openfold3`, `esmfold2`, `esmfold2_fast` (see [Choosing engines](fold.md#choosing-engines); ESMFold2 weights and caveats are covered in [ESMFold2](fold.md#esmfold2)) |
| `--af3_model_dir` | `models/alphafold3` | AlphaFold3 weights directory (`af3` only; see [AlphaFold3 weights](fold.md#alphafold3-weights)) |
| `--af3_paired_msa` / `--af3_templates` / `--af3_run_data_pipeline` | `true` / `inline` / `false` | How AF3 is fed and run. `-profile af3_germinal_parity` sets the combination Germinal uses (pairing off, AF3's own template search, data pipeline on) and needs `--af3_db_dir`; see [Germinal: external folding validation](germinal.md#external-folding-validation) |
| `--msa_method` | `jackhmmer_af2` | `jackhmmer_af2` or `mmseqs2_colabfold` |
| `--create_target_msa` | `false` | Build MSA for each target |
| `--create_binder_msa` | `false` | Build MSA for each binder (usually leave off for de novo binders). AF2/`af2_mono` always fold the binder chain query-only regardless of this flag (see `bin/fold/assemble_af2_multimer_msas.py`); if AF2/`af2_mono` are the only selected `--methods`, the pipeline warns and skips building the binder MSA entirely. |
| `--n_predictions` | `5` | Samples per complex per method. Matches the native default of the diffusion engines, and with `--af2_keep_models all` it is what AF2's single run already produces |
| `--skip_engens` | `true` | EnGens is off by default (would emit N × M reports) |
| `--consensus_metric` | `ipsae` | `ipsae` or `iptm`; which metric's per-tool z-scores are averaged into `consensus_z` |
| `--z_stat` | `max` | `max` or `mean`; the per-complex statistic over samples that gets standardised |
| `--z_scope` | `target` | `target` standardises within `(target, tool)`; `global` pools all targets into one distribution |
| `--min_pool` | `10` | Pools holding fewer complexes than this are flagged `z_pool_small` in the summary and warned about on stderr |

Every pair is a 2-chain complex, so `af2` and `af2_mono` under
`--msa_method jackhmmer_af2` both need a `uniprot/`-bearing `--af2_db_path`
(default `alphafold_20211129`; the [Fold](fold.md) workflow instead defaults to
`alphafold_20240229`, which is monomer-only). Target MSAs for AF2 are built once
(jackhmmer, or a ColabFold/mmseqs2 a3m converted into AF2's per-chain format)
and combined with the query-only binder chain into one AF2 multimer input per
pair.

Every method-specific flag documented under `--method fold --help` (`--af2_*`,
`--boltz_*`, `--rf3_*`, `--protenix_*`, `--af3_*`, `--openfold3_*`,
`--esmfold2_*`) applies here too — the same engines, run per target×binder pair instead of per input FASTA.
One worth knowing about: `--af2_keep_models` controls which of AF2's 5 models
per run are kept. It defaults to `all` here, so one AF2 run yields the 5
structures `--n_predictions` asks for at no extra cost (`--method fold` defaults
to `best`). For `af2_mono`, `best` may still be preferable — without an initial
guess its monomer models often fail to dock, and ranking is what separates the
good pose.

## Output

Default layout under `--outdir` (`params.json` and `logs/` sit at the outdir
root, alongside `fold_pulldown/`):

```
results/
├── params.json
├── logs/
└── fold_pulldown/
    ├── pairs.tsv                    # id, target, binder — one row per co-folded pair
    ├── msa/<msa_method>/             # shared target + binder MSAs
    ├── msa/paired/                   # per-chain MSA rendered into each engine's native format
    ├── <tool>/<target>_and_<binder>/ # per-tool, per-pair engine outputs (raw + confidence JSON)
    ├── <tool>/<tool>_fold_scores.tsv # per-tool score table
    ├── predictions/                  # flat gather: <tool>_<target>_and_<binder>_*.cif
    ├── fold_pulldown_scores.tsv       # master score table: one row per predicted structure
    ├── fold_pulldown_summary.tsv      # one row per (target, binder, tool)
    └── fold_pulldown_report.html      # Quarto overview
```

Each co-folded pair gets a complex id `<target>_and_<binder>` (from the
sanitised FASTA header ids), used for its `<tool>/<target>_and_<binder>/`
directory and as the stem of its files under `predictions/`. `pairs.tsv` maps
each id back to its `target` and `binder`, and is what a downstream join
against `fold_pulldown_scores.tsv` should use.

`msa/paired/` holds those per-chain renders (RF3 `TaxID=` a3m, Protenix
mnemonic-headers a3m, Boltz `key,sequence` CSV). Each chain is written on its
own; whether an engine then pairs them is covered in
[Taxonomic pairing](#taxonomic-pairing). AF3, OpenFold3 and ESMFold2 re-render
their pairing input inside the predict task, so those a3ms are not published
here.

## Interpreting scores

Different predictors have different absolute score scales. The summary table
provides **within-tool z-scores** (`iptm_z`, `ipsae_z`) so complexes can be
compared on a common footing inside each model, and a **`consensus_z`** (the mean
over tools of one metric's per-tool z for that complex) for a cross-model ranking.
Both z-columns are always written, whichever metric `consensus_z` is built from.
The two aggregations are separate steps: `--z_stat` picks the statistic over
*samples* within one tool (`max` by default), and `consensus_z` is always the
mean of the resulting z-scores across *tools*.

Three defaults decide that ranking, and each can be reverted:

- **`--consensus_metric ipsae`.** ipSAE ([Dunbrack
  2025](https://doi.org/10.1101/2025.02.10.637595)) normalises by interface size and
  isolates the interface from whole-complex confidence, which `iptm` and `ptm` mix
  together. Pass `--consensus_metric iptm` to rank on ipTM instead.
- **`--z_stat max`.** One diffusion sample per complex is noisy, so the statistic
  worth standardising is the best of the samples rather than their average. This is
  also what makes `--n_predictions > 1` valuable. Pass `--z_stat mean` for the
  average - but be aware `mean` can be heavily impacted by outliers with poor scores.
- **`--z_scope target`.** The default is to pool only within-target scores, since raw co-folding scores are typically not comparable across targets, and the distribution of between-target scores often varies more than within-target scores, so pooling several targets into one distribution makes a complex rank partly on which target it was paired with. If you'd like to pool all targets together anyway, use `--z_scope global`.

Each summary row records `z_basis` (metric, statistic and scope, e.g.
`ipsae_max/target`), `n_pool`, the number of complexes that actually contributed a
value to that z-score's `--consensus_metric` (not just the pool's row count —
`iptm_n_pool` / `ipsae_n_pool` give the count for each metric individually), and
`z_pool_small`, `True` where that pool was smaller than `--min_pool`.
Watch those two: with *k* complexes in a pool the largest possible absolute z is
(*k*−1)/√*k*, so a two-complex pool can only ever report ±0.707 and the z-score
carries the ordering and nothing else.

`consensus_z` falls back from `--consensus_metric` to the other metric only when
an entire `(target, tool)` pool lacks the primary metric (e.g. a tool that never
emits `ipsae`) — never per-complex, so one complex's z is never averaged from a
different metric than its pool-mates. `consensus_z_metric` records which metric
actually contributed each row (comma-joined if tools within the same pair
contributed via different metrics), and `n_tools` is the number of tools that
contributed a non-missing value to that complex's `consensus_z`. `iptm_sd` /
`ipsae_sd` are blank for a complex with only one replicate — a single point has no
spread to report.

For custom statistics (Mann–Whitney, mixed models, score calibration), use
`fold_pulldown_scores.tsv` / `fold_pulldown_summary.tsv` directly — the HTML
report stays deliberately simple.

## Related

- [Fold](fold.md) — fold individual FASTA complexes (no target × binder fan-out)
