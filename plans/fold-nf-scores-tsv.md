# fold.nf master scores TSV

## Goal

Produce a single master TSV (`<outdir>/fold/fold_scores.tsv`) with **one row per
generated structure** across all four engines (AF2, Boltz-2, RF3, Protenix),
plus per-tool combined TSVs. Compute ipSAE (via `bin/ipsae.py`) for AF2. Column
names are **normalized** across engines for equivalent scores; plddt is
**rescaled to 0–1** everywhere; asymmetric per-chain-pair scores are dropped.

## Structure (matches existing Boltz pattern)

```
per-structure parse (one process invocation per model/sample, emits a TSV row to stdout)
   -> collectFile -> per-tool combined TSV  <outdir>/fold/<tool>/<tool>_fold_scores.tsv
   -> FOLD_MERGE_SCORES (bin/merge_fold_scores.py) -> master <outdir>/fold/fold_scores.tsv
```

Boltz already does the first two (`FOLD_PARSE_BOLTZ_CONFIDENCE` +
`collectFile` -> `boltz_fold_scores.tsv`). Replicate for AF2/RF3/Protenix, then
add the master merge over all four per-tool TSVs.

**Normalization lives in the master merge** (`bin/merge_fold_scores.py`), one
place with a per-tool column map. Per-tool parsers stay simple (flatten native
JSON + provenance columns), like `bin/parse_boltz_confidence.py`.

## Canonical (master) schema, in order

Provenance: `tool`, `id`, `model`, `original_file`, `predictions_file`
Scores: `ranking_score`, `ptm`, `iptm`, `plddt` (0–1), `pae`, `pde`, `has_clash`
ipSAE block: `ipsae`, `ipsae_d0chn`, `ipsae_d0dom`, `pdockq`, `pdockq2`, `lis`

Blank/NaN where an engine does not report a metric.

- `predictions_file` = the renamed file in `<outdir>/fold/predictions/`.
- `original_file` = the engine-native structure filename (basename).
- `model` = the per-structure index (AF2 model N / Boltz model_N / Protenix
  sample_N / RF3 seed-S_sample-N), batch-qualified when namespaced.

## Per-engine source → canonical map (verified against example results)

Example results live under `examples/fold-multimer/results/fold/`.

### AF2 (`fold/af2/predictions/<id>/`)
Native, per kept model N (1..5 for keep=all; ranked_0's model for keep=best):
- `result_model_N_multimer_v3_pred_0.pkl` (needs numpy+pickle): `ptm`, `iptm`,
  `ranking_confidence`, `plddt` (per-res ndarray, 0–100 → mean/100).
- ipSAE: `bin/ipsae.py pae_model_N_..json relaxed_model_N_..pdb 10 10` → the
  `Type == 'min'` row of `*_ipsae.tsv`: `ipSAE`→ipsae, `ipSAE_d0chn`,
  `ipSAE_d0dom`, `pDockQ`→pdockq, `pDockQ2`→pdockq2, `LIS`→lis.
- Map: ranking_score←ranking_confidence, ptm←ptm, iptm←iptm,
  plddt←mean(plddt)/100. No pae/pde/has_clash.
- keep=best: only `ranked_0.{pdb,cif}` published; the parse process still
  receives the full `out/<id>/**` channel (all pkls/pae present in work dir).
  Map ranked_0 → `ranking_debug.json` `order[0]` model → its pkl/pae. model
  label = that model's N (or "ranked_0").

### Boltz (`fold/boltz/.../confidence_complex_model_N.json`)
The BOLTZ module already writes ipSAE into this JSON. Existing
`boltz_fold_scores.tsv` (native) is the tool-level table — LEAVE IT AS-IS
(`bin/parse_boltz_confidence.py` is shared with boltz_pulldown.nf; only ADD
optional `--original-file`/`--predictions-file` args, never change defaults).
- Map: ranking_score←confidence_score, ptm←ptm, iptm←iptm,
  plddt←complex_plddt (already 0–1), pde←complex_pde, ipsae←ipsae_min. No pae,
  no has_clash. DROP `chains_ptm`, `pair_chains_iptm` (asymmetric).

### Protenix (`fold/protenix/.../<id>_summary_confidence_sample_N.json`)
- Map: ranking_score←ranking_score, ptm←ptm, iptm←iptm, plddt←plddt/100 (0–100!),
  pde←gpde, has_clash←has_clash. No pae, no ipSAE.
- DROP `chain_*`, `chain_pair_*` arrays (asymmetric).

### RF3 (`fold/rf3/.../<id>_..._summary_confidences.json`)
- Map: ranking_score←ranking_score, ptm←ptm, iptm←iptm, plddt←overall_plddt
  (0–1), pae←overall_pae, pde←overall_pde, has_clash←has_clash. No ipSAE.
- DROP `chain_ptm`, `chain_pair_*` arrays (asymmetric).

## `predictions_file` — single source of truth

The renamed flat-gather name is currently computed in each module's SECOND
`publishDir { saveAs }`. Add `lib/FoldNaming.groovy` with a static
`predictionName(tool, meta, nativeBasename)` that reproduces those rules, use it
in each module's second saveAs AND pass its result to the parse process, so the
`predictions_file` column can never drift from the published name.

Current rules (namespaced = `meta.fold_namespaced`; msa bit = `_msa<tag>` when
`meta.msa_depth_tag`):
- af2:   `af2_<id>_run<af2_run><msa>_<bn>`  (bn = relaxed_model_N..cif / ranked_0.cif / unrelaxed_model_N..cif)
- boltz: namespaced? `boltz_batch<fold_batch><msa>_<bn>` : `boltz<msa>_<bn>`  (bn = <id>_model_N.cif)
- rf3:   namespaced? `rf3_batch<fold_batch><msa>_<bn>` : `rf3<msa>_<bn>`  (bn = <id>_seed-S_sample-N_model.cif)
- protenix: namespaced? `protenix_batch<fold_batch><msa>_<bn>` : `protenix<msa>_<bn>`  (bn = <id>_sample_N.cif)

Uniqueness is already guaranteed by tool prefix + model/sample index + batch
namespacing; the helper just centralizes it.

## Nextflow wiring

Container for parse/merge/ipsae CPU tasks: reuse
`ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6`
(has pandas+numpy; used by FOLD_PARSE_BOLTZ_CONFIDENCE). ipsae.py needs only
numpy.

- New modules under `modules/fold/<tool>/`:
  `fold_parse_af2_confidence.nf` (+ `fold_af2_ipsae.nf`),
  `fold_parse_rf3_confidence.nf`, `fold_parse_protenix_confidence.nf`.
- New `modules/fold/common/fold_merge_scores.nf` = `FOLD_MERGE_SCORES`
  (executor local, cpus 1).
- In `alphafold2.nf` / `rosettafold3_fold.nf` / `protenix_fold.nf`
  subworkflows: fan the predictions/confidence channel to one parse call per
  structure, `collectFile` -> `<outdir>/fold/<tool>/<tool>_fold_scores.tsv`,
  emit `tsv`.
- In `fold.nf`: collect the four per-tool `tsv` emits, feed FOLD_MERGE_SCORES ->
  `fold/fold_scores.tsv`. Skip tools not in `--methods`.
- Add cpus=1 local selectors for the new parse/merge processes to
  `nextflow.config` + `conf/platforms/m3.config` (mirror the existing
  `BOLTZ_FOLD:FOLD_PARSE_BOLTZ_CONFIDENCE` block; use fully-qualified
  `withName: '<SUBWF>:<PROC>'` names).

## Validation (no GPU needed)

1. Offline-test each parser against the committed example JSONs under
   `examples/fold-multimer/results/fold/` before wiring.
2. Scoped `-resume` run from `examples/fold-multimer/` (predictions are cached;
   only the new parse/ipsae/merge tasks run — CPU only). Requires
   `export NXF_VER=24.10.0`. Confirm `fold/fold_scores.tsv` has the canonical
   columns and one row per structure in `fold/predictions/` (16 for the current
   example: 5 af2 + 1 boltz + 5 protenix + 5 rf3, pre-`n_predictions=5`; a fresh
   run would give more).

## Docs

Document `fold/fold_scores.tsv` + the per-tool tables and the normalized schema
in `docs/docs/workflows/fold.md` (Outputs section).
