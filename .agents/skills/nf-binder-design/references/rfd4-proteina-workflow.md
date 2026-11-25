# RFD4-Proteina Workflow - Full Reference

## Overview

`--method rfd4_proteina` designs binders with [RFD4-Proteina](https://github.com/RosettaCommons/RFD4-Proteina), the flow-matching model from the Baker Lab and NVIDIA GenAIR.

**Pipeline**: `rfd4-generate` → (optional) structural filters → (optional) refold scoring → (optional) FoldSeek → `combined_scores.tsv`.

The difference from `rfd3` that matters for everything downstream: RFD4-Proteina co-designs sequence and structure in one pass, so there is **no inverse-folding step and no cofold of its own**. Each design comes out with both a backbone and a binder sequence, and nothing has scored the pair. Validation is entirely the refold's job.

Inline metrics (`ppi_refolding`, `hotspot_conditioning`) run by default via the OpenFold3 weights baked into the container; `--rfd4_metrics_tags false` switches them off (`~metrics`). They are ignored with `--rfd4_inference_experiment quickstart`. Their per-design values are merged into `combined_scores.tsv` prefixed `rfd4_`, with dots replaced by underscores (`rfd4_oracle_ppi_overall_iptm`, `rfd4_hotspot_obedience_hotspot_obedience`, `rfd4_ppi_pass`); the join is done in `bin/rfd4_proteina/collect_rfd4_designs.py`, inside the task, because rfd4-generate's example ids repeat across batches. Nothing filters on `rfd4_ppi_pass` -- `--rfd4_filters` only sees structural metrics.

## Inputs

Two mutually exclusive modes.

**Condition specification (`--rfd4_config`)** - a YAML or JSON file in RFD4-Proteina's [condition format](https://rfd4-proteina.org/latest/reference/yaml-json.html). This is the full conditioning vocabulary: hotspots, burial, hydrogen bonds, coordinate groups, pi-stacking, CATH labels, sequence probabilities.

```yaml
pdl1_binder:
  input: PDL1.cif                                   # relative to this file
  contig: "50-120, /0, A17-131"                     # binder chain, then target
  conditions:
    C_CRDGRP: {1: [{select: "sel('A/*/*/[N,CA,C,O]')"}]}   # target backbone, free pose
    C_SEQ: {true: [{select: "sel('A')"}]}                  # target sequence fixed
    C_HOT: {true: [{select: "sel('A/*/56')"}, {select: "sel('A/*/115')"}]}
```

Structure paths inside the spec resolve relative to the spec file's own directory. Do not pass `--input_pdb` as well.

**Parameter mode (`--input_pdb` + `--contigs`)** - generates a spec automatically with standard binder defaults (target backbone as a free-pose coordinate group, fixed target sequence, no cysteines in the binder, hotspots as given). The generated spec is written to `<outdir>/rfd4_proteina/input/config/<design_name>.json`, so it can be copied and hand-edited as the starting point for a condition-spec run.

The designed segment is named `BINDER` in the generated contig, so the cysteine exclusion can select it:

```yaml
  contig: "BINDER:50-120, /0, A17-131"
  conditions:
    C_SEQ_PROB_CYS: {0.0: [{select: "sel('*/BINDER')"}]}
```

`--contigs` may carry the name itself, with or without a space after the colon (`"BINDER: 50-120, /0, A17-131"`); the name given is kept and the cysteine exclusion selects it. A contig with more than one designed segment is left unnamed and gets no cysteine exclusion, since the name would be ambiguous. To allow cysteines, copy the generated spec, drop `C_SEQ_PROB_CYS`, and rerun with `--rfd4_config`.

**Residue numbering gotcha.** Residue numbers in `contig` and `conditions` are mmCIF `label_seq_id`, **not** author numbering. PyMOL shows `auth_seq_id` unless you `set cif_use_auth, off`. This is the most likely reason a hotspot lands on the wrong residue.

### Contig syntax

Comma-separated, and not the same syntax as `rfd` or `rfd3`:

| Token | Meaning |
|-------|---------|
| `/0` | Chainbreak: ends the current chain, starts a new one |
| `A17-131` | Preserved segment: kept from the input with its chain and residue numbering |
| `50-120` | Expandable segment: variable length, sampled by the model |
| `E1:2-20` | Named expandable segment, so `conditions:` can target it |

Tokens between chainbreaks make one output chain. A chain holding a preserved segment inherits that segment's chain ID; a chain with no preserved segment is the binder and takes the first unused letter. So `"50-120, /0, A17-131"` gives target `A` and binder `B`. Naming a segment does not change this: `BINDER:` is a segment name, not a chain ID, and the binder chain is still `B`. Override with `--rfd4_target_chain` / `--rfd4_binder_chain`. Every contig needs at least one preserved segment.

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--rfd4_config` | — | Condition spec (YAML or JSON). Mutually exclusive with `--input_pdb` |
| `--input_pdb` | — | Target structure, parameter mode |
| `--contigs` | — | Contig string, parameter mode |
| `--hotspot_res` | — | Hotspots for parameter mode, eg `A56,A115,A123` |
| `--design_name` | `rfd4` | Prefix for generated output files |
| `--rfd4_n_designs` | 1 | Total designs to generate |
| `--rfd4_designs_per_task` | 1 | Designs per Nextflow task; raising it amortises model loading |
| `--rfd4_batch_size` | 1 | Designs concurrent on the GPU within a task |
| `--rfd4_num_workers` | 0 | Dataloader workers. **Keep at 0** - see Troubleshooting |
| `--rfd4_seed` | 0 | Base seed; task seeds are base + task index |
| `--rfd4_nsteps` | 60 (model default) | Integration steps in the sampling schedule |
| `--rfd4_recycle_strategy` | `medium/front_loaded` | Trunk recycle schedule, eg `thinking_medium`, `high/midpeak` |
| `--rfd4_metrics_tags` | `ppi_refolding,hotspot_conditioning` | Inline metrics, comma-separated; `false` disables; ignored with `quickstart` |
| `--rfd4_inference_experiment` | `inference_on_user_inputs` | Preset; `quickstart` for fast test runs |
| `--rfd4_compile_lengths` | `true` | `false` runs the trunk eagerly, skipping `torch.compile` warm-up |
| `--rfd4_partial_trajectory_fraction` | — | Partial diffusion: run only this final fraction of `nsteps` |
| `--rfd4_inverse_folding` | `false` | Hold the input backbone fixed, resample sequence only |
| `--rfd4_filters` | — | Semicolon-separated filters on generated designs, eg `rg<25` |
| `--rfd4_extra_args` | — | Extra Hydra `key=value` overrides for `rfd4-generate` |
| `--rfd4_target_chain` / `--rfd4_binder_chain` | derived | Chain IDs in the output complex |

`rfd4-generate` is Hydra-driven, so `--rfd4_extra_args` takes `key=value` pairs: `--rfd4_extra_args "generation.args.guidance_w=1.2 precision=bf16-mixed"`.

## FoldSeek: two separate jobs

FoldSeek is wired in for two unrelated purposes, with separate flags, subworkflows and
output directories. They share only the binder-chain extraction step
(`FOLDSEEK_PREPARE_QUERIES`), which the clustering subworkflow includes under an alias so
both can run in the same execution.

| | Database search | Design clustering |
|--|--|--|
| Flag | `--do_foldseek_search` | `--do_foldseek_cluster` |
| Question | Has this fold been seen before? | How many distinct binders did this run make? |
| Compares against | CATH50/PDB/AFDB | The run's own designs |
| Subworkflow | `FOLDSEEK_SEARCH` | `FOLDSEEK_CLUSTER_DESIGNS` |
| Output | `foldseek/<database>/` | `foldseek_cluster/` |
| Scores columns | none | `foldseek_cluster_*` on `combined_scores.tsv` |

Clustering runs `foldseek easy-cluster` over the binder chains with the thresholds from
upstream's `scripts/launch_generate.sh` (`--alignment-type 1 --cov-mode 0 --min-seq-id 0.0
--tmscore-threshold 0.5`), then `bin/foldseek/cluster_table.py` turns the
representative/member table into one row per design keyed on `id`, which
`COMBINE_RFD4_SCORES` merges with `bin/merge_scores.py`. The conversion is a separate
process (`FOLDSEEK_CLUSTER_TABLE`) because the foldseek image has no python3.

`--foldseek_cluster_multimer` switches to `easy-multimercluster` over whole complexes,
which skips chain extraction. `easy-interfacecluster` is absent from foldseek release 10,
so it is not yet an option.

Columns added: `foldseek_cluster_rep`, `foldseek_cluster_id` (1 = largest cluster),
`foldseek_cluster_size`, `foldseek_cluster_is_rep`.

## Refold Scoring

The shared refold step, with the same `--refold_with` / `--refold_filter` / `--refold_sort_by` / `--refold_max` as every other workflow. See the refold scoring section of `SKILL.md`.

```bash
  --refold_with boltz,chai1,openfold3 \
  --refold_max 50 \
  --use_msa_server
```

**Selection is the one place rfd4_proteina differs.** There are no cofold metrics to gate on, so `--refold_filter` and `--refold_sort_by` can only read the structural metrics produced by `--rfd4_filters` (`rg`, plus anything in `bin/filters.d/`). With `--rfd4_filters` unset there is nothing to rank by, and `--refold_max` alone takes the first N designs in whatever order they finish. Write the sort with an `=`: `--refold_sort_by='rg'`.

Selection verdicts land in `rfd4_proteina/refold/refold_selection.tsv`.

See `fold-stack.md` for how `--n_predictions`, the per-engine batch sizes and the seeds interact.

RF3 used as a refold engine takes `--refold_rf3_*` names (`--refold_rf3_ckpt_path`, `--refold_rf3_num_steps`, `--refold_rf3_n_recycles`, `--refold_rf3_batch_size`, `--refold_rf3_early_stopping_plddt_threshold`, `--refold_rf3_seed`), kept distinct from the generation-time `--rf3_*` set.

## Outputs

| Path | Contents |
|------|----------|
| `rfd4_proteina/combined_scores.tsv` | Binder sequences and every computed metric |
| `rfd4_proteina/designs/` | Generated complexes, one mmCIF per design |
| `rfd4_proteina/designs/rfd4_designs_batch*.tsv` | Design id, source example and seed per batch |
| `rfd4_proteina/batches/batch<n>/` | Raw `rfd4-generate` run dir per batch: inference metrics, each sample's generated structure, inline refold and confidences |
| `rfd4_proteina/input/config/<design_name>.json` | Generated condition spec (parameter mode) |
| `rfd4_proteina/input/prepared/` | Input structure converted to mmCIF |
| `rfd4_proteina/binders.fasta` | Generated binder sequences |
| `rfd4_proteina/filter/` | Filter scores and accept/reject lists (with `--rfd4_filters`) |
| `rfd4_proteina/refold/` | Refold structures, per-engine scores, interface tables (with `--refold_with`) |
| `foldseek_cluster/` | `clusters.tsv`, `cluster_table.tsv`, representative sequences (with `--do_foldseek_cluster`) |

Structure filenames are `<design_name>_<uid>_b<task>_<example>_<replicate>.cif`.

## arm64 / GB10

The `rfd4-proteina` container is built for `linux/arm64` and is the only generator that runs on NVIDIA DGX Spark GB10. Other pipeline images are amd64, so on an arm64 host add `conf/platforms/arm64.config`, which remaps tasks to arm64 builds:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  -c conf/platforms/arm64.config \
  -c nextflow.docker.config \
  --method rfd4_proteina --rfd4_config pdl1_binder.yaml
```

Covered on arm64: nf-binder-design-utils, Boltz-2, ESMFold2, RF3, Protenix, BindCraft, FoldSeek. Limits:

- FoldSeek runs CPU-only; upstream GPU binaries are x86_64.
- Refold engines available: `boltz`, `esmfold2`, `esmfold2_fast`, `rf3`, `protenix`. AF2, AF3, Chai-1 and OpenFold3 have no arm64 images.
- RF3 on GB10 needs the `-ngc2606` image. `rc-foundry:0.2.0_arm64-weights` is built on CUDA 12.8, whose nvrtc cannot target `sm_121`, and PyTorch's jiterator compiles `erfinv` for the host arch, so RF3 dies with `nvrtc: error: invalid value for --gpu-architecture (-arch)`. The config already points at the CUDA 13.3 rebuild.
- Protenix on GB10 needs the dedicated arm64 image: its fused layer-norm extension is built for `sm_70`-`sm_100`, so on `sm_121` inference fails with `no kernel image is available for execution on the device`. The arm64 image adds `sm_120`/`sm_121`.

On amd64 clusters none of this applies.

## Troubleshooting

- **Hangs before sampling, GPU idle.** `rfd4-generate` deadlocks on its dataloader when multiprocessing workers are on. Keep `--rfd4_num_workers 0`.
- **First task very slow.** `torch.compile` warms the trunk once per execution. For quick tests, `--rfd4_compile_lengths false --rfd4_inference_experiment quickstart`.
- **Out of memory.** Drop `--rfd4_batch_size` to 1. Memory scales with total sequence length across target and binder.
- **Hotspot appears to be on the wrong residue.** Contig and condition numbering is `label_seq_id`, not author numbering.
