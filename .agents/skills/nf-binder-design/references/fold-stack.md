# Fold Stack: Seeds, Samples and Batch Sizes

The same fold stack runs behind `--method fold`, `--method fold_pulldown` and the `--refold_with` step of every design workflow. How many structures it produces per input, how they are split into jobs, and which seed each job gets, are controlled the same way in all three. The defaults are **not** the same, and the parameter *names* differ between the refold and the standalone methods. Both catch people out.

## The model

Three things, in this order:

1. `--n_predictions` - how many structures per input per engine, in total.
2. `--<engine>_batch_size` - how many of those go in one job. The stack splits `n_predictions` into `ceil(n_predictions / batch_size)` jobs.
3. The base seed - job *i* runs with `seed + i`.

So `--n_predictions 6 --boltz_batch_size 3` is 2 jobs, seeds `base` and `base+1`, 3 samples each.

**Samples in one job are not independent.** Engines that batch samples share one trunk and recycling pass across the batch, so N samples from one seed differ only in the diffusion noise drawn afterwards. N separate seeds re-run the whole prediction. More seeds is more sampling diversity; a bigger batch is cheaper. The trade-off is model loading, paid once per job - roughly 45 s per Boltz load on a GB10.

## Defaults differ between the three entry points

| | `--refold_with` | `--method fold_pulldown` | `--method fold` |
|---|---|---|---|
| `--n_predictions` | **5** (`conf/refold_fold_stack.config`) | **5** (`workflows/fold_pulldown.nf`) | **unset**: each engine falls back to its own default |
| `--<engine>_batch_size` | **1**, so 5 predictions is 5 jobs of 1 sample, each with its own seed | unset: one job of N samples | unset: one job of N samples |
| Base seeds | fixed (42, or 1 for AF3), so runs are reproducible | mostly unset: each engine draws its own random seed | mostly unset |

The refold deliberately uses five model loads to get five independent seeds, because its whole job is to be an independent check. `fold` and `fold_pulldown` favour throughput: one job of N samples, which share a trunk pass. With `--n_predictions` unset in `fold`, Boltz, RF3, Protenix, AF3, OpenFold3, Chai-1 and ESMFold2 each emit 5 diffusion samples from a single job (Boltz is raised from its native 1 for parity across engines), and AF2 does one run keeping per `--af2_keep_models` (default `best`, so one structure).

Set `--n_predictions N` to pin every diffusion engine to exactly N, split by that engine's `--*_batch_size`. Leave `--*_batch_size` unset and you get one job of N.

**AF2 does not follow the model above.** It produces its 5 trained models per execution regardless of `--n_predictions`; `--af2_keep_models` decides how many executions are needed (`best` keeps the top model per run, so N runs; `all` keeps all 5, so `ceil(N/5)` runs).

## Parameter names: the `refold_` prefix rule

`modules/local/fold/common/fold_params.nf`: **where a `refold_`-prefixed copy of an engine parameter is declared, that prefixed name is the only one the refold reads.** Where no prefixed copy is declared, the refold reads the plain name.

This exists because a workflow's parameters are one flat namespace and some plain names already mean something else at design time - `rfd3`'s `--rf3_batch_size` is MPNN designs per RF3 task, not samples per fold job. So:

| Engine | Refold reads | Note |
|--------|--------------|------|
| Boltz | `--refold_boltz_batch_size`, `--refold_boltz_seed` | plain `--boltz_*` has no effect on the refold |
| RF3 | `--refold_rf3_batch_size`, `--refold_rf3_seed`, `--refold_rf3_num_steps`, `--refold_rf3_n_recycles`, `--refold_rf3_early_stopping_plddt_threshold`, `--refold_rf3_ckpt_path` | plain `--rf3_*` is `rfd3`'s design-time RF3 step |
| ESMFold2 | `--refold_esmfold2_batch_size`, `--refold_esmfold2_seeds` | |
| Chai-1 | `--refold_chai1_batch_size`, `--refold_chai1_seeds` | |
| OpenFold3 | `--refold_openfold3_batch_size`, `--refold_openfold3_seeds` | |
| AF3 | `--refold_af3_batch_size`, but **`--af3_seeds`** unprefixed | only the batch size has a prefixed copy |
| Protenix | `--protenix_batch_size`, `--protenix_seeds` - **unprefixed** | no prefixed copy is declared |
| AF2 | `--af2_random_seed`, `--af2_keep_models`, ... - **unprefixed** | |

In `--method fold` and `--method fold_pulldown` the plain names are always the right ones; `refold_`-prefixed names do nothing there. These parameters cannot simply be moved into `conf/refold.config`, which every entrypoint includes, because `FoldParams` decides by asking whether the prefixed name is declared at all - `--method fold` would then start reading the refold's values instead of its own.

## Seeds

In `fold` / `fold_pulldown`, seeds are unset by default and each engine draws its own, so pin `--af2_random_seed` / `--boltz_seed` / `--rf3_seed` / `--protenix_seeds` / `--chai1_seeds` when you want reproducibility. Do not inject a fresh random seed on every invocation if you intend to `-resume`: the changed value re-runs the tasks.

AlphaFold3, OpenFold3 and ESMFold2 require a seed, so they have fixed base defaults (`--af3_seeds 1`, `--openfold3_seeds 42`, `--esmfold2_seeds 42`). `--protenix_seeds` and `--chai1_seeds` follow the same `seed + i` per-batch convention. Without `--chai1_seeds`, each Chai-1 job makes its own seed and records it in the output filenames.

## Grouping seeds into one job

`--fold_seed_list` (default `true`) hands the whole seed list to engines that can run several seeds per invocation - `esmfold2`, `esmfold2_fast`, `openfold3`, `protenix`, and `af3` via `--af3_seeds` - so they pay one model load instead of N. Boltz and RF3 run one job per seed either way. `--fold_seed_list false` disables grouping for every engine.

## Aggregating back down

Several samples per design reduce to one row per design in `combined_scores.tsv`:

- Confidence metrics get direct aggregates: `<tool>_iptm_max`, `<tool>_iptm_mean`, `<tool>_iptm_n`, and so on, over all samples.
- Structural and interface metrics follow `--refold_structure_stat` (default `selected`): all of them are taken from the single sample ranked best by `--consensus_metric` (default `ipsae`), so they describe one physically coherent structure rather than a mix of poses. `<tool>_selected_batch` and `<tool>_selected_model` name that sample. Setting `min`, `mean`, `max` or `first` instead applies that one operator across every structural column.

`refold_sample_scores.tsv` keeps the unaggregated per-sample rows.

Cross-engine consensus (`consensus_z`, `n_tools`) is computed as in `fold_pulldown` (`bin/fold_pulldown_summarise.py`), controlled by `--consensus_metric`, `--z_stat` (default `max`), `--z_scope` (default `target`) and `--min_pool` (default 10).

## Practical settings

- **Refold, default:** leave it alone. 5 seeds x 1 sample is the independent check the refold is for.
- **Refold, short on GPU time:** `--n_predictions 3`, or raise `--refold_boltz_batch_size` to trade seed diversity for fewer model loads.
- **`fold` / `fold_pulldown`, reproducible:** set `--n_predictions` explicitly and pin every base seed; do not leave them random.
- **`fold` / `fold_pulldown`, many inputs:** leave `--*_batch_size` unset for one job of N samples per input, which is the fewest model loads, and accept that the samples share a trunk pass.
