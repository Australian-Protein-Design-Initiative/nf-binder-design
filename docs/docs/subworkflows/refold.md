# Refold Scoring (`--refold_with`)

Refolding is a held-out second opinion on designs the pipeline has already cofolded once with its own predictor.

Every design workflow folds its designs before this point - AlphaFold2 initial guess in `--method rfd` and `--method rfd_partial`, RosettaFold3 in `--method rfd3`, and BindCraft's and Germinal's own AF2 passes. Refolding does not replace that step. It takes a chosen subset of the designs, folds each binder sequence against its target again with one or more *independent* predictors, compares the result with the design, and scores the rebuilt interface. The first cofold's scores are what the refold selects on.

It is available in all six design workflows: `rfd`, `rfd_partial`, `rfd3`, `rfd4_proteina`, `bindcraft` and `germinal`. Every workflow takes the same `--refold_*` parameters and writes the same engine-prefixed columns into its own `combined_scores.tsv`.

Refolding is off by default in every workflow, because it is extra GPU hours per design on top of the pipeline's own cofold.

Implemented by `subworkflows/local/refold.nf` and `subworkflows/local/refold_select.nf`; the shared defaults are in `conf/refold.config` and `conf/refold_fold_stack.config`.

## Choosing the engines

```bash
--refold_with boltz,chai1,openfold3
```

Any engine `--method fold` supports can go in the list: `af2`, `af2_mono`, `boltz`, `rf3`, `protenix`, `af3`, `openfold3`, `chai1`, `esmfold2`, `esmfold2_fast`.

## Choosing which designs are refolded

Three parameters, and they compose. All are optional; with none of them set, every design is refolded.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--refold_filter` | none | Threshold gate on the first cofold's scores, semicolon-separated, in the same syntax as `--rfd_filters`. For example `'ipsae>=0.5;pae_interaction<=10'` |
| `--refold_sort_by` | none | Rank the designs that pass the gate by one metric, `-` prefix for descending |
| `--refold_max` | none | Cap on how many designs are refolded |

- **`--refold_filter` streams.** The gate is decided per design, so a design is refolded as soon as its own scores exist, while the rest of the set is still being scored.
- **`--refold_max` does not stream.** A cap makes the refold wait for the whole design set, with or without a ranking: "the best N" is not knowable until the last design has been scored, and "the first N" needs a count kept across every design. That trade-off is why the gate and the cap are separate parameters.
- **`--refold_max` with `--refold_sort_by` takes the best N; without it, the first N to be scored**, in whatever order they finish. The pipeline warns when a cap is set without a ranking.
- **`--refold_sort_by` on its own selects nothing.** With no cap every design that passes the gate is refolded anyway, so the ranking only changes an order nothing downstream depends on. That case keeps streaming rather than stalling for a ranking that decides nothing, and the pipeline warns.

A filter naming a metric the score table does not have stops the run and lists the columns that are present, rather than silently rejecting every design.

!!! warning "Write `--refold_sort_by` with an `=`"
    Use `--refold_sort_by='-ipsae'`, not `--refold_sort_by '-ipsae'`. Without the `=`, Nextflow reads the leading `-` as the start of one of its own options.

    Metrics where low is good (`pae_interaction`, `pair_pae_min`, an RMSD) need no prefix.

### The selection report

Every design that reached the step is accounted for in `<refold_publish_dir>/refold_selection.tsv`, so a thin refold set can be traced to the threshold or the ranking that caused it.

| Column | Contents |
|--------|----------|
| `design_id` | The design |
| `refold_selection` | `selected` or `rejected` |
| `refold_rank` | Position in the ranking, when `--refold_max` is set (ranked by `--refold_sort_by`, or arrival order without it) |
| `refold_sort_metric` | The metric ranked on |
| `refold_sort_value` | That design's value for it |

## Which chain is the binder

`--refold_binder_chain` and `--refold_target_chain` only apply to `--method bindcraft` and `--method germinal`, which have no derivation of their own (defaults: binder `B`, target `A`). `rfd3` and `rfd4_proteina` derive the chains from the contigs, and `rfd` and `rfd_partial` know them from how AF2 initial guess writes its complexes.

## Output columns

Every engine produces the same columns in `combined_scores.tsv`, prefixed with the engine name:

| Column | What it says |
|--------|--------------|
| `<tool>_iptm_max`, `<tool>_ipsae_max`, `<tool>_plddt_max`, ... | How confident the refold is, per `--refold_metrics` |
| `<tool>_target_aligned_binder_rmsd_all` | Superimpose on the target, score the binder: did it come back where it was designed to sit |
| `<tool>_complex_rmsd_all` | The whole complex against the design |
| `<tool>_monomer_vs_design_rmsd_all` | The binder folded alone, against its designed conformation |
| `<tool>_monomer_vs_complex_rmsd_all` | The binder folded alone, against the binder in the refolded complex |
| `<tool>_interface_dG`, `<tool>_interface_dSASA`, ... | BindCraft interface metrics on the refolded complex |
| `<tool>_selected_batch`, `<tool>_selected_model` | Which sample the RMSD and interface columns describe |
| `consensus_z`, `n_tools` | Per design, not per engine: a cross-engine z-score, so a design several independent predictors agree on ranks above one only its designer liked |

`refold_sample_scores.tsv`, under the same directory, has one row per refolded sample (keyed by `id`, `tool`, `batch`, `model`) with nothing reduced, so any statistic over the samples can be computed from it afterwards.

## Reporting parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--refold_with` | none | Comma-separated engines |
| `--refold_metrics` | iptm/ipsae/plddt/ptm/pae | Summary columns carried into `combined_scores.tsv` |
| `--refold_structure_stat` | `selected` | How several samples per design reduce to one RMSD and one interface score: `selected`, `min`, `first`, `mean` or `max` |
| `--refold_binder_monomer` | `true` | Also fold each binder alone, for the monomer-vs-complex and monomer-vs-design comparisons |
| `--refold_monomer_single_sequence` | `true` | Fold that monomer from sequence alone, with no MSA search |
| `--refold_interface_scoring` | `true` | BindCraft interface metrics on the refolded complex |
| `--output_rmsd_aligned` | `false` | Keep the superimposed structures from the RMSD comparisons |
| `--n_predictions` | `5` | Structures per design per engine, run as that many seeds of one sample each |

With `--refold_structure_stat selected`, the RMSD and interface columns are all read from the single best-ranked sample (by `--consensus_metric`, ipSAE) rather than aggregated per metric. The structural metrics all describe one set of coordinates, so taking the lowest RMSD and the lowest `interface_dG` independently would give a row that no predicted structure matches; and the interface metrics have no shared direction to take a best in anyway, since `interface_dG` is better lower while `interface_dSASA` and `interface_sc` are better higher. `min`, `mean`, `max` and `first` apply one statistic to every structural column instead.

## How the refold folds

The refold drives the same fold stack as `--method fold` and `--method fold_pulldown`, so MSA and template choices come from the fold stack's own parameters: `--msa_method`, `--create_target_msa`, `--use_msa_server`, `--templates`, `--colabfold_envdb`, `--uniref30`. Every engine flag from `--method fold --help` (`--af2_*`, `--boltz_*`, `--protenix_*`, `--af3_*`, `--openfold3_*`, `--chai1_*`, `--esmfold2_*`) applies here too.

RosettaFold3 is the one exception. Its refold flags are `--refold_rf3_ckpt_path`, `--refold_rf3_num_steps`, `--refold_rf3_n_recycles`, `--refold_rf3_batch_size`, `--refold_rf3_early_stopping_plddt_threshold` and `--refold_rf3_seed`. The plain `--rf3_*` names belong to the design-time RF3 step in `--method rfd3`, where they mean different things, and a workflow's parameters are one flat namespace. Where the prefixed names are declared, the fold stack reads only those.

## Renamed and removed parameters

The whole `--full_refold_*` family has been renamed without exception: `--full_refold_with` is now `--refold_with`, `--full_refold_max` is `--refold_max`, `--full_refold_rf3_batch_size` is `--refold_rf3_batch_size`, and so on. `--full_refold_filter_sort` is now `--refold_sort_by`.

Deprecated, still accepted with a warning:

| Old name | Use instead |
|----------|-------------|
| `--refold_af2ig_filters` | `--refold_filter` |
| `--refold_oracles` | `--refold_with` |
| `--refold_oracles_max` | `--refold_max` |

Removed. Refolding takes its MSA and template choices from the fold stack now, so these have no equivalent of their own. A run that sets one stops with the replacement named:

| Removed | Use instead |
|---------|-------------|
| `--refold_create_target_msa` | `--create_target_msa` |
| `--refold_use_msa_server` | `--use_msa_server` |
| `--refold_target_templates` | `--templates` |
| `--refold_alignment` | `--create_target_msa` with `--msa_method mmseqs2_colabfold` |
| `--refold_target_fasta` | No replacement; the target sequence is read from the design structures |

`--boltz_refold_batch_size` is removed too. It batched several designs into one Boltz refold comparison task, and that Boltz-specific path no longer exists. The nearest equivalent is `--refold_boltz_batch_size`, which is samples per fold job rather than designs per task.

## See also

- [RFdiffusion Workflows](../workflows/rfdiffusion.md)
- [RFdiffusion3 Workflow](../workflows/rfd3.md)
- [RFD4-Proteina Workflow](../workflows/rfd4-proteina.md)
- [Fold Workflow](../workflows/fold.md) - the engine flags the refold shares
