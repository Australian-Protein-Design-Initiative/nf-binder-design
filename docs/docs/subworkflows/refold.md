# Refold Scoring (`--refold_with`)

Refolding provides independent structural validation for designs previously predicted
during generation.

Every design workflow evaluates its generated backbones before this stage: AlphaFold2
initial guess in `--method rfd` and `--method rfd_partial`, RosettaFold3 in
`--method rfd3`, and native AF2 prediction passes in BindCraft and Germinal. Refolding
does not replace those initial evaluations. Instead, it takes selected designs,
co-folds each binder sequence against the target using independent prediction engines,
compares predicted conformations to the designed pose, and re-scores the interface.
Initial prediction scores determine which candidates are selected for refolding.

Refold scoring is available across six design workflows: `rfd`, `rfd_partial`,
`rfd3`, `rfd4_proteina`, `bindcraft`, and `germinal`. Each workflow accepts the
same `--refold_*` parameters and appends engine-prefixed score columns to its
`combined_scores.tsv`.

Refolding is disabled by default because it requires additional GPU compute per
candidate design.

The subworkflow is implemented in `subworkflows/local/refold.nf` and
`subworkflows/local/refold_select.nf`, with shared defaults defined in
`conf/refold.config` and `conf/refold_fold_stack.config`.

## Prediction engines

Specify engines as a comma-separated list:

```bash
--refold_with boltz,chai1,openfold3
```

All engines supported by `--method fold` are available: `af2`, `af2_mono`, `boltz`,
`rf3`, `protenix`, `af3`, `openfold3`, `chai1`, `esmfold2`, and `esmfold2_fast`.

## Candidate selection

Three parameters govern candidate selection and can be combined. When none are set,
all designs are refolded:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--refold_filter` | none | Threshold filter applied to initial scores, semicolon-separated (e.g. `'ipsae>=0.5;pae_interaction<=10'`). Uses the same syntax as `--rfd_filters`. |
| `--refold_sort_by` | none | Metric used to rank candidates passing the filter; prefix with `-` for descending order. |
| `--refold_max` | none | Maximum number of designs to refold. |

Filtering and batching behaviour:

- Streaming evaluation with `--refold_filter`: filters evaluate per design. A candidate
  is submitted for refolding as soon as its initial scores are available, without waiting
  for the remaining set.
- Global evaluation with `--refold_max`: capping the number of refolded designs requires
  evaluating the full candidate set before dispatching refold jobs, because ranking
  requires global comparison.
- Ranking with `--refold_max` and `--refold_sort_by`: the pipeline selects the top N
  designs according to the specified metric. If `--refold_max` is specified without
  `--refold_sort_by`, designs are selected in completion order, and the pipeline emits
  a warning.
- Sorting without `--refold_max`: specifying `--refold_sort_by` without `--refold_max`
  does not filter designs. All candidates passing the filter are refolded, and the
  pipeline emits a warning.

If a filter references a column missing from the score table, the pipeline terminates
with an error listing available columns rather than excluding all designs without
explanation.

!!! warning "Specify `--refold_sort_by` with an `=`"
    Use `--refold_sort_by='-ipsae'`, not `--refold_sort_by '-ipsae'`. Without the `=`,
    Nextflow interprets the leading hyphen as a Nextflow CLI option.

    Metrics where lower values are preferred (`pae_interaction`, `pair_pae_min`, RMSD)
    require no prefix.

### Selection report

All candidate designs evaluated for refolding are logged in
`<refold_publish_dir>/refold_selection.tsv`:

| Column | Description |
|--------|-------------|
| `design_id` | Design identifier |
| `refold_selection` | `selected` or `rejected` |
| `refold_rank` | Rank position when `--refold_max` is set (ranked by `--refold_sort_by`, or arrival order without it) |
| `refold_sort_metric` | Metric used for ranking |
| `refold_sort_value` | Metric value for that design |

!!! warning "Modifying selection parameters on resume requires a new `--outdir`"
    `refold_scores.tsv`, `refold_rmsd.tsv`, and `refold_selection.tsv` are gathered
    using `collectFile(storeDir:)`. Rerunning with modified selection parameters into an
    existing output directory retains previously gathered TSV files, causing mismatched
    reports. When changing selection criteria, specify a new `--outdir`.

## Binder chain identification

`--refold_binder_chain` and `--refold_target_chain` apply specifically to
`--method bindcraft` and `--method germinal` (defaults: binder `B`, target `A`).
`rfd3` and `rfd4_proteina` infer chains from contig specifications, whereas `rfd`
and `rfd_partial` identify chains from AF2 initial guess complexes.

## Target sequence configuration

By default, the refold target sequence matches the target from design generation:
the target chain sequence is parsed from selected design structures, and the consensus
sequence across designs is used for refolding. When designs report divergent target
sequences, the pipeline logs a warning and uses the most frequent sequence.

`--refold_target_fasta` supplies an alternate target sequence:

```bash
--refold_target_fasta full_length_target.fasta
```

This option allows testing binders against full-length targets when initial designs
were generated against truncated domains, assessing whether binding geometry is
maintained in the complete protein context.

The FASTA file must contain **exactly one record**. Because each binder is folded
against the target and scores are recorded one row per design in `refold_scores.tsv`,
evaluating multiple target variants requires separate pipeline executions with
distinct `--outdir` locations.

The pipeline verifies that the target sequence from the design structures appears as
a subsequence within the provided FASTA. If not found, a warning is logged.

### Target structures and MSAs

Three fold-stack parameters configure target inputs:

| Parameter | Description |
|-----------|-------------|
| `--templates` | Template structures matched to folded chains by sequence. Providing a full-length target structure enables prediction engines to model regions omitted during initial design. Only the target chain is templated unless `--binder_templates` is set. |
| `--target_msa` | Precomputed A3M alignment for the target. The first record must match the target sequence. Mutually exclusive with `--create_target_msa`. |
| `--create_target_msa` | Generates a target MSA using `--msa_method` alongside `--use_msa_server`, `--colabfold_envdb`, or `--uniref30`. |

These parameters omit the `--refold_` prefix because fold stack operations during
design workflows pertain exclusively to refolding.

!!! warning "Full-complex RMSD columns require equal sequence lengths"
    When `--refold_target_fasta` provides a target longer than the design structure,
    whole-complex alignment metrics requiring identical residue counts
    (`<tool>_complex_rmsd_all` and `<tool>_complex_tm_score_all`) evaluate to `NaN`.
    Structurally aligned pruned variants and target-aligned binder RMSD
    (`<tool>_target_aligned_binder_rmsd_*`) remain valid because the binder chain
    length is unchanged. The DockQ columns likewise remain valid: residues are paired
    by sequence alignment and each pair contributes the atoms both structures carry,
    so a longer target simply contributes the part it shares with the design.

!!! note "DockQ's receptor and ligand are the target and the binder"
    `--refold_target_chain` is the receptor and `--refold_binder_chain` the ligand,
    where the reference DockQ implementation instead always takes the larger chain as
    the receptor. `fnat` and `iRMSD` are symmetric and agree either way; `LRMSD`, and
    so the DockQ score, differ when the binder is the larger of the two. Requires
    `nf-binder-design-utils:0.1.8` or later, which carries `peppr`.

## Output columns

Prediction engines generate standardized columns in `combined_scores.tsv`, prefixed
with the engine name:

| Column | Description |
|--------|-------------|
| `<tool>_iptm_max`, `<tool>_ipsae_max`, `<tool>_plddt_max`, ... | Prediction confidence metrics, per `--refold_metrics` |
| `<tool>_target_aligned_binder_rmsd_all` | Binder RMSD relative to the design pose after superimposing on the target |
| `<tool>_target_aligned_binder_dockq` | DockQ of the refolded complex against the design pose, 0 to 1. The usual reading is <0.23 incorrect, 0.23-0.49 acceptable, 0.49-0.80 medium, >=0.80 high |
| `<tool>_target_aligned_binder_dockq_fnat`, `_dockq_fnonnat`, `_dockq_irmsd`, `_dockq_lrmsd` | The components DockQ is built from: fraction of design contacts recovered, fraction of predicted contacts that are not in the design, interface backbone RMSD and ligand (binder) RMSD |
| `<tool>_target_aligned_binder_interface_lddt` | Interface lDDT of the refolded complex against the design |
| `<tool>_target_aligned_binder_clash_count` | Atom clashes within the refolded complex. A property of the prediction alone, not a comparison |
| `<tool>_complex_rmsd_all` | Whole-complex RMSD relative to the design |
| `<tool>_monomer_vs_design_rmsd_all` | RMSD of binder folded alone versus its designed conformation |
| `<tool>_monomer_vs_complex_rmsd_all` | RMSD of binder folded alone versus binder in the refolded complex |
| `<tool>_interface_dG`, `<tool>_interface_dSASA`, ... | BindCraft interface metrics evaluated on the refolded complex |
| `<tool>_selected_batch`, `<tool>_selected_model` | Sample identifier corresponding to the reported RMSD and interface metrics |
| `consensus_z`, `n_tools` | Cross-engine consensus z-score and count of contributing tools, prioritizing designs supported by multiple independent predictors |

`refold_sample_scores.tsv`, published in the same directory, records individual sample
predictions (keyed by `id`, `tool`, `batch`, `model`) without aggregation, enabling
custom post-processing.

## Reporting parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--refold_with` | none | Comma-separated list of prediction engines |
| `--refold_metrics` | iptm/ipsae/plddt/ptm/pae | Summary metrics included in `combined_scores.tsv` |
| `--refold_structure_stat` | `selected` | Method for aggregating structural metrics across samples: `selected`, `min`, `first`, `mean`, or `max` |
| `--refold_binder_monomer` | `true` | Folds each binder alone for monomer-versus-complex and monomer-versus-design RMSD comparisons |
| `--refold_monomer_single_sequence` | `true` | Predicts binder monomer structures from sequence alone without MSA generation |
| `--refold_interface_scoring` | `true` | Calculates BindCraft interface metrics on refolded complexes |
| `--refold_dockq` | `true` | Calculates DockQ, its components, interface lDDT and a clash count for the target-aligned binder comparison |
| `--output_rmsd_aligned` | `false` | Saves superimposed coordinate files from RMSD comparisons |
| `--n_predictions` | `5` | Generated structures per design per engine |

With `--refold_structure_stat selected`, RMSD and interface metrics are extracted from
the top-ranked sample (ranked by `--consensus_metric`, default ipSAE) rather than
aggregated across samples independently. Because structural metrics describe a specific
molecular conformation, taking independent minimums across different metrics would produce
inconsistent composite values. Setting `min`, `mean`, `max`, or `first` applies that
statistical operator across all structural columns.

## Fold stack configuration

Refolding uses the same underlying prediction pipeline as `--method fold` and
`--method fold_pulldown`. Alignment and template configurations derive from standard
fold stack parameters: `--msa_method`, `--create_target_msa`, `--target_msa`,
`--use_msa_server`, `--templates`, `--colabfold_envdb`, and `--uniref30`. Prediction
engine flags (`--af2_*`, `--boltz_*`, `--protenix_*`, `--af3_*`, `--openfold3_*`,
`--chai1_*`, `--esmfold2_*`) apply as documented under `--method fold`.

RosettaFold3 requires specific refold parameter prefixes to distinguish them from
generation-time parameters in `--method rfd3`: `--refold_rf3_ckpt_path`,
`--refold_rf3_num_steps`, `--refold_rf3_n_recycles`, `--refold_rf3_batch_size`,
`--refold_rf3_early_stopping_plddt_threshold`, and `--refold_rf3_seed`. When declared,
the fold stack reads these prefixed parameters.

## Renamed and deprecated parameters

Parameters previously prefixed with `--full_refold_*` have been renamed to
`--refold_*` (`--full_refold_with` is `--refold_with`, `--full_refold_max` is
`--refold_max`, and `--full_refold_filter_sort` is `--refold_sort_by`).

Deprecated parameters accepted with a warning:

| Deprecated | Replacement |
|------------|-------------|
| `--refold_af2ig_filters` | `--refold_filter` |
| `--refold_oracles` | `--refold_with` |
| `--refold_oracles_max` | `--refold_max` |

Removed parameters (alignment and template options now use fold stack parameters directly):

| Removed | Replacement |
|---------|-------------|
| `--refold_create_target_msa` | `--create_target_msa` |
| `--refold_use_msa_server` | `--use_msa_server` |
| `--refold_target_templates` | `--templates` |
| `--refold_alignment` | `--target_msa` |

`--refold_target_fasta` remains under its own name because target sequence substitution
is unique to refolding (see [Target sequence configuration](#target-sequence-configuration)).

`--boltz_refold_batch_size` has been removed. Use `--refold_boltz_batch_size` to
configure sample batching per job.

## See also

- [RFdiffusion Workflows](../workflows/rfdiffusion.md)
- [RFdiffusion3 Workflow](../workflows/rfd3.md)
- [RFD4-Proteina Workflow](../workflows/rfd4-proteina.md)
- [Fold Workflow](../workflows/fold.md): details on prediction engine flags
