# RFD4-Proteina Workflow

!!! warning "Under development"
    The RFD4-Proteina workflow is currently under development. It is not feature complete and is not yet well tested.

The `--method rfd4_proteina` workflow designs binders with
[RFD4-Proteina](https://github.com/RosettaCommons/RFD4-Proteina), the flow-matching biomolecular
design model from the Baker Lab and NVIDIA GenAIR.

RFD4-Proteina co-designs sequence and structure in a single pass. Unlike `--method rfd3`, there is
no separate inverse-folding step between generation and scoring: each generated design contains
both backbone coordinates and an assigned binder sequence.

## Workflow overview

1. Generation: `rfd4-generate` samples designs from a condition specification. Requested designs
   are partitioned across Nextflow tasks with independent random seeds to run in parallel across GPUs.
2. Filtering (optional): structural filters evaluated on the generated binder chain (for example, `rg<25`).
3. Refolding (optional): independent structure predictors re-fold the complex and the binder monomer,
   evaluating prediction confidence, RMSD, and interface metrics using the shared
   [refold scoring](../subworkflows/refold.md) subworkflow.
4. FoldSeek (optional): structural similarity search against CATH50 or other target databases.
5. Aggregation: outputs are merged into `combined_scores.tsv` containing binder sequences and all
   computed metrics.

## Inputs

### Condition specification

The specification follows RFD4-Proteina's YAML or JSON
[condition format](https://rfd4-proteina.org/latest/reference/yaml-json.html), supplied via
`--rfd4_config`. This interface provides access to the complete conditioning vocabulary, including
hotspots, burial, hydrogen bonds, coordinate groups, pi-stacking, CATH labels, and sequence
probabilities:

```yaml
pdl1_binder:
  input: PDL1.cif                                   # relative to this file
  contig: "50-120, /0, A17-131"                     # 50-120 residue binder chain, then target
  conditions:
    C_CRDGRP: {1: [{select: "sel('A/*/*/[N,CA,C,O]')"}]}   # target backbone geometry, free pose
    C_SEQ: {true: [{select: "sel('A')"}]}                  # target sequence is fixed
    C_HOT: {true: [{select: "sel('A/*/56')"}, {select: "sel('A/*/115')"}]}   # hotspots
```

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd4_proteina \
  --rfd4_config pdl1_binder.yaml \
  --rfd4_n_designs 20
```

Input structure paths are defined inside the specification file and resolve relative to its parent
directory; do not supply `--input_pdb` when using `--rfd4_config`.

> **Note:** Residue numbers in `contig` and `conditions` refer to mmCIF `label_seq_id` rather than
> author numbering. PyMOL displays `auth_seq_id` unless configured with `set cif_use_auth, off`.

### Parameter generation mode

For standard binder design against a single target structure, specifying `--input_pdb` with `--contigs`
generates a specification file automatically. This uses standard protein binder defaults: target
backbone geometry as a free-pose coordinate group, fixed target sequence, cysteines excluded from
the binder, and specified hotspots:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd4_proteina \
  --input_pdb PDL1.cif \
  --contigs "50-120, /0, A17-131" \
  --hotspot_res "A56,A115,A123" \
  --rfd4_n_designs 20
```

The designed segment is named `BINDER` in the generated contig (`"BINDER:50-120, /0, A17-131"`),
which lets the cysteine exclusion select it as `C_SEQ_PROB_CYS: {0.0: [{select: "sel('*/BINDER')"}]}`.
`--contigs` may name the designed segment itself, with or without a space after the colon
(`"BINDER: 50-120, /0, A17-131"`); the name is kept as given and the cysteine exclusion selects it.
A contig with more than one designed segment is left unnamed and gets no cysteine exclusion.

The generated specification is saved to `<outdir>/rfd4_proteina/input/config/<design_name>.json` and can be copied
and modified for custom runs. To allow cysteines, copy that file, drop the `C_SEQ_PROB_CYS`
condition, and pass it with `--rfd4_config`.

### Contig syntax

RFD4-Proteina contig strings use a comma-separated format:

| Token | Meaning |
|-------|---------|
| `/0` | Chainbreak: terminates current chain and starts a new chain |
| `A17-131` | Preserved segment: retained from input structure with original chain and residue numbering |
| `50-120` | Expandable segment: variable-length segment sampled by the model |
| `E1:2-20` | Named expandable segment: can be targeted by `conditions:` |

Tokens between chainbreaks define a single output chain. A chain containing a preserved segment
inherits that segment's chain ID; a chain without preserved segments is designated as the binder
and receives the first available unused letter. For `"50-120, /0, A17-131"`, this produces target `A`
and binder `B`. Custom chain IDs can be set using `--rfd4_target_chain` and `--rfd4_binder_chain`.

Each contig must include at least one preserved segment.

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--rfd4_config` | | Condition specification (YAML or JSON). Mutually exclusive with `--input_pdb` |
| `--input_pdb` | | Target structure path for parameter mode |
| `--contigs` | | Contig definition string for parameter mode |
| `--hotspot_res` | | Hotspot residues (e.g. `A56,A115,A123`) for parameter mode |
| `--design_name` | `rfd4` | Prefix for generated output files |
| `--rfd4_n_designs` | `1` | Total designs to generate |
| `--rfd4_designs_per_task` | `1` | Designs generated per Nextflow task; increasing amortises model loading time |
| `--rfd4_batch_size` | `1` | Designs processed concurrently on the GPU within a single task |
| `--rfd4_num_workers` | `0` | Dataloader worker processes (keep at 0 to avoid dataloader deadlocks) |
| `--rfd4_seed` | `0` | Base seed; task seeds derive from base seed plus task index |
| `--rfd4_nsteps` | model default (60) | Number of integration steps in sampling schedule |
| `--rfd4_recycle_strategy` | `medium/front_loaded` | Trunk recycle schedule (e.g. `thinking_medium`, `high/midpeak`) |
| `--rfd4_metrics_tags` | `ppi_refolding,hotspot_conditioning` | Inline RFD4-Proteina metrics (comma-separated); `false` disables them. Ignored with `quickstart` |
| `--rfd4_inference_experiment` | `inference_on_user_inputs` | Experiment configuration preset (`quickstart` for fast test runs) |
| `--rfd4_compile_lengths` | `true` | Set to `false` to execute trunk eagerly without `torch.compile` warm-up |
| `--rfd4_partial_trajectory_fraction` | | Partial diffusion: executes only the specified final fraction of `nsteps` |
| `--rfd4_inverse_folding` | `false` | Keeps input backbone fixed and resamples sequence only |
| `--rfd4_filters` | | Semicolon-separated filters applied to generated designs (e.g. `rg<25`) |
| `--rfd4_extra_args` | | Additional Hydra `key=value` overrides passed to `rfd4-generate` |
| `--rfd4_target_chain` | derived | Target chain ID in output complex |
| `--rfd4_binder_chain` | derived | Binder chain ID in output complex |

`rfd4-generate` uses [Hydra](https://hydra.cc), so `--rfd4_extra_args` accepts `key=value` pairs
(for example, `--rfd4_extra_args "generation.args.guidance_w=1.2 precision=bf16-mixed"`).

### Refold validation

RFD4-Proteina's own inline metrics run by default (`--rfd4_metrics_tags
ppi_refolding,hotspot_conditioning`), using the OpenFold3 weights baked into the container; set
`--rfd4_metrics_tags false` to disable them.

Their per-design values are merged into `combined_scores.tsv`, prefixed `rfd4_` with dots
replaced by underscores: `rfd4_oracle_ppi_overall_iptm`, `rfd4_oracle_ppi_plddt_coi`,
`rfd4_hotspot_obedience_hotspot_obedience`, `rfd4_refolded_backbone_rmsd_rmsd`, `rfd4_ppi_pass`
and so on. The Hydra configuration columns that `per_example_metrics.csv` repeats on every row
are dropped. `rfd4_oracle_ppi_refolded_structure_path` and `rfd4_oracle_ppi_confidence_json_path`
are rewritten to point at the published copy under `batches/batch<n>/`, so they resolve relative
to `<outdir>/rfd4_proteina/`. The oracle used for `ppi_refolding` is OpenFold3, which is why its
own file names say `af3_ppi`.

Designs are additionally
validated downstream: each generated binder sequence is re-folded against its target structure using
independent prediction engines to score interface metrics and pose accuracy.

See [Refold Scoring](../subworkflows/refold.md) for shared options:

```bash
  --refold_with boltz,chai1,openfold3 \
  --refold_max 50 \
  --use_msa_server
```

Candidate selection uses `--refold_filter`, `--refold_sort_by`, and `--refold_max`. Because
RFD4-Proteina produces no intermediate co-folding metrics, filters evaluate against structural
metrics computed by `--rfd4_filters` (`rg` and scripts in `bin/filters.d/`). When `--rfd4_filters`
is unset, use `--refold_max` alone to select the first N designs. Format `--refold_sort_by` with
an equals sign (e.g. `--refold_sort_by='rg'`) to prevent Nextflow from interpreting leading hyphens
as CLI flags. Selection decisions are recorded in `refold/refold_selection.tsv`.

Supported engines include `af2`, `af2_mono`, `boltz`, `rf3`, `protenix`, `af3`, `openfold3`,
`chai1`, `esmfold2`, and `esmfold2_fast`. Each engine appends standard columns to `combined_scores.tsv`:

| Column | Description |
|--------|-------------|
| `<tool>_iptm_max`, `<tool>_ipsae_max`, `<tool>_plddt_max`, ... | Prediction confidence metrics, per `--refold_metrics` |
| `<tool>_target_aligned_binder_rmsd_all` | Binder RMSD relative to the designed pose after superimposing on the target |
| `<tool>_complex_rmsd_all` | Whole-complex RMSD relative to the design |
| `<tool>_monomer_vs_design_rmsd_all` | Binder monomer RMSD against designed conformation |
| `<tool>_monomer_vs_complex_rmsd_all` | Binder monomer RMSD against binder in refolded complex |
| `<tool>_interface_dG`, `<tool>_interface_dSASA`, ... | BindCraft interface metrics evaluated on refolded complex |
| `<tool>_selected_batch`, `<tool>_selected_model` | Sample identifier for reported structural and interface metrics |
| `consensus_z`, `n_tools` | Cross-engine consensus z-score and count of contributing tools |

#### Sample generation and aggregation

By default, refolding runs `--n_predictions 5` structures per design per engine as 5 separate
single-sample jobs with distinct seeds (42, 43, ...) rather than one job with 5 samples. Engines that
batch samples internally share trunk computation across samples, meaning multiple samples from one seed
differ only in diffusion noise, whereas multiple seeds evaluate independent trajectories. Running
separate seeds incurs additional model loading overhead (roughly 45 seconds per load for Boltz on a GB10).

Batching behaviour is controlled per engine via `--refold_<engine>_batch_size` (and `--protenix_batch_size`).
For example, `--n_predictions 6 --refold_boltz_batch_size 3` runs 2 seeds with 3 samples each.

Engines that support multi-seed execution within a single task (`esmfold2`, `esmfold2_fast`, `openfold3`,
`protenix`, and `af3` via `--af3_seeds`) execute all seeds in one invocation to minimize loading
overhead. Boltz and RF3 execute one task per seed. Disable grouping across all engines with
`--fold_seed_list false`.

Metrics in `combined_scores.tsv` are aggregated across samples:

| Group | Aggregation method | Description |
|-------|--------------------|-------------|
| Confidence (`_max`, `_mean`, `n`) | Direct summary aggregates | `iptm_max` reports maximum sample score; `iptm_mean` reports average |
| Structural and interface metrics | `--refold_structure_stat` (default: `selected`) | Extracted from the top-ranked sample based on `--consensus_metric` (default: ipSAE) |

Structural and interface values are taken from a single representative sample so that reported
metrics correspond to a physically consistent coordinate set. Selecting on the same metric used for
confidence ranking ensures `<tool>_ipsae_max` and `<tool>_complex_rmsd_all` describe the identical
predicted structure. `<tool>_selected_batch` and `<tool>_selected_model` identify the source sample.

`refold_sample_scores.tsv` records unaggregated values for all samples. Setting
`--refold_structure_stat min|mean|max|first` applies that single operator across all structural columns.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--refold_with` | none | Comma-separated list of prediction engines |
| `--refold_filter` | none | Filter expression applied to design scores (requires `--rfd4_filters`) |
| `--refold_sort_by` | none | Metric used to rank candidates before applying `--refold_max` |
| `--refold_max` | all | Maximum designs to refold |
| `--n_predictions` | `5` | Generated structures per design per engine |
| `--fold_seed_list` | `true` | Passes full seed list to engines supporting multi-seed invocation |
| `--refold_metrics` | iptm/ipsae/plddt/ptm/pae | Summary metrics included in `combined_scores.tsv` |
| `--refold_structure_stat` | `selected` | Structural metric aggregation method |
| `--refold_binder_monomer` | `true` | Generates isolated binder predictions |
| `--refold_monomer_single_sequence` | `true` | Folds isolated binders without MSA generation |
| `--refold_interface_scoring` | `true` | Calculates BindCraft interface metrics on refolded complexes |

Fold stack options match those in `--method fold` and `--method fold_pulldown` (`--n_predictions`,
`--msa_method`, `--create_target_msa`, `--use_msa_server`, `--templates`).

RosettaFold3 parameters for refolding use dedicated prefixes to avoid collision with generation-time
parameters: `--refold_rf3_ckpt_path`, `--refold_rf3_num_steps`, `--refold_rf3_n_recycles`,
`--refold_rf3_batch_size`, `--refold_rf3_early_stopping_plddt_threshold`, and `--refold_rf3_seed`.

The subworkflow is implemented in `subworkflows/local/refold.nf`, accepting input tuples of
`tuple(meta, structure, scores)`. Shared defaults reside in `conf/refold.config` and
`conf/refold_fold_stack.config`.

FoldSeek structural similarity search is enabled with `--do_foldseek` (see the
[FoldSeek subworkflow](../subworkflows/foldseek.md)).

## Outputs

| Path | Contents |
|------|----------|
| `<outdir>/rfd4_proteina/input/config/<design_name>.json` | Generated condition spec (parameter mode) |
| `<outdir>/rfd4_proteina/input/prepared/` | Input structure converted to mmCIF |
| `<outdir>/rfd4_proteina/designs/` | Generated complex structures (one mmCIF per design) |
| `<outdir>/rfd4_proteina/designs/rfd4_designs_batch*.tsv` | Design identifier, source example, and seed per batch |
| `<outdir>/rfd4_proteina/batches/batch<n>/` | Raw `rfd4-generate` run directory for each batch: inference metrics, and each sample's generated structure with the inline refold and its confidences |
| `<outdir>/rfd4_proteina/combined_scores.tsv` | Consolidated table of binder sequences and all computed metrics |
| `<outdir>/rfd4_proteina/binders.fasta` | Generated binder sequences in FASTA format |
| `<outdir>/rfd4_proteina/filter/` | Filter scores and accepted/rejected design lists (with `--rfd4_filters`) |
| `<outdir>/rfd4_proteina/refold/` | Refolded structures, per-engine scores, and interface tables (with `--refold_with`) |
| `<outdir>/rfd4_proteina/refold/refold_sample_scores.tsv` | Unaggregated per-sample scores for all refold evaluations |
| `<outdir>/rfd4_proteina/refold/refold_selection.tsv` | Candidate selection records including verdicts, ranks, and metric values |

Output structure filenames follow the pattern `<design_name>_<uid>_b<task>_<example>_<replicate>.cif`.

## Running on arm64

The `rfd4-proteina` container is built for `linux/arm64` (compatible with NVIDIA DGX Spark GB10 systems).

Because other pipeline images target `linux/amd64`, executing them on an arm64 host requires
`conf/platforms/arm64.config`, which maps tasks to arm64-compatible container builds:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  -c conf/platforms/arm64.config \
  -c nextflow.docker.config \
  --method rfd4_proteina --rfd4_config pdl1_binder.yaml
```

The `conf/platforms/arm64.config` file supports `nf-binder-design-utils`, Boltz-2, ESMFold2, RF3,
Protenix, BindCraft, and FoldSeek. Current arm64 platform limitations:

- FoldSeek runs in CPU mode because upstream GPU binaries are compiled for x86_64 only.
- Supported refold engines on arm64: `boltz`, `esmfold2`, `esmfold2_fast`, `rf3`, and `protenix`.
  AF2, AF3, Chai-1, and OpenFold3 do not currently provide arm64 container images.

!!! note "RF3 on GB10 requires the `-ngc2606` container image"
    `rc-foundry:0.2.0_arm64-weights` is compiled against CUDA 12.8, whose nvrtc compiler does not
    target `sm_121`. PyTorch's jiterator compiles `erfinv` for the host architecture, causing RF3 to
    fail with `nvrtc: error: invalid value for --gpu-architecture (-arch)`. The configuration points
    to `0.2.0_arm64-weights-ngc2606`, built on CUDA 13.3, resolving this issue.

!!! note "Protenix on GB10 requires the dedicated arm64 image"
    Protenix uses a fused layer-norm CUDA extension built for `sm_70` through `sm_100`. On GB10 systems
    (`sm_121`), inference fails at layer normalization with `no kernel image is available for execution on the device`.
    The arm64 container image compiles support for `sm_120` and `sm_121`.

On amd64 clusters, all pipeline steps execute using default configurations without these options.

## Troubleshooting

- Execution hangs before sampling with GPU idle: `rfd4-generate` can deadlock on its dataloader
  when multiprocessing workers are active. Keep `--rfd4_num_workers 0`.
- Slow initial task execution: `torch.compile` warms up the model trunk during the first run.
  This compilation occurs once per execution; for rapid testing, specify `--rfd4_compile_lengths false`
  and `--rfd4_inference_experiment quickstart`.
- Out of memory: decrease `--rfd4_batch_size` to 1. Memory consumption scales with total sequence
  length across both target and binder chains.
