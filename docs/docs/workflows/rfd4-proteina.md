# RFD4-Proteina Workflow

The `--method rfd4_proteina` workflow designs binders with
[RFD4-Proteina](https://github.com/RosettaCommons/RFD4-Proteina), the flow-matching biomolecular
design model from the Baker Lab and NVIDIA GenAIR.

RFD4-Proteina co-designs sequence and structure in one pass. Unlike `--method rfd3`, there is no
inverse-folding step between generation and scoring: each design comes out of the model as a
complex that already carries a binder sequence.

## Overview

1. **Generate**: `rfd4-generate` samples designs from a condition spec. The requested designs are
   split across Nextflow tasks, each with its own seed, so they run in parallel across GPUs.
2. **Filter** (optional): structural filters on the generated binder chain, for example `rg<25`.
3. **Refold** (optional): one or more independent predictors fold the complex and the binder alone,
   giving confidence, RMSD and interface metrics. This is the same
   [refold scoring](../subworkflows/refold.md) every design workflow runs.
4. **FoldSeek** (optional): structural similarity search against CATH50 or another database.
5. **Combine**: one `combined_scores.tsv` with the binder sequences and every metric computed.

## Inputs

### Condition spec (the default interface)

The spec is RFD4-Proteina's own YAML or JSON
[condition specification](https://rfd4-proteina.org/latest/reference/yaml-json.html), passed with
`--rfd4_config`. It is the baseline interface here because it gives access to the full
conditioning vocabulary: hotspots, burial, hydrogen bonds, coordinate groups, pi-stacking, CATH
labels and sequence probabilities.

```yaml
pdl1_binder:
  input: PDL1.cif                                   # relative to this file
  contig: "50-120, /0, A17-131"                     # a new 50-120 residue chain, then the target
  conditions:
    C_CRDGRP: {1: [{select: "sel('A/*/*/[N,CA,C,O]')"}]}   # target backbone geometry, free pose
    C_SEQ: {true: [{select: "sel('A')"}]}                  # target sequence is known
    C_HOT: {true: [{select: "sel('A/*/56')"}, {select: "sel('A/*/115')"}]}   # hotspots
```

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd4_proteina \
  --rfd4_config pdl1_binder.yaml \
  --rfd4_n_designs 20
```

Input structure paths are read from the spec, so do not also pass `--input_pdb`. Paths inside the
spec resolve against the spec file's own directory.

> **Note:** residue numbers in a `contig` and in `conditions` refer to the mmCIF `label_seq_id`,
> not the author numbering. PyMOL shows `auth_seq_id` unless you `set cif_use_auth, off`.

### Params mode

For a plain binder design against one target, `--input_pdb` with `--contigs` generates a spec for
you, using the upstream default recipe for protein binders: target backbone geometry as one
free-pose coordinate group, the target sequence held fixed, and hotspots where given.

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd4_proteina \
  --input_pdb PDL1.cif \
  --contigs "50-120, /0, A17-131" \
  --hotspot_res "A56,A115,A123" \
  --rfd4_n_designs 20
```

The generated spec is published to `<outdir>/rfd4_proteina/<design_name>.json`, so it is a
reasonable starting point to copy and extend as a `--rfd4_config` of your own.

### Contig syntax

RFD4-Proteina's contigs differ from RFdiffusion's. One comma-separated string:

| Token | Meaning |
|-------|---------|
| `/0` | Chainbreak; ends the current output chain and starts a new one |
| `A17-131` | Preserved segment, kept from the input by its original chain and residue numbering |
| `50-120` | Expandable segment of variable length, which the model samples |
| `E1:2-20` | Named expandable segment, so `conditions:` can select it |

A run of tokens between chainbreaks is one output chain. A run containing a preserved segment
takes that segment's chain id; a run with none is the designed binder and gets the first unused
letter. For `"50-120, /0, A17-131"` that means target `A` and binder `B`, which is what the
workflow assumes unless you set `--rfd4_target_chain` / `--rfd4_binder_chain`.

A contig must contain at least one preserved segment.

## Key parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--rfd4_config` | | Condition spec (YAML or JSON). Mutually exclusive with `--input_pdb` |
| `--input_pdb` | | Target structure, params mode |
| `--contigs` | | Contig string, params mode |
| `--hotspot_res` | | Hotspot residues, eg `A56,A115,A123`, params mode |
| `--design_name` | `rfd4` | Prefix for output file names |
| `--rfd4_n_designs` | `1` | Total designs to generate |
| `--rfd4_designs_per_task` | `1` | Designs per Nextflow task; raise it to amortise model loading |
| `--rfd4_batch_size` | `1` | Designs processed in parallel on the GPU within a task |
| `--rfd4_num_workers` | `0` | Dataloader workers; leave at 0 (see Troubleshooting) |
| `--rfd4_seed` | `0` | Base seed; each task uses this plus its index |
| `--rfd4_nsteps` | model default (60) | Integration steps in the sampling schedule |
| `--rfd4_recycle_strategy` | model default | Trunk recycle schedule, eg `thinking_medium`, `high/midpeak` |
| `--rfd4_inference_experiment` | `inference_on_user_inputs` | Experiment preset; `quickstart` for a fast smoke test |
| `--rfd4_compile_lengths` | `true` | Set `false` to run the trunk eagerly and skip the `torch.compile` warm-up |
| `--rfd4_partial_trajectory_fraction` | | Partial diffusion: run only this final fraction of `nsteps` |
| `--rfd4_inverse_folding` | `false` | Hold the input backbone fixed and resample sequence only |
| `--rfd4_filters` | | Semicolon-separated filters on generated designs, eg `rg<25` |
| `--rfd4_extra_args` | | Additional Hydra `key=value` overrides passed to `rfd4-generate` |
| `--rfd4_target_chain` | derived | Target chain id in the generated complex |
| `--rfd4_binder_chain` | derived | Binder chain id in the generated complex |

`rfd4-generate` is a [Hydra](https://hydra.cc) application, so `--rfd4_extra_args` takes
`key=value`, never `--key value`. Anything in the
[CLI options reference](https://rfd4-proteina.org/latest/reference/cli-options.html) can go there,
for example `--rfd4_extra_args "generation.args.guidance_w=1.2 precision=bf16-mixed"`.

### Refold scoring

Inline scoring is switched off (`~metrics`), because RFD4-Proteina's own refolding metrics need an
AlphaFold 3 / OpenFold3 environment that is not in the container. Designs are scored downstream
instead, by folding them again: each design's binder sequence is folded against its target, the
result is compared with the design, and the rebuilt interface is scored.

This is the same refold step every design workflow runs; see [Refold Scoring](../subworkflows/refold.md) for the shared parameters. What follows is how it behaves here.

`--refold_with` names the engines:

```bash
  --refold_with boltz,chai1,openfold3 \
  --refold_max 50 \
  --use_msa_server
```

`--refold_filter`, `--refold_sort_by` and `--refold_max` choose which designs go through. Unlike the other design workflows, this one has no cofold of its own, so the scores those two select on are the structural metrics `--rfd4_filters` computes (`rg` and the rest of `bin/filters.d/`), and they exist only when `--rfd4_filters` is set. Without it, use `--refold_max` on its own, which takes the first N designs to be generated. Write `--refold_sort_by` with an `=`, as `--refold_sort_by='rg'`, otherwise Nextflow reads a leading `-` as one of its own options. Every design that reached the step is accounted for in `refold/refold_selection.tsv`.

Any engine `--method fold` supports can go in the list: `af2`, `af2_mono`, `boltz`, `rf3`,
`protenix`, `af3`, `openfold3`, `chai1`, `esmfold2`, `esmfold2_fast`. Every one of them produces
the same columns in `combined_scores.tsv`, prefixed with the engine name:

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

#### Samples per design

The refold runs `--n_predictions` structures per design per engine, 5 by default, as 5 fold jobs
of one sample each with seeds 42, 43, ... rather than one job of 5 samples. The engines that batch
samples share a single trunk and recycling pass across them, so five samples from one seed differ
only in the diffusion noise, while five seeds re-run the whole prediction. The disadvantage is five
model loads rather than one - on a GB10, Boltz pays roughly 45 s of fixed cost per load.

The split is set per engine by `--refold_<engine>_batch_size 1` (and `--protenix_batch_size`,
which this workflow declares itself). Setting one of them higher gives that engine a mixture:
`--n_predictions 6 --refold_boltz_batch_size 3` is 2 seeds x 3 samples.

Where an engine takes every seed in one invocation, it is given them all rather than run once per
seed, so the model is loaded once instead of five times. That applies to `esmfold2`,
`esmfold2_fast`, `openfold3` and `protenix`; AF3 does the same under its own `--af3_seeds` comma
list. Boltz and RF3 are single-seed per process, so they still run one job per seed.
`--fold_seed_list false` turns the grouping off everywhere.

Everything in `combined_scores.tsv` is then reduced across those samples, in two ways:

| Group | Reduced by | Meaning |
|-------|-----------|---------|
| Confidence (`_max`, `_mean`, `n`) | the summary's own aggregates | `iptm_max` is the best of the samples, `iptm_mean` their average |
| RMSD and interface scores | `--refold_structure_stat`, `selected` by default | read from the best-ranked sample, by `--consensus_metric` (ipSAE) |

The structural metrics are read from one sample rather than aggregated per metric, because they
all describe one set of coordinates. Taking the lowest RMSD and the lowest `interface_dG`
independently would give a row that no predicted structure matches, and the interface metrics have
no shared direction to take a best in anyway - `interface_dG` is better lower, `interface_dSASA`
and `interface_sc` better higher. Selecting on the same quantity the confidence ranking uses means
`<tool>_ipsae_max` and `<tool>_complex_rmsd_all` describe the same structure. `<tool>_selected_batch` and
`<tool>_selected_model` record which sample that was. Both are needed: an engine run once per seed
writes one sample per batch, while an engine handed the whole seed list returns every sample from
batch 1, and only the model label tells those apart.

`refold_sample_scores.tsv` has every sample's own values, unreduced, so any statistic over them can
be computed from it afterwards. `--refold_structure_stat min|mean|max|first` overrides the
selection with one statistic applied to every structural column, which suits RMSD and does not suit
the interface metrics.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--refold_with` | none | Comma-separated engines |
| `--refold_filter` | none | Threshold gate on the design scores, eg `'rg<25'` (needs `--rfd4_filters`) |
| `--refold_sort_by` | none | Rank by one metric before `--refold_max` takes the best; `-` prefix for descending |
| `--refold_max` | all | Maximum designs to refold |
| `--n_predictions` | `5` | Structures per design per engine, as that many seeds of one sample |
| `--fold_seed_list` | `true` | Give an engine every seed in one invocation where it supports that |
| `--refold_metrics` | iptm/ipsae/plddt/ptm/pae | Summary columns carried into `combined_scores.tsv` |
| `--refold_structure_stat` | `selected` | How several samples per design reduce to one RMSD and one interface score |
| `--refold_binder_monomer` | `true` | Also fold each binder alone |
| `--refold_monomer_single_sequence` | `true` | Fold that monomer from sequence alone |
| `--refold_interface_scoring` | `true` | BindCraft interface metrics on the refolded complex |

How the refold folds is governed by the fold-stack parameters - `--n_predictions`,
`--msa_method`, `--create_target_msa`, `--use_msa_server`, `--templates` and each engine's own
flags - the same ones `--method fold` and `--method fold_pulldown` use.

RF3 is the one exception to "each engine's own flags": here they are `--refold_rf3_ckpt_path`,
`--refold_rf3_num_steps`, `--refold_rf3_n_recycles`, `--refold_rf3_batch_size`,
`--refold_rf3_early_stopping_plddt_threshold` and `--refold_rf3_seed`. The plain
`--rf3_*` names belong to the design-time RF3 step in `--method rfd3`, where they mean different
things, and a workflow's parameters are one flat namespace. Where the prefixed names are declared
the fold stack reads only those, so passing `--rf3_num_steps` here has no effect.

Implemented by `subworkflows/local/refold.nf`, which holds nothing rfd4-specific. It takes designs
as `tuple(meta, structure, scores)`, and every design workflow uses it: `rfd`, `rfd_partial`,
`rfd3`, `rfd4_proteina`, `bindcraft` and `germinal`. Nextflow does not apply `params` declarations
made in an included script, so the shared defaults live in `conf/refold.config` and
`conf/refold_fold_stack.config` rather than in each workflow; a workflow declares only its own
`--refold_publish_dir`.

!!! note "Renamed and removed parameters"
    The whole `--full_refold_*` family is now `--refold_*` (`--full_refold_with` is `--refold_with`,
    `--full_refold_max` is `--refold_max`, `--full_refold_filter_sort` is `--refold_sort_by`).
    `--refold_oracles` and `--refold_oracles_max` were second parameters for the same idea, needed
    only because the refold used to be a Boltz-specific process that could not take another engine;
    they are accepted for now, with a warning, and folded into `--refold_with` and `--refold_max`.
    The MSA and template parameters that belonged to the Boltz path -
    `--refold_create_target_msa`, `--refold_use_msa_server`, `--refold_alignment` and
    `--refold_target_templates` - are gone; a run that sets one stops with the fold-stack
    replacement named. `--refold_target_fasta` is not one of them: it kept its name and its
    meaning, and refolds the selected designs against a target other than the one they were
    built against.

FoldSeek search is enabled with `--do_foldseek`. See the
[FoldSeek subworkflow](../subworkflows/foldseek.md) for its options.

## Outputs

| Path | Contents |
|------|----------|
| `<outdir>/rfd4_proteina/designs/` | Generated complexes, one CIF per design |
| `<outdir>/rfd4_proteina/designs/rfd4_designs_batch*.tsv` | Design id, source example and seed per batch |
| `<outdir>/rfd4_proteina/combined_scores.tsv` | Every design with its sequence and metrics |
| `<outdir>/rfd4_proteina/binders.fasta` | Binder sequences |
| `<outdir>/rfd4_proteina/filter/` | Filter scores and accepted/rejected designs, with `--rfd4_filters` |
| `<outdir>/rfd4_proteina/refold/` | Refolded structures, per-engine scores, RMSD and interface tables, with `--refold_with` |
| `<outdir>/rfd4_proteina/refold/refold_sample_scores.tsv` | One row per refolded sample, nothing reduced: `id`, `tool`, `batch`, `model` and every confidence, RMSD and interface score that sample has |
| `<outdir>/rfd4_proteina/refold/refold_selection.tsv` | Every design considered for refolding, with its selection verdict, rank and sort value |

Design file names are `<design_name>_<uid>_b<task>_<example>_<replicate>.cif`. The task index keeps
ids unique across the parallel tasks.

## Running on arm64

The `rfd4-proteina` container is built for `linux/arm64` (DGX Spark GB10 and similar), since
upstream does not publish an arm64 image.

The other images this workflow uses are published for `linux/amd64`, and on an arm64 host their
tasks fail immediately with `exec format error`. Add `conf/platforms/arm64.config`, which points
them at arm64 builds of the same versions, alongside a runtime config:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  -c conf/platforms/arm64.config \
  -c nextflow.docker.config \
  --method rfd4_proteina --rfd4_config pdl1_binder.yaml
```

`nextflow.docker.config` from the example directory runs the containers under Docker rather than
Apptainer, which is usually what is available on such a workstation.

That config covers `nf-binder-design-utils` (which runs `PREPARE_RFD4_INPUT`,
`RFD4_FILTER_DESIGNS` and `COMBINE_RFD4_SCORES`), Boltz-2, ESMFold2, RF3, Protenix, BindCraft and
FoldSeek, so generation, filtering, refolding, FoldSeek search and scoring all run. Two limitations
remain:

- FoldSeek runs on CPU, because upstream's GPU build is x86_64 only. The searches here are CPU
  searches, so this is a difference in speed rather than in what works.
- `--refold_with` accepts engines that have no arm64 image. `boltz`, `esmfold2`,
  `esmfold2_fast`, `rf3` and `protenix` run on this platform; AF2, AF3, Chai-1 and OpenFold3 do
  not.

!!! note "RF3 on a GB10 needs the -ngc2606 image"
    `rc-foundry:0.2.0_arm64-weights` is built on CUDA 12.8, whose nvrtc cannot target `sm_121`.
    Most kernels JIT from that image's `compute_120` PTX, but PyTorch's jiterator compiles
    `erfinv` for the exact device architecture, so RF3 dies with `nvrtc: error: invalid value for
    --gpu-architecture (-arch)` once the weights have loaded. The config points at
    `0.2.0_arm64-weights-ngc2606`, the same release on CUDA 13.3, which does not have this
    problem.

!!! note "Protenix on a GB10 needs the arm64 image specifically"
    Protenix's fused layer-norm is a CUDA extension compiled when the image is built, and upstream
    builds it for `sm_70` to `sm_100` only. On a GB10 (`sm_121`) inference reaches the first
    layer-norm call and dies with `no kernel image is available for execution on the device`. The
    arm64 image adds `sm_120` and `sm_121` to that list and rebuilds the extension.

On an amd64 cluster every step runs as normal without that config.

## Troubleshooting

**The run hangs before the first sampling step, with the GPU idle.** `rfd4-generate` can deadlock
on its dataloader when it uses worker processes. `--rfd4_num_workers` defaults to `0` for that
reason; if you raised it, lower it back.

**The first task is slow.** `torch.compile` warms up the trunk on the first run, which takes
several minutes. It is worth paying once for a real run; for a smoke test add
`--rfd4_compile_lengths false` and `--rfd4_inference_experiment quickstart`.

**Out of memory.** Lower `--rfd4_batch_size` to 1 and confirm the longest input fits before
raising it. Memory use scales with the token count of the whole complex, target included.
