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
3. **Refold** (optional): Boltz-2 refolds the complex and the binder alone, giving confidence and
   RMSD metrics. This reuses the same scoring machinery as the RFdiffusion workflows.
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

### Refolding and search

Inline scoring is switched off (`~metrics`), because RFD4-Proteina's own refolding metrics need an
AlphaFold 3 / OpenFold3 environment that is not in the container. Scoring is done downstream with
Boltz-2 instead, using the same options as the other workflows:

```bash
  --full_refold_with boltz \
  --full_refold_max 50 \
  --full_refold_use_msa_server true
```

FoldSeek search is enabled with `--do_foldseek`. See the
[FoldSeek subworkflow](../subworkflows/foldseek.md) for its options.

### Scoring with other folding oracles (prototype)

`--full_refold_with boltz` is tied to Boltz-2, because its RMSD, ipSAE and BindCraft metrics come
out of Boltz-specific processes. `--refold_oracles` is the oracle-agnostic alternative: it takes
each design's binder sequence, co-folds it against the target with whichever predictors you name,
and adds tool-prefixed columns to `combined_scores.tsv`.

```bash
  --refold_oracles boltz,chai1,openfold3 \
  --refold_oracles_max 50 \
  --n_predictions 3
```

Any engine `--method fold_pulldown` supports can go in the list: `af2`, `af2_mono`, `boltz`,
`rf3`, `protenix`, `af3`, `openfold3`, `chai1`, `esmfold2`, `esmfold2_fast`. The columns added are
`<tool>_iptm_max`, `<tool>_ipsae_max` and so on, plus one `consensus_z` and `n_tools` per design:
a cross-tool z-score, so a design several independent oracles agree on ranks above one that only
its designer liked.

This path is governed by the fold-family parameters (`--n_predictions`, `--msa_method`,
`--create_target_msa`, `--use_msa_server`, and each engine's own flags), not by the
`--full_refold_*` ones, which stay with the Boltz path. The two are independent, and a run may use
both.

The trade-off against `--full_refold_with boltz`: breadth across oracles and a consensus score,
but no RMSD-to-design and no binder-monomer metrics, because the oracles are given the binder
sequence rather than the design structure. For triage you usually want both - Boltz-2 for how
closely the refold reproduces the designed pose, the oracle consensus for how many independent
models believe the interface at all.

Implemented by `subworkflows/local/refold_oracles.nf`, which holds nothing rfd4-specific. It takes
designs as `tuple(meta, structure)`, so the other design workflows can adopt it unchanged. Nextflow
does not apply `params` declarations made in an included script, so the adopting workflow also has
to copy the parameter block marked "Oracle refold scoring" from `workflows/rfd4_proteina.nf`.

## Outputs

| Path | Contents |
|------|----------|
| `<outdir>/rfd4_proteina/designs/` | Generated complexes, one CIF per design |
| `<outdir>/rfd4_proteina/designs/rfd4_designs_batch*.tsv` | Design id, source example and seed per batch |
| `<outdir>/rfd4_proteina/combined_scores.tsv` | Every design with its sequence and metrics |
| `<outdir>/rfd4_proteina/binders.fasta` | Binder sequences |
| `<outdir>/rfd4_proteina/filter/` | Filter scores and accepted/rejected designs, with `--rfd4_filters` |
| `<outdir>/boltz_refold/` | Boltz-2 refolds and scores, with `--full_refold_with boltz` |
| `<outdir>/rfd4_proteina/refold_oracles/` | Per-oracle predictions and scores, with `--refold_oracles` |

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
`RFD4_FILTER_DESIGNS` and `COMBINE_RFD4_SCORES`), Boltz-2 and FoldSeek, so generation, filtering,
Boltz-2 refolding, FoldSeek search and scoring all run. Two limitations remain:

- FoldSeek runs on CPU, because upstream's GPU build is x86_64 only. The searches here are CPU
  searches, so this is a difference in speed rather than in what works.
- `--refold_oracles` accepts engines that have no arm64 image. `boltz` is the one that runs on
  this platform; AF2, AF3, Protenix, Chai-1, RosettaFold3, OpenFold3 and ESMFold2 do not.

On an amd64 cluster every step runs as normal without that config, except that there is no amd64
`rfd4-proteina` image built yet.

## Troubleshooting

**The run hangs before the first sampling step, with the GPU idle.** `rfd4-generate` can deadlock
on its dataloader when it uses worker processes. `--rfd4_num_workers` defaults to `0` for that
reason; if you raised it, lower it back.

**The first task is slow.** `torch.compile` warms up the trunk on the first run, which takes
several minutes. It is worth paying once for a real run; for a smoke test add
`--rfd4_compile_lengths false` and `--rfd4_inference_experiment quickstart`.

**Out of memory.** Lower `--rfd4_batch_size` to 1 and confirm the longest input fits before
raising it. Memory use scales with the token count of the whole complex, target included.
