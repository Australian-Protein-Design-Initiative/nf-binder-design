# BindCraft2 Workflow

!!! warning "Under development"
    The BindCraft2 workflow is currently under development. It is not feature complete and is not yet well tested.

Parallel [BindCraft2](https://github.com/PacesaLab/BindCraft2) execution for protein binder design across multiple GPUs.

!!! warning "BindCraft2 cannot be offered as a hosted service"
    BindCraft2 is under a [source-available, hosting-restricted license](https://github.com/PacesaLab/BindCraft2/blob/main/LICENSE). It is free to use for any purpose, including commercial use, on your own or your organisation's infrastructure, and you may share the designs it produces. However, you may not provide BindCraft2, or a service whose functionality derives from it, to third parties as a hosted, cloud, API, web application, workflow platform or SaaS offering without a separate commercial license from the copyright holders.

    This restriction carries through to `nf-binder-design`: the `--method bc2` workflow must not be deployed as a hosted service (or as an invocable tool or workflow step in a hosted platform) for third parties.

## Overview

The `--method bc2` workflow (alias `--method bindcraft2`) runs BindCraft2 campaigns in parallel batches, then gathers them with `bindcraft campaign_output` and renders an HTML report.

BindCraft2 is a rewrite of [BindCraft](bindcraft.md), not a version bump. It runs as a separate `--method bc2` workflow instead of a flag on `--method bindcraft`, adding:

- One campaign JSON instead of three settings files; set any option from the CLI via `--bc2_set`.
- Twelve binder modalities (`VHH`, `scFv`, `Fab`, `ARP`, `cyclic_peptide`, `multidomain`, `induced_fit`, ...) and nine composable design properties (`--bc2_humanize`, `--bc2_protease_stable`, ...).
- Multi-target campaigns, off-target detargeting, and PDB, mmCIF or FASTA inputs.
- mmCIF structure output and `i_pDAE` candidate ranking.

## Differences from standalone BindCraft2

Standalone BindCraft2 runs on a single GPU until it produces `N` accepted designs. This pipeline:

- Runs a **fixed trajectory budget** (`--bc2_n_traj`), split into batches of `--bc2_batch_size`
- Exits once the budget finishes, keeping runtime predictable
- Parallelises batches across available GPUs
- Assigns a distinct `campaign_seed` per batch to prevent duplicate trajectories

To target a specific count of accepted designs, start with a small budget (`--bc2_n_traj 100`) to measure the acceptance rate, then scale `--bc2_n_traj`.

### Adaptive tuning under batching

BindCraft2 includes two adaptive mechanisms that update settings mid-campaign based on the current campaign's trajectory table:

| Mechanism | Engages after |
|---|---|
| `autotune` | every 10 trajectories in one campaign |
| `desperation` | 750 fruitless trajectories in one campaign, then one escalation rung per 50 |

Each batch executes `--bc2_batch_size` trajectories, so both mechanisms remain inactive at the default batch size of 1. That suits fixed-budget, reproducible comparisons, though long standalone runs on difficult targets behave differently once adaptive tuning activates. The pipeline logs inactive mechanisms with the batch plan. Increase `--bc2_batch_size` if you want `autotune` to activate.

## Command-line Options

See available options with `--method bc2 --help`:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method bc2 --help
```

## Example Usage

### Params mode

Generate the campaign JSON from flat parameters, using the same parameter names as `--method bindcraft`:

```bash
#!/bin/bash

nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method bc2 \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56,A125" \
  --binder_length_range "60-150" \
  --bc2_n_traj 8 \
  --bc2_batch_size 2 \
  -profile local \
  -resume
```

### Config mode

Pass a BindCraft2 campaign JSON directly. Nextflow detects referenced files (`targets[].target_path` and `binder_scaffold`), stages them into each task, and rewrites their paths so relative entries resolve from any launch directory:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method bc2 \
  --bc2_config campaign.json \
  --bc2_modality VHH \
  --bc2_humanize \
  --bc2_n_traj 100 \
  --bc2_batch_size 5 \
  -profile local
```

### Built-in target presets

BindCraft2 includes preconfigured targets inside the container. Setting `"target"` in the campaign JSON loads the target structure, chains and hotspots from the container without staging external files:

```json
{
  "campaign_name": "hPDL1_VHH",
  "target": "hPDL1",
  "modality": "VHH",
  "humanize": true
}
```

Available presets: `hPDL1`, `mPDL1`, `hPD1`, `hIL2R`, `hIL7RA` and `dynorphin_a`. Run `bindcraft design --list-targets` inside the container to list current presets.

Modality and design properties are standard campaign JSON keys, so you can set them in JSON instead of passing `--bc2_modality` or `--bc2_*` CLI flags.

For running on SLURM, use `-profile slurm` with `--slurm_account=xxYY` and the appropriate `-c` site-specific configuration file (refer to and adapt from the [rfdiffusion example](rfdiffusion.md#parallel-tasks-on-an-hpc-cluster)).

## Key Parameters

| Flag | Description |
|------|-------------|
| `--bc2_config` | BindCraft2 campaign settings JSON. Mutually exclusive with `--input_pdb` |
| `--input_pdb` | Target PDB, mmCIF or FASTA; a directory or quoted glob runs each structure. Mutually exclusive with `--bc2_config` |
| `--target_chains` | Target chain IDs, comma-separated (params mode) |
| `--hotspot_res` | Hotspot residues, e.g. `A54,A56` or `B12-16` (params mode) |
| `--hotspot_subsample` | Fraction of each target's hotspots to keep per batch (both modes) |
| `--binder_length_range` | `60-150` range, `80` exact, or `60,80,100` discrete choices (params mode) |
| `--bc2_n_traj` | Total trajectories to run |
| `--bc2_batch_size` | Trajectories per parallel task |
| `--bc2_seed` | Base `campaign_seed`. Batch seed is `bc2_seed + seedIndex` (`seedIndex` is the batch index in config mode, or `pdbIndex * nBatches + batchIndex` for several structures in params mode) |
| `--gpu_devices` | Specify multiple GPUs, e.g. `--gpu_devices=0,1` (use only for `-profile local`) |
| `--do_foldseek` | Run a FoldSeek structural novelty search over the ranked designs |

## Hotspot subsampling

`--hotspot_subsample` samples `ceil(N * fraction)` hotspots per target (minimum 1) independently for each batch, so batches explore different parts of the epitope.

```bash
--hotspot_res "A54,A56,A66,A115" --hotspot_subsample 0.5
```

With these flags, each batch selects two of the four hotspots using its `campaign_seed`. Sampling uses `--bc2_seed` for reproducibility; `--method bindcraft` seeds hotspot subsampling from the system clock instead.

Subsampling behaviour:

- Works in both params and config mode, since subsampling runs while staging the campaign JSON into each task.
- Treats each comma-separated entry as one hotspot, keeping or dropping ranges like `B12-16` as a unit.
- Samples each target independently in multi-target campaigns.
- Applies only to hotspots defined under `targets[]` in the campaign JSON. Hotspots from a built-in `"target"` preset are merged by BindCraft2 after staging and are not subsampled; the pipeline logs a warning if you combine `--hotspot_subsample` with a preset target.

## Passing BindCraft2 Settings

### `--bc2_set`

BindCraft2 exposes around 235 settings. Pass any of them via `--bc2_set` as `KEY=VALUE` pairs forwarded to `bindcraft design --set`:

```bash
--bc2_set 'binder_lengths=[70,90];filters.i_pTM.threshold=0.8;design_models=3'
```

- Separate assignments with `;` instead of `,`, because JSON values often contain commas.
- Dotted keys update nested fields: `filters.i_pTM.threshold=0.8` changes only that threshold.
- Values parse as JSON (falling back to plain strings), so `[70,90]` becomes a list and `humanize=true` becomes a boolean.
- Full setting list: [`settings/core/reference.json`](https://github.com/PacesaLab/BindCraft2/blob/main/settings/core/reference.json).

The pipeline manages `project_folder`, `campaign_seed` and `max_trajectories` per batch and rejects them in `--bc2_set`. Set `--bc2_seed`, `--bc2_n_traj` and `--bc2_batch_size` instead.

!!! warning "Paths in `--bc2_set` are not staged"
    Files referenced via `--bc2_set` are not staged; BindCraft2 resolves those paths relative to the task work directory. Define file inputs in the campaign JSON so they are detected and staged.

### Modality and core presets

| Flag | Description |
|------|-------------|
| `--bc2_modality` | Binder format, e.g. `VHH`, or several combined: `VHH,induced_fit`. One of `binder`, `large_binder`, `peptide`, `cyclic_peptide`, `homo_oligomer`, `multidomain`, `VHH`, `scFv`, `Fab`, `ARP`, `induced_fit`, `fold_switch` |
| `--bc2_core` | Core profile. `--bc2_core benchmark` pins the seed and disables `autotune` and `desperation`, for reproducible comparisons |
| `--bc2_metadata` | JSON file of fields to record with the campaign |

### Design properties

Composable boolean flags, each corresponding to a `settings/property/*.json` upstream:

| Flag | Description |
|------|-------------|
| `--bc2_humanize` | Human-like sequence preferences and an MHC-II anchor-score objective |
| `--bc2_bigbang` | Start gradient stages from the coordinates on hand |
| `--bc2_disulfide_staple` | Allow cysteine, require a detected disulfide |
| `--bc2_forced_targeting` | Concentrate the paratope on declared hotspots |
| `--bc2_initial_guess` | Re-predict the candidate from the trajectory pose |
| `--bc2_mixed_topology` | At most 50% helix, at least 20% beta sheet |
| `--bc2_protease_stable` | Penalise cleavage motifs and exposed loops and termini |
| `--bc2_termini_accessible` | Point both chain ends away from the target |
| `--bc2_termini_together` | Keep the N and C termini within 10 A |

## Reproducibility

Each batch uses its own `campaign_seed`, computed as `--bc2_seed + seedIndex` (batch index in config mode; `pdbIndex * nBatches + batchIndex` when params mode folds several structures). Keep `--bc2_seed`, `--bc2_n_traj` and `--bc2_batch_size` constant to reproduce a run; changing `--bc2_batch_size` alters batch seeds and resulting designs.

Pass `--bc2_core benchmark` to also fix BindCraft2's internal seed and disable adaptive tuning.

## Output Structure

Results are saved to `--outdir` in the `bc2/` subdirectory:

```
── bc2
│   ├── batches
│   │   ├── PDL1_0
│   │   │   ├── 1_Trajectories
│   │   │   │   └── !_Trajectories.csv
│   │   │   ├── 2_Refolded
│   │   │   │   └── !_Refolded.csv
│   │   │   ├── 3_Ranked
│   │   │   │   ├── !_Ranked.csv
│   │   │   │   └── PDL1_l82_s0a1b2c3_1.cif
│   │   │   ├── campaign_metadata.json
│   │   │   ├── summary.csv
│   │   │   ├── bc2.log
│   │   │   └── bc2_campaign.staged.json
│   │   └── PDL1_1
│   │       └── ...
│   ├── ranked
│   │   └── PDL1_l82_s0a1b2c3_1.cif
│   ├── bc2_report.html
│   ├── summary.csv
│   ├── ranked.csv
│   ├── refolded.csv
│   └── trajectories.csv
└── logs
    ├── dag_20250725_084959.html
    ├── report_20250725_084959.html
    ├── timeline_20250725_084959.html
    └── trace_20250725_084959.txt
```

- `bc2/ranked/`: flat directory of accepted mmCIF structures from all batches.
- `summary.csv`: combined campaign summary generated by `bindcraft campaign_output`.
- `ranked.csv`, `refolded.csv` and `trajectories.csv`: per-stage tables concatenated across batches with an added `campaign` column (`campaign_output` does not merge these).
- `bc2_report.html`: summary HTML report.
- `${outdir}/foldseek/{database}/`: FoldSeek output directory when `--do_foldseek` is enabled (outside `bc2/`). See [FoldSeek output format](../subworkflows/foldseek.md#output-format).

## Structural novelty search (FoldSeek)

Enable `--do_foldseek` to run a FoldSeek structural similarity search on ranked designs (using the same flags as the other workflows):

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method bc2 \
  --input_pdb 'input/PDL1.pdb' \
  --target_chains "A" \
  --do_foldseek \
  --foldseek_database CATH50 \
  -profile local
```

High TM-score hits indicate similarity to known folds; designs without close hits are structurally novel. FoldSeek outputs go to `${outdir}/foldseek/{database}/` (outside `bc2/`) and are not merged into `ranked.csv` or `bc2_report.html`.

FoldSeek runs after `bindcraft campaign_output` on the mmCIF files in `bc2/ranked/`, covering only designs retained after re-ranking. If no designs pass filtering, the pipeline skips FoldSeek without error.

Set `--foldseek_use_webserver true` to query the FoldSeek web API instead of a local database. See the [FoldSeek subworkflow docs](../subworkflows/foldseek.md#command-line-options) for all `--foldseek_*` options (database, search mode, output options, CATH annotation).

### Binder chain for FoldSeek

FoldSeek queries only the binder chains, excluding the target. BindCraft2 writes target chains first and binder chains last. The pipeline reads `binder_chains` from `campaign_metadata.json` and extracts the trailing binder chains.

By contrast, the default heuristic selects the chain with the fewest atoms, which fails for two BindCraft2 use cases:

- Multi-chain binders (`--bc2_modality Fab`, `scFv`, `multidomain`, or symmetric oligomers with `copies > 1`), where the heuristic would extract only one chain.
- Targets smaller than the binder (for example cropped or peptide targets), where the heuristic would extract the target chain.

If `campaign_metadata.json` is unreadable, the pipeline logs a warning and falls back to the smallest-chain heuristic.

## Unwrapped BindCraft2 subcommands

This pipeline invokes two BindCraft2 subcommands: `design` (per batch) and `campaign_output` (final aggregation). You can configure all `bindcraft design` options via the campaign JSON, `--bc2_set`, `--bc2_modality`, `--bc2_core`, `--bc2_metadata`, or the `--bc2_*` property flags.

Post-processing subcommands operate on completed campaigns and are not wrapped. Execute them directly on a batch directory in `results/bc2/batches/` using the BindCraft2 container:

| Subcommand | Purpose |
|---|---|
| `bindcraft rank <campaign>` | Re-rank a campaign on a different metric |
| `bindcraft filter <campaign> --where METRIC>=VALUE` | Re-filter accepted designs without re-running |
| `bindcraft score <design.cif>` | Score one structure |
| `bindcraft archive` / `unarchive` | Pack or unpack trajectory folders for copying |

## Post-completion crashes

BindCraft2 sometimes exits 139 (SIGSEGV) or 134 (SIGABRT) while JAX/CUDA is shutting down, after it has already printed `campaign done:` and written the batch. The task treats that as success when `${batch_id}/summary.csv` and `${batch_id}/campaign_metadata.json` are present and non-empty, and logs a warning. A crash before that point, or any other non-zero status, still fails the task.

## Examples

The [examples/](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/tree/main/examples) directory contains complete working examples:

- `examples/pdl1-bc2`: PD-L1 binder design from a local target PDB, in both params and config mode
- `examples/pdl1-bc2-mini`: a four-trajectory smoke test, single- or dual-GPU
- `examples/bindcraft2-vhh`: humanised VHH (nanobody) design against PD-L1, including a config-mode run using the built-in `hPDL1` target preset

See the [examples/README.md](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/blob/main/examples/README.md) for details.
