# BindCraft2 Workflow Reference

## Overview

The `--method bc2` workflow (alias `--method bindcraft2`) runs [BindCraft2](https://github.com/PacesaLab/BindCraft2) campaigns in parallel batches, then gathers them with `bindcraft campaign_output`.

BindCraft2 (BC2) is a rewrite of BindCraft, not a version bump. It uses a different CLI, settings format and output layout, so `--method bindcraft` and `--method bc2` are separate workflows and only share common input flags.

Additions over `--method bindcraft`:

- One campaign JSON instead of three settings files, with all ~235 settings configurable via `--bc2_set`.
- Twelve binder modalities: `binder`, `large_binder`, `peptide`, `cyclic_peptide`, `homo_oligomer`, `multidomain`, `VHH`, `scFv`, `Fab`, `ARP`, `induced_fit`, `fold_switch`.
- Nine composable design properties (humanisation, protease stability, termini constraints, ...).
- Multi-target campaigns and off-target detargeting.
- PDB, mmCIF or FASTA targets; structures are written as **mmCIF**, not PDB.
- Candidates ranked by `i_pDAE`.

Like `--method bindcraft`, `--method bc2` runs a **fixed trajectory budget** and exits instead of looping until N designs pass filtering. Start with a small budget to measure the acceptance rate, then scale `--bc2_n_traj`.

## Input modes

Mutually exclusive; specify one.

### Config mode (`--bc2_config campaign.json`)

Provide a BindCraft2 campaign settings JSON. Only `targets[].target_path` and `binder_scaffold` are file paths; the pipeline resolves both relative to the campaign JSON directory, stages them into every task, and rewrites the JSON entries to basenames. Relative paths work from any launch directory.

```json
{
  "campaign_name": "PDL1_bc2",
  "targets": [
    {"name": "PDL1", "target_path": "input/PDL1.pdb", "chains": "A", "hotspots": "A56"}
  ],
  "binder_lengths": [55, 120]
}
```

Built-in presets (such as `{"target": "hPDL1"}`) resolve their files inside the container and require no file staging.

### Params mode (`--input_pdb`)

Generates `campaign.json` from flat CLI flags using the same input parameter names as `--method bindcraft`.

## Parameters

Use `--method bc2 --help` for the definitive list.

### Input / Output

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--bc2_config` | Campaign settings JSON. Mutually exclusive with `--input_pdb` | - |
| `--input_pdb` | Target PDB/mmCIF/FASTA; a directory or quoted glob runs each structure | - |
| `--outdir` | Output directory | `results` |

### Params mode only

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--target_chains` | Target chain IDs, comma-separated | - |
| `--hotspot_res` | Hotspots, e.g. `A54,A56` or `B12-16`. Chain prefix optional, ranges allowed | none |
| `--binder_length_range` | `60-150` range, `80` exact, or `60,80,100` discrete choices | `60-150` |
| `--design_name` | `campaign_name` for the generated config | `bc2_design` |

BC2 hotspot syntax is more permissive than BindCraft v1's `{chain}{resnum}`: ranges and unprefixed residues are accepted.

### Both modes

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--hotspot_subsample` | Fraction of each target's hotspots to keep per batch, `ceil(N*fraction)` with a minimum of 1 | `1.0` |
| `--do_foldseek_search` | FoldSeek structural novelty search over the ranked designs | `false` |

`--hotspot_subsample` runs while staging the campaign JSON into each task, so it applies in both params and config mode. Each batch samples hotspots using its `campaign_seed` instead of the system clock. Built-in `"target"` presets are unaffected because BindCraft2 merges preset hotspots after staging.

`--do_foldseek_search` queries only binder chains. BindCraft2 writes target chains first and binder chains last. The pipeline reads `binder_chains` from `campaign_metadata.json` and extracts the trailing binder chains. The default smallest-chain heuristic fails for multi-chain binders (Fab, scFv, multidomain, symmetric oligomer) or targets smaller than the binder. FoldSeek writes outputs to `${outdir}/foldseek/{database}/` without merging them into `ranked.csv` or the HTML report. All other `--foldseek_*` options match the other workflows.

### Batching

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--bc2_n_traj` | Total trajectories to run | `10` |
| `--bc2_batch_size` | Trajectories per parallel task | `1` |
| `--bc2_seed` | Base `campaign_seed`. Batch seed is `bc2_seed + seedIndex`, where `seedIndex` is the batch index in config mode, or `pdbIndex * nBatches + batchIndex` when several structures are folded in params mode | `0` |

### BindCraft2 settings

| Parameter | Description |
|-----------|-------------|
| `--bc2_set` | Any setting as `KEY=VALUE`, `;`-separated |
| `--bc2_modality` | Binder format, e.g. `VHH`, or combined: `VHH,induced_fit` |
| `--bc2_core` | Core profile, e.g. `benchmark` |
| `--bc2_metadata` | JSON file of fields to record with the campaign |
| `--bc2_extra_args` | Extra arguments appended to `bindcraft design` |

### Design properties (boolean flags)

`--bc2_humanize`, `--bc2_bigbang`, `--bc2_disulfide_staple`, `--bc2_forced_targeting`, `--bc2_initial_guess`, `--bc2_mixed_topology`, `--bc2_protease_stable`, `--bc2_termini_accessible`, `--bc2_termini_together`.

Each corresponds to a `settings/property/*.json` upstream and is composable with the others and with `--bc2_modality`.

## `--bc2_set`

Pass any BindCraft2 setting via `--bc2_set`, which forwards `KEY=VALUE` pairs to `bindcraft design --set`:

```bash
--bc2_set 'binder_lengths=[70,90];filters.i_pTM.threshold=0.8;design_models=3'
```

- Use `;` between assignments because JSON values often contain commas.
- CLI runs accept a single semicolon-separated string (repeated `--bc2_set` flags overwrite earlier ones); `-params-file` and config files also accept a YAML/Groovy list.
- Dotted keys update nested fields (`filters.i_pTM.threshold=0.8` leaves sibling filter settings intact).
- Values parse as JSON, falling back to plain strings.
- Full setting list: [`settings/core/reference.json`](https://github.com/PacesaLab/BindCraft2/blob/main/settings/core/reference.json).

The driver rejects `project_folder`, `campaign_seed` and `max_trajectories` before launching GPU tasks, pointing to `--bc2_seed`, `--bc2_n_traj` and `--bc2_batch_size`.

Files passed via `--bc2_set` are not staged; BindCraft2 resolves those paths relative to the task work directory instead of the config directory. Put file inputs in the campaign JSON.

## Seeding and reproducibility

BindCraft2 derives each trajectory RNG from `fold_in(PRNGKey(campaign_seed), trajectory_number)`, numbering trajectories from 1 in every batch. To avoid duplicate outputs across batches, each batch receives `--bc2_seed + seedIndex` (`seedIndex` is the batch index in config mode, or `pdbIndex * nBatches + batchIndex` when params mode folds several structures).

Reproducing a run requires keeping `--bc2_seed`, `--bc2_n_traj` and `--bc2_batch_size` constant. Changing `--bc2_batch_size` alters batch seeds and resulting designs.

Verify per-batch seeds after a run:

```bash
grep -h campaign_seed results/bc2/batches/*/campaign_metadata.json
```

Duplicate seeds across batches or repeated filenames in `3_Ranked/` indicate a seeding misconfiguration.

Pass `--bc2_core benchmark` to fix BindCraft2's internal seed and disable `autotune` and `desperation`.

## Adaptive tuning under batching

BindCraft2 has two adaptive mechanisms that inspect the trajectory table inside a single campaign's `project_folder`:

| Mechanism | Engages after | Tunable |
|---|---|---|
| `autotune` | every 10 trajectories in one campaign | No, a module constant |
| `desperation` | 750 fruitless trajectories, then one escalation rung per 50 | Yes, `desperation_trajectories` |

At the default `--bc2_batch_size` of 1, both mechanisms stay inactive. Even `--bc2_batch_size 50` stays below the 750-trajectory desperation threshold. The pipeline logs inactive mechanisms with the batch plan.

Keeping settings static across trajectories matches `--bc2_core benchmark` and supports reproducible fixed-budget comparisons. On difficult targets, however, long standalone BindCraft2 runs enable `initial_guess`, `target_flexibility`, multimer validation and extra recycles after 750 failed trajectories.

Increase `--bc2_batch_size` to let adaptive tuning activate, or apply specific desperation settings directly via `--bc2_set`:

```bash
--bc2_set 'validation_model=multimer;initial_guess=true;target_flexibility=0.5'
```

## No accepted-design quota

There is no `--bc2_n_designs` flag. Active jobs cannot stop early once a global design count is met, so a quota across M batches makes M separate campaigns each wait for their own full quota. On a hard target without a trajectory cap, a batch times out before publishing any output.

If you want an early stop within each batch, pass `--bc2_set 'number_of_final_designs=3'`.

## Output Structure

```
results/bc2/
├── batches/
│   ├── PDL1_0/
│   │   ├── 1_Trajectories/!_Trajectories.csv
│   │   ├── 2_Refolded/!_Refolded.csv
│   │   ├── 3_Ranked/!_Ranked.csv  + *.cif
│   │   ├── campaign_metadata.json
│   │   ├── summary.csv
│   │   ├── bc2.log
│   │   └── bc2_campaign.staged.json   # the path-rewritten config BC2 was given
│   └── PDL1_1/ ...
├── ranked/                        # accepted mmCIFs from every batch, flattened
├── summary.csv                    # from `bindcraft campaign_output`
├── ranked.csv                     # stage tables concatenated across batches,
├── refolded.csv                   #   each with a `campaign` column
├── trajectories.csv
└── bc2_report.html
```

`bindcraft campaign_output` updates per-batch rankings and writes the combined `summary.csv`. The pipeline concatenates the stage CSVs across batches separately.

## Processes

| Process | Role |
|---|---|
| `BC2_GENERATE_CONFIG` | Params mode only; builds `campaign.json`. Runs locally |
| `BINDCRAFT2` | One GPU task per batch; runs `bindcraft design` |
| `BC2_CAMPAIGN_OUTPUT` | Stages every batch under one parent, runs `bindcraft campaign_output`, concatenates stage tables |
| `BC2_REPORTING` | Renders `assets/bc2_reporting.qmd` to `bc2_report.html` |

## Containers and weights

`ghcr.io/australian-protein-design-initiative/containers/bindcraft2:v1.0.1`, mirrored for Monash as `ghcr.io-australian-protein-design-initiative-containers-bindcraft2-v1.0.1.img`.

AlphaFold2 weights (5.2 GB) are bundled in the container image at `/opt/bindcraft/bindcraft/weights/alphafold`, avoiding runtime downloads. Use `--bc2_af2_params` only with custom images built without bundled weights.

## Examples

- `examples/pdl1-bc2/`: params mode, config mode, a humanised VHH campaign, and a dual-GPU run
- `examples/pdl1-bc2-mini/`: four-trajectory smoke test (single- or dual-GPU)
- `examples/bindcraft2-vhh/`: humanised VHH against PD-L1, including config mode with the built-in `hPDL1` target preset
