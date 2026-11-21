# VHH (nanobody) binder design with BindCraft2

These scripts design humanised single-domain antibody (VHH) binders against human PD-L1 using `--method bc2 --bc2_modality VHH`. They share the target from `examples/pdl1-bc2` for direct comparison with de novo campaigns.

| Script | What it does |
|---|---|
| `run.sh` | Params mode on a single GPU: `--input_pdb` plus `--bc2_modality VHH --bc2_humanize` |
| `run-dual-gpu.sh` | The same campaign spread across two local GPUs |
| `run-config.sh` | Config mode with BindCraft2's built-in `hPDL1` target (no structure file to stage) |

## VHH modality settings

`--bc2_modality VHH` switches BindCraft2 from de novo backbone generation to fold conditioning on the built-in VHH scaffold inside the container (`settings/modality/VHH.json`). The container image includes the scaffold and CDR mutation mask, so no external scaffold file is needed.

Key effects of the VHH preset:

- `--binder_length_range` is ignored: BindCraft2 gets the binder length from the scaffold and only varies CDR regions defined in the mask. Params mode still writes `binder_lengths` to the generated JSON, and BindCraft2 ignores it.
- Each trajectory takes longer than de novo design because the preset sets `design_models` to 3 and `validation_model` to `multimer`.
- Acceptance rates differ from de novo runs at the same trajectory count. Four trajectories only tests execution; increase `--bc2_n_traj` for production runs.

`--bc2_humanize` adds human-like sequence preferences and an MHC-II anchor-score objective. Like the other eight `--bc2_*` property flags, it combines with `--bc2_modality`.

## Built-in target presets

`run-config.sh` sets `"target": "hPDL1"` in `campaign.json` to load the target structure, chains and hotspots directly from the container:

```json
{
  "campaign_name": "hPDL1_VHH",
  "target": "hPDL1",
  "modality": "VHH",
  "humanize": true
}
```

Because modality and design properties are standard JSON keys, a campaign file supports the same configuration as CLI flags. Available built-in targets: `hPDL1`, `mPDL1`, `hPD1`, `hIL2R`, `hIL7RA` and `dynorphin_a`.

## Output

Outputs are written to `results/bc2/` with the same layout as `examples/pdl1-bc2`:

- `batches/PDL1_0/` ...: per-batch BindCraft2 campaign folders
- `ranked/`: flat folder of accepted mmCIF structures from all batches
- `summary.csv`, `ranked.csv`, `refolded.csv`, `trajectories.csv`
- `bc2_report.html`: HTML summary report

See the [BindCraft2 workflow docs](https://australian-protein-design-initiative.github.io/nf-binder-design/workflows/bc2/)
for the full option list.
