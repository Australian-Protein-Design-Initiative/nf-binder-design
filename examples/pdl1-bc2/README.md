# PD-L1 binder design with BindCraft2

These scripts run small [BindCraft2](https://github.com/PacesaLab/BindCraft2) campaigns against human PD-L1 with `--method bc2`. They use the same target as `examples/pdl1-bindcraft` for direct comparison between the two methods.

| Script | What it does |
|---|---|
| `run.sh` | Params mode: the campaign JSON is generated from `--input_pdb`, `--target_chains` and `--hotspot_res` |
| `run-config.sh` | Config mode: `campaign.json` drives the run, with one filter threshold adjusted via `--bc2_set` |
| `run-vhh.sh` | A humanised VHH (nanobody) campaign, via `--bc2_modality VHH --bc2_humanize` |
| `run-dual-gpu.sh` | The params-mode run spread across two local GPUs |

Setting `--bc2_n_traj 4` with `--bc2_batch_size 1` launches four parallel tasks with distinct `campaign_seed` values. Four trajectories only tests execution and may yield zero accepted designs; increase `--bc2_n_traj` for production runs.

## Output

Outputs are written to `results/bc2/`:

- `batches/PDL1_0/` ...: per-batch BindCraft2 campaign folders (`1_Trajectories/`, `2_Refolded/`, `3_Ranked/`, `summary.csv`, `campaign_metadata.json`)
- `ranked/`: flat folder of accepted mmCIF structures from all batches
- `summary.csv`: combined summary from `bindcraft campaign_output`
- `ranked.csv`, `refolded.csv`, `trajectories.csv`: per-stage tables concatenated across batches with a `campaign` column
- `bc2_report.html`: HTML summary report

See the [BindCraft2 workflow docs](https://australian-protein-design-initiative.github.io/nf-binder-design/workflows/bc2/)
for the full option list.
