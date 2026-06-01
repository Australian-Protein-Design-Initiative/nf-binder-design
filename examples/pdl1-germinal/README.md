PD-L1 VHH (nanobody) design with [Germinal](https://github.com/SantiagoMille/germinal) and Protenix structure prediction.

This example uses a combined Hydra config (`configs/pdl1_vhh.yaml`) for a PDL1 nanobody design with AbLang and Protenix. The target structure is in `pdbs/pdl1.pdb`; the nanobody scaffold (`nb.pdb`) is supplied from the Germinal container if not present locally.

Trajectory count is controlled by Nextflow (`--germinal_n_traj`), not the values in the YAML config file.

To run locally on a single GPU:

```bash
./run.sh
```

To run locally on a dual-GPU machine:

```bash
./run-dual-gpu.sh
```

To run on M3 BDI (SLURM):

```bash
./run-m3-bdi.sh
```

The BDI GPU partition (`bdi`/`bdiq`) requires a SLURM account with QOS access. The script defaults to `yt41`; override with `SLURM_ACCOUNT=your_account ./run-m3-bdi.sh`.

Results are written under `results/germinal/`, with merged CSVs at the top level (e.g. `all_trajectories.csv`, `accepted_designs.csv`, `failure_counts.csv`), the resolved Hydra config in `results/germinal/config/final_config.yaml`, and accepted structures in `results/germinal/accepted/structures/`.
