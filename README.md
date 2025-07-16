# nf-binder-design

Nextflow pipelines for de novo protein binder design.

- RFDiffusion -> ProteinMPNN -> AlphaFold2 initial guess
- RFDiffusion Partial Diffusion
- BindCraft (in parallel across multiple GPUs)
- "Boltz Pulldown" (an AlphaPulldown-like protocol using Boltz-2)

> Note well: Components of these workflows use RFDiffusion and BindCraft, which depend on PyRosetta/Rosetta, which is free for non-commercial use, however commercial use requires a paid license agreement with University of Washington: https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md and https://rosettacommons.org/software/licensing-faq/

## Setup

Install [Nextflow](https://www.nextflow.io/docs/latest/install.html).

Clone the git repository:

```bash
git clone https://github.com/Australian-Protein-Design-Initiative/nf-binder-design
```

## Binder design with RFDiffusion

Example:

```bash
OUTDIR=results
mkdir -p $OUTDIR/logs

nextflow run main.nf \
    --input_pdb target.pdb \
    --outdir $OUTDIR \
    --contigs "[A371-508/A753-883/A946-1118/A1135-1153/0 70-100]" \
    --hotspot_res "[A473,A995,A411,A421]" \
    --rfd_n_designs=10 \
    --rfdiffusion_batch_size 1 \
    -with-report $OUTDIR/logs/report_$(date +%Y%m%d_%H%M%S).html \
    -with-trace $OUTDIR/logs/trace_$(date +%Y%m%d_%H%M%S).txt \
    -resume \
    -profile local
```

> If you are working on a specific cluster like M3 or MLeRP, you should omit `-profile local` and add the `-c` flag pointing to the specific platform config, eg `-c conf/platforms/m3.config` for M3.

A more complex example, as a wrapper script for the M3 HPC cluster, using a the site-specific config (`-c`), a specific RFDiffusion model (`--rfd_model_path`), a radius of gyration filter on the generated RFDiffusion backbones (`--rfd_filters`), custom ProteinMPNN weights (`--pmpnn_weigths`) and radius of gyration potentials (`--rfd_extra_args`):

```bash
#!/bin/bash
# CHANGE THIS - this is the path where your git clone of this repo is
WF_PATH="/some/path/to/nf-binder-design"

mkdir -p results/logs
DATESTAMP=$(date +%Y%m%d_%H%M%S)

# CHANGE THIS to a path in scratch or scratch2 to act as the cache directory for apptainer
# Containers will be automatically downloaded to this path.
# You can add it to ~/.bashrc if you prefer
export NXF_APPTAINER_CACHEDIR=/some/path/to/scratch2/apptainer_cache
export NXF_APPTAINER_TMPDIR=/some/path/to/scratch2/apptainer_cache/tmp

# CHANGE the --slurm_account to match the project ID you wish to run SLURM jobs under
nextflow \
-c ${WF_PATH}/conf/platforms/m3.config run \
${WF_PATH}/main.nf  \
--slurm_account=ab12 \
--input_pdb 'input/target_cropped.pdb' \
--design_name my-binder \
--outdir results \
--contigs "[B346-521/B601-696/B786-856/0 70-130]" \
--hotspot_res "[B472,B476,B484,B488]" \
--rfd_n_designs=1000 \
--rfd_batch_size=5 \
--rfd_filters="rg<20" \
--pmpnn_seqs_per_struct=2 \
--pmpnn_weigths="/models/HyperMPNN/retrained_models/v48_020_epoch300_hyper.pt" \
--rfd_model_path="/models/rfdiffusion/Complex_beta_ckpt.pt" \
--rfd_extra_args='potentials.guiding_potentials=[\"type:binder_ROG,weight:7,min_dist:10\"] potentials.guide_decay="quadratic"' \
-resume \
-with-report results/logs/report_${DATESTAMP}.html \
-with-trace results/logs/trace_${DATESTAMP}.txt
```

### Summarize af2_initial_guess scores

This happens by default when the pipeline successfully completes, however you may want to run it mid-run to monitor how things are going:

```bash
OUTDIR=results

bin/af2_combine_scores.py -o $OUTDIR/combined_scores.tsv -p $OUTDIR/af2_results
```

## Partial diffusion on binder designs

> NOTE: It seems with output from previous designs that the binder is always named chain A, and your other chains are named B, C, etc - irrespective of the chain ID in the original target PDB file. Residue numbering is 1 to N, sequential irrespective of gaps in the chain, rather than original target chain numbering.

```bash
OUTDIR=results
mkdir -p $OUTDIR/logs

# Generate 10 partial designs for each binder, in batches of 5
# Note the 'single quotes' around the '*.pdb' glob pattern !
nextflow run partial.nf  \
    --input_pdb 'my_designs/*.pdb' \
    --rfd_n_partial_per_binder=10 \
    --rfd_batch_size=5 \
    --hotspot_res "[A473,A995,A411,A421]" \
    --rfd_partial_T=2,5,10,20 \
    -with-report $OUTDIR/logs/report_$(date +%Y%m%d_%H%M%S).html \
    -with-trace $OUTDIR/logs/trace_$(date +%Y%m%d_%H%M%S).txt \
    -profile local
```

## Design filter plugin system

The pipeline supports a plugin system for calculating and filtering on custom metrics for designs. This is currently controlled by the `--rfd_filters` parameter.

Filters are simple Python scripts located in the `bin/filters.d/` directory. Any `*.py` file in this directory will be automatically discovered.

Currently, only a radius of gyration (`rg`) filter is implemented in `bin/filters.d/rg.py`, which can be used as a template for creating new filters.

### Plugin API

A filter plugin must implement two functions:

1.  `register_metrics() -> list[str]`
    This function should return a list of the metric names (as strings) that the plugin can calculate. For example: `return ["rg", "my_custom_score"]`.

2.  `calculate_metrics(pdb_files: list[str], binder_chains: list[str]) -> pd.DataFrame`
    This function takes a list of PDB file paths and a list of binder chain IDs. It should perform its calculations and return a `pandas.DataFrame` with the following structure:
    *   The **index** of the DataFrame must be the design ID (i.e., the PDB filename without the `.pdb` extension).
    *   The **columns** must correspond to the metric names returned by `register_metrics()`.

The main `bin/filter_designs.py` script will call these plugins as needed based on the filter expressions provided to the pipeline (e.g., `--rfd_filters "rg<20"`).

## Binder design with BindCraft

The `bindcraft.nf` helps run [BindCraft](https://github.com/martinpacesa/BindCraft) trajectories in parallel across multiple GPUs.
This is particularly well suited for running BindCraft on an HPC cluster, or a workstation with multiple GPUs.

Unlike the default BindCraft configuration which runs for an indeterminate amount of time until a number of accepted designs are found,
this pipeline will run a fixed number of trajectories `--bindcraft_n_designs` and stop.

Example:

```bash
DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run bindcraft.nf  \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56" \
  --binder_length_range "55-120" \
  --bindcraft_n_designs 2 \
  --bindcraft_batch_size 1 \
  --bindcraft_advanced_settings_preset "default_4stage_multimer" \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
```

If you have multiple GPUs per compute node, you can specify them with the `--gpu_devices` flag, eg `--gpu_devices=0,1`.

Results are saved to the `--outdir` directory, in the `bindcraft` subdirectory, with CSV outputs from each batch combined into single tables, eg `bindcraft/final_design_stats.csv`, eg:

```
── bindcraft
│   ├── accepted
│   │   └── results
│   │       └── Accepted
│   │           ├── bindcraft_design_1_l57_s942028_mpnn6_model1.pdb
│   │           └── bindcraft_design_1_l57_s942028_mpnn8_model2.pdb
│   ├── batches
│   │   ├── 0
│   │   │   └── results
│   │   │       ├── failure_csv.csv
│   │   │       ├── final_design_stats.csv
│   │   │       ├── mpnn_design_stats.csv
│   │   │       ├── Trajectory
│   │   │       └── trajectory_stats.csv
│   │   └── 1
│   │       └── results
│   │           ├── Accepted
│   │           ├── failure_csv.csv
│   │           ├── final_design_stats.csv
│   │           ├── MPNN
│   │           ├── mpnn_design_stats.csv
│   │           ├── Rejected
│   │           ├── Trajectory
│   │           └── trajectory_stats.csv
│   ├── failure_stats_summed.csv
│   ├── final_design_stats.csv
│   ├── mpnn_design_stats.csv
│   └── trajectory_stats.csv
└── logs
    ├── report_20250725_084959.html
    ├── trace_20250725_084959.txt
```

## Boltz pulldown

An [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown)-like protocol, using [Boltz](https://github.com/jwohlwend/boltz). This is essentially running multimer predictions for all sequences in the target set (`--targets targets.fasta`) against all sequences in the binder set (`--binders binders.fasta`).

By default, no multiple sequence alignments are used. If you set the `--create_target_msa=true` and `--create_binder_msa=true` flags, target and binder multiple sequence alignments will be generated, respectively.

For _de novo_ designed binders a typical pattern might be to use `--create_target_msa=true` to allow the target to use an MSA to improve prediction accuracy, but not attempt to find homologs for the _de novo_ designed partner.

Specifying `--use_msa_server` will use the remote ColabFold mmseq2 server to generate multiple sequence alignments.

If generating the MSAs locally, indexed mmseqs2 databases can be downloaded and generated with the ColabFold [setup_databases.sh](https://github.com/sokrypton/ColabFold/blob/main/setup_databases.sh) script. You'll need to specify the `--uniref30` and `--colabfold_envdb` paths to point to these databases.

You can add additional args to the Boltz command line via `ext.args` for the `BOLTZ` process in `nextflow.config`, eg:

```
process {
    withName: BOLTZ {
        accelerator = 1
        time = 2.hours
        memory = '8g'
        cpus = 2
        
        // if using CPU only
        ext.args = "--accelerator cpu"
    }
}
```

By default, `boltz_pulldown.nf` will output to `results/boltz_pulldown`, which includes the Boltz outputs with scores and predicted structures, and a summary table `boltz_pulldown.tsv`. It also outputs `boltz_pulldown_report.html` with summary statistics and plots of the binder/target ipTM scores.

## License

MIT

> Note that some software dependencies of the pipeline are under less permissive licenses - in particular, RFDiffusion and BindCraft use [Rosetta/PyRosetta](https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md) which is **only free for Non-Commercial use**.
