# nf-binder-design

[![DOI](https://zenodo.org/badge/921968505.svg)](https://doi.org/10.5281/zenodo.16809704) | [![Documentation](https://img.shields.io/badge/docs-online-blue)](https://australian-protein-design-initiative.github.io/nf-binder-design/)

Nextflow pipelines for de novo protein binder design.

![RFdiffusion workflow](docs/docs/images/rfd-workflow.png)

- RFdiffusion → ProteinMPNN → AlphaFold2 initial guess → Boltz-2 refolding
- RFdiffusion Partial Diffusion → Boltz-2 refolding
- BindCraft (in parallel across multiple GPUs)
- "Boltz Pulldown" (an AlphaPulldown-like protocol using Boltz-2)

> ⚠️ Note: Components of these workflows use RFdiffusion and BindCraft, which depend on PyRosetta/Rosetta, which is free for non-commercial use. Commercial use requires a paid license agreement with University of Washington: https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md and https://rosettacommons.org/software/licensing-faq/

----

Documentation: https://australian-protein-design-initiative.github.io/nf-binder-design/

- [nf-binder-design](#nf-binder-design)
  - [Setup](#setup)
  - [Examples](#examples)
    - [Commandline options](#commandline-options)
  - [Binder design with RFdiffusion](#binder-design-with-rfdiffusion)
    - [Single node or local workstation](#single-node-or-local-workstation)
    - [Parallel on an HPC cluster](#parallel-on-an-hpc-cluster)
  - [Partial diffusion on binder designs](#partial-diffusion-on-binder-designs)
  - [Binder design with BindCraft](#binder-design-with-bindcraft)
  - [Boltz pulldown](#boltz-pulldown)
  - [Utility scripts](#utility-scripts)
  - [Design filter plugin system](#design-filter-plugin-system)
    - [Plugin API](#plugin-api)
  - [License](#license)


## Setup

Install [Nextflow](https://www.nextflow.io/docs/latest/install.html).

Clone the git repository:

```bash
git clone https://github.com/Australian-Protein-Design-Initiative/nf-binder-design
```

## Examples

See the [examples](examples/) directory for examples.

### Commandline options

For any of the workflows, you can see the commandline options with `--help`, eg:

```bash
nextflow run main.nf --help
```

Any of the `--params` commandline options can alternatively be defined in a `params.json` file and passed to the workflow with `-params-file params.json`.

## Binder design with RFdiffusion

### Single node or local workstation

Simple example (single 'local' compute node):

```bash
OUTDIR=results
mkdir -p $OUTDIR/logs

nextflow run main.nf \
    --input_pdb target.pdb \
    --outdir $OUTDIR \
    --contigs "[A371-508/A753-883/A946-1118/A1135-1153/0 70-100]" \
    --hotspot_res "A473,A995,A411,A421" \
    --rfd_n_designs=10 \
    --rfdiffusion_batch_size 1 \
    -with-report $OUTDIR/logs/report_$(date +%Y%m%d_%H%M%S).html \
    -with-trace $OUTDIR/logs/trace_$(date +%Y%m%d_%H%M%S).txt \
    -resume \
    -profile local
```

> If you are working on a specific HPC cluster like M3 or MLeRP, you should omit `-profile local` and add the `-c` flag pointing to the specific platform config, eg `-c conf/platforms/m3.config` for M3.

### Parallel on an HPC cluster

A more complex example, as a wrapper script for the M3 HPC cluster, using a the site-specific config (`-c`), a specific RFdiffusion model (`--rfd_model_path`), a radius of gyration filter on the generated RFdiffusion backbones (`--rfd_filters`), custom ProteinMPNN weights (`--pmpnn_weigths`) and radius of gyration potentials (`--rfd_extra_args`):

```bash
#!/bin/bash
# CHANGE THIS - this is the path where your git clone of this repo is
WF_PATH="/some/path/to/nf-binder-design"

mkdir -p results/logs
DATESTAMP=$(date +%Y%m%d_%H%M%S)

# Ensure our tmp directory is in a location with enough space
export TMPDIR=$(realpath ./tmp)
export NXF_TEMP=$TMPDIR
mkdir -p $TMPDIR

# CHANGE THIS to a path in scratch or scratch2 to act as the cache directory for apptainer
# Containers will be automatically downloaded to this path.
# You can add it to ~/.bashrc if you prefer
export NXF_APPTAINER_CACHEDIR=/some/path/to/scratch2/apptainer_cache
export NXF_APPTAINER_TMPDIR=$TMPDIR

# There's a module for Nextflow on M3
module load nextflow/24.04.3 || true

# CHANGE the --slurm_account to match the project ID you wish to run SLURM jobs under
nextflow \
-c ${WF_PATH}/conf/platforms/m3.config run \
${WF_PATH}/main.nf  \
--slurm_account=ab12 \
--input_pdb 'input/target_cropped.pdb' \
--design_name my-binder \
--outdir results \
--contigs "[B346-521/B601-696/B786-856/0 70-130]" \
--hotspot_res "B472,B476,B484,B488" \
--rfd_n_designs=1000 \
--rfd_batch_size=5 \
--rfd_filters="rg<20" \
--pmpnn_seqs_per_struct=2 \
--pmpnn_relax_cycles=1 \
--pmpnn_weigths="/models/HyperMPNN/retrained_models/v48_020_epoch300_hyper.pt" \
--rfd_model_path="/models/rfdiffusion/Complex_beta_ckpt.pt" \
--rfd_extra_args='potentials.guiding_potentials=[\"type:binder_ROG,weight:7,min_dist:10\"] potentials.guide_decay="quadratic"' \
-resume \
-with-report results/logs/report_${DATESTAMP}.html \
-with-trace results/logs/trace_${DATESTAMP}.txt
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
    --hotspot_res "A473,A995,A411,A421" \
    --rfd_partial_T=2,5,10,20 \
    -with-report $OUTDIR/logs/report_$(date +%Y%m%d_%H%M%S).html \
    -with-trace $OUTDIR/logs/trace_$(date +%Y%m%d_%H%M%S).txt \
    -profile local
```

## Binder design with BindCraft

![BindCraft workflow](docs/docs/images/bindcraft-workflow.png)

The `bindcraft.nf` helps run [BindCraft](https://github.com/martinpacesa/BindCraft) trajectories in parallel across multiple GPUs.
This is particularly well suited for running BindCraft on an HPC cluster, or a workstation with multiple GPUs.

Unlike the default BindCraft configuration which runs for an indeterminate amount of time until a number of accepted designs are found,
this pipeline will run a fixed number of trajectories `--bindcraft_n_traj` and stop.

Example:

```bash
DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run bindcraft.nf  \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56,A125" \
  --hotspot_subsample 0.5 \
  --binder_length_range "55-120" \
  --bindcraft_n_traj 2 \
  --bindcraft_batch_size 1 \
  --bindcraft_advanced_settings_preset "default_4stage_multimer" \
  --bindcraft_filters_preset "default_filters" \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
```

`--bindcraft_advanced_settings_preset` and `--bindcraft_filters_preset` are are those available in the BindCraft [settings_advanced](https://github.com/martinpacesa/BindCraft/tree/main/settings_advanced) and [settings_filters](https://github.com/martinpacesa/BindCraft/tree/main/settings_filters) directories (without the .json extension).

`--hotspot_subsample` randomly takes this random proportion of the hotspot residues for each design, allowing the impact of hotspot selection to be explored in a single run.

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
│   ├── bindcraft_report.html
│   ├── failure_csv.csv
│   ├── final_design_stats.csv
│   ├── mpnn_design_stats.csv
│   └── trajectory_stats.csv
└── logs
    ├── report_20250725_084959.html
    ├── trace_20250725_084959.txt
```

A report summarizing the results is generated in `bindcraft_report.html`.

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

## Utility scripts

The `bin/` directory contains utility scripts, most of which are used internally by the workflows, but can also be run as standalone scripts. Run `uv run bin/somescript.py --help` for usage _(use [uv](https://docs.astral.sh/uv/getting-started/installation/) to run the scripts and dependencies with be automatically dealt with)_

The `bin/af2_combine_scores.py` can be useful for the RFdiffusion pipelines to monitor mid-run how things are going:

```bash
OUTDIR=results

uv run bin/af2_combine_scores.py -o $OUTDIR/combined_scores.tsv -p $OUTDIR/af2_results
```

## Design filter plugin system

The `main.nf` and `partial.nf` pipelines support a plugin system for calculating and filtering on custom metrics for designs. This is currently controlled by the `--rfd_filters` parameter.

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

## License

MIT

> Note that some software dependencies of the pipeline are under less permissive licenses - in particular, RFdiffusion and BindCraft use [Rosetta/PyRosetta](https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md) which is **only free for Non-Commercial use**.
