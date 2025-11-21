# BoltzGen Workflow

![BoltzGen workflow](../images/boltzgen-workflow.png)

The `boltzgen.nf` workflow automates the design of binders using the [BoltzGen](https://github.com/HannesStark/boltzgen) generative model.
The key advantage of this workflow is improved parallelization across multiple GPUs and better resource allocation on HPC clusters.

It supports the `protein-anything`, `peptide-anything`, `protein_small-molecule` and `nanobody-anything` protocols.

## Overview

The BoltzGen workflow performs the following steps, as per `boltzgen run`:

-  **Design**: Generates binder backbone structures based on a YAML configuration.
-  **Inverse Folding**: Generates sequences for the designed backbones.
-  **Folding**: Re-folds the designed binder in complex with the target.
-  **Design Folding**: Re-folds the designed binder alone.
-  **Affinity**: (Optional -  `protein_small-molecule` protocol only) Calculates affinity scores for protein-small molecule complexes.
- (**Merge batches**: nf-binder-design `bolzgen.nf` merges independently generated batches of designs, ready for analysis and filtering)
-  **Analysis & Filtering**: Aggregates scores and filters designs based on user-defined criteria.

## General Information

### Command-line Options

You can see available options with `--help`:

```bash
nextflow run boltzgen.nf --help
```

### YAML Configuration

The BoltzGen workflow is controlled by a YAML configuration file passed via `--config_yaml`. This file defines the input structure, constraints, and design objectives. Refer to the [BoltzGen documentation](https://github.com/HannesStark/boltzgen?tab=readme-ov-file#how-to-make-a-design-specification-yaml) on how to prepare your configuration file.

> **Note:** ⚠️ BoltzGen expects residues specified as indices starting at 1, irresepctive of the numbering in your PDB/mmCIF file. You may find it easier to use `bin/renumber_chains.py` to renumber your input PDB to be sequentially numbered starting at 1 to simplify choosing hotspot residues in ChimeraX/Pymol/Molstar etc.

## Usage

### Basic Execution

To run the workflow with a configuration file:

```bash
nextflow run boltzgen.nf --config_yaml config/my_design.yaml
```

### Key Parameters

- `--config_yaml`: **(Required)** Path to the BoltzGen YAML configuration file.
- `--outdir`: Output directory for results (default: `results`).
- `--design_name`: Name of the design, used for output file prefixes. Defaults to the basename of the config file.
- `--protocol`: Protocol type. Options:
    - `protein-anything` (default)
    - `peptide-anything`
    - `protein-small_molecule`
    - `nanobody-anything`
- `--num_designs`: Total number of designs to generate (default: 100).
- `--inverse_fold_num_sequences`: Number of sequences to generate per backbone during the inverse folding step  (default: 1).
- `--batch_size`: Number of designs to process per batch (default: 10).
- `--budget`: Final number of designs to keep after diversity optimization and filtering (default: 10).
- `--devices`: Number of GPU devices to use (default: 1).
- `--num_workers`: Number of DataLoader workers.

## Key Outputs

The `results` directory (or specified `--outdir`) will contain:

- `params.json`: A record of the parameters used for the run.
- `boltzgen/batches`: Each independent design batch, in folders for each step (`design`, `inverse_folding`, `folding`, `design_folding`)
- `bolzgen/merged`: All batches merged, after analysis and filtering. Final designs are in `final_ranked_designs`

`boltzgen/merged` should be equivalent to the output of a non-Nextflow execution of `boltzgen run` - see the [BoltzGen docs for details on the pipeline output](https://github.com/HannesStark/boltzgen?tab=readme-ov-file#pipeline-output).

## Examples

The [examples/](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/tree/main/examples) directory contains complete working examples for BoltzGen workflows:

- `examples/boltzgen-protein`: Protein binder design (protocol `protein-anything`)
- `examples/boltzgen-cyclic-peptide`: Cyclic peptide binder design (protocol `peptide-anything`)
- `examples/boltzgen-small-molecule-binder`: Small molecule binder design (protocol `protein-small_molecule`)

See the [examples/README.md](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/blob/main/examples/README.md) for details.
