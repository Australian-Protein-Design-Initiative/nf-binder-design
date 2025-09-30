# Boltz Pulldown

An experimental [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown)-like protocol using [Boltz](https://github.com/jwohlwend/boltz) for multimer predictions.

## Overview

Boltz Pulldown runs multimer predictions from sequence for all target-binder pairs, allowing you to assess binding interactions between multiple potential interaciton partners.

## Command-line Options

See available options with `--help`:

```bash
nextflow run boltz_pulldown.nf --help
```

## Multiple Sequence Alignments (MSAs)

By default, **no MSAs are used**. To enable:

- `--create_target_msa=true`: Generate MSAs for targets
- `--create_binder_msa=true`: Generate MSAs for binders

> ⚠️ For _de novo_ designed binders, use `--create_target_msa=true` to improve _target_ prediction accuracy, but skip generating _binder_ MSAs (`--create_binder_msa=false`, the default) as these are less informative since no homologs exist.

## MSA Generation Options

### Option 1: Remote Server

```bash
--use_msa_server
```

Uses the remote ColabFold mmseqs2 server.

### Option 2: Local Databases

Download and set up databases with the ColabFold [setup_databases.sh](https://github.com/sokrypton/ColabFold/blob/main/setup_databases.sh) script, then specify paths:

```bash
--uniref30 /path/to/uniref30
--colabfold_envdb /path/to/colabfold_envdb
```

## Custom Boltz Arguments

Add extra Boltz command-line arguments via `ext.args` in `nextflow.config`:

```groovy
process {
    withName: BOLTZ {
        accelerator = 1
        time = 2.hours
        memory = '8g'
        cpus = 2
        
        // Example: CPU-only mode
        ext.args = "--accelerator cpu"
    }
}
```

## Output

Default output directory: `results/boltz_pulldown/`

Includes:

- Boltz predictions with scores and structures
- `boltz_pulldown.tsv` - summary table
- `boltz_pulldown_report.html` - report with statistics and ipTM score plots
