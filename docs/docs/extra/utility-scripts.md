# Utility Scripts

The `bin/` directory contains utility scripts that are used internally by the workflows but can also be run as standalone tools.

## Running Scripts

Run with [uv](https://docs.astral.sh/uv/getting-started/installation/) for automatic dependency management - each script has a `--help` option to show usage information:

```bash
uv run bin/somescript.py --help
```

## Available Scripts

### af2_combine_scores.py

Combines AlphaFold2 scores from multiple predictions into a single summary table.
This can be useful to monitor mid-run progress for RFdiffusion pipelines (it's run automatically at the end of the pipeline).

**Usage:**
```bash
OUTDIR=results
uv run bin/af2_combine_scores.py -o $OUTDIR/combined_scores.tsv -p $OUTDIR/af2_results
head $OUTDIR/combined_scores.tsv
```

### calculate_shape_scores.py

Calculates shape-based scoring metrics for protein designs (eg radius of gyration, etc.).

### create_bindcraft_settings.py

Generates BindCraft configuration files.

### create_boltz_yaml.py

Creates YAML configuration files for Boltz predictions.

### filter_designs.py

Main script for the design filter plugin system. Automatically discovers and calls filter plugins based on filter expressions. Required filter plugins in `filters.d/`. Could be used for post-pipeline filtering, but is largely intended for internal use.

### get_contigs.py

Extracts 'contig' information from protein structures in RFdiffusion syntax - useful for determining the contig ranges from a hand-cropped structure.

### merge_scores.py

Merges scoring tables from multiple sources.

### pdb_to_fasta.py

Extracts the FASTA sequence of chains in a PDB files.

### renumber_chains.py

Renumbers a chain in a PDB file.

### trim_to_contigs.py

Trims protein structures to specified contig regions, using the RFdiffusion contig syntax.

### bindcraft_scoring.py

Design scoring code extracted from BindCraft - outputs a subset of BindCraft's scores for any set of designs. Required PyRosetta - should be run using the BindCraft container.


## Filter Plugins

Custom filter plugins are located in `bin/filters.d/`. Any `*.py` file in this directory will be automatically discovered.

### Available Filters

- **rg** (radius of gyration) - in `bin/filters.d/rg.py`

### Creating Custom Filters

Create a new `.py` file in `bin/filters.d/` implementing two functions:

#### 1. `register_metrics() -> list[str]`

Returns list of metric names:

```python
def register_metrics() -> list[str]:
    return ["rg", "my_custom_score"]
```

#### 2. `calculate_metrics(pdb_files: list[str], binder_chains: list[str]) -> pd.DataFrame`

Calculates metrics and returns a DataFrame:

```python
def calculate_metrics(pdb_files: list[str], binder_chains: list[str]) -> pd.DataFrame:
    # Perform calculations
    # Return DataFrame with:
    #   - Index: design ID (PDB filename without .pdb)
    #   - Columns: metric names from register_metrics()
    return results_df
```
