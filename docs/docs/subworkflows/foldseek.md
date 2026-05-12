# FoldSeek Structural Search

Run [FoldSeek](https://github.com/steineggerlab/foldseek) structural similarity searches against protein structure databases directly within the pipeline.

Useful for identifying structurally similar proteins to your designed binders and annotating designs with known structural folds.

## Quick Start

```bash
# Search designs against CATH50 (database auto-downloaded)
nextflow run main.nf --method foldseek \
    --input_pdbs 'results/designs/*.pdb'

# Search mmCIF designs (plain or gzipped)
nextflow run main.nf --method foldseek \
    --input_pdbs 'results/designs/*.cif.gz'

# Search with all hits, HTML report, and gzip output
nextflow run main.nf --method foldseek \
    --input_pdbs 'results/designs/*.pdb' \
    --foldseek_maxaccept 0 \
    --outdir results/foldseek_run
```

## Usage as a Subworkflow

The `FOLDSEEK_SEARCH` subworkflow can be imported into other workflows to add structural search as a step:

```groovy
include { FOLDSEEK_SEARCH } from '../subworkflows/local/foldseek_search'

workflow {
    ch_pdbs = Channel.fromPath('designs/*.pdb')
    ch_meta = Channel.of([id: 'my_search'])

    FOLDSEEK_SEARCH(ch_pdbs, ch_meta)

    FOLDSEEK_SEARCH.out.tsv.subscribe { f ->
        log.info "Results: ${f}"
    }
}
```

### Subworkflow Outputs

| Channel | Description |
|---------|-------------|
| `tsv` | TSV alignment results (or `.tsv.gz` if gzip enabled) |
| `html` | Pretty-printed HTML alignment report |
| `all_output` | All output files collected |
| `db_name_out` | Name of the database searched |

## Command-line Options

```bash
nextflow run main.nf --method foldseek --help
```

### Required

| Flag | Description |
|------|-------------|
| `--input_pdbs` | Glob pattern or path to input structure files (PDB `.pdb`, mmCIF `.cif`; plain or `.gz` compressed) |

### Database

| Flag | Default | Description |
|------|---------|-------------|
| `--foldseek_database` | `CATH50` | Database to search against (see table below) |
| `--foldseek_databases_path` | *(auto-download)* | Path to local databases directory. Each database in a subdirectory: `{path}/{db_name}/`. If the database doesn't exist there, it is downloaded and published to this path for reuse. If unset, the database is downloaded to the Nextflow work directory (cached across runs). |
| `--foldseek_use_webserver` | `false` | Use the FoldSeek web API instead of local search |

### Search

| Flag | Default | Description |
|------|---------|-------------|
| `--foldseek_mode` | `3diaa` | FoldSeek search mode |
| `--foldseek_maxaccept` | `1` | Max accepted alignments per query. `1` = top hit only, `0` = unlimited |

### Output

| Flag | Default | Description |
|------|---------|-------------|
| `--foldseek_gzip_output` | `false` | Gzip the TSV output (`.tsv.gz`) |
| `--foldseek_include_html_output` | `false` | Include HTML alignment report (off by default — HTML can be 100s of MB for large result sets) |

### CATH Annotation

When using a database whose name starts with `CATH` (e.g. `CATH50`), FoldSeek results are **automatically** annotated with CATH hierarchy descriptions. No flag needed — it just works.

By default (`--foldseek_maxaccept 1`), only the top FoldSeek hit per design is returned, so the CATH annotation describes the most structurally similar known fold. This is useful for quickly surveying what structural class your designs resemble.

| Flag | Default | Description |
|------|---------|-------------|
| `--foldseek_cath_names_path` | *(auto-download)* | Path to a pre-downloaded `cath-names.txt` file. If unset, the file is downloaded from the CATH FTP server and cached in `--foldseek_databases_path`. |

## Available Databases

Databases available via `foldseek databases`:

| Name | Source | Remote API | Description |
|------|--------|------------|-------------|
| `CATH50` | [CATH](https://www.cath.info) | `cath50` | CATH domain structures at 50% sequence identity clustering. Good balance of coverage and speed. **Default.** |
| `PDB` | [RCSB PDB](https://www.rcsb.org) | `pdb100` (complex-aware) | Full Protein Data Bank (local); PDB100 (complex-aware, incl. multimers) via remote. |
| `Alphafold/UniProt` | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | — | All AlphaFold predicted structures. Very large. Only available locally. |
| `Alphafold/UniProt50` | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | `afdb50` | AlphaFold structures at 50% clustering. Reduced redundancy. |
| `Alphafold/UniProt50-minimal` | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | — | Minimal version. Only available locally. |
| `Alphafold/Proteome` | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | `afdb-proteome` | AlphaFold proteome-level structures. |
| `Alphafold/Swiss-Prot` | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | `afdb-swissprot` | AlphaFold structures for Swiss-Prot reviewed proteins. |
| `ESMAtlas30` | [ESM Atlas](https://esmatlas.com) | `mgnify_esm30` | ESM Metagenomic Atlas at 30% clustering. Metagenomic predicted structures. |
| `BFMD` | [BFMD](https://foldseek.steineggerlab.workers.dev/bfmd.version) | `bfmd` | BFMD database. |
| `BFVD` | [BFVD](https://bfvd.steineggerlab.workers.dev) | `BFVD` | BFVD database. |
| `ProstT5` | [HuggingFace](https://huggingface.co/Rostlab/ProstT5) | — | Protein language model for 3Di prediction. Not a search database — only available locally. |

> **Tip:** For binder design validation, `CATH50` is a good default — it covers known structural domains and is fast. Use `PDB` for comprehensive searches against all experimentally determined structures.

## Database Caching

FoldSeek databases can be large (CATH50 ~1GB, PDB ~10s of GB). Two caching strategies are available:

### Auto-download (default)

When `--foldseek_databases_path` is not set, the database is downloaded once per pipeline run and cached in the Nextflow work directory. Subsequent runs with `-resume` reuse the cached download.

### Persistent local path

Set `--foldseek_databases_path` to a directory on shared storage:

```bash
# First run: downloads and publishes to /data/foldseek/
nextflow run main.nf --method foldseek \
    --input_pdbs 'designs/*.pdb' \
    --foldseek_databases_path /data/foldseek

# Subsequent runs: reuses existing /data/foldseek/CATH50/
nextflow run main.nf --method foldseek \
    --input_pdbs 'designs/*.pdb' \
    --foldseek_databases_path /data/foldseek
```

The database directory is published via `publishDir`, ensuring compatibility with distributed execution backends (SLURM, AWS Batch, etc.) where the Nextflow work directory is shared but the launch filesystem is not.

## Output Format

### TSV (`foldseek_results.tsv`)

Tab-separated with column headers (`--format-mode 4`):

| Column | Description |
|--------|-------------|
| `query` | Query structure name |
| `target` | Target structure name |
| `fident` | Fractional sequence identity |
| `alnlen` | Alignment length |
| `mismatch` | Number of mismatches |
| `qstart` | Query alignment start |
| `qend` | Query alignment end |
| `qlen` | Query sequence length |
| `tstart` | Target alignment start |
| `tend` | Target alignment end |
| `tlen` | Target sequence length |
| `evalue` | E-value |
| `bits` | Bit score |
| `qtmscore` | Query TM-score |
| `ttmscore` | Target TM-score |
| `alntmscore` | Alignment TM-score |
| `rmsd` | RMSD of aligned residues |
| `lddt` | lDDT score |
| `prob` | Probability |

### HTML (`foldseek_results.html`)

Pretty-printed interactive alignment report (`--format-mode 3`). Can be very large for big result sets — set `--foldseek_include_html_output true` to enable.

### Annotated TSV (`foldseek_results_annotated.tsv`)

When `--foldseek_database CATH50` (or any database name starting with `CATH`) is used with local search, an additional annotated TSV is produced automatically with these extra columns:

| Column | Description |
|--------|-------------|
| `cath_class` | CATH Class level (e.g. "Mainly Alpha", "Mainly Beta", "Alpha Beta") |
| `cath_architecture` | CATH Architecture level (e.g. "Orthogonal Bundle", "Sandwich", "Alpha Horseshoe") |
| `cath_topology` | CATH Topology level (e.g. "Immunoglobulin-like", "Rossmann fold") |
| `cath_homologous_superfamily` | CATH Homologous Superfamily level (e.g. "Immunoglobulins", "UBA domain") |
| `cath_code` | Full CATH code extracted from the FoldSeek target (e.g. "1.25.40.10") |

The CATH code is extracted from the FoldSeek target name, which for CATH50 has the format:
`af_{ACCESSION}_{START}_{END}_{CATHCODE}` (e.g. `af_A0A096MJB5_3_128_1.25.40.10`).

For non-CATH databases, annotation columns are left empty (no error).

## Examples

### Top hit per design against CATH50 (with automatic annotation)

```bash
nextflow run main.nf --method foldseek \
    --input_pdbs 'results/designs/*.pdb' \
    --outdir results/foldseek_cath50
```

This produces both `foldseek_results.tsv` (raw) and `foldseek_results_annotated.tsv` (with CATH columns). Annotation is automatic for any database whose name starts with `CATH`.

### All hits against PDB with gzip

```bash
nextflow run main.nf --method foldseek \
    --input_pdbs 'results/designs/*.pdb' \
    --foldseek_database PDB \
    --foldseek_maxaccept 0 \
    --foldseek_gzip_output true \
    --foldseek_include_html_output false \
    --outdir results/foldseek_pdb
```

### Using a pre-existing database on shared storage

```bash
nextflow run main.nf --method foldseek \
    --input_pdbs 'results/designs/*.pdb' \
    --foldseek_databases_path /shared/foldseek_dbs \
    --foldseek_database CATH50 \
    --outdir results/foldseek
```

### Remote search via FoldSeek web API

```bash
nextflow run main.nf --method foldseek \
    --input_pdbs 'results/designs/*.pdb' \
    --foldseek_use_webserver true \
    --foldseek_database pdb100 \
    --outdir results/foldseek_remote
```
