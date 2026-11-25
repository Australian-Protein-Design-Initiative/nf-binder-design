# FoldSeek Structural Search

Run [FoldSeek](https://github.com/steineggerlab/foldseek) structural similarity searches against protein structure databases directly within the pipeline.

Useful for identifying structurally similar proteins to your designed binders and annotating designs with known structural folds.

!!! note "FoldSeek is used for two different jobs"

    **Database search** (`--do_foldseek_search`, or `--method foldseek`) compares each design
    against a database such as CATH50, answering "has this fold been seen before?". That is
    what this page covers.

    **Design clustering** (`--do_foldseek_cluster`) compares the designs against each other
    with no database involved, answering "how many structurally distinct binders did this
    run produce?". See [Design clustering](#design-clustering) below.

    They have separate flag namespaces and separate output directories, and either, both or
    neither can be enabled.

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
| `--foldseek_include_html_output` | `false` | Include HTML alignment report (off by default: HTML reports can reach hundreds of megabytes for large result sets) |

### CATH Annotation

When using a database whose name starts with `CATH` (e.g. `CATH50`), FoldSeek results are **automatically** annotated with CATH hierarchy descriptions without requiring additional flags.

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
| `Alphafold/UniProt` | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | None | All AlphaFold predicted structures. Very large. Only available locally. |
| `Alphafold/UniProt50` | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | `afdb50` | AlphaFold structures at 50% clustering. Reduced redundancy. |
| `Alphafold/UniProt50-minimal` | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | None | Minimal version. Only available locally. |
| `Alphafold/Proteome` | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | `afdb-proteome` | AlphaFold proteome-level structures. |
| `Alphafold/Swiss-Prot` | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | `afdb-swissprot` | AlphaFold structures for Swiss-Prot reviewed proteins. |
| `ESMAtlas30` | [ESM Atlas](https://esmatlas.com) | `mgnify_esm30` | ESM Metagenomic Atlas at 30% clustering. Metagenomic predicted structures. |
| `BFMD` | [BFMD](https://foldseek.steineggerlab.workers.dev/bfmd.version) | `bfmd` | BFMD database. |
| `BFVD` | [BFVD](https://bfvd.steineggerlab.workers.dev) | `BFVD` | BFVD database. |
| `ProstT5` | [HuggingFace](https://huggingface.co/Rostlab/ProstT5) | None | Protein language model for 3Di prediction. Not a search database: only available locally. |

> **Tip:** For binder design validation, `CATH50` is a suitable default: it covers known structural domains and provides fast execution. Use `PDB` for comprehensive searches against all experimentally determined structures.

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

Interactive alignment report (`--format-mode 3`). Can be very large for big result sets: set `--foldseek_include_html_output true` to enable.

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

## Design clustering

Enabled with `--do_foldseek_cluster`, this clusters a run's own designs against each other
rather than searching a database, to show how many structurally distinct binders the run
produced. Supported by the `rfd4_proteina` workflow.

```bash
nextflow run main.nf --method rfd4_proteina \
    --input_pdb target.pdb --contigs "50-120, /0, A18-132" \
    --do_foldseek_cluster
```

By default it runs `foldseek easy-cluster` over the extracted binder chains, using the
same thresholds as RFD4-Proteina's own `scripts/launch_generate.sh`.

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--do_foldseek_cluster` | `false` | Enable all-vs-all clustering of the designs |
| `--foldseek_cluster_max` | all | Maximum designs to cluster |
| `--foldseek_cluster_tmscore` | `0.5` | TM-score cluster boundary; lower values merge more designs |
| `--foldseek_cluster_coverage` | `0.0` | Coverage threshold (`-c`) |
| `--foldseek_cluster_cov_mode` | `0` | Coverage mode (`--cov-mode`) |
| `--foldseek_cluster_min_seq_id` | `0.0` | Sequence identity threshold; `0` keeps clustering purely structural |
| `--foldseek_cluster_alignment_type` | `1` | Alignment type; `1` is TM-align |
| `--foldseek_cluster_multimer` | `false` | Cluster whole complexes with `easy-multimercluster` |
| `--foldseek_cluster_multimer_tm` | `0.7` | `--multimer-tm-threshold` |
| `--foldseek_cluster_chain_tm` | `0.7` | `--chain-tm-threshold` |
| `--foldseek_cluster_interface_lddt` | `0.3` | `--interface-lddt-threshold` |
| `--foldseek_cluster_extra_args` | — | Extra arguments passed to `foldseek` |

`--foldseek_cluster_multimer` clusters whole complexes instead of binder chains, requiring
every chain to match. For a set of designs that share one target this is much stricter than
clustering the binders alone.

`easy-interfacecluster` is not yet available. FoldSeek release 10 (January 2025) is still
the most recent release and does not include it; interface clustering was added to the
development branch in February 2026, so it needs a build from master rather than a newer
release or container.

### Outputs

| Path | Contents |
|------|----------|
| `<outdir>/foldseek_cluster/clusters.tsv` | FoldSeek's representative/member table |
| `<outdir>/foldseek_cluster/cluster_table.tsv` | One row per design, keyed on `id` |
| `<outdir>/foldseek_cluster/cluster_rep_seq.fasta` | Representative sequences |
| `<outdir>/foldseek_cluster/cluster_report.tsv` | Multimer report (with `--foldseek_cluster_multimer`) |

`cluster_table.tsv` is merged into `combined_scores.tsv`, adding:

| Column | Meaning |
|--------|---------|
| `foldseek_cluster_rep` | Id of the cluster's representative design |
| `foldseek_cluster_id` | Cluster number; `1` is the most populated cluster |
| `foldseek_cluster_size` | Number of designs in this cluster |
| `foldseek_cluster_is_rep` | `1` for the representative, `0` otherwise |

Filtering a ranked table to `foldseek_cluster_is_rep == 1` gives one design per distinct
fold, which is usually what an ordering by score alone fails to provide.

### Usage as a subworkflow

```groovy
include { FOLDSEEK_CLUSTER_DESIGNS } from '../subworkflows/local/foldseek_cluster'

workflow {
    FOLDSEEK_CLUSTER_DESIGNS(Channel.fromPath('designs/*.cif'))
}
```

It emits `cluster_tsv`, `cluster_table` and `rep_seq`. The binder-chain extraction process
is included under an alias, so a workflow that also runs `FOLDSEEK_SEARCH` can use both in
one execution.
