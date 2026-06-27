# FoldSeek Workflow — Full Reference

## Overview

FoldSeek searches designed binder structures against protein structure databases to identify known folds and annotate designs with structural similarity.

Two usage modes:

1. **Inline flag** — add `--do_foldseek` to `rfd`, `rfd3`, `bindcraft`, or `boltzgen` to search outputs automatically
2. **Standalone workflow** — `--method foldseek` to search any PDB/mmCIF files

Default database: `CATH50` (CATH domains at 50% clustering). CATH annotation columns are added automatically for CATH databases.

## Inline: `--do_foldseek`

Add to any supported design workflow:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --input_pdb 'input/target.pdb' \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs 10 \
  --do_foldseek \
  -profile local -resume
```

Works with `--method rfd`, `rfd3`, `bindcraft`, and `boltzgen`. FoldSeek runs on the final designed binder chains after the main pipeline completes.

## Standalone: `--method foldseek`

Search existing structure files without running a design workflow:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method foldseek \
  --input_pdbs 'results/designs/*.pdb' \
  --outdir results/foldseek \
  -profile local -resume
```

Supports PDB (`.pdb`) and mmCIF (`.cif`, `.cif.gz`).

### All hits with gzip output

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method foldseek \
  --input_pdbs 'results/designs/*.pdb' \
  --foldseek_database PDB \
  --foldseek_maxaccept 0 \
  --foldseek_gzip_output true \
  --outdir results/foldseek_pdb \
  -profile local -resume
```

## Key Parameters

Use `--method foldseek --help` for the full list.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_pdbs` | Required | Glob or path to input structures |
| `--foldseek_database` | `CATH50` | Database name (see table below) |
| `--foldseek_databases_path` | Auto-download | Persistent local database directory |
| `--foldseek_maxaccept` | `1` | Max alignments per query (`0` = unlimited) |
| `--foldseek_use_webserver` | `false` | Use FoldSeek web API instead of local search |
| `--foldseek_gzip_output` | `false` | Gzip TSV output |
| `--foldseek_include_html_output` | `false` | HTML report (can be very large) |

## Databases

| Name | Description |
|------|-------------|
| `CATH50` | CATH domains, 50% clustering — **default**, fast, good for fold classification |
| `PDB` | Full Protein Data Bank |
| `Alphafold/UniProt50` | AlphaFold structures, 50% clustering |

For binder validation, `CATH50` is usually sufficient. Use `PDB` for comprehensive experimental structure searches.

## Database Caching

Set `--foldseek_databases_path` on shared storage to avoid re-downloading (~1 GB for CATH50):

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method foldseek \
  --input_pdbs 'results/designs/*.pdb' \
  --foldseek_databases_path /scratch/shared/foldseek_dbs \
  --outdir results/foldseek \
  -profile slurm -resume
```

## Outputs

| File | Description |
|------|-------------|
| `foldseek_results.tsv` | Raw alignment results |
| `foldseek_results_annotated.tsv` | With CATH hierarchy columns (CATH databases only) |
| `foldseek_results.html` | Interactive report (if `--foldseek_include_html_output true`) |

Annotated columns for CATH: `cath_class`, `cath_architecture`, `cath_topology`, `cath_homologous_superfamily`, `cath_code`.
