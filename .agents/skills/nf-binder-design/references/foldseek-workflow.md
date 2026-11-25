# FoldSeek Workflow — Full Reference

## Overview

FoldSeek is wired in for two unrelated purposes. Everything on this page except the
[Design clustering](#design-clustering) section is about the first.

| | Database search | Design clustering |
|--|--|--|
| Flag | `--do_foldseek_search`, `--method foldseek` | `--do_foldseek_cluster` |
| Question | Has this fold been seen before? | How many distinct binders did this run make? |
| Compares against | CATH50 / PDB / AFDB | The run's own designs |
| Subworkflow | `FOLDSEEK_SEARCH` | `FOLDSEEK_CLUSTER_DESIGNS` |
| Output | `foldseek/<database>/` | `foldseek_cluster/` |

The two have separate flag namespaces (`--foldseek_*` versus `--foldseek_cluster_*`), and
either, both or neither can be enabled in one run.

FoldSeek searches designed binder structures against protein structure databases to identify known folds and annotate designs with structural similarity.

Two usage modes:

1. **Inline flag** — add `--do_foldseek_search` to `rfd`, `rfd3`, `bindcraft`, or `boltzgen` to search outputs automatically
2. **Standalone workflow** — `--method foldseek` to search any PDB/mmCIF files

Default database: `CATH50` (CATH domains at 50% clustering). CATH annotation columns are added automatically for CATH databases.

## Inline: `--do_foldseek_search`

Add to any supported design workflow:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --input_pdb 'input/target.pdb' \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs 10 \
  --do_foldseek_search \
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

## Design clustering

`--do_foldseek_cluster` clusters a run's designs against each other, with no database, to
measure how many structurally distinct binders it produced. Wired into `rfd4_proteina`;
the `FOLDSEEK_CLUSTER_DESIGNS` subworkflow itself is method-agnostic.

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd4_proteina \
  --input_pdb 'input/PDL1.pdb' \
  --contigs "50-120, /0, A18-132" \
  --do_foldseek_cluster \
  --outdir results -profile local -resume
```

Runs `foldseek easy-cluster` over the extracted binder chains with the thresholds upstream
RFD4-Proteina uses in `scripts/launch_generate.sh`: `--alignment-type 1 --cov-mode 0
--min-seq-id 0.0 --tmscore-threshold 0.5`.

| Flag | Default | Description |
|------|---------|-------------|
| `--do_foldseek_cluster` | false | Enable clustering |
| `--foldseek_cluster_max` | all | Cap on designs clustered |
| `--foldseek_cluster_tmscore` | 0.5 | TM-score cluster boundary; lower merges more |
| `--foldseek_cluster_coverage` | 0.0 | `-c` |
| `--foldseek_cluster_cov_mode` | 0 | `--cov-mode` |
| `--foldseek_cluster_min_seq_id` | 0.0 | 0 keeps it purely structural |
| `--foldseek_cluster_alignment_type` | 1 | 1 is TM-align |
| `--foldseek_cluster_multimer` | false | `easy-multimercluster` over whole complexes |
| `--foldseek_cluster_multimer_tm` | 0.7 | `--multimer-tm-threshold` |
| `--foldseek_cluster_chain_tm` | 0.7 | `--chain-tm-threshold` |
| `--foldseek_cluster_interface_lddt` | 0.3 | `--interface-lddt-threshold` |
| `--foldseek_cluster_extra_args` | — | Passthrough |

### Outputs

`foldseek_cluster/clusters.tsv` (foldseek's representative/member table),
`cluster_table.tsv` (one row per design), `cluster_rep_seq.fasta`, and
`cluster_report.tsv` in multimer mode.

`cluster_table.tsv` is merged into `combined_scores.tsv`, adding `foldseek_cluster_rep`,
`foldseek_cluster_id` (1 = most populated), `foldseek_cluster_size` and
`foldseek_cluster_is_rep`. Filtering a ranked table to `foldseek_cluster_is_rep == 1`
gives one design per distinct fold.

### Implementation notes

- `FOLDSEEK_CLUSTER` runs foldseek; `FOLDSEEK_CLUSTER_TABLE` converts its output with
  `bin/foldseek/cluster_table.py`. Two processes because the foldseek image has no
  `python3` - running the script inside it fails with `python3: command not found`.
- The subworkflow includes `FOLDSEEK_PREPARE_QUERIES` under an alias, so a workflow running
  both the search and clustering does not reuse one process instance.
- Multimer mode skips chain extraction, passing the complexes in whole.
- `easy-interfacecluster` is **not** in foldseek release 10 (`941cd33`), in either the
  biocontainer or the rfd4-proteina image, and there is no newer container to switch to:
  release 10 (January 2025) is still the latest tagged release, and
  `quay.io/biocontainers/foldseek` has no tag beyond `10.941cd33`.

  `easy-interfacecluster` and `easy-interfacesearch` were added to master in February 2026
  (`data/easyinterfacecluster.sh`, registered in `src/FoldseekBase.cpp`) and are unreleased.
  Getting them needs a build from master - either upstream's continuous static build at
  `https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz`, which tracks master, or an APDI
  image built from master alongside the existing `foldseek:10-941cd33-arm64`.

  Its output names match `easy-cluster` (`<prefix>_cluster.tsv`, `<prefix>_rep_seq.fasta`),
  so `FOLDSEEK_CLUSTER` would need only a new mode branch, not new output handling.
  `--foldseek_cluster_multimer` is today's approximation.
