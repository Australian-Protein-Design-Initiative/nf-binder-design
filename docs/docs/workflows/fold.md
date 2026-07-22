# Fold Workflow (`fold.nf`)

Standalone multi-method structure prediction for monomer FASTA inputs. Predicts
structures with any combination of AlphaFold2, Boltz-2, RosettaFold3 and
Protenix, sharing one MSA-generation stage, then (by default) clusters the
ensemble with EnGens.

> **Phase 1: monomer inputs only.** A FASTA with more than one record fails
> fast. Multimer / paired-MSA support is planned for a later phase. See
> [`examples/fold/`](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/tree/main/examples/fold).

## Overview

`fold.nf` is a **standalone** entry point (not `--method` under `main.nf`):

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design/fold.nf \
  --input 'input/*.fasta' \
  --outdir results \
  --methods af2,boltz,rf3,protenix \
  --msa_method jackhmmer_af2 \
  -profile local
```

From a git clone:

```bash
nextflow run /path/to/nf-binder-design/fold.nf \
  --input UL119_domain.fasta \
  --outdir results \
  --methods af2,boltz,rf3,protenix \
  --msa_method jackhmmer_af2 \
  -profile slurm,m3
```

Each input FASTA is one prediction unit (single sequence). Shared MSAs feed all
selected predictors; optional MSA subsample and EnGens clustering produce a
conformational ensemble from the combined predictions.

## Command-line Options

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design/fold.nf --help
```

### Key Parameters

| Flag | Description |
|------|-------------|
| `--input` | Single FASTA, glob, or directory of FASTA files (required) |
| `--outdir` | Output directory (default: `results`) |
| `--methods` | Comma-separated: `af2`, `boltz`, `rf3`, `protenix` (default: `af2`) |
| `--msa_method` | `jackhmmer_af2` (default) or `mmseqs2_colabfold` |
| `--n_predictions` | Total structures per input, per method (split by method batch size) |
| `--msa_subsample` | Off by default; `true` (default depth list) or a custom `max_seq:max_extra_seq` list. Depths with `max_seq >=` MSA size are skipped |
| `--msa_subsample_include_full` | Keep one full-MSA job when subsampling (default: `true`) |
| `--skip_engens` | Skip post-prediction EnGens clustering |
| `--engens_clustering` | `hdbscan` (default), `gmm`, `km`, or comma-separated |

Method-specific flags (`--af2_*`, `--boltz_*`, `--rf3_*`, `--protenix_*`) are
documented in `--help`. Seeds are unset by default so each engine draws its own
random seed (pin `--af2_random_seed` / `--boltz_seed` / `--rf3_seed` /
`--protenix_seeds` for reproducibility; do not inject a fresh random seed on
every CLI invocation if you want `-resume` to cache).

## Multiple Sequence Alignments (MSAs)

All selected predictors share one MSA stage controlled by `--msa_method`.

### Option 1: AlphaFold jackhmmer / HHblits (`jackhmmer_af2`)

Default route. Runs AlphaFold's MSA pipeline against local genetic databases
under `--af2_db_path`, then converts the resulting MSAs to an a3m for
Boltz / RF3 / Protenix. AF2 predict reuses the same features (including
templates) unless MSA subsample rebuilds features from a shallow a3m.

```bash
--msa_method jackhmmer_af2 \
--af2_db_path /path/to/alphafold_dbs \
--af2_db_preset full_dbs
```

Databases used under `--af2_db_path` (paths are fixed relative to that root):

| Database | Role |
|----------|------|
| UniRef90 | jackhmmer |
| MGnify | jackhmmer |
| BFD + UniRef30 | HHblits (`full_dbs`) |
| PDB70 + PDB mmCIF | templates |

Individual DB paths are **not** separate CLI knobs — set the root with
`--af2_db_path` (or `params.af2_db_path` in a config). See
[Setting up databases](#setting-up-databases).

### Option 2: ColabFold MMseqs2 (`mmseqs2_colabfold`)

#### Remote server

```bash
--msa_method mmseqs2_colabfold --use_remote_server true
```

Queries the public ColabFold MMseqs2 API (no local DBs). Suitable for small /
occasional runs; for heavy use prefer local databases.

#### Local databases

```bash
--msa_method mmseqs2_colabfold \
--uniref30 /path/to/colabfold_dbs/uniref30 \
--colabfold_envdb /path/to/colabfold_dbs/colabfold_envdb
```

`--uniref30` must contain `uniref30_*` MMseqs2 DB files;
`--colabfold_envdb` must contain `colabfold_envdb*` files (layout produced by
[`scripts/download_colabfold_dbs.sh`](#colabfold-mmseqs2-databases)).

> There is no site-wide default ColabFold DB path on M3 yet — use
> `--use_remote_server true` or install local DBs yourself.

### MSA subsample (optional)

CF-random-style shallow random MSAs per predict task:

```bash
--msa_subsample true
# or a custom list, e.g. --msa_subsample '1:2,8:16,64:128'
```

`--msa_subsample true` uses depths `1:2,2:4,...,64:128`.
`--msa_subsample_include_full` (default `true`) also keeps one full-MSA job
(AF2 retains jackhmmer templates on that path). Depths with `max_seq` greater
than or equal to the MSA sequence count are skipped (they would only shuffle
the full MSA). For each depth job (including full), `fold/msa_ids/` lists
`header_line<TAB>id` where `header_line` is the 0-based file line of that
sequence's `>` header in the original a3m (disambiguates duplicate accessions).

## Example Usage

Minimal AF2-only run with jackhmmer MSAs:

```bash
nextflow run /path/to/nf-binder-design/fold.nf \
  --input input/pdl1.fasta \
  --outdir results \
  --methods af2 \
  --msa_method jackhmmer_af2 \
  --af2_db_path /mnt/datasets/alphafold/alphafold_20240229 \
  -profile local
```

Multi-method ensemble (25 structures per method) with ColabFold remote MSA:

```bash
nextflow run /path/to/nf-binder-design/fold.nf \
  --input UL119_domain.fasta \
  --outdir results \
  --methods af2,boltz,rf3,protenix \
  --msa_method mmseqs2_colabfold \
  --use_remote_server true \
  --n_predictions 25 \
  --af2_keep_models all \
  -profile slurm,m3
```

See `examples/fold/run-m3.sh` and `examples/fold/run-local.sh` for complete
HPC / workstation wrappers (including Apptainer bind mounts for AF2 DBs).

## EnGens clustering

After prediction, EnGens runs by default (UMAP + HDBSCAN) and writes
`results/engens/<id>/clusters.html` plus representative conformations.

| Flag | Description |
|------|-------------|
| `--skip_engens` | Skip clustering |
| `--engens_clustering` | `hdbscan` (default), `gmm`, `km`, or comma-separated |
| `--engens_min_structures` | Minimum structures before clustering (default: 3) |
| `--engens_max_clusters` | Upper bound for auto cluster-count search |

To cluster an existing folder of `.cif` / `.pdb` without re-folding, use the
standalone `engens.nf` workflow:

```bash
nextflow run /path/to/nf-binder-design/engens.nf \
  --input results/fold/predictions/ \
  --id UL119_domain \
  --outdir results \
  -profile slurm,m3
```

## Output

Default layout under `--outdir`:

```
results/
├── fold/
│   ├── msa/<msa_method>/     # shared MSAs + a3m (jackhmmer_af2 or mmseqs2_colabfold)
│   ├── af2/msas/             # AF2-only features.pkl (not under fold/msa/)
│   ├── af2/ … boltz/ … rf3/ … protenix/
│   ├── predictions/          # flat gather: af2_*, boltz_*, rf3_*, protenix_* mmCIF
│   ├── msa_ids/              # when --msa_subsample: header_line<TAB>id (0-based '>' line)
│   ├── params.json
│   └── logs/
└── engens/<id>/              # clusters.html + representative conformations (HDBSCAN by default)
```

## Setting up databases

You need local databases only for:

- `--msa_method jackhmmer_af2` (AlphaFold genetic DBs), and/or
- `--msa_method mmseqs2_colabfold` **without** `--use_remote_server true`
  (ColabFold MMseqs2 DBs).

Model weights for Boltz / RF3 / Protenix are typically baked into the pipeline
containers; AF2 params are downloaded with the AlphaFold DB tree (`params/`).

Helper scripts live in the repo [`scripts/`](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/tree/main/scripts)
directory.

### AlphaFold genetic databases

Official source: [google-deepmind/alphafold](https://github.com/google-deepmind/alphafold)
([genetic databases](https://github.com/google-deepmind/alphafold#genetic-databases)).

**Requirements:** `aria2c`, `rsync`, `git`. Full databases are ~556 GB download
and ~2.62 TB unzipped (SSD recommended).

```bash
# From a clone of nf-binder-design:
./scripts/download_alphafold_dbs.sh /data/alphafold_dbs
# or reduced set:
./scripts/download_alphafold_dbs.sh /data/alphafold_dbs reduced_dbs
```

This clones DeepMind's repo shallowly and runs their
`scripts/download_all_data.sh`, which fetches BFD (or small BFD), MGnify,
PDB70, PDB mmCIF, UniRef30, UniRef90, UniProt, PDB seqres, and model params.

Equivalent manual invocation:

```bash
git clone --depth 1 https://github.com/google-deepmind/alphafold.git
bash alphafold/scripts/download_all_data.sh /data/alphafold_dbs full_dbs
```

Point the fold workflow at the download root:

```bash
--msa_method jackhmmer_af2 \
--af2_db_path /data/alphafold_dbs \
--af2_db_preset full_dbs
```

Expected layout (abbreviated):

```
$AF2_DB_PATH/
  bfd/                 # full_dbs only
  small_bfd/           # reduced_dbs only
  mgnify/mgy_clusters_2022_05.fa
  params/
  pdb70/
  pdb_mmcif/
  uniref30/UniRef30_2021_03*
  uniref90/uniref90.fasta
  uniprot/             # multimer (future)
  pdb_seqres/          # multimer (future)
```

`fold.nf` currently hardcodes the relative paths used by the DeepMind download
scripts (e.g. `mgnify/mgy_clusters_2022_05.fa`, `uniref30/UniRef30_2021_03`).
The M3 default is `/mnt/datasets/alphafold/alphafold_20240229` (group
`alphafold`); bind-mount it in Apptainer `runOptions` when using `-profile m3`
(see `examples/fold/nextflow.m3.config`).

Ensure the tree is readable by compute jobs (`chmod -R a+rX` if needed).

### ColabFold MMseqs2 databases

Official downloads and setup notes: [colabfold.mmseqs.com](https://colabfold.mmseqs.com/)
and ColabFold's
[`setup_databases.sh`](https://github.com/sokrypton/ColabFold/blob/main/setup_databases.sh).

**Requirements:** `mmseqs` in `PATH`, plus `aria2c` or `curl`/`wget`. Indexing
is memory-heavy (ColabFold documents on the order of hundreds of GB RAM for
full indexes / single-query search with indexes preloaded).

```bash
# Install MMseqs2 first, then:
./scripts/download_colabfold_dbs.sh /data/colabfold_dbs
```

The helper fetches upstream `setup_databases.sh`, downloads UniRef30 +
ColabFold env DB (prebuilt expandable-profile archives by default), runs
`mmseqs createindex` unless `MMSEQS_NO_INDEX=1`, and organises outputs into:

```
/data/colabfold_dbs/
  uniref30/           # pass to --uniref30
  colabfold_envdb/    # pass to --colabfold_envdb
```

Useful environment overrides:

| Variable | Effect |
|----------|--------|
| `SKIP_TEMPLATES=1` | Skip PDB mmCIF / Foldseek template downloads |
| `MMSEQS_NO_INDEX=1` | Skip `createindex` (smaller disk; slower search) |
| `DOWNLOADS_ONLY=1` | Download archives only |
| `GPU=1` | GPU-capable indexes (needs GPU-enabled MMseqs2) |
| `UNIREF30DB` / `CFDB` | Archive stems (defaults: `uniref30_2302`, `colabfold_envdb_202108`) |

Equivalent manual setup:

```bash
wget https://raw.githubusercontent.com/sokrypton/ColabFold/main/setup_databases.sh
chmod +x setup_databases.sh
./setup_databases.sh /data/colabfold_dbs
# Then point --uniref30 / --colabfold_envdb at dirs containing
# uniref30_* and colabfold_envdb* MMseqs2 files (or use the helper).
```

Databases provided by ColabFold (see [colabfold.mmseqs.com](https://colabfold.mmseqs.com/)):

1. **UniRef30** — 30% identity clustered UniRef100
2. **ColabFold env DB** — environmental sequences (BFD/MGnify-derived plus
   metagenomic sources); alternatively BFD/MGnify-only archives are listed on
   the download page
3. Optional template DBs (PDB100, etc.)

Use with fold:

```bash
--msa_method mmseqs2_colabfold \
--uniref30 /data/colabfold_dbs/uniref30 \
--colabfold_envdb /data/colabfold_dbs/colabfold_envdb
```

Bind-mount `/data/colabfold_dbs` (or the paths you pass) into Apptainer when
running under Singularity/Apptainer profiles.

### Which MSA route should I use?

| Situation | Recommendation |
|-----------|----------------|
| Site already has AF2 DBs (e.g. M3 `/mnt/datasets/alphafold/...`) | `--msa_method jackhmmer_af2 --af2_db_path …` |
| No local DBs, small number of sequences | `--msa_method mmseqs2_colabfold --use_remote_server true` |
| Heavy ColabFold-style search on your cluster | Install local ColabFold DBs with `scripts/download_colabfold_dbs.sh` |

## Related

- Example run directory: [`examples/fold/`](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/tree/main/examples/fold)
- Boltz Pulldown also accepts `--uniref30` / `--colabfold_envdb` for local MSAs
  ([Boltz Pulldown](boltz-pulldown.md))
- Standalone EnGens: `engens.nf`
