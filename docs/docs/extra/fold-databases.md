# Fold: Setting up databases

Database setup for the [Fold](../workflows/fold.md) and [Fold Pulldown](../workflows/fold-pulldown.md)
workflows' `--msa_method` options.

You need local databases only for:

- `--msa_method jackhmmer_af2` (AlphaFold genetic DBs), and/or
- `--msa_method mmseqs2_colabfold` **without** `--use_remote_server true`
  (ColabFold MMseqs2 DBs).

Model weights for Boltz / RF3 / Protenix / OpenFold3 are baked into the pipeline
containers; AF2 params are downloaded with the AlphaFold DB tree (`params/`),
though the container also bundles a copy (see [AF2 multimer needs the 2021 DB
snapshot](../workflows/fold.md#af2-multimer-needs-the-2021-db-snapshot)).

Helper scripts live in the repo [`scripts/`](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/tree/main/scripts)
directory.

## AlphaFold genetic databases

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
  uniprot/             # AF2 multimer (2021 snapshot only)
  pdb_seqres/          # AF2 multimer (2021 snapshot only)
```

The fold workflow defaults the per-DB relative paths to this DeepMind
download-script layout (e.g. `mgnify/mgy_clusters_2022_05.fa`,
`uniref30/UniRef30_2021_03`); `--af2_uniref30_subpath`, `--af2_mgnify_subpath`,
`--af2_uniprot_subpath`, `--af2_pdb_seqres_subpath` and `--af2_pdb70_subpath`
override an individual one when a snapshot uses different filenames — needed
for AF2 multimer against the 2021 snapshot, whose HHblits DB is `uniclust30`.

On M3, `-profile m3` defaults to `alphafold_20211129`, which has `uniprot/` +
`pdb_seqres/` for AF2 multimer, together with its sub-paths (`uniclust30`,
`mgy_clusters_2018_12.fa`). `alphafold_20240229` is monomer-only and names
those two DBs differently, so pointing `--af2_db_path` at it under
`-profile m3` also needs `--af2_uniref30_subpath uniref30/UniRef30_2021_03
--af2_mgnify_subpath mgnify/mgy_clusters_2022_05.fa` (as
`examples/fold/nextflow.m3.config` does). Both live under `/mnt/datasets/alphafold/` (group
`alphafold`); bind-mount that path in Apptainer `runOptions` when using
`-profile m3` (see `examples/fold/nextflow.m3.config` /
`examples/fold-multimer/nextflow.m3.config`).

Ensure the tree is readable by compute jobs (`chmod -R a+rX` if needed).

## ColabFold MMseqs2 databases

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

## Which MSA route should I use?

| Situation | Recommendation |
|-----------|----------------|
| Site already has AF2 DBs (e.g. M3 `/mnt/datasets/alphafold/...`) | `--msa_method jackhmmer_af2 --af2_db_path …` |
| No local DBs, small number of sequences | `--msa_method mmseqs2_colabfold --use_remote_server true` |
| Heavy ColabFold-style search on your cluster | Install local ColabFold DBs with `scripts/download_colabfold_dbs.sh` |

## Related

- [Fold](../workflows/fold.md), [Fold Pulldown](../workflows/fold-pulldown.md)
- Boltz Pulldown also accepts `--uniref30` / `--colabfold_envdb` for local MSAs
  ([Boltz Pulldown](../workflows/boltz-pulldown.md))
