# Database download helpers

Scripts to install local MSA databases used by `fold.nf` (and Boltz Pulldown
local ColabFold search).

| Script | Purpose |
|--------|---------|
| `download_alphafold_dbs.sh` | AlphaFold genetic DBs via [DeepMind alphafold](https://github.com/google-deepmind/alphafold) `download_all_data.sh` |
| `download_colabfold_dbs.sh` | ColabFold MMseqs2 DBs via [setup_databases.sh](https://github.com/sokrypton/ColabFold/blob/main/setup_databases.sh) (includes `mmseqs createindex`) |

Full usage, layouts, and fold.nf flags: [docs — Fold / Setting up databases](../docs/docs/workflows/fold.md#setting-up-databases).

```bash
./scripts/download_alphafold_dbs.sh /data/alphafold_dbs
./scripts/download_colabfold_dbs.sh /data/colabfold_dbs
```
