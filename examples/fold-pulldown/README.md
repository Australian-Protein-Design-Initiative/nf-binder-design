# Fold Pulldown — Mosaic Multispecifics × PD-L1 / IL-7Ra

Co-fold [Mosaic Multispecifics](https://proteinbase.com/collections/mosaic-multispecifics)
miniprotein binders (Escalante Bio) against **PD-L1** and **IL-7Ra**, using
`--method fold_pulldown` and ColabFold remote MSAs for targets.

See the [Fold Pulldown docs](../../docs/docs/workflows/fold-pulldown.md).

## Inputs

| File | Description |
|------|-------------|
| `input/binders.fasta` | 2 Mosaic binder sequences — a quick smoke-test subset |
| `input/binders.all.fasta` | All 10 Mosaic binder sequences (IDs from Proteinbase; sequences from the [Proteinbase master table](https://storage.proteinbase.com/proteinbase_all_data_28_01_2026.csv)) |
| `input/targets.fasta` | `PDL1` (115 aa from `examples/pdl1-rfd3/input/PDL1.pdb`) and `IL7RA` (198 aa from `IL7RA.cif`, converted to PDB) |
| `input/PDL1.pdb` / `input/IL7RA.pdb` | Source structures (pipeline uses FASTA; PDBs are for reference) |

The 10 binder IDs in `binders.all.fasta`: `radiant-bat-ruby`,
`bright-panther-frost`, `golden-cat-lava`, `solid-ram-wave`,
`bright-otter-reed`, `scarlet-kiwi-thorn`, `small-goat-ice`,
`quick-falcon-stone`, `brisk-tiger-oak`, `misty-cat-stone`.

`run-m3.sh` / `run-local.sh` default to `input/binders.fasta` (2 binders × 2
targets = 4 complexes). Pass `--binders input/binders.all.fasta` for the full
10-binder × 2-target (20-complex) pulldown:

```bash
./run-m3.sh --binders input/binders.all.fasta
```

With `--n_predictions 5` and five engines (see `run-m3.sh`), the full pulldown
predicts 500 structures (plus MSA jobs).

## MSA notes

- `--msa_method mmseqs2_colabfold --use_remote_server true --create_target_msa true`
  builds target MSAs via the ColabFold API (no local ColabFold DBs required).
- Binder MSAs are off (`--create_binder_msa false`): these are designed
  miniproteins with no useful homologs.
- **AF2** reuses the ColabFold target a3m (materialised into AF2's per-chain MSA
  files for the target chain, plus a multimer `features.pkl`); the binder chain
  stays query-only regardless of `--create_binder_msa`. Boltz / RF3 / Protenix /
  OpenFold3 use the same ColabFold a3ms directly.

To use native jackhmmer MSAs for AF2 (and for the shared a3m route) instead:

```bash
./run-m3.sh --msa_method jackhmmer_af2 --use_remote_server false
```

## Running

M3 / SLURM:

```bash
./run-m3.sh
```

Local GPU workstation:

```bash
./run-local.sh
```

Subset of engines or fewer samples:

```bash
./run-m3.sh --methods boltz,rf3 --n_predictions 1
```

## Outputs

Under `results/fold_pulldown/`:

- `pairs.tsv` — id (`<target>_and_<binder>`), target, binder
- `fold_pulldown_scores.tsv` — one row per predicted structure
- `fold_pulldown_summary.tsv` — per-(target, binder, tool) aggregates + z-scores
- `fold_pulldown_report.html` — Quarto overview
- Per-engine predictions / MSAs under the same publish tree
  (`<tool>/<target>_and_<binder>/`, `predictions/`, `msa/`)

`results/params.json` and `results/logs/` sit at the outdir root, not under
`results/fold_pulldown/`.

## AlphaFold3 + Protenix variant

`run-af3-m3.sh` / `run-af3-local.sh` run only `--methods af3,protenix` into
`results-af3/`. AlphaFold3 weights are not bundled: download them first with
`../../models/download_af3_weights.sh` (after reading the terms it prints), or
point `AF3_MODEL_DIR` at an existing weights directory, eg:

```bash
AF3_MODEL_DIR=/path/to/af3_weights ./run-af3-m3.sh
```
