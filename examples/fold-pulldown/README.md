# Fold Pulldown — Mosaic Multispecifics × PD-L1 / IL-7Ra

Co-fold the 10 [Mosaic Multispecifics](https://proteinbase.com/collections/mosaic-multispecifics)
miniprotein binders (Escalante Bio) against **PD-L1** and **IL-7Ra**, using all
fold engines (AF2, Boltz-2, RF3, Protenix) and ColabFold remote MSAs for targets.

See the [Fold Pulldown docs](../../docs/docs/workflows/fold-pulldown.md).

## Inputs

| File | Description |
|------|-------------|
| `input/binders.fasta` | 10 Mosaic binder sequences (IDs from Proteinbase; sequences from the [Proteinbase master table](https://storage.proteinbase.com/proteinbase_all_data_28_01_2026.csv)) |
| `input/targets.fasta` | `PDL1` (115 aa from `examples/pdl1-rfd3/input/PDL1.pdb`) and `IL7RA` (198 aa from `IL7RA.cif`, converted to PDB) |
| `input/PDL1.pdb` / `input/IL7RA.pdb` | Source structures (pipeline uses FASTA; PDBs are for reference) |

Binder IDs: `radiant-bat-ruby`, `bright-panther-frost`, `golden-cat-lava`,
`solid-ram-wave`, `bright-otter-reed`, `scarlet-kiwi-thorn`, `small-goat-ice`,
`quick-falcon-stone`, `brisk-tiger-oak`, `misty-cat-stone`.

This is a **20-complex** pulldown (2 targets × 10 binders). With
`--n_predictions 5` and four engines, expect on the order of ~400 predicted
structures (plus MSA jobs).

## MSA notes

- `--msa_method mmseqs2_colabfold --use_remote_server true --create_target_msa true`
  builds target MSAs via the ColabFold API (no local ColabFold DBs required).
- Binder MSAs are off (`--create_binder_msa false`): these are designed
  miniproteins with no useful homologs.
- **AF2** reuses the ColabFold target a3m (materialised into AF2’s per-chain MSA
  files for chain A, plus a multimer `features.pkl`); the binder chain stays
  query-only. Boltz / RF3 / Protenix use the same ColabFold a3ms directly.

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

- `fold_pulldown_scores.tsv` — one row per predicted structure
- `fold_pulldown_summary.tsv` — per-(target, binder, tool) aggregates + z-scores
- `fold_pulldown_report.html` — Quarto overview
- Per-engine predictions / MSAs under the same publish tree
