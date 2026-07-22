Multi-method structure folding with the standalone `fold.nf` workflow.

**Phase 1 of `plans/fold-nf-multi-method-folding.md`: monomer inputs only.**
A FASTA file with more than one record (multimer) fails fast with a clear
error; multimer support (paired MSAs per method) is a later phase.

This example folds human PD-L1 (`input/pdl1.fasta`, reused from
`examples/pdl1-rfd/input/full/3BIK_B.fasta`) as a single-chain monomer with any
combination of AlphaFold2, Boltz-2, RosettaFold3 and Protenix, sharing one
MSA-generation stage selected by `--msa_method`.

`input/complex.fasta` is a 2-record placeholder kept only to exercise (and
demonstrate) the Phase 1 monomer-only guard - running `fold.nf` against it
fails fast with a "multimer folding not yet supported" error; it is not a
real multimer test case yet.

## Param harmonization

| Shared | AF2 | Boltz | RF3 | Protenix |
|---|---|---|---|---|
| `--input`, `--outdir` | `--af2_db_path` | `--use_msa_server` (Boltz's own MSA server, independent of `--msa_method`) | `--rf3_ckpt_path` | `--protenix_seeds` (`--seeds`) |
| `--methods` (`af2,boltz,rf3,protenix`) | `--af2_model_preset` | `--templates` | `--rf3_num_steps` | `--protenix_cycle` (`--cycle`) |
| `--msa_method` (`jackhmmer_af2`\|`mmseqs2_colabfold`) | `--af2_db_preset` | `--boltz_recycling` (`--recycling_steps`) | `--rf3_n_recycles` | `--protenix_step` (`--step`) |
| `--use_remote_server`, `--uniref30`, `--colabfold_envdb` (ColabFold local/remote MSA) | `--af2_max_template_date` | `--boltz_batch_size` (`--diffusion_samples` per job) | `--rf3_batch_size` (`diffusion_batch_size` per job) | `--protenix_batch_size` (`--sample` per job) |
| `--n_predictions` (total structures; split by method batch size) | `--af2_random_seed` | `--boltz_sampling_steps` (`--sampling_steps`) | `--rf3_early_stopping_plddt_threshold` | `--protenix_model_name`, `--protenix_use_msa` |
| `--gpu_devices` | `--af2_keep_models`, `--af2_no_relax` | `--boltz_seed` | `--rf3_seed` | `--protenix_seeds` |
| EnGens (default on): `--skip_engens`, `--engens_clustering` (`gmm`\|`km`), `--engens_min_structures`, `--engens_max_clusters` | | | | |

Method-namespaced params (`--af2_*`, `--boltz_*`,
`--rf3_*`, `--protenix_*`) are kept explicit rather than unified, since
recycles/samples mean different things per engine (see `fold.nf --help`).

Protenix's per-chain a3m contract (confirmed against
`protenix:v2.0.0-weights` on 2026-07-17) is the `unpairedMsaPath` field on a
`proteinChain` input entry - the shared `--msa_method` a3m is fed there
directly, same monomer contract as Boltz's `msa:`/RF3's `msa_path`. Weights
(`protenix_base_default_v1.0.0.pt` etc) are baked into the container under
`/models/protenix/checkpoint`, so no download happens at predict time.

## Multimer / paired-MSA notes (deferred to Phase 2)

- AF2 multimer needs the 2021 snapshot (`alphafold_20211129`, has
  `uniprot/`+`pdb_seqres/`) for internal species pairing; the current
  monomer-only default (`alphafold_20240229`) cannot be used for multimer.
- Boltz/RF3 multimer will default to `--use_msa_server` (MSA-server-side
  pairing) rather than local paired ColabFold search, at least initially.
- Protenix multimer needs each chain's `pairedMsaPath` set alongside
  `unpairedMsaPath` (see `protenix/data/msa/msa_featurizer.py`); deferred
  alongside the other engines.
- See `plans/fold-nf-multi-method-folding.md` §4 for the full strategy.

## Running

To run on M3 (submits each stage via SLURM):

```bash
./run-m3.sh
```

To run locally (e.g. on a GPU workstation):

```bash
./run-local.sh
```

Both default to `--methods af2,boltz,rf3,protenix --msa_method jackhmmer_af2`.
Pass `--methods` to select a subset - e.g. AlphaFold2 only (the native
jackhmmer/hhblits MSA + GPU predict route):

```bash
./run-m3.sh --methods af2
```

To try the ColabFold MSA route instead (bridging the resulting a3m into AF2 too):

```bash
./run-m3.sh --msa_method mmseqs2_colabfold --use_remote_server true
```

(No M3 default local ColabFold DB path exists yet - see
`plans/fold-nf-multi-method-folding.md` - so `--use_remote_server true` is
required unless you have your own `colabfold_search`-format `--uniref30`/
`--colabfold_envdb` DBs.)

## Outputs

- Shared MSAs under `results/fold/msa/<msa_method>/` (`jackhmmer_af2` or
  `mmseqs2_colabfold`), including the a3m derived for Boltz/RF3/Protenix.
- AF2-only `features.pkl` under `results/fold/af2/msas/` (not under `fold/msa/`).
- Per-method predictions under `results/fold/` and a flat gather of mmCIF
  structures in `results/fold/predictions/` (tool-prefixed filenames: `af2_`,
  `boltz_`, `rf3_`, `protenix_`).
- EnGens (unless `--skip_engens`): `results/engens/<id>/clusters.html` plus
  representative structures under
  `results/engens/<id>/clustering/<featurizer>/<gmm|km>/conformations/`
  (e.g. `residue_mindist`, `backbone_torsions`, `backbone_torsions-residue_mindist`).
- `results/fold/params.json` written on completion.

To re-cluster an existing folder of structures without re-running prediction:

```bash
nextflow run engens.nf --input results/fold/predictions/ --id UL119_domain \
    --outdir results -profile slurm,m3
```
