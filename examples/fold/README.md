Multi-method structure folding with the standalone `fold.nf` workflow.

See the [Fold workflow docs](../../docs/docs/workflows/fold.md) for MSA options,
EnGens, outputs, and AlphaFold / ColabFold database setup
(`scripts/download_alphafold_dbs.sh`, `scripts/download_colabfold_dbs.sh`).

Monomer **and multimer** inputs are supported. A single-record FASTA folds as a
monomer; a multi-record FASTA folds as one complex (chains A, B, C, …) with a
taxonomically-paired MSA per engine (see the [multimer docs](../../docs/docs/workflows/fold.md#multimer-complexes)).

`input/pdl1.fasta` folds human PD-L1 as a single-chain monomer. For a multimer
(protein complex) run, see the sibling [`examples/fold-multimer`](../fold-multimer)
example (2-chain PD-L1 homodimer with per-engine paired MSAs).

## Param harmonization

| Shared | AF2 | Boltz | RF3 | Protenix |
|---|---|---|---|---|
| `--input`, `--outdir` | `--af2_db_path` | `--use_msa_server` (Boltz's own MSA server, independent of `--msa_method`) | `--rf3_ckpt_path` | `--protenix_seeds` (`--seeds`) |
| `--methods` (`af2,boltz,rf3,protenix`) | `--af2_model_preset` | `--templates` | `--rf3_num_steps` | `--protenix_cycle` (`--cycle`) |
| `--msa_method` (`jackhmmer_af2`\|`mmseqs2_colabfold`) | `--af2_db_preset` | `--boltz_recycling` (`--recycling_steps`) | `--rf3_n_recycles` | `--protenix_step` (`--step`) |
| `--use_remote_server`, `--uniref30`, `--colabfold_envdb` (ColabFold local/remote MSA) | `--af2_max_template_date` | `--boltz_batch_size` (`--diffusion_samples` per job) | `--rf3_batch_size` (`diffusion_batch_size` per job) | `--protenix_batch_size` (`--sample` per job) |
| `--n_predictions` (total structures; split by method batch size) | `--af2_random_seed` | `--boltz_sampling_steps` (`--sampling_steps`) | `--rf3_early_stopping_plddt_threshold` | `--protenix_model_name`, `--protenix_use_msa` |
| `--msa_subsample` (`false`\|`true`\|`max:extra,...`), `--msa_subsample_include_full` | `--af2_keep_models`, `--af2_no_relax` | `--boltz_seed` | `--rf3_seed` | `--protenix_seeds` |
| `--gpu_devices` | | | | |
| EnGens (default on): `--skip_engens`, `--engens_clustering` (`hdbscan`\|`gmm`\|`km`), `--engens_min_structures`, `--engens_max_clusters` | | | | |

Method-namespaced params (`--af2_*`, `--boltz_*`,
`--rf3_*`, `--protenix_*`) are kept explicit rather than unified, since
recycles/samples mean different things per engine (see `fold.nf --help`).

Protenix's per-chain a3m contract (confirmed against
`protenix:v2.0.0-weights` on 2026-07-17) is the `unpairedMsaPath` field on a
`proteinChain` input entry (monomer), plus `pairedMsaPath` for multimer - same
per-chain contract as Boltz's `msa:`/RF3's `msa_path`. Weights
(`protenix_base_default_v1.0.0.pt` etc) are baked into the container under
`/models/protenix/checkpoint`, so no download happens at predict time.

## Multimer / paired MSAs

- One canonical taxonomy parse (`bin/msa_taxonomy.py`) renders each engine's
  native paired format: RF3 `TaxID=` a3m, Protenix species-mnemonic
  paired/unpaired a3m, Boltz `key,sequence` CSV. AF2 uses its own native
  multimer pipeline. See `plans/fold-nf-multimer-paired-msa.md`.
- Use `--msa_method jackhmmer_af2`: only its rich headers carry the taxonomy
  pairing needs. ColabFold headers are taxonomy-less, so ColabFold multimer runs
  unpaired — use `--use_msa_server true` (Boltz) instead.
- AF2 multimer needs the 2021 snapshot (`alphafold_20211129`, has
  `uniprot/`+`pdb_seqres/`); the monomer-only default (`alphafold_20240229`)
  cannot be used (fold.nf fails fast). Its HHblits DB is `uniclust30`, so set
  `--af2_uniref30_subpath` (see [`examples/fold-multimer`](../fold-multimer)).

## Running

To run on M3 (submits each stage via SLURM):

```bash
./run-m3.sh
```

To run locally (e.g. on a GPU workstation):

```bash
./run-local.sh
```

To fold a multimer (protein complex) instead, see the sibling
[`examples/fold-multimer`](../fold-multimer) example.

Both default to `--methods af2,boltz,rf3,protenix --msa_method jackhmmer_af2`.
Pass `--methods` to select a subset - e.g. AlphaFold2 only (the native
jackhmmer/hhblits MSA + GPU predict route):

```bash
./run-m3.sh --methods af2
```

Optional CF-random-style MSA subsample (shallow random MSA depths per predict
task; `--msa_subsample true` uses depths `1:2,2:4,...,64:128`, or pass a custom
`max_seq:max_extra_seq` list). Depths with `max_seq >=` MSA size are skipped;
`--msa_subsample_include_full` (default) keeps one full-MSA job. Full-MSA jobs
keep AF2 templates; shallow AF2 jobs rebuild `features.pkl` from the subsampled
a3m without templates:

```bash
./run-m3.sh --methods af2,boltz --msa_subsample true --n_predictions 5
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
- With `--msa_subsample`, sequence ID lists for each depth job under
  `results/fold/msa_ids/` (`header_line<TAB>id`; 0-based `>` line in the
  original a3m; filenames include method, batch/run and `_msa<depth>`).
  Depths with `max_seq >=` MSA size are skipped.
- EnGens (unless `--skip_engens`): `results/engens/<id>/clusters.html` plus
  representative structures under
  `results/engens/<id>/clustering/<featurizer>/<gmm|km|hdbscan>/conformations/`
  (e.g. `residue_mindist`, `backbone_torsions`, `backbone_torsions-residue_mindist`).
- `results/fold/params.json` written on completion.

To re-cluster an existing folder of structures without re-running prediction:

```bash
nextflow run engens.nf --input results/fold/predictions/ --id UL119_domain \
    --outdir results -profile slurm,m3
```
