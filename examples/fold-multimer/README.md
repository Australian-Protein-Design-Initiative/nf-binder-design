Multimer (protein complex) folding with `--method fold`.

See the [Fold workflow docs](../../docs/docs/workflows/fold.md#multimer-complexes)
for the full multimer / paired-MSA strategy. This example is the multimer
counterpart of [`examples/fold`](../fold) (which folds monomers).

`input/complex.fasta` is a 2-record FASTA that folds as one complex: each record
is a chain (A, B in file order). Here it is the **human PD-L1 / PD-1 ectodomain
complex** — chain A human PD-L1 (UniProt Q9NZQ7, the ectodomain from PDB 3BIK)
and chain B human PD-1 (UniProt Q15116, the IgV ectodomain, residues 21–147).
Swap in any multi-record FASTA (up to 26 chains) to fold a different complex; a
homo-oligomer is expressed as repeated identical records.

`input/complex.human-mouse.fasta` is an alternate input using the original 3BIK
pairing — human PD-L1 + **mouse** PD-1 (UniProt Q02242) — for comparison against
the human/human default (`input/complex.human-human.fasta` is an identical copy
of `input/complex.fasta`).

## How multimer pairing works here

One FASTA → per-chain MSA search → one canonical taxonomy parse
(`bin/fold/msa_taxonomy.py`) renders each engine's native paired-MSA format:

| Engine | Pairing key | Fed |
|--------|-------------|-----|
| AF2 | native multimer pipeline (species pairing) | whole complex + `--model_preset=multimer`, 2021 DB snapshot |
| RF3 | numeric `TaxID=` | per-chain a3m with `TaxID=` headers |
| Protenix | species mnemonic (`_HUMAN`, `_9BETA`) | per-chain `pairedMsaPath` + `unpairedMsaPath` |
| Boltz-2 | taxid `key` | per-chain `key,sequence` CSV |
| OpenFold3 | species mnemonic | per-chain `colabfold_main.a3m` + `uniprot_hits.a3m` (`tr\|ACC\|ACC_SPECIES/1-N` headers, used for pairing only) |
| ESMFold2 | `key=<taxid>` (done by ESMFold2 itself) | per-chain a3m with `key=<taxid>` headers; rows with no taxonomy are kept |

RF3 / Protenix / Boltz-2's rendered per-chain files are published under
`results/fold/msa/paired/`. AF3, OpenFold3 and ESMFold2 (not shown above; see
the [Fold docs](../../docs/docs/workflows/fold.md#paired-msas-how-each-engine-differs))
render their own pairing input inline in their predict/input-prep tasks
instead.

**`--msa_method jackhmmer_af2` is required for paired multimers** — only its
rich UniProt/UniRef headers carry taxonomy. ColabFold headers are taxonomy-less,
so a ColabFold multimer folds unpaired; for that route use `--use_msa_server true`
(Boltz fetches + pairs its own MSA) and drop `af2` from `--methods`.

## AF2 databases (2021 snapshot)

AF2 multimer needs the 2021 snapshot (`alphafold_20211129`), which ships
`uniprot/` + `pdb_seqres/`; the default `alphafold_20240229` is monomer-only and
the fold workflow fails fast if `af2` is requested for a multimer against it.
The snapshot's DB filenames differ from the 20240229 defaults, so
`nextflow.m3.config` overrides `--af2_uniref30_subpath` (uniclust30),
`--af2_mgnify_subpath` (2018_12), `--af2_uniprot_subpath` and
`--af2_pdb_seqres_subpath`.

`nextflow.m3.config` also overrides `--af2_data_dir` to a host params
directory, so both DB snapshots stay bind-mounted on M3 for `run_alphafold.py`
to read weights from; this is not required in general — the `alphafold2`
container already bundles **multimer_v3** weights alongside the monomer/ptm
ones at its default `--af2_data_dir` (`/app/alphafold`).

## Running

On M3 (submits each stage via SLURM):

```bash
./run-m3.sh
```

Locally (e.g. a GPU workstation with the 2021 DB mounted):

```bash
./run-local.sh
```

Pass `--methods` to select a subset, e.g. skip AF2 (and its 2021-DB dependency)
and let Boltz pair its own MSA:

```bash
./run-m3.sh --methods boltz,rf3,protenix,openfold3
# or the MSA-server route for Boltz:
./run-m3.sh --methods boltz --use_msa_server true
```

## Outputs

Same layout as `examples/fold` (`results/params.json` and `results/logs/` at
the outdir root; per-method predictions under `results/fold/`, a flat mmCIF
gather in `results/fold/predictions/`), plus the per-chain paired MSAs under
`results/fold/msa/paired/`.

## AlphaFold3 + Protenix variant

`run-af3-m3.sh` / `run-af3-local.sh` run only `--methods af3,protenix` into
`results-af3/`. AlphaFold3 weights are not bundled: download them first with
`../../models/download_af3_weights.sh` (after reading the terms it prints), or
point `AF3_MODEL_DIR` at an existing weights directory, eg:

```bash
AF3_MODEL_DIR=/path/to/af3_weights ./run-af3-m3.sh
```
