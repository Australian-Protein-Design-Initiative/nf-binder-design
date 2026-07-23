Multimer (protein complex) folding with the standalone `fold.nf` workflow.

See the [Fold workflow docs](../../docs/docs/workflows/fold.md#multimer-complexes)
for the full multimer / paired-MSA strategy. This example is the multimer
counterpart of [`examples/fold`](../fold) (which folds monomers).

`input/complex.fasta` is a 2-record FASTA that folds as one complex: each record
is a chain (A, B in file order). Here it is the **human PD-L1 / PD-1 ectodomain
complex** — chain A human PD-L1 (UniProt Q9NZQ7, the ectodomain from PDB 3BIK)
and chain B human PD-1 (UniProt Q15116, the IgV ectodomain, residues 21–147).
Swap in any multi-record FASTA (up to 26 chains) to fold a different complex; a
homo-oligomer is expressed as repeated identical records.

Two alternate inputs are kept alongside it: `input/complex.human-mouse.fasta`
(the original 3BIK pairing — human PD-L1 + **mouse** PD-1, UniProt Q02242) and
`input/complex.human-human.fasta` (a copy of the default). Results from the
human/mouse fold are archived under `results.human-mouse/`.

## How multimer pairing works here

One FASTA → per-chain MSA search → one canonical taxonomy parse
(`bin/msa_taxonomy.py`) renders each engine's native paired-MSA format:

| Engine | Pairing key | Fed |
|--------|-------------|-----|
| AF2 | native multimer pipeline (species pairing) | whole complex + `--model_preset=multimer`, 2021 DB snapshot |
| RF3 | numeric `TaxID=` | per-chain a3m with `TaxID=` headers |
| Protenix | species mnemonic (`_HUMAN`, `_9BETA`) | per-chain `pairedMsaPath` + `unpairedMsaPath` |
| Boltz-2 | taxid `key` | per-chain `key,sequence` CSV |

The rendered per-chain files are published under `results/fold/msa/paired/`.

**`--msa_method jackhmmer_af2` is required for paired multimers** — only its
rich UniProt/UniRef headers carry taxonomy. ColabFold headers are taxonomy-less,
so a ColabFold multimer folds unpaired; for that route use `--use_msa_server true`
(Boltz fetches + pairs its own MSA) and drop `af2` from `--methods`.

## AF2 databases (2021 snapshot)

AF2 multimer needs the 2021 snapshot (`alphafold_20211129`), which ships
`uniprot/` + `pdb_seqres/`; the default `alphafold_20240229` is monomer-only and
`fold.nf` fails fast if `af2` is requested for a multimer against it. The
snapshot's DB filenames differ from the 20240229 defaults, so
`nextflow.m3.config` overrides `--af2_uniref30_subpath` (uniclust30),
`--af2_mgnify_subpath` (2018_12), `--af2_uniprot_subpath` and
`--af2_pdb_seqres_subpath` (all verified against the on-disk layout).

One extra wrinkle: the container's `run_alphafold.py` loads **multimer_v3**
weights, which the 2021 snapshot lacks (it only ships v1 multimer params). Since
`--data_dir` only locates `params/`, the config reads weights from the 20240229
snapshot (`--af2_data_dir`) while keeping the genetic DBs on the 2021 one.

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
./run-m3.sh --methods boltz,rf3,protenix
# or the MSA-server route for Boltz:
./run-m3.sh --methods boltz --use_msa_server true
```

## Outputs

Same layout as `examples/fold` (per-method predictions under `results/fold/`, a
flat mmCIF gather in `results/fold/predictions/`, `results/fold/params.json`),
plus the per-chain paired MSAs under `results/fold/msa/paired/`.
