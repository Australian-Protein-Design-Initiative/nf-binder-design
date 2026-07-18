# `fold.nf` — multi-method structure folding workflow

## STATUS

**Phase 1 (monomer) — DONE & verified on M3 (2026-07-17).** `fold.nf --methods af2,boltz,rf3
--msa_method jackhmmer_af2` folded PD-L1 end-to-end: AF2 pLDDT 91.9–92.8 (reproduces the
standalone af2.nf result), Boltz confidence 0.87 / pTM 0.73 / pLDDT 0.91, RF3 pLDDT 0.834 /
pTM 0.76. The shared jackhmmer MSA + `af2_msas_to_a3m` converter fed Boltz and RF3 cleanly.
Output layout harmonized to `results/{af2,boltz,rf3}/`. The critical a3m→AF2 `features.pkl`
bridge assumption was verified against the container (`run_alphafold.py` loads a pre-written
`features.pkl` under `--use_precomputed_msas`). `ipsae.py` degrades gracefully on monomers.

**Still open (deferred):**
- **ColabFold route (`--msa_method mmseqs2_colabfold`) NOT yet run** — no local ColabFold DBs
  found on M3 (`--uniref30`/`--colabfold_envdb`); needs those or `--use_remote_server true`
  (compute-node internet). This is the route that exercises the a3m→AF2 bridge dynamically.
- **Phase 2: multimer + paired MSAs** (§4). Generators are written N-chain-generic but gated
  off by fold.nf's monomer-only guard. **Cross-engine species-pairing is the core Phase 2
  problem** (monomer needs none): AF2-multimer pairs internally via the `uniprot` DB (2021
  snapshot only); Boltz/RF3 need paired a3m blocks; **Protenix** needs species/taxonomy that
  `colabfold_search` omits — use Protenix's bundled `scripts/colabfold_msa.py` (in the
  container) which adds pseudo-taxonomy IDs to the paired MSA, per
  https://github.com/bytedance/Protenix/blob/main/docs/colabfold_compatible_msa.md
  (feeds `pairedMsaPath`+`unpairedMsaPath` per chain). Monomer Protenix uses just
  `unpairedMsaPath` (a single a3m) — already wired, no pairing fixup needed.
- **Phase 3: protenix** — container `ghcr.io/.../protenix:v2.0.0-weights` (pending mirror).
- **Robustness:** run the Nextflow driver as its own `sbatch` job (login-node/session
  processes get reaped mid-run — cost us two killed runs, recovered via `-resume`).

## Phase 3 spec — add `protenix` (monomer)

Add Protenix as a 4th `--methods` engine, mirroring the RF3/Boltz folding subworkflows.
Monomer-only (consistent with Phase 1); multimer is Phase 2 for all engines.

**Container:** module `container` directive = `oras://ghcr.io/australian-protein-design-initiative/containers/protenix:v2.0.0-weights`
(weights baked in, like rf3's rc-foundry). **M3 pulls fail from ghcr.io** — add a
`withName: PROTENIX_FOLD` override in `conf/platforms/monash_containers.config` →
`${monash_container_base}/ghcr.io-australian-protein-design-initiative-containers-protenix-v2.0.0-weights.img`
(already mirrored/confirmed present).

**Probe FIRST (read-only, in the mirrored .img — DO NOT pull from ghcr.io):**
1. `protenix --help` / `protenix predict --help` (or whatever the entrypoint is) — the predict
   subcommand, input flag (JSON? `--input`), `--out_dir`, seeds, `--n_sample`/diffusion samples,
   `--n_step`/`--cycle`/recycles, model/checkpoint flag (weights baked in — confirm the path or
   that it's the default), and GPU behavior.
2. **Input JSON schema** — Protenix uses an AF3-style JSON (list of sequences with
   `proteinChain`/`sequence`, etc.). Determine the exact schema for a single protein chain.
3. **Precomputed MSA** — CRITICAL for consistency with boltz/rf3: how to feed our shared
   ColabFold/jackhmmer-derived **a3m** so Protenix does NOT run its own MSA search (e.g. a
   `precomputed_msa_dir`, an `msa` field in the JSON, or a `--use_msa`/`--msa_dir` flag).
   Protenix expects paired+unpaired a3m for multimers, but for monomer a single a3m/`pairing`
   is enough — find the monomer contract. If Protenix can't take our a3m cleanly, document it
   and fall back to letting Protenix build its own MSA (note the inconsistency).
4. **Output layout** — cif/pdb + confidence/summary JSON; where written under `--out_dir`.

**Files (mirror rf3):**
- `subworkflows/local/protenix_fold.nf` (`PROTENIX_FOLD`): `take tuple(meta, fasta, a3m)` →
  `GENERATE_PROTENIX_INPUT` → `PROTENIX_FOLD` process; `emit predictions, confidence_json`.
- `modules/local/protenix/protenix_fold.nf` (GPU, same GPU-guard idiom as `rf3_fold.nf`,
  publishDir `${params.outdir}/protenix/<meta.id>/` via the recursive-glob + `saveAs` strip
  pattern — see rf3_fold.nf; watch for the same double-nesting trap).
- `modules/local/protenix/generate_protenix_input.nf` + `bin/make_protenix_input.py`
  (FASTA + optional a3m → Protenix input JSON; N-chain-generic but monomer-gated by fold.nf).
- **fold.nf:** add `'protenix'` to `VALID_METHODS`, `params.protenix_*` defaults (seeds,
  n_sample, n_step/cycle as probed), a `FOLD_MSA.out.for_protenix` emit (same
  `tuple(meta,fasta,a3m)` as `for_rf3`), and `if ('protenix' in methods) PROTENIX_FOLD(...)`.
- **Config:** `withName: PROTENIX_FOLD` GPU resources in `nextflow.config` + `conf/platforms/m3.config`
  (`--gres=gpu:1 --partition=gpu`), container mirror override in `monash_containers.config`.
- **examples/fold:** update README + run scripts to include `protenix` in the `--methods` example.

**Verify on M3:** `--methods protenix --msa_method jackhmmer_af2` on `input/pdl1.fasta` (reuse
the cached jackhmmer MSA via `-resume` if possible). Confirm a sane monomer pLDDT/confidence and
`results/protenix/pdl1/` layout consistent with af2/boltz/rf3. Then a combined
`--methods af2,boltz,rf3,protenix` smoke test.

**Note:** a local ColabFold `mmseqs2_colabfold` test run is IN PROGRESS writing to
`examples/fold/results-colabfold/` (driver + SLURM jobs owned by the coordinator). Do NOT kill
SLURM jobs, touch `results-colabfold/`, or launch a conflicting `fold.nf` run in
`examples/fold/` while it's live — coordinate via the coordinator.

## Context

Extend the standalone AF2 work (`af2.nf`) into a general **folding dispatcher**: one
entrypoint `fold.nf` that predicts structures for FASTA inputs with any combination of
`--methods af2,boltz,rf3` (protenix later), sharing a single MSA-generation stage with a
selectable `--msa_method`. This deliberately overlaps `boltz_pulldown.nf`; it is a
**standalone prototype** now and will be folded into a renamed/expanded `boltz_pulldown`
inside `main.nf` **later (not in this work)**. No `main.nf` edits.

### Decisions locked in with the user
- **Standalone `fold.nf`** (root-level `workflow {}`, like `af2.nf`/`boltzgen_filter.nf`).
  Reusable per-method subworkflows under `subworkflows/local/` so the future
  `boltz_pulldown` merge is a plain `include`.
- **`af2.nf` is refactored** to call a new `ALPHAFOLD2` subworkflow; the existing
  `examples/af2/` entrypoint keeps working unchanged.
- **RF3: build a NEW standalone folding module**, independent of the binder-design-coupled
  `ROSETTAFOLD3` (which is hard-wired to a 2-body target+binder complex + a binder
  postprocess). Must fold a monomer OR a multi-chain multimer. Harmonization with the rfd3
  pipeline's RF3 is explicitly deferred to the future `boltz_pulldown` merge.
- **`--msa_method jackhmmer_af2 | mmseqs2_colabfold`** as a shared knob. ColabFold MSAs
  must be able to feed **all three** engines, including AF2 — so we build an a3m→AF2
  `msas/` bridge (approved, despite being custom).
- **Inputs: monomer + protein multimer.** No ligands this round. Multimers must use
  **correct species-paired MSAs** (see §4) — this is a first-class requirement, not an
  afterthought.

### Ground-truth findings from repo/M3 investigation (do not re-derive)
- **AF2 (existing, keep):** `modules/local/af2/{alphafold2_jackhmmer_msa,alphafold2}.nf`,
  container `.../alphafold_cuda12_upstream-...sif`, `python /app/alphafold/run_alphafold.py`.
- **Boltz:** `modules/local/common/boltz.nf` (`BOLTZ`, container
  `ghcr.io/.../boltz:v2.2.1-2`) consumes a **YAML** built by `bin/create_boltz_yaml.py`;
  Boltz reads a per-chain **a3m** directly (`msa: <path>`). There's an **unwired**
  `modules/local/common/create_boltz_yaml.nf` with a `CREATE_BOLTZ_YAML_MONOMER` variant —
  best template for a generic single-target fold. Boltz tuning currently rides only on
  `task.ext.args` (no first-class params).
- **RF3 (reference only):** `modules/local/rfd3/rosettafold3.nf` (`ROSETTAFOLD3`, container
  `oras://ghcr.io/.../rc-foundry:0.2.0-weights`), command `rf3 fold inputs=<json> out_dir=...
  ckpt_path= num_steps= n_recycles= diffusion_batch_size= early_stopping_plddt_threshold=`.
  Input JSON is a component list (per-chain `{seq, chain_id, msa_path?}`) built by
  `bin/rfd3/make_rf3_input_spec.py` (only emits target+binder). RF3 reads an a3m per chain.
- **ColabFold MSA:** `modules/local/common/mmseqs_colabfoldsearch.nf`
  (`MMSEQS_COLABFOLDSEARCH`, container `quay.io/nf-core/proteinfold_colabfold:1.1.1`),
  emits `**.a3m`. Local mode uses `colabfold_search` against `--uniref30`+`--colabfold_envdb`
  (both default `false`, no M3 default path); remote mode uses `bin/colabfold_remote_msa.py`.
  **Currently produces UNPAIRED per-chain a3m only** (no `--use-env-pairing`/pairing).
- **No a3m→AF2 `msas/` converter exists anywhere** — must be built.
- **AF2 DB on M3:** default `params.alphafold2_db_path=/mnt/datasets/alphafold/alphafold_20240229`
  is **monomer-only** (no `uniprot/`, no `pdb_seqres/`). The **2021 snapshot**
  `/mnt/datasets/alphafold/alphafold_20211129` HAS `uniprot/` (102G) + `pdb_seqres/` and is a
  complete multimer set (bfd/mgnify/pdb70/uniref90 symlinked; note it ships `uniclust30`, not
  `uniref30`). **AF2 multimer is only possible against the 2021 snapshot.**
- `paramsToMap`+`onComplete` params.json helper is copy-pasted in `main.nf` and `af2.nf`;
  factor into `lib/`. No `--methods`-style multi-value parsing exists; closest idiom is
  `params.x.split(',').collect{ it.trim() }` (bindcraft).

---

## Architecture

```
fold.nf                              # standalone entrypoint: parse --methods/--msa_method,
                                     #   build input channel, call FOLD_MSA once, dispatch
lib/ParamsHelper.groovy              # shared paramsToMap (de-dupe from main.nf/af2.nf)
subworkflows/local/
  fold_msa.nf        # FOLD_MSA: (meta,fasta) + msa_method -> { af2_msas, a3m } (per engine need)
  alphafold2.nf      # ALPHAFOLD2 subworkflow (wraps existing 2 modules) — af2.nf reuses this
  boltz_fold.nf      # BOLTZ_FOLD: (meta,fasta,a3m?) -> Boltz prediction + confidence
  rosettafold3_fold.nf # ROSETTAFOLD3_FOLD: (meta,fasta,a3m?) -> RF3 prediction + confidence
modules/local/af2/
  colabfold_a3m_to_af2_msas.nf   # NEW: bridge a3m -> AF2 precomputed msas/ layout
  af2_msas_to_a3m.nf             # NEW: derive a single a3m from AF2 jackhmmer msas/ (for boltz/rf3)
modules/local/boltz/
  fold_create_boltz_yaml.nf      # generic FASTA->YAML (adapt the unwired MONOMER module)
modules/local/rf3/
  rf3_fold.nf                    # NEW generic RF3 fold process (monomer/multimer, no binder coupling)
  generate_rf3_fold_input.nf     # NEW generic FASTA(+a3m) -> rf3 input JSON
bin/
  colabfold_a3m_to_af2_msas.py   # NEW converter (see §4)
  af2_msas_to_a3m.py             # NEW converter (see §4)
  make_rf3_fold_spec.py          # NEW generic RF3 JSON generator (monomer/multimer)
examples/fold/                   # run-m3.sh, run-local.sh, nextflow.m3.config, input/, README
```

Reuse `MMSEQS_COLABFOLDSEARCH` and the existing AF2 MSA process as-is (no edits to their
current callers). Keep `bin/create_boltz_yaml.py` (extend only if needed).

---

## §1 Input & meta contract (harmonized)

`fold.nf` input follows **`af2.nf`'s convention** (one FASTA file = one prediction unit;
multiple records in a file = multimer chains). Accept a single file, a glob, or a directory:

```groovy
def p = params.input
ch_input = ( file(p).isDirectory() ? Channel.fromPath("${p}/*.{fasta,fa,faa}")
                                    : Channel.fromPath(p) )
    .map { f -> [ [id: f.baseName, n_chains: countRecords(f)], f ] }
```

`meta.id = f.baseName` (matches AF2's output-subdir contract). `meta.n_chains` (parsed via
`splitFasta`) drives monomer-vs-multimer branching (model preset, pairing). Every subworkflow
`take`s `tuple(meta, fasta)` (+ MSA) and `emit`s a harmonized `tuple(meta, path(prediction_dir))`.

---

## §2 `fold.nf` dispatch

```groovy
params.methods    = 'af2'                 // comma list: af2,boltz,rf3
params.msa_method = 'jackhmmer_af2'        // jackhmmer_af2 | mmseqs2_colabfold
def methods = params.methods.split(',').collect { it.trim().toLowerCase() }
def valid = ['af2','boltz','rf3']
methods.each { assert it in valid : "Unknown method '${it}' (valid: ${valid})" }

FOLD_MSA(ch_input, params.msa_method)      // one MSA stage, shared across methods
if ('af2'   in methods) ALPHAFOLD2(FOLD_MSA.out.af2_msas)
if ('boltz' in methods) BOLTZ_FOLD(FOLD_MSA.out.for_boltz)
if ('rf3'   in methods) ROSETTAFOLD3_FOLD(FOLD_MSA.out.for_rf3)
```

Help block + `onComplete` params.json (shared helper). Outputs published per method under
`${params.outdir}/{af2,boltz,rf3}/...` (each subworkflow owns its publishDir, mirroring the
existing `af2/{msas,predictions}` layout).

---

## §3 MSA subworkflow (`FOLD_MSA`)

`FOLD_MSA` centralizes MSA generation and adapts each engine's required format. It computes
only what the selected methods need (wire emits lazily / guard by `params.methods`).

| `--msa_method` | native form | AF2 needs `msas/` dir | Boltz/RF3 need a3m |
|---|---|---|---|
| `jackhmmer_af2` | AF2 `msas/` (per-source .sto/.a3m/.hhr) | ✓ native (ALPHAFOLD2_JACKHMMER_MSA) | convert via `af2_msas_to_a3m` |
| `mmseqs2_colabfold` | ColabFold `.a3m` | bridge via `colabfold_a3m_to_af2_msas` | ✓ native |

**Emits:** `af2_msas = tuple(meta, fasta, msas_dir)`, `for_boltz = tuple(meta, fasta, a3m)`,
`for_rf3 = tuple(meta, fasta, a3m)`.

**Two new converters (both containerized; run in the AF2 or colabfold container as
appropriate):**
- `colabfold_a3m_to_af2_msas.py` — a3m → a per-target `msas/` dir AF2's
  `--use_precomputed_msas` accepts. AF2 monomer reads `uniref90_hits.sto`, `mgnify_hits.sto`,
  `bfd_uniref_hits.a3m`/`bfd_uniref_hits.a3m`, `pdb_hits.hhr`. Strategy: write the ColabFold a3m
  as `bfd_uniref_hits.a3m` (a3m is accepted directly) and provide empty/degenerate `.sto`/`.hhr`
  for the other sources so AF2's parser is satisfied. **VERIFY** which files this AF2 build's
  monomer feature pipeline strictly requires (probe `run_alphafold.py`/`pipeline.py`), and that
  templates can be disabled cleanly (empty `pdb_hits.hhr`). Document quality caveat: ColabFold a3m
  ≠ the full jackhmmer(uniref90+mgnify)+hhblits(bfd) blend.
- `af2_msas_to_a3m.py` — AF2 `msas/` → single a3m for Boltz/RF3. Simplest robust route: take
  `bfd_uniref_hits.a3m` (already a3m) as the base and merge `uniref90_hits.sto`/`mgnify_hits.sto`
  (reformat sto→a3m, e.g. the container's `reformat.pl` or Biopython) into one a3m with the query
  first. Keep the query-first invariant Boltz/RF3 expect.

---

## §4 Paired MSAs for multimers (the risk area — explicit strategy)

The current jackhmmer (monomer preset) and ColabFold (unpaired) paths do **not** produce
species-paired multimer MSAs. Per engine:

- **AF2 multimer (jackhmmer_af2):** set `--model_preset=multimer`; AlphaFold's multimer data
  pipeline does species pairing **internally** using the `uniprot` all-seqs MSA. Requires the
  **2021 multimer snapshot** (`alphafold_20211129`: has `uniprot/`+`pdb_seqres/`). Add
  `--uniprot_database_path` + `--pdb_seqres_database_path` flags **only** in multimer mode; keep
  `pdb70`/template flags for monomer. **This is the recommended, correct paired-MSA route for AF2
  multimers.** Auto-select the DB snapshot by `meta.n_chains` OR require the user to point
  `--alphafold2_db_path` at a multimer-complete snapshot (validate presence of `uniprot/` when a
  multimer input is seen; fail fast with a clear message otherwise).
- **Boltz/RF3 multimer:** two options — (a) **MSA server** (`--use_msa_server`): Boltz/ColabFold
  server does species pairing automatically (simplest correct route; recommended default for
  multimers); (b) **local paired ColabFold**: enable `colabfold_search` pairing
  (`--use-env-pairing 1` / pair mode) to emit paired+unpaired a3m, then route per-chain. Boltz
  reads per-chain a3m; RF3 takes `msa_path` per chain. **VERIFY** how Boltz-2 v2.2.1 and this RF3
  build expect paired blocks to be delimited before committing to local pairing.

### §4a Shared paired-MSA builder (preferred Phase 2 direction, from user)

Rather than per-engine pairing hacks, build **one engine-agnostic paired-MSA step** and fan its
output to every engine. Pipeline:
1. `colabfold_search` (already have `MMSEQS_COLABFOLDSEARCH`) → per-chain **unpaired** a3m.
2. **Add species/taxonomy** that `colabfold_search` omits — reuse Protenix's bundled
   `scripts/colabfold_msa.py` (in the protenix container; adds pseudo-taxonomy IDs), or apply the
   same header-regex species extraction directly.
3. **Pair by species** — port/adapt the AlphaFold-Multimer-style algorithm Protenix implements in
   `protenix/data/msa/msa_utils.py`: `MSAPairingEngine.pair_chains_by_species()` +
   `get_species_ids()` (UniProt `(?:tr|sp)|...|..._<SPECIES>` and `UniRef100_<acc>_<species>`
   regexes → group rows by species → stack one representative per chain into paired rows). The
   core routine is reusable numpy/dict logic (depends on Protenix constants + the feature-dict
   shape, not its serialization) — reimplement at the **a3m level** in a standalone
   `bin/pair_msa.py` (per-chain a3m in → paired + unpaired a3m out) so it has no Protenix runtime
   dependency and can feed all engines.
4. **Fan out** the paired/unpaired a3m to each engine's native slot: Protenix
   `pairedMsaPath`+`unpairedMsaPath` per chain; Boltz per-chain `msa:`; RF3 `msa_path` per chain;
   (optionally) AF2-multimer via the ColabFold→AF2 features bridge extended to paired features.

This makes local (offline) multimer folding consistent across engines and removes the reliance on
each tool's own MSA server. AF2-multimer can still use its native `uniprot`-based internal
pairing (2021 snapshot) as the reference/cross-check. Candidate: a `PAIR_MSA` process/subworkflow
consumed by `FOLD_MSA` when `meta.n_chains > 1`.

**Recommended phasing:** implement **monomer end-to-end first** (all 3 methods, both MSA
methods, AF2 bridge) — DONE — then multimer via the §4a shared builder. Multimer AF2 can also use
the 2021 snapshot as a correctness reference; `log()` clearly when a multimer falls back to
unpaired MSAs so it's never silent.

---

## §5 Per-method subworkflows

### `ALPHAFOLD2` (refactor of `af2.nf` wiring)
- `take: ch(meta, fasta, msas_dir)`; call existing `ALPHAFOLD2` GPU module with precomputed MSAs
  (already built). `emit: predictions`. `af2.nf` becomes: build input → `FOLD_MSA` → `ALPHAFOLD2`.
- Multimer: pass `--model_preset=multimer` + multimer DB flags (§4). Keep monomer path identical
  to today (verified working, pLDDT ~91–93).

### `BOLTZ_FOLD`
- `take: ch(meta, fasta, a3m)`. Generate YAML from FASTA via a generic
  `FOLD_CREATE_BOLTZ_YAML` (adapt the unwired `CREATE_BOLTZ_YAML_MONOMER`; for multimer emit one
  `protein:` entry per record with `id: [A,B,...]` and each chain's `msa:`), then call the existing
  `BOLTZ` process (`step_name='boltz'`). Parse confidence via `bin/parse_boltz_confidence.py`.
  Publish `${outdir}/boltz/...`.
- New first-class params → `task.ext.args` (via a `withName: '.*:BOLTZ'` selector or an `ext.args`
  string built in the subworkflow): `--boltz_recycling`, `--boltz_diffusion_samples`,
  `--boltz_sampling_steps` (map to boltz's `--recycling_steps`/`--diffusion_samples`/
  `--sampling_steps`). Confirm exact flag names against `boltz predict --help` in v2.2.1-2.

### `ROSETTAFOLD3_FOLD` (new, decoupled)
- `take: ch(meta, fasta, a3m)`. New `GENERATE_RF3_FOLD_INPUT` runs `bin/make_rf3_fold_spec.py`
  to build the RF3 JSON generically: one `{seq, chain_id}` component per FASTA record
  (chains A,B,C…), attaching `msa_path` per chain when an a3m is provided; **no** target/binder
  roles, **no** template, **no** binder postprocess. Then new `RF3_FOLD` process runs
  `rf3 fold inputs=<json> out_dir=output ...` (same container/GPU-guard idiom as `ROSETTAFOLD3`,
  minus the postprocess). Emit + publish `${outdir}/rf3/...` (raw rf3 output +
  `*_summary_confidences.json`). Params: `--rf3_num_steps`, `--rf3_n_recycles`,
  `--rf3_diffusion_batch_size`, `--rf3_ckpt_path`, `--rf3_early_stopping_plddt_threshold`
  (reuse the rfd3 defaults: 50 / 10 / 5 / the ckpt path / 0.5).

---

## §6 CLI param harmonization

Shared (apply to all methods): `--input`, `--outdir`, `--methods`, `--msa_method`,
`--use_msa_server`, `--gpu_devices`/`--require_gpu` (existing). MSA DBs:
`--alphafold2_db_path` (jackhmmer/AF2), `--uniref30` + `--colabfold_envdb` (ColabFold);
document that multimer AF2 needs a multimer-complete `--alphafold2_db_path`.

Method-namespaced tuning (kept explicit rather than falsely unified, since recycles/samples
mean different things per engine): `--af2_*` (model_preset auto by n_chains, models_to_relax,
random_seed, num_predictions_per_model), `--boltz_*` (recycling, diffusion_samples,
sampling_steps), `--rf3_*` (num_steps, n_recycles, diffusion_batch_size, ckpt_path,
early_stopping_plddt_threshold). Rename current `--alphafold2_*` params? Keep as-is for
`af2.nf` back-compat; `fold.nf` can accept both. Provide a harmonization table in the help
block and README.

---

## §7 Config additions

- `nextflow.config`: `withName:` resources for new processes (RF3_FOLD → accelerator=1 like
  ROSETTAFOLD3; converters → local/light; FOLD_CREATE_BOLTZ_YAML → local). Reuse existing
  BOLTZ / MMSEQS_COLABFOLDSEARCH / ALPHAFOLD2* selectors.
- `conf/platforms/m3.config`: partitions for RF3_FOLD (`--gres=gpu:1 --partition=gpu`),
  converters on `comp`. MMSEQS already configured.
- `conf/platforms/monash_containers.config`: container overrides for new GPU processes
  (RF3_FOLD → the rc-foundry Monash mirror `.img`, as `ROSETTAFOLD3` already has).
- Example `examples/fold/nextflow.m3.config`: AF2 DB binds (as af2 example) **plus** ColabFold
  DB paths (`--uniref30`, `--colabfold_envdb`) and a `-B` for the 2021 multimer snapshot when used.

---

## §8 Example `examples/fold/`
- `input/pdl1.fasta` (monomer, reuse), `input/complex.fasta` (2-record multimer test).
- `run-m3.sh` (SLURM) and `run-local.sh` mirroring the af2 example (keep `sg alphafold`).
  Demonstrate `--methods af2,boltz,rf3 --msa_method mmseqs2_colabfold`.
- `README.md` with the param harmonization table and the multimer/paired-MSA notes.

---

## §9 Verification (end-to-end, incremental)
0. **Probes first (read-only):** (a) `boltz predict --help` in `boltz:v2.2.1-2` → confirm
   `--recycling_steps`/`--diffusion_samples`/`--sampling_steps` names. (b) In the AF2 container,
   confirm exactly which `msas/` files the monomer pipeline requires under
   `--use_precomputed_msas` (for the bridge) and that empty template hits are accepted.
   (c) `rf3 fold` accepts a generic multi-component JSON with `msa_path` and no template.
1. **Monomer, jackhmmer_af2:** `--methods af2,boltz,rf3` on `pdl1.fasta`. AF2 matches the known
   good result (pLDDT ~91–93). Boltz/RF3 produce structures + confidences. Verify the
   `af2_msas_to_a3m` a3m feeds Boltz/RF3.
2. **Monomer, mmseqs2_colabfold:** same input; verify the a3m→AF2 bridge yields a sane AF2
   prediction (compare pLDDT to route 1), and Boltz/RF3 consume the a3m natively.
3. **Multimer:** `complex.fasta`. AF2 multimer against the 2021 snapshot (paired MSA internal);
   Boltz/RF3 via `--use_msa_server`. Sanity-check interface confidence (ipTM/PAE).
4. `nextflow config` resolves selectors; help block prints on no `--input`; `params.json` written.

---

## §10 Risks / open items
1. **a3m→AF2 bridge fidelity** — vanilla AF2's precomputed-MSA path is picky; the bridge may
   need the exact per-source filenames and non-empty query. Highest-risk new code — probe first.
2. **AF2 multimer DB** — only the 2021 snapshot works; it ships `uniclust30` (not `uniref30`) and
   older uniref90/mgnify/bfd. Confirm the container's `run_alphafold.py` multimer flags accept it.
3. **Local paired ColabFold MSA** — deferred in favor of `--use_msa_server` for multimers; if
   local pairing is needed, verify per-engine paired-block format.
4. **RF3 generic JSON** — `make_rf3_fold_spec.py` is new; confirm chain_id handling and that
   omitting templates/roles is valid for `rf3 fold`.
5. **Scope** — implement monomer fully first; multimer second. `log()` any unpaired-MSA fallback.
6. **Not touching `main.nf`/`boltz_pulldown.nf`** this round (future merge).
```
