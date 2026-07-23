# `fold.nf` — multimer support with taxonomically-paired MSAs

> **Implementation status (2026-07-23): DONE (pending cluster verification).**
> Implemented: `bin/msa_taxonomy.py` (canonical parse + RF3/Protenix/Boltz
> renderers) with pytest at `tests/bin/test_msa_taxonomy.py` (11 tests, green);
> `bin/split_complex_fasta.py`; per-chain MSA args in `make_rf3_fold_spec.py` /
> `make_protenix_input.py` + new `make_boltz_complex_yaml.py`; `fold.nf` guard
> removal + real `n_chains` + multimer validation (≤26 chains, non-empty seqs,
> af2-needs-jackhmmer, af2-needs-uniprot/, colabfold-unpaired warning,
> subsample-monomer-only); `FOLD_MSA` multimer branch (split → per-chain search
> → `ANNOTATE_MSA` → grouped per-tool bundles) with the monomer path kept
> byte-identical; `SPLIT_COMPLEX_FASTA` + `ANNOTATE_MSA` + `*_COMPLEX` generate
> processes; dual-input routing in the boltz/rf3/protenix subworkflows; AF2
> multimer DB-flag plumbing (multimer preset + uniprot/pdb_seqres, parametrised
> uniref30/uniclust30 subpath) in both AF2 modules; docs + `examples/fold/`
> multimer input/run script/config. Verified: pytest, Nextflow compile + `-preview`
> DAG build, pure-python end-to-end data-flow (split→render→generators), and all
> fail-fast/warn guards. **Not yet verified (needs the cluster — §6.1):** GPU
> end-to-end folds, actual per-engine pairing depth, and the 2021-snapshot
> uniclust30 sub-path + multimer DB vintage (§8 risk #2).

## Goal

Lift `fold.nf`'s Phase-1 monomer-only restriction so a multi-record FASTA folds
as a **protein complex** across all four engines (AF2, Boltz-2, RF3, Protenix).
The hard part is not the plumbing — the per-tool input generators are already
written N-chain-generic (see §2) — it is producing **correct species-/taxonomy-
paired MSAs** and presenting them in each tool's native format (§3, §4). No
ligands/nucleic acids this round: protein complexes only.

**Prior-art anchor.** This supersedes/executes Phase 2 of
[`fold-nf-multi-method-folding.md`](fold-nf-multi-method-folding.md) §4/§4a,
which already recorded the ground-truth investigation (AF2 2021 DB snapshot has
`uniprot/`; ColabFold module emits unpaired only; Protenix ships a pairing
script; `MSAPairingEngine.pair_chains_by_species()` is the reference algorithm).
Read that §4/§4a before starting — do not re-derive those findings.

EnGens on complexes is **not** a priority (§7): feed the complex `.cif`/`.pdb`
straight in, accept whatever it does, apply only minimal fixes if it crashes.

---

## §0 Probe findings (2026-07-23) — ground truth, do not re-derive

All §6.0 probes were run read-only against the cached SIFs. Results:

- **(a) AF2 multimer.** `run_alphafold.py` exposes `--model_preset=multimer`,
  `--num_multimer_predictions_per_model`, `--uniprot_database_path`,
  `--pdb_seqres_database_path`, `--uniref30_database_path`. The 2021 snapshot
  `/mnt/datasets/alphafold/alphafold_20211129` has **real** `uniprot/` + `pdb_seqres/`
  dirs (mode 750, group `alphafold`) and `params/`; `bfd`, `mgnify`, `pdb70`,
  `uniref90` are symlinks into `/mnt/reference/alphafold/alphafold_20210726`; it
  ships **`uniclust30`, not `uniref30`**. The default `alphafold_20240229` snapshot
  has neither `uniprot/` nor `pdb_seqres/` → monomer-only, confirmed.
- **(b) Boltz-2.** MSA input is dispatched **by file extension** in
  `boltz/main.py:615` (`.a3m → parse_a3m`, `.csv → parse_csv`). `parse_csv` requires
  columns **exactly `["key", "sequence"]`**; the `key` **is** the taxonomy id used
  for cross-chain pairing (empty/`nan`/`""` → `-1` = unpaired). `parse_a3m` only
  extracts taxonomy for `>UniRef100…` headers **and needs a taxonomy DB (Redis)**
  to map UniRef→taxid — not available offline. **⇒ For offline paired Boltz input,
  emit a per-chain `.csv` with `key=<taxid>`,`sequence`.** (`--use_msa_server` is
  the alternative: Boltz fetches+pairs itself.)
- **(c) RF3 (rc-foundry `rf3 fold`).** Inference reads **one `msa_path` a3m per
  chain component** (`rf3/utils/inference.py` → `build_msa_paths_by_chain_id_from_
  component_list`, in `atomworks/io/tools/inference.py`). `atomworks` `parse_a3m`
  extracts taxonomy via **`re.search(r"TaxID=(\d+)", header)`** and pairs across
  chains internally by matching tax id (`atomworks/ml/transforms/msa/
  _msa_pairing_utils.py`: "Joins MSAs by matching sequences with the same taxonomic
  ID"). **⇒ RF3 needs per-chain a3m whose hit headers contain `TaxID=<n>`.** It
  pairs itself; we do not pre-pair for RF3.
- **(d) Protenix.** Per-chain **`unpairedMsaPath` + `pairedMsaPath`** fields
  (`protenix/data/msa/msa_featurizer.py:613-639`; the older `precomputed_msa_dir`
  with `pairing.a3m`/`non_pairing.a3m` still works but warns). Pairing fires when
  `need_pairing = len(unique_prot_seqs) > 1` and runs
  `MSAPairingEngine.pair_chains_by_species()`, extracting species with
  `get_species_ids()` →
  `_UNIPROT_REGEX = (?:tr|sp)\|[A-Z0-9]{6,10}(?:_\d+)?\|[A-Z0-9]{1,10}_(?P<SpeciesId>[A-Z0-9]{1,5})`
  **or** `_UNIREF_REGEX = ^UniRef100_[^_]+_([^_/]+)`. **⇒ Protenix pairs on the
  species *mnemonic* (`_HUMAN`, `_9BETA`), NOT the numeric `TaxID=`.** It needs
  per-chain a3m headers in `tr|…|…_SPECIES` or compact `UniRef100_ACC_SPECIES` form.
- **(e) ColabFold.** `run_mmseqs2(..., use_pairing=True, pairing_strategy="greedy"|
  "complete")` posts to the **`ticket/pair`** endpoint (vs `ticket/msa`); local
  `colabfold_search` pairs via `--use-env-pairing 1 --db4 <pairing_db>
  --pairing_strategy N`. Our `bin/colabfold_remote_msa.py` currently hits only
  `ticket/msa` (unpaired).

**The decisive finding.** A real ColabFold/jackhmmer a3m header mix (from the
UL119 run) looks like:
```
>tr|Q8QRZ0|Q8QRZ0_9BETA Membrane glycoprotein UL119 OS=Panine betaherpesvirus 2 ...
>UniRef100_A0A7T7DGJ2 Ba152/Ba151 n=1 Tax=Baboon cytomegalovirus TaxID=120505 RepID=A0A7T7DGJ2_9BETA
```
RF3 parses `TaxID=120505` (matches the `UniRef100…` lines, **not** the `tr|` lines);
Protenix parses the species mnemonic `9BETA` (matches the `tr|` lines via UNIPROT
regex, **not** the space-separated `UniRef100 …` lines, whose compact
`UniRef100_ACC_SPECIES` form is absent). **⇒ No single a3m header form feeds both
RF3 and Protenix at full taxonomy coverage.** This forces a **shared-core + per-tool
renderer** design (§4), not one universal a3m.

---

## §1 What "multimer" means here + the meta contract change

- **One FASTA file = one complex.** Each record = one chain, in file order →
  chain IDs `A, B, C, …`. This is already the documented convention
  (`fold.nf` header, `make_rf3_fold_spec.py`, `make_protenix_input.py`).
- Homo-oligomers: represented as repeated identical records (simplest, keeps the
  generic path). A `count:`/copy-number shorthand is a **later** nicety, not this
  round — document that N copies = N records for now.
- `meta.n_chains` is already computed in `fold.nf` (line 358-360) and stashed on
  the meta map (currently hard-set to `1` at line 372). The only change needed:
  set `n_chains` to the real count and **remove the fail-fast guard**
  (lines 354-367). Everything downstream can branch on `meta.n_chains > 1`.

**Sub-task:** replace the guard with validation that is still strict about the
cases we genuinely can't do yet (e.g. > 26 chains → error; a record whose
sequence is empty → error).

---

## §2 What is already multimer-ready (do NOT rewrite)

Confirmed by reading the generators — each loops over all FASTA records and only
*gates* multi-chain via `len(sequences) == 1` on the **MSA** attachment, not the
sequence emission:

| File | Multimer-ready? | The one gated line |
|------|-----------------|--------------------|
| `bin/make_rf3_fold_spec.py` | ✅ emits one `{seq, chain_id}` per record | `if a3m_path and len(sequences) == 1:` (line 75) — only attaches MSA for monomer |
| `bin/make_protenix_input.py` | ✅ emits one `proteinChain` per record | `if a3m_path is not None and len(sequences) == 1:` (line 76) |
| `bin/create_boltz_yaml.py` | ⚠️ has target+binder complex mode but keyed on `--target/--binder`, not an N-record FASTA | needs a generic N-chain path (§4.2) |

So RF3 and Protenix generators need **only** a per-chain MSA argument added
(a list of a3m paths, or paired+unpaired paths for Protenix). Boltz needs a new
generic multi-chain YAML path. **The sequence-emission logic is done.**

The **RF3 / Protenix / Boltz predict processes** and their subworkflows do not
care about chain count — they consume whatever JSON/YAML the generator writes.
The publishDir/aggregation logic is per-model, not per-chain, so it is unaffected.

---

## §3 Per-tool multimer configuration & how each consumes a paired MSA

This is the crux the user flagged: **each tool presents multimers differently and
wants paired MSAs in a different format.** Ground-truth for the uncertain rows
must be re-probed in-container before coding (probes listed in §6.0).

### AF2 (AlphaFold2) — the odd one out
- **Different model.** Multimer needs `--model_preset=multimer` (the current
  default is `monomer_ptm`). This loads entirely different network weights
  (`params_model_*_multimer_v3`) and a **different data pipeline** that does
  species pairing **internally**.
- **Different databases.** Multimer pairing reads the `uniprot/` all-seqs DB +
  `pdb_seqres/`. The M3 default snapshot `alphafold_20240229` is **monomer-only**
  (no `uniprot/`). The **2021 snapshot `/mnt/datasets/alphafold/alphafold_20211129`**
  is the complete multimer set (`uniprot` 102 G, `pdb_seqres`; bfd/mgnify/pdb70/
  uniref90 symlinked; ships `uniclust30` not `uniref30`). **AF2 multimer is only
  possible against the 2021 snapshot.**
- **Recommended route:** for AF2 multimer, use `--msa_method jackhmmer_af2`
  pointed at the 2021 snapshot; AF2's own multimer pipeline builds the paired MSA
  correctly. This is the *reference-quality* paired MSA and needs no bespoke
  pairing code from us. Do **not** try to bridge our ColabFold a3m into AF2
  multimer features this round (paired-feature injection is fragile — defer).
- **`--num_multimer_predictions_per_model`** applies here (currently a param,
  `params.af2_num_predictions_per_model`); wire it in multimer mode.
- **fold.nf behaviour:** when `meta.n_chains > 1` and `af2` is requested,
  auto-select `--model_preset=multimer` and require `af2_db_path` to contain
  `uniprot/`; fail fast with a clear message pointing at the 2021 snapshot if not.

### Boltz-2
- **Presentation:** one `protein:` entry per chain in the YAML, each with its own
  `id:` and `msa:` path. Multiple copies of one chain → `id: [A, B]` list form.
- **Pairing — CONFIRMED (§0b):** dispatched by extension (`boltz/main.py:615`).
  Feed a per-chain **`.csv`** with columns exactly `key,sequence`, `key=<taxid>`;
  Boltz pairs rows across chains sharing a key (`""`/`nan` → unpaired). `.a3m`
  pairing is offline-infeasible (needs a Redis taxonomy DB). Alternative:
  `--use_msa_server` (Boltz fetches + pairs itself; already `params.use_msa_server`).

### RF3 (rc-foundry `rf3 fold`)
- **Presentation:** one `{seq, chain_id, msa_path?}` component per chain (already
  emitted). Homo-oligomer copies = repeated components with distinct `chain_id`.
- **Pairing — CONFIRMED (§0c):** RF3 reads one **a3m per chain** and pairs
  internally by matching `TaxID=(\d+)` parsed from hit headers (atomworks
  `parse_a3m` + `_msa_pairing_utils`). **⇒ render a per-chain a3m whose hit headers
  carry `TaxID=<n>`.** No pre-pairing needed.

### Protenix
- **Presentation:** one `proteinChain` per chain (already emitted).
- **Pairing — CONFIRMED (§0d):** each `proteinChain` takes explicit
  **`pairedMsaPath` + `unpairedMsaPath`** a3m files; Protenix pairs internally by
  **species mnemonic** (`_HUMAN`, `_9BETA`) via `_UNIPROT_REGEX`/`_UNIREF_REGEX` —
  **not** numeric `TaxID=`. **⇒ render the `pairedMsaPath` a3m with
  `tr|ACC|NAME_<SPECIES>` or compact `UniRef100_ACC_<SPECIES>` headers**;
  `unpairedMsaPath` = the plain merged a3m.

**Summary of divergence (all confirmed):** AF2 pairs internally (multimer model +
2021 `uniprot` DB); Boltz wants a `key,sequence` CSV keyed on taxid; RF3 wants an
a3m with `TaxID=` headers; Protenix wants an a3m with species-*mnemonic* headers.
Because RF3 and Protenix key on *different* header fields, the shared builder must
carry a canonical `(sequence, taxid, species_mnemonic)` record per hit and **render
each tool's format** (§4) — a single a3m cannot serve both at full coverage.

---

## §4 Paired-MSA generation strategy

### §4.1 The shared builder: one taxonomy-annotated per-chain MSA, rendered per tool

Because RF3 pairs on numeric `TaxID=` and Protenix on the species *mnemonic* (§0),
and each engine **pairs itself** from per-chain files, the builder does **not**
pre-pair into register. Instead it produces, per chain, a **canonical annotated
alignment** and then renders each tool's required file. Gated on `meta.n_chains > 1`.

**Step 1 — per-chain MSA search (already done).** Reuse the existing route
(`jackhmmer_af2` or `mmseqs2_colabfold`, remote or local). For a multi-record
FASTA the search runs **per chain** (one query each), producing one unpaired a3m
per chain — the same a3m fold.nf already builds for monomers, just N of them.

**Step 2 — canonical taxonomy annotation (`bin/msa_taxonomy.py`, NEW).** Parse
each hit header once and resolve a canonical record `(sequence, tax_id,
species_mnemonic)`, drawing from whatever the header carries (verified present in
real a3m headers, §0): `TaxID=(\d+)`; the `sp|/tr|…_<SPECIES>` mnemonic; the
`UniRef100_ACC … RepID=ACC_<SPECIES>` mnemonic; or `OS=`/`Tax=<name>` mapped to a
taxid via a small bundled name→taxid table when the numeric id is absent. Hits
with no resolvable taxonomy are kept as unpaired-only rows.

**Step 3 — render per tool (from the same canonical records):**
- **RF3 a3m** — rewrite each hit header to include `TaxID=<tax_id>` (drop rows with
  no tax_id from the *paired* portion; keep them in the unpaired tail). RF3's
  atomworks parser then pairs across chains by tax_id.
- **Protenix `pairedMsaPath` a3m** — rewrite headers to compact
  `>UniRef100_<acc>_<species_mnemonic>` (or `>tr|<acc>|X_<species_mnemonic>`) so
  `get_species_ids()` resolves them; `unpairedMsaPath` = the plain merged a3m.
- **Boltz `.csv`** — emit `key,sequence` with `key=<tax_id>` (blank when unknown →
  Boltz treats as unpaired). This is the cleanest offline pairing for Boltz.
- **AF2** — no rendering; uses native jackhmmer-multimer (§3), 2021 `uniprot` DB.

Keeping one canonical parse and N renderers (rather than N bespoke parsers) means
the taxonomy logic lives in exactly one tested place.

**Optional depth booster (follow-up, not v1): true pairing search.** ColabFold can
return deeper, mmseqs-paired hits via the **`ticket/pair`** endpoint
(`use_pairing=True`, `pairing_strategy` greedy/complete) — extend
`bin/colabfold_remote_msa.py` to add it; locally, `colabfold_search --use-env-pairing 1
--db4 <pairing_db> --pairing_strategy N` (needs the pairing DB, absent on M3 today).
This only *increases* the paired-row pool; the render/annotate step (§4.1 step 2-3)
is unchanged. Ship v1 with header-derived taxonomy from the normal search; add the
pair endpoint if paired depth proves too shallow. Whatever runs, `log()` the paired
row count per chain so a thin/absent pairing is never silent.

### §4.2 Wiring into `FOLD_MSA`

`FOLD_MSA` currently emits `for_boltz = for_rf3 = for_protenix = ch_a3m` (one
unpaired a3m). Changes:

- Keep the **monomer** emits exactly as today (single a3m, no annotation) so the
  verified monomer path is byte-for-byte unchanged.
- For `meta.n_chains > 1`: run the per-chain search (step 1), then a `PAIR_MSA`
  (a.k.a. `ANNOTATE_MSA`) process wrapping `bin/msa_taxonomy.py` (step 2-3),
  emitting per-chain rendered files. Emit tuples like
  `tuple(meta, fasta, [rf3_a3m…], [protenix_paired_a3m…], [protenix_unpaired_a3m…], [boltz_csv…])`
  (or a per-chain map keyed by chain id — cleaner than positional lists).
- **AF2** ignores these; its multimer path pulls the jackhmmer-multimer MSA.

### §4.3 Generator edits (small, additive)

- `bin/make_protenix_input.py`: accept per-chain `--paired-a3m` + `--unpaired-a3m`
  (chain-id-keyed, or lists in record order); drop the `len(sequences) == 1` gate;
  set both `pairedMsaPath`/`unpairedMsaPath` per chain.
- `bin/make_rf3_fold_spec.py`: accept per-chain `--a3m` (chain-id-keyed list); drop
  the `len(sequences) == 1` gate; attach the (TaxID-annotated) `msa_path` per
  component.
- Boltz: add a generic N-chain YAML path (new `bin/make_boltz_complex_yaml.py`,
  cleaner than overloading the target/binder `bin/create_boltz_yaml.py`) — one
  `protein:` entry per record, each `msa:` pointing at that chain's **`.csv`** (or
  omit `msa:` when `--use_msa_server`). Leave `boltz_pulldown`'s current usage
  untouched.
- `modules/fold/{boltz,rf3,protenix}` input generators: stage and pass the
  per-chain rendered MSA files.

---

## §5 File-by-file change list

**Edited**
- `fold.nf` — remove monomer guard (§1); set real `meta.n_chains`; multimer
  validation (≤26 chains, non-empty seqs); auto `af2_model_preset=multimer` +
  `uniprot/` presence check when `n_chains > 1`; help-text + docstring updates.
- `subworkflows/local/fold_msa.nf` — add `PAIR_MSA` (multimer); emit per-chain
  paired/unpaired bundles; keep monomer emits identical.
- `subworkflows/local/alphafold2.nf` — multimer preset + 2021-snapshot DB flags
  (`--uniprot_database_path`, `--pdb_seqres_database_path`); pass
  `--num_multimer_predictions_per_model`.
- `modules/fold/af2/alphafold2*.nf` — multimer DB flags + preset plumbing.
- `bin/make_protenix_input.py`, `bin/make_rf3_fold_spec.py` — per-chain MSA args,
  drop monomer gate (§4.3).
- `bin/colabfold_remote_msa.py` — per-chain search for multi-record FASTA; (later)
  `ticket/pair` depth booster (§4.1 optional).
- `modules/local/common/mmseqs_colabfoldsearch.nf` — (later) local `--use-env-pairing`
  depth booster, gated on multimer + pairing DB.
- `modules/fold/{boltz,rf3,protenix}/*generate*/*yaml*.nf` — stage per-chain MSAs.
- `docs/docs/workflows/fold.md` — multimer section (per-tool notes, DB
  requirements, pairing, homo-oligomer = repeated records).
- `examples/fold/` — add `input/complex.fasta` (2-chain test) + a multimer run
  script; note the 2021-snapshot bind mount for AF2 multimer.

**New**
- `bin/msa_taxonomy.py` — canonical header→(seq, tax_id, species_mnemonic) parser +
  per-tool renderers (RF3 `TaxID=` a3m, Protenix mnemonic a3m, Boltz `key,sequence`
  CSV). The one place taxonomy logic lives.
- `PAIR_MSA`/`ANNOTATE_MSA` process/subworkflow wrapping `bin/msa_taxonomy.py`.
- `bin/make_boltz_complex_yaml.py` — generic N-chain Boltz YAML.

---

## §6 Verification

### §6.0 Probes — DONE (2026-07-23)

All source/DB probes are complete; results in **§0**. Remaining verification is the
empirical end-to-end below (§6.1) — chiefly that header-derived taxonomy yields
non-trivial paired depth per chain, and that the 2021 snapshot drives AF2 multimer.
One still-empirical unknown carried into §6.1: whether the `uniclust30`-vintage
HHblits DB passed via `--uniref30_database_path` is accepted by this container's
multimer pipeline (format is compatible; vintage is old).

### §6.1 End-to-end (incremental, on a small 2-chain complex)
1. **Monomer regression** — rerun a monomer input; byte-for-byte-equivalent
   behaviour to today (guard removal must not perturb the monomer path). Confirm
   `-resume` still caches.
2. **AF2 multimer** — `complex.fasta`, `--methods af2 --msa_method jackhmmer_af2`
   against the 2021 snapshot. Sanity-check interface pTM/ipTM and PAE; confirm the
   multimer model actually loaded.
3. **Protenix multimer** — paired+unpaired a3m from `PAIR_MSA`; check interface
   confidence and that Protenix skipped its own MSA search.
4. **Boltz multimer** — per-chain `key,sequence` CSV; confirm the interface
   confidence reflects pairing (compare vs `--use_msa_server` as a cross-check).
5. **RF3 multimer** — per-chain `TaxID=`-annotated a3m; confirm chains pair
   (interface ipTM/PAE non-degenerate).
6. **Combined** `--methods af2,boltz,rf3,protenix` smoke test on `complex.fasta`.
7. **EnGens on complex** (§7) — non-blocking; note pass/fail.

---

## §7 EnGens on complexes (explicitly low priority)

Feed the complex `.cif`/`.pdb` straight into the existing `ENGENS_CLUSTER`
subworkflow and see what happens. EnGens' featurization may assume a single chain
(e.g. residue indexing, contact maps). **Do not invest in multimer-aware EnGens
this round.** If it crashes, apply the smallest possible fix (e.g. concatenate
chains / restrict to CA of all chains) or, failing that, auto-skip EnGens for
`n_chains > 1` with a `log.warn`. Correctness of the clustering on complexes is
out of scope.

---

## §8 Risks / open items
1. **Header taxonomy coverage** (now the top risk) — pairing quality depends on how
   many hit headers carry resolvable taxonomy. Real headers are a *mix* of
   `TaxID=`-bearing and mnemonic-bearing lines (§0), so RF3 and Protenix will each
   pair a *subset*. `bin/msa_taxonomy.py` must squeeze every cue (TaxID=, mnemonic,
   `OS=`/`Tax=` name→taxid table) and `log()` per-chain paired depth. If depth is
   too shallow, add the `ticket/pair` / `--use-env-pairing` depth booster (§4.1).
2. **AF2 2021 snapshot compatibility** — `uniclust30` (not `uniref30`), older
   uniref90/mgnify/bfd. DBs exist and flags accept the paths (§0a); the remaining
   unknown is whether the multimer pipeline runs cleanly on the older
   HHblits/uniclust30 vintage — verify in §6.1 step 2.
3. **`msa_taxonomy.py` correctness** — the single source of truth for pairing;
   unit-test the header regexes against real jackhmmer *and* ColabFold a3m samples
   (they differ), and validate that each engine actually pairs (non-degenerate
   interface score) rather than silently running unpaired.
4. **Homo-oligomer ergonomics** — repeated records is clunky; a `count:` shorthand
   is deferred. Document clearly.
5. **EnGens** — may not handle multi-chain; low priority per user (§7).
6. **Scope creep** — no ligands/nucleic acids; protein complexes only this round.
7. **Not touching `main.nf`/`boltz_pulldown.nf`** (future merge, per prior plan).
