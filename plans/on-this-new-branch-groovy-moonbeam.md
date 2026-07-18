# AlphaFold2 (`af2.nf`) workflow

## Context

We want a first-class, standalone AlphaFold2 prediction workflow in nf-binder-design
so we can fold arbitrary sequences (monomers and multimers) with the Monash custom
AlphaFold2 CUDA-12 container. Today AF2 only appears as the `AF2_INITIAL_GUESS` scoring
module buried inside the RFdiffusion workflow — there is no way to just predict a
structure from a FASTA. The two prediction stages (jackhmmer MSA generation, then GPU
structure prediction) are built as **reusable, independent processes** so a future
`boltz_pulldown` / `boltz_refolding` subworkflow can adopt AF2 as an alternative
folding/scoring engine alongside Boltz-2.

The reference for how the container is actually driven on M3 lives in `./tmp/`
(`01_msa_cpu_simple.sbatch`, `02_predict_cuda12.sbatch`, `run_alphafold_cuda12`) and in
`/scratch2/yt41/andrewpe/tmp/af2-test/NOTES.md`, which documents the hard-won gotchas
(required `--use_gpu_relax`, all DB flags required on both stages, the `alphafold` unix
group / NFS-16-group problem, and L40S/CUDA-12 compatibility).

### Decisions locked in with the user
- **Standalone entrypoint.** A root-level `af2.nf` with its own `workflow {}` (like
  `boltzgen_filter.nf`), run via `nextflow run af2.nf`. **No `main.nf` edits.** The two
  processes live in reusable module files under `modules/local/af2/` so pipeline
  integration later is a plain `include`.
- **Only natively-supported CLI knobs.** Wire exactly the flags `run_alphafold.py`
  actually exposes (verified via `--helpfull`, step 0). Do **not** build the
  `CUSTOM_CONFIG_PATH` / patched-`config.py` fallback for `--num_recycle` / model-subset.
  If those flags exist natively in this custom build, wire them; if not, omit those knobs
  and document them as unsupported for now.

---

## Step 0 — Verification probes the coder MUST run first (read-only)

The container command and the exact flag set are the only real unknowns. Resolve them
before writing process bodies. Requires the `alphafold` group + the cached `.sif`
(pull once from the URL below). Run inside an M3 interactive session:

```bash
SIF=<cached path to alphafold_cuda12_upstream-c77e5d2_custom-57618c5.sif>

# (a) how AF is invoked under `apptainer exec` (runscript is bypassed by Nextflow)
apptainer inspect --runscript "$SIF"
apptainer exec "$SIF" bash -lc 'ls -l /opt/alphafold/run_alphafold.py; python -c "import alphafold; print(alphafold.__file__)"'

# (b) enumerate the REAL flags — decides which knobs we can wire
apptainer exec "$SIF" python /opt/alphafold/run_alphafold.py --helpfull 2>&1 \
  | grep -iE "recycle|models_to_relax|model_preset|random_seed|num_multimer|generate_msas|precomputed|models_to_use|model_names|num_predict|use_gpu_relax"
```

Expected fallback command line (high confidence, given the wrapper's `--workdir /opt/alphafold/`):
`python /opt/alphafold/run_alphafold.py <flags>`. Do **not** `cd /opt/alphafold`; invoking
by absolute path puts the package on `sys.path` and data files load relative to the
install. Keep the Nextflow task dir as cwd and pass `--output_dir` as an absolute path.

---

## File layout

**New files**
- `af2.nf` — root-level standalone workflow: params defaults, help block, input channel, `MSA → predict` wiring, `params.json` `onComplete` (copy the `paramsToMap` helper from `boltzgen_filter.nf`).
- `modules/local/af2/alphafold2_jackhmmer_msa.nf` — `ALPHAFOLD2_JACKHMMER_MSA` (CPU).
- `modules/local/af2/alphafold2.nf` — `ALPHAFOLD2` (GPU) + precomputed-MSA staging.
- `examples/af2/input/pdl1.fasta`
- `examples/af2/run-m3.sh`
- `examples/af2/nextflow.m3.config`
- `examples/af2/README.md` (brief, matching other examples)

**Edited files**
- `nextflow.config` — add generic resource `withName:` selectors for the two processes (§Config).
- `conf/platforms/m3.config` — add `withName:` `clusterOptions`/partition blocks (§Config).
- `conf/platforms/monash_containers.config` — **optional**: add `withName` overrides. The
  process `container` directive is already the Monash mirror URL, so an override is
  redundant; add it only for the "keep selectors in sync" convention noted in that file's header.

---

## Processes (`modules/local/af2/`)

Both processes hardcode the container directive to the mirror URL (there is no ghcr
equivalent to fall back to):

```groovy
container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/alphafold_cuda12_upstream-c77e5d2_custom-57618c5.sif'
```

**Shared DB-flag helper** (define once at top of `af2.nf`; all flags are required on
*both* stages or arg-parsing fails — see NOTES.md). References `params.alphafold2_db_path`
(default `/mnt/datasets/alphafold/alphafold_20240229`):

```
--data_dir=${d}
--uniref90_database_path=${d}/uniref90/uniref90.fasta
--mgnify_database_path=${d}/mgnify/mgy_clusters_2022_05.fa
--bfd_database_path=${d}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
--uniref30_database_path=${d}/uniref30/UniRef30_2021_03
--pdb70_database_path=${d}/pdb70/pdb70
--template_mmcif_dir=${d}/pdb_mmcif/mmcif_files
--obsolete_pdbs_path=${d}/pdb_mmcif/obsolete.dat
--max_template_date=${params.alphafold2_max_template_date}
--db_preset=${params.alphafold2_db_preset}
--model_preset=${params.alphafold2_model_preset}
```
> Note for multimer: this DB layout (`pdb70/`, no `pdb_seqres`/`uniprot`) is
> monomer-oriented. If `--alphafold2_model_preset multimer` is needed, verify the multimer
> DBs exist and swap to `--pdb_seqres_database_path` + `--uniprot_database_path`; otherwise
> document multimer as unsupported with this DB snapshot.

### `ALPHAFOLD2_JACKHMMER_MSA` (CPU)
- `tag "${meta.id}"`; `publishDir "${params.outdir}/af2/msas", pattern: "${meta.id}/msas/**", mode: 'copy'`.
- **input:** `tuple val(meta), path(fasta)`
- **output:** `tuple val(meta), path(fasta), path("${meta.id}"), emit: msa` (the dir AF writes, named by FASTA stem == `meta.id`).
- **script:** run with `--output_dir=$PWD --generate_msas_only=true --use_precomputed_msas=false --use_gpu_relax=false --models_to_relax=none` + DB flags. AF writes `$PWD/${meta.id}/msas/` and `features.pkl`.

### `ALPHAFOLD2` (GPU)
- `tag "${meta.id}"`; `publishDir "${params.outdir}/af2/predictions", pattern: "${meta.id}/**", mode: 'copy'`.
- **input:** `tuple val(meta), path(fasta), path(msa_dir)` — `msa_dir` is staged with name `== meta.id`.
- **output:** `tuple val(meta), path("out/${meta.id}"), emit: predictions`.
- **GPU guard preamble:** copy the idiom from `modules/local/rfd/af2_initial_guess.nf:19-24`
  (`find_available_gpu.py` + `CUDA_VISIBLE_DEVICES` when `params.gpu_devices` set) plus the
  `params.require_gpu` / `nvidia-smi` fail-fast check from `rfdiffusion.nf`.
- **MSA staging (load-bearing):**
  ```bash
  mkdir -p out
  cp -rL "${msa_dir}" "out/${meta.id}"   # -L dereferences symlink-staged files -> real, writable
  python /opt/alphafold/run_alphafold.py \
      --fasta_paths=${fasta} --output_dir=$PWD/out \
      --use_precomputed_msas=true --generate_msas_only=false --use_gpu_relax=true \
      <knob flags> <DB flags>
  ```
  `--use_precomputed_msas=true` reuses the raw `msas/*` files at `out/${meta.id}/msas/`;
  the FASTA stem (`meta.id`) makes the subdir name line up exactly. `cp -rL` is required
  because `stageInMode='symlink'` stages symlinks — AF must rewrite `features.pkl` and
  prediction outputs into real local files for `-resume` and clean `publishDir`.

**Knob → flag mapping (natively-supported only; confirm each in step 0):**
- `--alphafold2_random_seed` (if set) → `--random_seed=<n>`
- `--alphafold2_model_preset` → `--model_preset=` (in DB helper)
- `--alphafold2_models_to_relax` → `--models_to_relax=` (all|best|none)
- `--alphafold2_num_predictions_per_model` → `--num_multimer_predictions_per_model=` **(multimer only)**
- `--num_recycle` / model-subset: wire **only if** step 0 shows native flags; otherwise omit.

> "Structures per FASTA" nuance to document: monomer presets emit one structure per model
> (5 total for the standard preset — model-subset would be needed to reduce, currently
> deferred); `--num_multimer_predictions_per_model` only applies in multimer mode.

---

## `af2.nf` workflow wiring

Params defaults at top (natively-supported knobs only):
```
params.input = false
params.outdir = 'results'
params.alphafold2_db_path = '/mnt/datasets/alphafold/alphafold_20240229'
params.alphafold2_model_preset = 'monomer_ptm'      // monomer|monomer_ptm|multimer
params.alphafold2_db_preset = 'full_dbs'
params.alphafold2_max_template_date = '2024-01-01'
params.alphafold2_random_seed = false               // set an int to fix
params.alphafold2_num_predictions_per_model = 1      // multimer only
params.alphafold2_models_to_relax = 'best'           // all|best|none
// require_gpu / gpu_devices / gpu_allocation_detect_process_regex already defined in nextflow.config
```

Input channel — **each FASTA file = one prediction unit** (the *opposite* of
`boltz_pulldown.nf`'s per-record `.splitFasta`; multiple records in one file are multimer
chains). Accept a single file, a glob, or a directory:
```groovy
def p = params.input
ch_input = ( file(p).isDirectory() ? Channel.fromPath("${p}/*.{fasta,fa,faa}")
                                    : Channel.fromPath(p) )
    .map { f -> [ [id: f.baseName], f ] }
```
Directory idiom parallels `boltzgen_filter.nf:117`. `meta.id = f.baseName` guarantees the
MSA/predict subdir-name contract.

Wiring + help block (mirror `boltz_pulldown.nf:126-151` for the "no `--input` → print help
and `exit(1)`" pattern):
```groovy
ALPHAFOLD2_JACKHMMER_MSA(ch_input)
ALPHAFOLD2(ALPHAFOLD2_JACKHMMER_MSA.out.msa)
```

---

## Config additions

**`nextflow.config`** — generic resources (mirror `AF2_INITIAL_GUESS` at lines 146-151):
```groovy
withName: ALPHAFOLD2_JACKHMMER_MSA { cpus = 8; memory = '64.GB'; time = 8.hours }
withName: ALPHAFOLD2              { accelerator = 1; cpus = 4; memory = '64.GB'; time = 2.hours }
```

**`conf/platforms/m3.config`** — partitions/clusterOptions (mirror lines 145-190). MSA on
`comp` (do **not** use `--qos=shortq`, it caps at 30 min); predict on `gpu` with **no L40S
exclusion** (this CUDA-12 container requires L40S-class GPUs):
```groovy
withName: ALPHAFOLD2_JACKHMMER_MSA {
    clusterOptions = '--partition=comp' + (params.slurm_account ? " --account=${params.slurm_account} " : '')
    time = 8.hours; memory = '64g'; cpus = 8
}
withName: ALPHAFOLD2 {
    accelerator = 1
    clusterOptions = '--gres=gpu:1 --partition=gpu' + (params.slurm_account ? " --account=${params.slurm_account} " : '')
    time = 2.hours; memory = '64g'; cpus = 4
}
```

**DB bind mounts** — the AF DB lives under `/mnt/datasets/...`, outside `autoMounts`. Put
this in the **example's** `nextflow.m3.config` (keeps the global config lean) and bind the
param path itself so the mount can never drift from the DB flags. `runOptions` is a full
string replacement, so restate the existing m3.config options:
```groovy
params.alphafold2_db_path = '/mnt/datasets/alphafold/alphafold_20240229'
apptainer {
    runOptions = "--nv --cleanenv --env PYTHONNOUSERSITE=1 -B \$HOME -B /scratch2 -B /fs04 -B ${params.alphafold2_db_path} -B /mnt/reference/alphafold"
}
```

---

## Example `examples/af2/`

- **`input/pdl1.fasta`** — human PD-L1. Reuse the sequence already in the repo at
  `examples/pdl1-rfd/input/full/3BIK_B.fasta` (single-record monomer test).
- **`run-m3.sh`** — re-exec under the `alphafold` group, then launch Nextflow (SLURM
  account auto-detect copied verbatim from `examples/pdl1-rfd/run-m3-full.sh`):
  ```bash
  #!/bin/bash
  set -euo pipefail
  # AF2 DBs at /mnt/datasets/alphafold are group=alphafold, mode 750. Re-exec under that
  # group so the sbatch-submitted jobs inherit the GID (see propagation note below).
  if [ "$(id -gn)" != "alphafold" ]; then exec sg alphafold -c "$0 $*"; fi

  PIPELINE_DIR=../..
  DATESTAMP=$(date +%Y%m%d_%H%M%S)
  DEFAULT_SLURM_ACCOUNT=$(sacctmgr --parsable2 show user -s ${USER} | tail -1 | cut -f 2 -d \|)

  nextflow run ${PIPELINE_DIR}/af2.nf \
    -c nextflow.m3.config \
    --slurm_account ${DEFAULT_SLURM_ACCOUNT} \
    --input 'input/pdl1.fasta' \
    --outdir results \
    --alphafold2_model_preset monomer_ptm \
    -profile slurm,m3 -resume \
    -with-report results/logs/report_${DATESTAMP}.html \
    -with-trace results/logs/trace_${DATESTAMP}.txt
  ```
- **`nextflow.m3.config`** — sets `params.alphafold2_db_path` and the DB bind mounts (above).

### `alphafold` group propagation (why `sg` in run-m3.sh)
`sg alphafold -c` gives the shell primary-GID `alphafold`; Nextflow's child `sbatch`
inherits it, and `slurmstepd` `setgid()`s the job to that GID on the compute node, so the
`apptainer exec` reading `/mnt/datasets/alphafold` (group `alphafold`, mode 750) succeeds.
The user is also already a supplementary member of `alphafold`, which often suffices via
`initgroups`, but `sg` makes it site-independent. **Verify once** with a probe job before a
full run: `sbatch --partition=comp --wrap='id; apptainer exec -B /mnt/datasets $SIF ls /mnt/datasets/alphafold/alphafold_20240229/uniref90/'`
and confirm `id` shows `gid=…(alphafold)` and the `ls` works. Fallback if propagation
fails: add `--gid=alphafold` to the two AF processes' `clusterOptions`.

---

## Verification (end-to-end)

1. **Step 0 probes** (above) — confirm the invocation command and finalize the flag set.
2. **Group probe job** — confirm GID propagation into the container (above).
3. **Syntax:** `nextflow run af2.nf` with no `--input` prints the help block and exits 1.
4. **Config resolves:** `nextflow config -profile slurm,m3 af2.nf` shows the AF2 `withName`
   selectors, resources, and the DB bind mounts in `apptainer.runOptions`.
5. **Full PD-L1 run** from `examples/af2/`: `./run-m3.sh`. Expect:
   - `ALPHAFOLD2_JACKHMMER_MSA` completes on `comp`, producing `results/af2/msas/pdl1/msas/*`.
   - `ALPHAFOLD2` runs on a GPU node (incl. L40S), reuses the MSAs, and writes
     `results/af2/predictions/pdl1/` with `ranked_0.pdb`, `ranking_debug.json`,
     `result_model_*.pkl`, and (with `monomer_ptm`) PAE/pLDDT. Sanity-check the pLDDT range.
   - `results/params.json` is written on completion.

## Risks / open items
1. **Container command** (step 0a) — `python /opt/alphafold/run_alphafold.py` is the
   high-confidence fallback; confirm no `cd` needed.
2. **`--num_recycle` / model-subset** (step 0b) — wire only if native flags exist; else omit
   (no config.py fallback, per decision).
3. **Multimer DBs** — the given snapshot is monomer-oriented; verify before advertising multimer.
4. **`cp -rL` staging** — the load-bearing assumption; validate AF finds `out/${meta.id}/msas/`
   end-to-end on the PD-L1 run.
5. **GID propagation** — verify with the probe job before trusting a full run.
