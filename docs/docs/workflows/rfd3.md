## RFDiffusion3 Workflow (`--method rfd3`)

RFdiffusion3-based workflow for _de novo_ protein binder design using the [RostettaCommons/foundry](https://github.com/RosettaCommons/foundry) RFDiffusion3 (`rfd3`), MPNN (`mpnn` for ProteinMPNN, SolubleMPNN or LigandMPNN), and RosettaFold3 (`rf3`).

See also: [Official RFDiffusion3 documentation](https://rosettacommons.github.io/foundry/models/rfd3/protein_binder_design.html).

### Overview

The RFdiffusion3 workflow currently includes:

- **`--method rfd3`**: Complete protein binder design pipeline  
  (RFDiffusion3 → ProteinMPNN → RosettaFold3 structure prediction).

> Partial diffusion refinement (as in `--method rfd_partial`) is **not yet implemented** for RFDiffusion3. This may be added in a future version.

### General Information

#### Command-line Options

You can see the available options for the workflow with:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd3 --help
```

#### Parameters File

As with other workflows, parameter command-line options (those prefixed with `--`) can also be defined in a `params.json` file:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd3 \
  -params-file params.json
```

Example `params.json`:

```json
{
  "contigs": "A17-131,/0,50-120",
  "hotspot_res": "A56,A115,A123",
  "rfd3_n_designs": 10
}
```

Parameter names typically mirror the equivalent options in the underlying tools, prefixed with `rfd3_` for RFdiffusion3 and `mpnn_` for ProteinMPNN.

### Key Inputs and Modes

The RFDiffusion3 workflow supports **two ways** of specifying the model configuration:

1. **Config mode (JSON/YAML file)** – recommended for advanced use  
2. **Params mode (command-line options)** – mostly compatible with the RFdiffusion v1 (`rfd`) workflow

#### 1. Config Mode (`--rfd3_config`)

In this mode you provide a full RFDiffusion3 input file (JSON or YAML) that follows the Foundry documentation  
[`RFdiffusion3 — Protein binder design examples`](https://rosettacommons.github.io/foundry/models/rfd3/protein_binder_design.html).

Example (single-target JSON, adapted from the Foundry docs):

```json
{
  "pdl1": {
      "dialect": 2,
      "infer_ori_strategy": "hotspots",
      "input": "input/PDL1.pdb",
      "contig": "50-120,/0,A18-132",
      "select_hotspots": {
          "A56": "CG,OH"
      },
      "is_non_loopy": true
  }
}
```

You can then run:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd3 \
  --input_pdb input/pdl1.pdb \
  --rfd3_config pdl1_rfd3.json \
  --outdir results
```

In this mode, the workflow passes your config directly to RFDiffusion3; some sampler options defined outside the config file (e.g. `inference_sampler.step_scale`, `inference_sampler.gamma_0`) are provided as commandline options (eg `--rfd3_step_scale`, `--rfd3_gamma_0`), the rest can be passed directly via `--rfd3_extra_args`.

#### 2. Params Mode (`--contigs`, `--hotspot_res`, etc.)

In params mode, the workflow builds a minimal RFDiffusion3 JSON config for you internally. This is convenient for quick runs and for users familiar with the RFdiffusion v1 (`rfd`) workflow.

Minimal example:

```bash
OUTDIR=results
mkdir -p "$OUTDIR/logs"

nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd3 \
  --input_pdb input/pdl1.pdb \
  --outdir "$OUTDIR" \
  --design_name pdl1_rfd3 \
  --contigs "A17-131,/0,50-120" \
  --hotspot_res "A56,A115,A123" \
  --rfd3_n_designs 10 \
  --rfd3_batch_size 1 \
  -with-report "$OUTDIR/logs/report_$(date +%Y%m%d_%H%M%S).html" \
  -with-trace "$OUTDIR/logs/trace_$(date +%Y%m%d_%H%M%S).txt" \
  -profile local
```

Key behaviour in params mode:

- **Contigs (`--contigs`)**  
  - Accepts both RFdiffusion **v1-style** contigs (e.g. `[A18-132/0 65-120]`) and RFDiffusion3 **v3-style** contigs (e.g. `A18-132,/0,65-120`).  
  - The workflow automatically normalises v1-style strings to v3 format before building the config.  
  - **RFDiffusion3 v3-style contigs are recommended** for new projects.

- **Hotspots (`--hotspot_res`)**  
  - Accepts a comma-separated list of residue IDs, e.g. `--hotspot_res "A56,A115,A123"`.  
  - The internal config generator converts these to `select_hotspots` entries with `"ALL"` atoms for each residue (all atoms of each residue are targeted).  
  - In many cases it is preferable to target **specific atoms** (e.g. `CG,OH` for a tyrosine). To do this with RFD3, you should instead use a full JSON/YAML config with an explicit `select_hotspots` dictionary, as described in the official RFD3 documentation:  
    [`RFdiffusion3 — Protein binder design examples`](https://rosettacommons.github.io/foundry/models/rfd3/protein_binder_design.html).

### Key Parameters

Below is a summary of the most important parameters for `--method rfd3`.

#### Core RFdiffusion3 Parameters

- **`--input_pdb`**  
  Target PDB/CIF file for binder design (cropped to the chains and residues of interest).

- **`--rfd3_config`**  
  Path to a RFDiffusion3 JSON/YAML config file (config mode).  
  When set, this overrides params-mode config generation.

- **`--contigs`** (params mode)  
  Contig definition for RFDiffusion3, specifying the binder length range and target residues.  
  - Example (v3 style): `--contigs "A17-131,/0,50-120"`.  
  - Example (v1 style): `--contigs "[A18-132/0 65-120]"` (auto-translated to v3).

- **`--hotspot_res`** (params mode)  
  Hotspot residues as a comma-separated list, e.g. `--hotspot_res "A56,A115,A123"`.  
  This is converted to residue-level `select_hotspots` with `"ALL"` atoms per residue. For atom-level hotspot control, use a JSON/YAML config with `select_hotspots`.

- **`--rfd3_n_designs`**  
  Total number of RFDiffusion3 backbones to generate.

- **`--rfd3_batch_size`**  
  Number of designs per RFD3 batch (`diffusion_batch_size`).

- **`--rfd3_step_scale`**  
  Sets `inference_sampler.step_scale`. Default is `3`, following the Foundry recommendation for improved PPI designability.

- **`--rfd3_gamma_0`**  
  Sets `inference_sampler.gamma_0`. Default is `0.2`, again matching the Foundry recommendation.

- **`--rfd3_allow_loopy`**  
  When `false` (default), the generated config enables `is_non_loopy: true`, encouraging more structured binders with fewer loops.  
  When `true`, `is_non_loopy` conditioning is disabled and the model is allowed to make loopier designs.

- **`--rfd3_extra_args`**  
  Additional CLI arguments passed through to the RFDiffusion3 CLI (e.g. extra sampler overrides).

- **`--rfd3_filters`**  
  Semicolon-separated list of filters applied to RFD3 backbone structures **immediately after** RFDiffusion3, before MPNN. Uses the same filter system as the RFdiffusion v1 workflow (e.g. `bin/filter_designs.py` and `bin/filters.d/`). Example: `--rfd3_filters "rg<25"` to keep only binders with a radius of gyration (Rg) below 25 Å.

#### ProteinMPNN Parameters

ProteinMPNN is run on each backbone produced by RFD3. The workflow supports both **modern `mpnn_` options** and a set of **backwards-compatible `pmpnn_` options**:

> For consistency with the previous RFdiffusion v1 (`rfd`) workflow, some legacy `pmpnn_` options are supported and applied as the equivalent `mpnn_` option. **We encourage you to use only the `--mpnn_` options**, since the backward compatible `--pmpnn_` options in RFD3 may be removed in a later version.

Common options:

- **`--mpnn_model_type`**  
  Model type for ProteinMPNN (e.g. `protein_mpnn`). **Default:** `protein_mpnn`.

- **`--mpnn_legacy_weights`**  
  Whether to use the legacy ProteinMPNN weights format. **Default:** `true`.

- **`--mpnn_designed_chains`**  
  Chain IDs to redesign. **Default:** `A`.

- **`--mpnn_batch_size`** (or the legacy **`--pmpnn_seqs_per_struct`**)  
  Number of sequences to sample per backbone. **Default:** 1 sequence per backbone (via `--pmpnn_seqs_per_struct=1` when `--mpnn_batch_size` is not set).

- **`--mpnn_temperature`** (or the legacy **`--pmpnn_temperature`**)  
  Sampling temperature for ProteinMPNN. **Default:** `0.1`.

- **`--mpnn_structure_noise`** (or the legacy **`--pmpnn_augment_eps`**)  
  Amount of structural noise used for augmentation. **Default:** `0` (no structure noise).

- **`--mpnn_omit`** (or the legacy **`--pmpnn_omit_aas`**)  
  Residues to omit from the sequence design; 1-letter codes (e.g. `CX`) are automatically mapped to appropriate 3-letter lists for MPNN. **Default:** `CX` (omit cysteine and unknown).

- **`--mpnn_checkpoint_path`** (or the legacy **`--pmpnn_weights`**)  
  Path to custom ProteinMPNN weights; falls back to the default container weights if not set. **Default:** use the built-in ProteinMPNN checkpoint in the container.

#### RosettaFold3 Parameters

- **`--rf3_ckpt_path`**  
  Path to the RF3 checkpoint used for structure prediction of designed binders (and complexes). A sensible default path inside the container is set in the workflow, but it can be overridden if needed.

#### Refolding Parameters

The workflow supports an optional secondary refolding step (e.g. using Boltz-2) to validate the ROSETTAFOLD3 predicted structures, optionally using a non-truncated target sequence, with an optional target MSA and/or target template structure.

- **`--refold_with`**  
  Tool to use for refolding (e.g. `boltz`). **Default:** disabled.

- **`--refold_max`**  
  Maximum number of top scoring ROSETTAFOLD3 designs to pass to the refolding step.

- **`--refold_filter_sort`**  
  The metric used to sort the ROSETTAFOLD3 outputs before prioritizing them for refolding via `--refold_max`.  
  Use a `-` prefix to sort in descending order (e.g. `-plddt`).  
  **Default:** `pair_pae_min` (ascending).

#### RMSD (design vs refold)

After RosettaFold3, the workflow compares the C-alpha RMSD of each design (MPNN output) to the refolded structure (RF3 output).

### Key Outputs

By default outputs are written under `results/rfd3/` (or the directory specified by `--outdir`). The key outputs of the `--method rfd3` workflow include:

- **RFDiffusion3 backbones**  
  `rfdiffusion3/output/` - contains the CIF files and JSON-format metrics for each generated backbone design.

- **ProteinMPNN sequences and structures**  
  - `mpnn/output/` - contained the designed sequences for each backbone as FASTA files and structure files (CIF)
 with these sidechains.

- **RosettaFold3 predictions**  
  - `rosettafold3/output/` - contains the RosettaFold3-predicted structures for each designed binder+target.

- **RMSD (design vs refold)**  
  - `rfd3/rosettafold3/rmsd/` - TSV files comparing each design (MPNN CIF) to its refolded structure (RF3 CIF):  
    `rmsd_target_aligned_binder.tsv`, `rmsd_complex.tsv`, `rmsd_binder_aligned_binder.tsv`, `rmsd_target_aligned_target.tsv`.  
  Per-design TSVs are also published in the same directory.

- **`combined_scores.tsv`**  
  A single TSV under `rfd3/` merging RFDiffusion3 and RosettaFold3 per-design scores (e.g. `ranking_score`, `iptm`, `plddt`, `pair_pae_min`, `ptm_binder`), sorted by `pair_pae_min`. The RMSD values (initial design vs refold) are also included as `refold_rmsd_*` columns.

The exact directory layout follows the `modules/local/rfd3/` processes (`rfdiffusion3`, `mpnn`, `rosettafold3`); you can examine those modules or a completed run for the precise structure.

### Examples

For complete, runnable examples using RFDiffusion3, see:

- `examples/pdl1-rfd3/` – protein binder design with `rfd3`, both cli params and JSON config examples

For context on the original RFdiffusion v1 workflows (and how RFD3 compares conceptually), see `docs/docs/workflows/rfdiffusion.md`.

