# RFDiffusion3 Workflow — Full Reference

## Overview

The `--method rfd3` workflow runs the newer RFDiffusion3 model for _de novo_ protein binder design.
**Pipeline**: RFDiffusion3 → MPNN (ProteinMPNN/SolubleMPNN/LigandMPNN) → RosettaFold3 (rf3) structure prediction → (optional) Boltz-2 refolding.

> **Note**: Partial diffusion refinement (`rfd_partial` style) is not yet implemented for RFDiffusion3.

## Execution Modes

There are two ways to configure the `rfd3` workflow:

### 1. Params Mode (Recommended for simple runs)
Uses CLI parameters to build a minimal RFD3 JSON config.

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd3 \
  --input_pdb input/target.pdb \
  --outdir results \
  --contigs "A17-131,/0,50-120" \
  --hotspot_res "A56,A115,A123" \
  --rfd3_n_designs 10 \
  -profile local -resume
```

**Contig Syntax in rfd3**:
- RFD3 uses v3-style contigs: `A17-131,/0,50-120` (comma-separated, target first or binder first).
- The workflow auto-converts v1-style contigs (e.g. `[A18-132/0 65-120]`) to v3-style.
- **Chain IDs in RFD3**: RFD3 assigns chains A, B... based on contig polymer order. E.g., for `C18-132,/0,50-120`, the target becomes A and the binder becomes B. The workflow infers the binder chain automatically with `--mpnn_designed_chains='auto'`.

**Hotspots in Params Mode**:
- `--hotspot_res "A56,A115"` targets *all atoms* of those residues.
- `--rfd3_hotspot_subsample 0.66` (or `--hotspot_subsample`): Subsamples hotspots randomly per batch.

### 2. Config Mode (For advanced runs)

Provide a full RFDiffusion3 JSON or YAML config via `--rfd3_config`. **Use config mode when you need atom-level hotspots, partial flexibility in the target structure, or other [InputSelection](https://rosettacommons.github.io/foundry/models/rfd3/input.html#the-inputselection-mini-language) fields** that params mode cannot express.

See the [PPI design tutorial](https://rosettacommons.github.io/foundry/models/rfd3/tutorials/ppi_design_tutorial.html) and [input specification](https://rosettacommons.github.io/foundry/models/rfd3/input.html) for full RFD3 option reference.

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd3 \
  --input_pdb input/target.pdb \
  --rfd3_config input/rfd3_config.json \
  --outdir results \
  -profile local -resume
```

#### Atom-level hotspots (`select_hotspots`)

In params mode, `--hotspot_res "A56,A115"` targets **all atoms** of each residue. For protein–protein interface design, RFD3 supports **atom-level hotspots** in the config: name both the residue and the specific atoms that should lie within ~4.5 Å of the designed binder ([PPI tutorial](https://rosettacommons.github.io/foundry/models/rfd3/tutorials/ppi_design_tutorial.html)).

Hotspot residues must appear in the `contig` string so RFD3 is aware of them. Pair with `infer_ori_strategy: "hotspots"` to place the ORI token near the hotspot patch.

```json
{
  "target_binder": {
    "dialect": 2,
    "input": "input/target.pdb",
    "contig": "40-120,/0,E6-155",
    "length": "190-270",
    "infer_ori_strategy": "hotspots",
    "is_non_loopy": true,
    "select_hotspots": {
      "E64": "CD2,CZ",
      "E88": "CG,CZ",
      "E96": "CD1,CZ"
    }
  }
}
```

Dictionary values use the [InputSelection mini-language](https://rosettacommons.github.io/foundry/models/rfd3/input.html#the-inputselection-mini-language): comma-separated atom names per residue key (`E64`, `E88`, …).

#### Target flexibility with fixed sequence (`select_fixed_atoms`)

Allow **structural flexibility** on parts of the input target while keeping sequence fixed — useful when the binding site may adopt alternate conformations ([PPI tutorial — Other Useful Settings](https://rosettacommons.github.io/foundry/models/rfd3/tutorials/ppi_design_tutorial.html)).

```json
{
  "target_binder": {
    "dialect": 2,
    "input": "input/target.pdb",
    "contig": "40-120,/0,E6-155",
    "infer_ori_strategy": "hotspots",
    "is_non_loopy": true,
    "select_hotspots": {
      "E64": "CD2,CZ",
      "E88": "CG,CZ"
    },
    "select_fixed_atoms": {
      "E25": [],
      "E26": "BKBN",
      "E27": "CA,CB,OG"
    }
  }
}
```

| Value | Meaning |
|-------|---------|
| `[]` (empty list) | All atoms in that residue are **flexible** (coordinates can move) |
| `"BKBN"` | Backbone fixed (`N,CA,C,O`); side-chain atoms can move |
| `"CA,CB,OG"` | Only named atoms fixed; other atoms in the residue can move |
| `"ALL"` | All atoms fixed |
| `"TIP"` | Common tip atom for the residue type |

Other InputSelection shorthands (`ALL`, `TIP`, contig strings) are documented in the [foundry input spec](https://rosettacommons.github.io/foundry/models/rfd3/input.html#select-fixed-atoms).

#### Minimal config example (params-equivalent)

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

See `examples/pdl1-rfd3/pdl1_rfd3.json` for a repo example.

## Parameters

Use `--method rfd3 --help` for the definitive list.

### RFDiffusion3 Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--rfd3_n_designs` | Total backbones to generate | Required in Params mode |
| `--rfd3_batch_size` | Designs per RFD3 batch (larger batch = faster overall run, but more VRAM) | 1 |
| `--rfd3_step_scale` | `inference_sampler.step_scale` | `3` (improves PPI designability) |
| `--rfd3_gamma_0` | `inference_sampler.gamma_0` | `0.2` (improves PPI designability) |
| `--rfd3_is_non_loopy` | `true` = more all-helical binders, `false` = more loops/beta sheets | Unset |
| `--rfd3_filters` | Filter expression applied immediately after RFD3 (e.g. `"rg<25"`) | None |
| `--rfd3_extra_args` | Additional arguments passed to RFD3 CLI | None |

### MPNN Parameters
**Important**: Use the `--mpnn_` prefix for RFD3 (legacy `--pmpnn_` is supported but discouraged).

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--mpnn_preset` | Checkpoint preset: `vanilla`, `soluble`, or `hyper` (HyperMPNN) | `vanilla` (via default pt) |
| `--mpnn_weights_noise` | Training noise tier (`005`, `010`, `020`, `030`) | `020` (0.2 Å) |
| `--mpnn_batch_size` | Sequences sampled per backbone | 1 |
| `--mpnn_temperature` | Sampling temperature | 0.1 |
| `--mpnn_structure_noise` | Noise added to coordinates at inference (in Å) | 0 |
| `--mpnn_omit` | Amino acids to omit | `CX` |

### RosettaFold3 (rf3) Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--rf3_num_steps` | Diffusion steps for RF3 | 50 (faster than RF3 default of 200) |
| `--rf3_n_recycles` | Structure module recycles | 10 |
| `--rf3_batch_size` | Structures folded per `rf3 fold` invocation | 1 |
| `--rf3_early_stopping_plddt_threshold` | Exit early if pLDDT is below this after 1st recycle | 0.5 |
| `--rf3_create_target_msa` | Generate target MSA via ColabFold | `false` |
| `--rf3_use_msa_server` | Use public API for MSA (needs `--rf3_create_target_msa=true`) | `false` |
| `--rf3_target_template` | Alternative PDB to use as the RF3 target template | None |

### Boltz-2 Full Refolding (Optional)
To validate RF3 predictions with Boltz-2:
```bash
--full_refold_with boltz
--full_refold_max 100
--full_refold_filter_sort "-plddt"
```

## Outputs

Default output directory is `results/rfd3/` (or your `--outdir`).

- `rfdiffusion3/output/` — Backbones (CIF files and JSON metrics).
- `mpnn/output/` — Designed sequences and sidechain-packed structures.
- `rosettafold3/output/` — RF3 predicted structures.
- `rfd3/rosettafold3/rmsd/` — TSV files comparing design (MPNN) to refold (RF3).
- `rfd3/combined_scores.tsv` — The merged summary of RFD3, RF3, and Boltz (if enabled) metrics, sorted by `pair_pae_min`.

## Advanced: Boltz Refolding with MSAs and Templates

When designing against a cropped target, validate binders against the full-length target using Boltz-2 full refold:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd3 \
  --rfd3_config input/rfd3_config.json \
  --outdir results \
  --rfd3_n_designs 100 \
  --rfd3_filters "rg<=20" \
  --rf3_use_msa_server true \
  --rf3_create_target_msa true \
  --full_refold_with 'boltz' \
  --full_refold_filter_sort 'pair_pae' \
  --full_refold_max 10 \
  --full_refold_create_target_msa true \
  --full_refold_use_msa_server true \
  --full_refold_target_templates 'input/templates/' \
  --output_rmsd_aligned true \
  -profile slurm -resume
```

See `examples/pdl1-rfd3/` for a complete working example.
