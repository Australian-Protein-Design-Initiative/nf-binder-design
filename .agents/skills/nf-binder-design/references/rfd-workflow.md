# RFdiffusion Workflows — Full Reference

## Table of Contents

- [Method: rfd (Full Pipeline)](#method-rfd-full-pipeline)
- [Method: rfd_partial (Partial Diffusion)](#method-rfd_partial-partial-diffusion)
- [Contig Syntax](#contig-syntax)
- [Hotspot Residues](#hotspot-residues)
- [Design Filters](#design-filters)
- [Refold Scoring](#refold-scoring)
- [RMSD Comparisons](#rmsd-comparisons)
- [Custom Filter Plugins](#custom-filter-plugins)
- [Complete Examples](#complete-examples)

---

## Method: rfd (Full Pipeline)

**Pipeline**: RFdiffusion → ProteinMPNN → AlphaFold2 (initial guess) → (optional) refold scoring

### VRAM Requirements

- **~400–500 total residues** (target + max binder length) is the ceiling for 12 GB VRAM at `--rfd_batch_size 1`. Larger targets need more VRAM or further cropping.
- AF2 initial guess and refold scoring are also VRAM-intensive; they run one at a time per GPU regardless of batch size.


### All Parameters

Use `--method rfd --help` for the definitive list. Key parameters:

#### Input / Output
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input_pdb` | Target PDB file(s). Use single-quoted glob for multiple: `'input/*.pdb'` | Required |
| `--outdir` | Output directory | `results` |
| `--design_name` | Name prefix for outputs | Derived from input |

#### RFdiffusion
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--contigs` | Contig definition string (see syntax below) | Required |
| `--hotspot_res` | Hotspot residues, comma-separated (e.g. `"A56,A125"`) | Optional |
| `--rfd_n_designs` | Total number of backbone designs | Required |
| `--rfd_batch_size` | Designs per parallel batch | 1 |
| `--rfd_model_path` | Path to RFdiffusion model checkpoint (inside container) | Default model |
| `--rfd_extra_args` | Extra arguments passed directly to RFdiffusion | None |
| `--rfd_filters` | Filter expression for backbones (e.g. `"rg<=20"`) | None |

#### ProteinMPNN
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--pmpnn_seqs_per_struct` | Sequences to generate per backbone | 1 |
| `--pmpnn_relax_cycles` | FastRelax cycles | 0 |
| `--pmpnn_weights` | Custom ProteinMPNN weights path (e.g. HyperMPNN) | Default |

#### AF2 Initial Guess
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--af2ig_recycle` | Number of AF2 recycles | 1 |

#### Refold Scoring (optional)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--refold_with` | Comma-separated engines to refold with: `af2`, `af2_mono`, `boltz`, `rf3`, `protenix`, `af3`, `openfold3`, `chai1`, `esmfold2`, `esmfold2_fast` | None (no refolding) |
| `--refold_filter` | Threshold gate on the AF2IG scores, semicolon-separated | None |
| `--refold_sort_by` | Rank by one metric before `--refold_max` takes the best; `-` prefix for descending. Write it with an `=` | None |
| `--refold_max` | Maximum designs to refold. With `--refold_sort_by` the best N, without it the first N to be scored | No limit |
| `--use_msa_server` | Use ColabFold MMseqs2 server for MSA | `false` |
| `--templates` | Known structures, matched to the folded chains by sequence alignment | None |
| `--output_rmsd_aligned` | Output RMSD-aligned structures | `false` |

`--refold_with` turns refolding on; the filter and the cap only select which designs go through it. Refolding is a held-out second opinion **in addition to** AF2 initial guess, not a replacement. See the Refold Scoring section below.

#### SLURM
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--slurm_account` | SLURM account/project ID | None |

---

## Method: rfd_partial (Partial Diffusion)

Refines existing binder designs by applying partial diffusion — denoising from a partially noised version of the input structure.

### Parameters specific to rfd_partial

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input_pdb` | Existing binder-target complex PDBs (glob in single quotes) | Required |
| `--contigs` | Contig definition (same syntax as rfd) | Required |
| `--hotspot_res` | Hotspot residues | Optional |
| `--rfd_n_partial_per_binder` | Number of partial designs per input binder | 10 |
| `--rfd_batch_size` | Designs per batch | 1 |
| `--rfd_partial_T` | Comma-separated noise timesteps (e.g. `"5,15"` or `"2,5,10,20"`) | Required |
| `--rfd_filters` | Filter expression for generated backbones | None |

All `--refold_*` parameters from the `rfd` method are also available.

### Important: Chain ID Renumbering

When applying partial diffusion to designs output from `--method rfd`:
- The **binder is always chain A**, regardless of original target chain IDs
- Other chains are named B, C, etc.
- **Residue numbering is sequential 1 to N**
- Adjust your `--hotspot_res` accordingly

---

## Contig Syntax

Contigs define which parts of the target to use and how long the designed binder should be, using RFdiffusion syntax.

### Format: `"[<target_ranges>/0 <binder_length_range>]"`

- **Target ranges**: `ChainStart-End` segments separated by `/`
- **`/0`**: Variable-length linker (zero-length gap between target and binder)
- **Binder length**: Single number or range like `65-120`

### Examples

```
"[A18-132/0 65-120]"
```
Use residues 18–132 of chain A as target, design a binder of 65–120 residues.

```
"[A371-508/A753-883/A946-1118/A1135-1153/0 70-100]"
```
Use multiple segments of chain A (non-contiguous binding surface), design a binder of 70–100 residues.

```
"[B346-521/B601-696/B786-856/0 70-130]"
```
Use segments from chain B.

### Determining Contigs

- Use `bin/get_contigs.py` to extract contigs from a hand-cropped PDB structure
- Use `bin/trim_to_contigs.py` to trim a structure to specific contig regions
- Visually inspect your target in ChimeraX/PyMOL to identify the binding surface residues

---

## Hotspot Residues

Hotspot residues bias the diffusion process toward specific target residues. Format: `ChainResidue`, comma-separated.

```
--hotspot_res "A56"
--hotspot_res "A56,A125"
--hotspot_res "B472,B476,B484,B488"
```

These correspond to residue numbers in the input PDB.

---

## Design Filters

### Backbone Filters (`--rfd_filters`)

Filter RFdiffusion backbone designs before passing to ProteinMPNN and AF2 initial guess:

```bash
--rfd_filters="rg<=20"
```

Available metrics for backbone filtering: `rg` (radius of gyration), plus any custom filter plugins.

### AF2 Initial Guess Score Filters (`--refold_filter`)

Select AF2 initial guess designs for refolding. Semicolon-separated conditions:

```bash
--refold_filter="pae_interaction<=10;plddt_binder>=80"
```

(`--refold_af2ig_filters` is the old name. It still works and warns.)

Available AF2IG scores:
- `pae_interaction` — PAE between binder and target (lower is better)
- `plddt_binder` — pLDDT of the binder (higher is better)
- `plddt_target`, `plddt_total`
- `pae_binder`, `pae_target`
- `binder_aligned_rmsd`, `target_aligned_rmsd`

Shape/size scores (also available):
- `rg` (radius of gyration), `dmax`, `asphericity`, `approx_rh`

### Filter Syntax

```
metric<=value      # less than or equal
metric>=value      # greater than or equal
metric<value       # less than
metric>value       # greater than
```

Multiple filters separated by `;` (AND logic).

---

## Refold Scoring

When `--refold_with` names at least one engine, the selected designs are re-predicted with it:
- Complex prediction (binder + target)
- Monomer prediction (binder alone), with `--refold_binder_monomer` (default `true`)

This provides an orthogonal validation: if an independent predictor reproduces a similar complex structure, the design is more likely to work experimentally. It runs in addition to AF2 initial guess, which still scores every design.

### Which designs are refolded

Three parameters that compose, all optional; with none set, every design is refolded.

- `--refold_filter 'pae_interaction<=10;plddt_binder>=80'` - threshold gate on the AF2IG scores. Decided per design, so refolding starts while the rest of the set is still being scored.
- `--refold_sort_by='-ipsae'` - rank by one metric, `-` prefix for descending. Cannot stream; the refold waits for every design to be scored. **Write it with an `=`**, otherwise Nextflow parses the leading `-` as its own option.
- `--refold_max 100` - cap. With `--refold_sort_by` the best N, without it the first N to be scored.

`<outdir>/rfd/refold/refold_selection.tsv` records every design's verdict, rank and sort value.

### Refolding Tips

- `--use_msa_server=true` uses the public ColabFold MMseqs2 server for MSA generation (recommended for targets with known homologues)
- `--templates 'input/templates/'` passes known structures, matched to the folded chains by sequence alignment, which helps when the design input was a cropped domain
- Name more than one engine (`--refold_with boltz,chai1`) to get a cross-engine `consensus_z` column

### Renamed and removed

The `--full_refold_*` family is now `--refold_*` throughout. `--refold_af2ig_filters` (use `--refold_filter`), `--refold_oracles` (use `--refold_with`) and `--refold_oracles_max` (use `--refold_max`) still work and warn. These are **removed**, and a run that sets one stops with the replacement named: `--refold_create_target_msa` (use `--create_target_msa`), `--refold_use_msa_server` (use `--use_msa_server`), `--refold_target_templates` (use `--templates`), `--refold_alignment` (use `--create_target_msa` with `--msa_method mmseqs2_colabfold`), `--refold_target_fasta` (no replacement; the target sequence is read from the design structures). `--boltz_refold_batch_size` is also gone.

---

## RMSD Comparisons

When refold scoring is enabled, several C-alpha RMSD comparisons are calculated and merged into `combined_scores.tsv` as engine-prefixed columns (`<tool>` is the engine name, eg `boltz`):

| Column | Superimpose On | Measure | Interpretation |
|--------|----------------|---------|----------------|
| `<tool>_target_aligned_binder_rmsd_all` | Target | Binder | Binding pose deviation after refolding |
| `<tool>_complex_rmsd_all` | Both chains | Both chains | Overall structural agreement, design vs refold |
| `<tool>_monomer_vs_design_rmsd_all` | Binder | Binder | Binder fold change, monomer vs the design |
| `<tool>_monomer_vs_complex_rmsd_all` | Binder | Binder | Binder fold change, monomer vs the refolded complex |

**Key quality metrics**:
- `<tool>_target_aligned_binder_rmsd_all` < ~3.5 Å → binder maintains its binding pose
- `<tool>_monomer_vs_complex_rmsd_all` < ~3.5 Å → binder structure is stable whether bound or unbound

The confidence columns alongside them are `<tool>_iptm_max`, `<tool>_ipsae_max`, `<tool>_plddt_max` and the rest of `--refold_metrics`, plus `<tool>_interface_dG` and the other BindCraft interface metrics on the refolded complex, and `consensus_z` across engines.

---

## Custom Filter Plugins

Create a `.py` file in `bin/filters.d/` with two functions:

```python
def register_metrics() -> list[str]:
    return ["my_metric"]

def calculate_metrics(pdb_files: list[str], binder_chains: list[str]) -> pd.DataFrame:
    # Return DataFrame indexed by design ID, columns = metric names
    ...
```

The `filter_designs.py` script automatically discovers plugins in `bin/filters.d/`.

---

## Complete Examples

### Cropped Target with Full-Length Refolding

When designing against a cropped target, validate the binder against the **full-length** target to ensure it does not clash with removed domains. `--templates` takes the full-length structures, matched to the folded chains by sequence alignment:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --input_pdb 'input/cropped_target.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs 100 \
  --rfd_filters "rg<=20" \
  --refold_with boltz \
  --refold_filter "pae_interaction<=10;plddt_binder>=80" \
  --templates 'input/templates/' \
  --use_msa_server true \
  -profile local -resume
```

### Minimal Local Run (RFdiffusion)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs=4 \
  -profile local \
  -resume
```

### Full Pipeline with Refold Scoring (Local)

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs=4 \
  --rfd_batch_size=1 \
  --pmpnn_seqs_per_struct=2 \
  --pmpnn_relax_cycles=3 \
  --rfd_filters="rg<=20" \
  --af2ig_recycle=3 \
  --refold_with=boltz \
  --refold_filter="pae_interaction<=10;plddt_binder>=80" \
  --refold_sort_by='pae_interaction' \
  --refold_max=100 \
  --use_msa_server=true \
  --templates='input/full/' \
  --output_rmsd_aligned=true \
  -profile local \
  -resume
```

### HPC Cluster with SLURM

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --slurm_account=ab12 \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs=100 \
  --rfd_batch_size=5 \
  --pmpnn_seqs_per_struct=2 \
  --pmpnn_relax_cycles=3 \
  --rfd_filters="rg<=20" \
  --af2ig_recycle=3 \
  --refold_with=boltz \
  --refold_filter="pae_interaction<=10;plddt_binder>=80" \
  --refold_sort_by='pae_interaction' \
  --refold_max=100 \
  --use_msa_server=true \
  -profile slurm \
  -resume
```

Use `-profile slurm,<site>` for site-specific configs (see `conf/platforms/`).

### Advanced: Custom Model and Potentials

```bash
--rfd_model_path="/models/rfdiffusion/Complex_beta_ckpt.pt" \
--rfd_extra_args='potentials.guiding_potentials=["type:binder_ROG,weight:7,min_dist:10"] potentials.guide_decay="quadratic"' \
--pmpnn_weights="/models/HyperMPNN/retrained_models/v48_020_epoch300_hyper.pt"
```

Model paths refer to paths **inside the container**.

### Partial Diffusion

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd_partial \
  --input_pdb 'input/*.pdb' \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_partial_per_binder=1 \
  --rfd_batch_size=1 \
  --rfd_partial_T "5,15" \
  --pmpnn_seqs_per_struct=1 \
  --pmpnn_relax_cycles=3 \
  --rfd_filters="rg<=20" \
  --af2ig_recycle=3 \
  --refold_with=boltz \
  --refold_filter="pae_interaction<=10;plddt_binder>=80" \
  --refold_sort_by='pae_interaction' \
  --refold_max=100 \
  --use_msa_server=true \
  -profile local \
  -resume
```
