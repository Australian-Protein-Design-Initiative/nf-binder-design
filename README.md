# nf-binder-design

> RFDiffusion -> ProteinMPNN -> AlphaFold2 initial guess

## Setup

Install [Nextflow](https://www.nextflow.io/docs/latest/install.html).

Clone the git repository:

```bash
git clone https://github.com/Australian-Protein-Design-Initiative/nf-binder-design
```

> We use Apptainer containers by default. If the https://github.com/Australian-Protein-Design-Initiative/containers repo is not yet public, you'll need to provide your own containers or software installation.

## Binder design

Example:

```bash
OUTDIR=results
mkdir -p $OUTDIR/logs

nextflow run main.nf \
    --input_pdb target.pdb \
    --outdir $OUTDIR \
    --contigs "[A371-508/A753-883/A946-1118/A1135-1153/0 70-100]" \
    --hotspot_res "[A473,A995,A411,A421]" \
    --rfd_n_designs=10 \
    --rfdiffusion_batch_size 1 \
    -with-report $OUTDIR/logs/report_$(date +%Y%m%d_%H%M%S).html \
    -with-trace $OUTDIR/logs/trace_$(date +%Y%m%d_%H%M%S).txt \
    -resume \
    -profile local
```

> If you are working on a specific cluster like M3 or MLeRP, you should omit `-profile local` and add the `-c` flag pointing to the specific platform config, eg `-c conf/platforms/m3.config` for M3.

A more complex example, as a wrapper script for M3, using a the site-specific config for M3 (`-c`), a specific RFDiffusion model (``), a custom ProteinMPNN weights (`--pmpnn_weigths`) and radius of gyration potentials (`--rfd_extra_args`):

```bash
#!/bin/bash
WF_PATH="/some/path/to/nf-binder-design" # change this

mkdir -p results/logs
DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow \
-c ${WF_PATH}/conf/platforms/m3.config run \
${WF_PATH}/main.nf  \
--slurm_account=ab12 \
--input_pdb 'input/target_cropped.pdb' \
--design_name my-binder \
--outdir results \
--contigs "[B346-521/B601-696/B786-856/0 70-130]" \
--hotspot_res "[B472,B476,B484,B488]" \
--rfd_n_designs=1000 \
--rfd_batch_size=5 \
--pmpnn_seqs_per_struct=2 \
--pmpnn_weigths="/models/HyperMPNN/retrained_models/v48_020_epoch300_hyper.pt" \
--rfd_model_path="/models/rfdiffusion/Complex_beta_ckpt.pt" \
--rfd_extra_args='potentials.guiding_potentials=[\"type:binder_ROG,weight:7,min_dist:10\"] potentials.guide_decay="quadratic"' \
-resume \
-with-report results/logs/report_${DATESTAMP}.html \
-with-trace results/logs/trace_${DATESTAMP}.txt
```

> Note: Ensure you set --slurm_account to your M3 account/project ID.

### Summarize af2_initial_guess scores

This happens by default when the pipeline successfully completes, however you may want to run it mid-run to monitor how things are going:

```bash
OUTDIR=results

bin/af2_combine_scores.py -o $OUTDIR/combined_scores.tsv -p $OUTDIR/af2_results

## Partial diffusion on binder designs

> NOTE: It seems with output from previous designs that the binder is always named chain A, and your other chains are named B, C, etc - irrespective of the chain ID in the original target PDB file. Residue numbering is 1 to N, sequential irrespective of gaps in the chain, rather than original target chain numbering.

```bash
OUTDIR=results
mkdir -p $OUTDIR/logs

# Generate 10 partial designs for each binder, in batches of 5
# Note the 'single quotes' around the '*.pdb' glob pattern !
nextflow run partial.nf  \
    --input_pdb 'my_designs/*.pdb' \
    --rfd_n_partial_per_binder=10 \
    --rfd_batch_size=5 \
    --hotspot_res "[A473,A995,A411,A421]" \
    --rfd_partial_T=2,5,10,20 \
    -with-report $OUTDIR/logs/report_$(date +%Y%m%d_%H%M%S).html \
    -with-trace $OUTDIR/logs/trace_$(date +%Y%m%d_%H%M%S).txt \
    -profile local
```

## License

MIT

> Note that some software dependencies of the pipeline are under less permissive licenses - in particular, components that use [Rosetta/PyRosetta](https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md) are only free for Non-Commercial use.