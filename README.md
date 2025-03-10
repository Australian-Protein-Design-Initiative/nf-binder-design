# nf-binder-design

> RFDiffusion -> ProteinMPNN -> AlphaFold2 initial guess

## Setup

Install [Nextflow](https://www.nextflow.io/docs/latest/install.html).

Clone the git repository:

```bash
git clone https://github.com/Australian-Protein-Design-Initiative/nf-binder-design
```

> We use Apptainer containers by default. If the https://github.com/Australian-Protein-Design-Initiative/containers repo is not yet public, you'll need to authenticate - see section below on authenticating with the Github Package Registry.

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

### Summarize af2_initial_guess scores

This happens by default when the pipeline successfully completes, however you may want to run it mid-run to monitor how things are going:

```bash
OUTDIR=results

bin/af2_combine_scores.py -o $OUTDIR/combined_scores.tsv -p $OUTDIR/af2_results
```

## Partial diffusion on binder designs

> NOTE: It seems with output from previous designs that the binder is always named chain A, and your other chains are named B, C, etc - irrespective of the chain ID in the original target PDB file. Residue numbering is 1 to N, sequential irrespective of gaps in the chain, rather than original target chain numbering.

```bash
# Generate 10 partial designs for each binder, in batches of 5
# Note the 'single quotes' around the '*.pdb' glob pattern !
nextflow run partial.nf  \
    --input_pdb 'my_designs/*.pdb' \
    --rfd_n_partial_per_binder=10 \
    --rfd_batch_size=5 \
    -profile local

# --require_gpu=false \
# --binder_chain='A' \
# --target_contigs='B90-550' \
```

## Authenticating with a private Github container registry

Go to https://github.com/settings/tokens and generate a secret Personal Access Token (classic) with `read:packages` permissions. 

Authenticate using you Github username and that token as the password:

```
# Old apptainer, pre-1.3
apptainer remote login -u mygithubusername docker://ghcr.io/australian-protein-design-initiative/containers

# New apptainer, >=1.3
# apptainer registry login -u mygithubusername docker://ghcr.io/australian-protein-design-initiative/containers
```