# nf-binder-design

> RFDiffusion -> ProteinMPNN -> AlphaFold2 initial guess

## Binder design

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

Summarize AlphaFold2 scores:

```bash
bin/af2_combine_scores.py -o $OUTDIR/combined_scores.tsv -p $OUTDIR/af2_results
```

## Partial diffusion on binder designs

> NOTE: It seems with output from previous designs that the binder is always named chain A, and your other chains are named B, C, etc - irrespective of the chain ID in the original target PDB file. Residue numbering is 1 to N, sequential irrespective of gaps in the chain, rather than original target chain numbering.

```bash
nextflow run partial.nf \
  --input_pdb 'my_designs/*.pdb' \
  -profile local 

#  --require_gpu=false \
#  --binder_chain='A' \
#  --target_contigs='B90-550' \
```

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
