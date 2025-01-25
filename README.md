# nf-binder-design

> RFDiffusion -> ProteinMPNN -> AlphaFold2 initial guess

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
