AlphaFold2 structure prediction with the standalone `af2.nf` workflow.

This example folds human PD-L1 (`input/pdl1.fasta`, reused from `examples/pdl1-rfd/input/full/3BIK_B.fasta`)
as a single-chain monomer with the Monash custom AlphaFold2 CUDA-12 container.
Each FASTA file passed to `af2.nf` is one prediction unit (multiple records in
a single file are treated as multimer chains).

To run on M3 (submits the CPU MSA and GPU prediction stages via SLURM):

```bash
./run-m3.sh
```

To run locally (both stages as local processes, e.g. on a GPU workstation):

```bash
./run-local.sh
```

Results are written under `results/af2/msas/pdl1/` (raw MSAs, from the CPU
jackhmmer/hhblits stage) and `results/af2/predictions/pdl1/` (GPU structure
prediction: `ranked_0.pdb`, `ranking_debug.json`, `result_model_*.pkl`, and
with `monomer_ptm`, per-residue pLDDT / predicted aligned error).
`results/params.json` is written on completion.
