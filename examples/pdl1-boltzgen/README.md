# BoltzGen PDL1 binder design with FoldSeek

Designs protein binders against PDL1 (chain A, residues 18–132) with hotspot at residue 56 (native numbering).
Note that because BoltzGen uses 1-based residue index numbering rather than native PDB numbering, we use `PDL1.renum.pdb` and the hotspot residue number is 39, so that the output complex numbering matches the input structure numbering.
After filtering, runs FoldSeek structural search against CATH50 on the final design chains. The CATH50 foldseek database is downloaded automatically and saved to `databases/foldseek/CATH50/` for reuse across runs.

To run locally on a single GPU:

```bash
./run-local.sh
```
