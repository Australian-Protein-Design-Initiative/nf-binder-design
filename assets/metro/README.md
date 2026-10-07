# Workflow metro maps

Source `.mmd` files for the metro-map diagrams at the top of the workflow docs,
rendered with [nf-metro](https://github.com/seqeralabs/nf-metro):

| Source | Rendered (docs/docs/images/) | Used in |
|--------|------------------------------|---------|
| `fold_metro_map.mmd` | `fold_metro_map.{svg,png}` | `docs/docs/workflows/fold.md` |
| `fold_pulldown_metro_map.mmd` | `fold_pulldown_metro_map.{svg,png}` | `docs/docs/workflows/fold-pulldown.md` |

Each fold engine is one line. `af2` covers `af2`/`af2_mono` and `esmfold2` covers
`esmfold2`/`esmfold2_fast`. EnGens is drawn after scoring for readability,
although it clusters the predicted structures rather than the score table.

The line order (rf3, protenix, esmfold2, af2, then the template-using af3, boltz,
openfold3) is deliberate: with the template-using lines split up, nf-metro 2.1.0
aborts with overlapping bypass curves around "Template matching".

## Regenerate

From the repository root:

```bash
pip install 'nf-metro==2.1.0'

for n in fold fold_pulldown; do
    nf-metro render "assets/metro/${n}_metro_map.mmd" --mode light \
        -o "docs/docs/images/${n}_metro_map.svg"
    nf-metro render "assets/metro/${n}_metro_map.mmd" --mode light --raster-width 2265 \
        -o "docs/docs/images/${n}_metro_map.png"
    sed -i -e '$a\' "docs/docs/images/${n}_metro_map.svg"
done
```

Use `nf-metro validate <file>.mmd` to check a map, and `--debug` to see ports and
grid lines while iterating.
