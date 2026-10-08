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

Every engine except ESMFold2 passes through "Template matching". The line order
(esmfold2 first, then rf3, protenix, af2, af3, boltz, openfold3) keeps ESMFold2 on
the top trunk; other orders crowd or overlap the bypass around that station, and
nf-metro 2.1.0 aborts on some of them.

Optional steps are marked two ways: the station label is wrapped in parentheses,
and `%%metro marker: <station> | square, open` draws it as a sharp-cornered square
instead of the usual rounded pill, with a matching `%%metro marker_legend:` row
explaining the glyph. `(Taxonomy pairing)` in the pulldown map uses this -- it only
runs for the engines that need paired MSAs, and only when `--create_target_msa` /
`--create_binder_msa` are enabled (both default to `false`).

Marker shapes are `circle`, `square` or `pill`; the fill is `open`, `solid` or a
literal colour. On a multi-line station `circle` and the default glyph are the same
rounded pill, and `open` differs from `solid` by too little to see, so `square` is
the one that actually reads as different here.

## Rendering

`hook.py` is registered under `hooks:` in `mkdocs.yml`, so `mkdocs build` and
`mkdocs serve` re-render any `.mmd` whose outputs are older than it. Editing a
source and rebuilding the docs is enough; the rendered SVG and PNG stay committed
so GitHub can show them when the workflow pages are read in the repository.

nf-metro is pinned in `docs/requirements.txt` alongside mkdocs:

```bash
pip install -r docs/requirements.txt
cd docs && mkdocs serve
```

Without nf-metro installed the hook logs a warning and the build falls back to the
committed images, so the docs still build.

To render by hand, or to force a re-render when the mtimes say otherwise:

```bash
for n in fold fold_pulldown; do
    nf-metro render "docs/metro/${n}_metro_map.mmd" --mode light \
        -o "docs/docs/images/${n}_metro_map.svg"
    nf-metro render "docs/metro/${n}_metro_map.mmd" --mode light --raster-width 2265 \
        -o "docs/docs/images/${n}_metro_map.png"
    sed -i -e '$a\' "docs/docs/images/${n}_metro_map.svg"
done
```

Use `nf-metro validate <file>.mmd` to check a map, and `--debug` to see ports and
grid lines while iterating.
