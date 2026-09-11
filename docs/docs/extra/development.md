# Development

Notes for contributors and developers working on the pipeline.

## Running the tests

The pipeline ships with a small [nf-test](https://www.nf-test.com/) suite under `tests/` for
process- and workflow-level unit tests.

Install nf-test (one-off):

```bash
curl -fsSL https://code.askimed.com/install/nf-test | bash
# move the resulting `nf-test` binary onto your PATH, eg:
mv nf-test $HOME/.local/bin/
```

Run all tests from the repository root:

```bash
nf-test test
```

Run a single test file, or filter by tag (CPU-only vs GPU, etc.):

```bash
nf-test test tests/modules/local/common/unique_id.nf.test
nf-test test --tag unique_id
# RFdiffusion process tests (GPU, Apptainer, ~1–2 min)
nf-test test --tag rfd
```

See the [Testing](../testing.md) page and the upstream
[nf-test docs](https://www.nf-test.com/docs/getting-started/) for layout, writing tests,
Apptainer cache notes, and more detail.

## Linting

Lint `.nf` and `config` files with `nextflow lint` before opening a PR:

```bash
# Lint the pipeline source (use a Nextflow with the linter, e.g. 25.04.x)
NXF_VER=25.04.7 nextflow lint -o concise main.nf workflows subworkflows modules
```

`nextflow lint` catches genuine issues such as unused variables and parameters
(prefix unused closure parameters with `_` to silence those warnings).

It parses with the strict (v2) syntax, so it also reports constructs that this
pipeline deliberately keeps for Nextflow 24.04.3 compatibility — the conditional
`include` statements in `main.nf`, the `@Field`/`import` in
`modules/local/rfd3/rfd3_utils.nf`, and multi-name function includes (reported as
`... is not defined`). These are expected; **do not** "fix" them in ways that
break 24.04.3. See [Nextflow version compatibility](../setup.md#nextflow-version-compatibility).

Avoid `nextflow lint -format` on Nextflow 25.10+, as it rewrites code to the
strict syntax; if you need to reformat, run it with `NXF_VER=25.04.7` so the
result stays compatible with parsers prior to 25.10.

`nextflow lint` always parses with the strict syntax; `NXF_SYNTAX_PARSER=v1` has
no effect on it (that variable only changes `nextflow run`).

## Releasing

The version number appears in three files, and there is no single source of truth --
bump all three together:

| File | Field |
| --- | --- |
| `nextflow.config` | `manifest.version` |
| `CITATION.cff` | `version`, `date-released`, and the `repository-code` tree URL |
| `CHANGELOG.md` | a new `## [X.Y.Z] - YYYY-MM-DD` heading below `## [Unreleased]` |

`docs/docs/changelog.md` is a symlink to the top-level `CHANGELOG.md`, so it updates
itself. Leave `doi:` in `CITATION.cff` alone -- it is the Zenodo *concept* DOI that
resolves to all versions; the per-version DOI is minted by Zenodo when the GitHub
release is published.

To cut a release:

```bash
# 1. Condense the [Unreleased] entries, then close the section as the new version
#    and bump nextflow.config + CITATION.cff.
$EDITOR CHANGELOG.md nextflow.config CITATION.cff

# 2. Commit the bump
git commit -am "Prepare 0.3.1 release: <summary>"

# 3. Tag with UNPREFIXED semver -- 0.3.1, not v0.3.1
git tag -a 0.3.1 -m "nf-binder-design 0.3.1"

# 4. Keep main and develop at the same commit, then push both plus the tag
git push origin develop
git push origin main
git push origin 0.3.1

# 5. Publish the release, using the new CHANGELOG section as the body
gh release create 0.3.1 --title "0.3.1" --notes-file <(...)
```

Tags must be unprefixed semver: `.github/workflows/docs.yml` builds the versioned
documentation from tags matching `v*` or `[0-9]+.[0-9]+.[0-9]+`, and every existing tag
in the repository is unprefixed (`0.3.0`, `0.2.0`, ...). Pushing the tag is what
publishes the immutable `/X.Y.Z/` docs; pushes to `main` and `develop` update the
floating docs for those branches.

Steps 4 and 5 above trigger the `docs` workflow up to three times at once, and each run
deploys to the same `gh-pages` branch. A `docs-deploy` concurrency group serialises them
and a rejected push is retried, so they queue rather than clobber each other -- but they
run one at a time, so give them a few minutes before checking that `/X.Y.Z/` is live.


## License

The `nf-binder-design` pipeline code is licensed under the MIT License.

> Note that some software dependencies of the pipeline are under less permissive licenses - in particular, RFdiffusion and BindCraft use [Rosetta/PyRosetta](https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md) which is **only free for Non-Commercial use**.

See also the [license section on the home page](../index.md#license) for commercial-use restrictions and citation guidance.
