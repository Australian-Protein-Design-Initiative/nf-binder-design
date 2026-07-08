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

## License

The `nf-binder-design` pipeline code is licensed under the MIT License.

> Note that some software dependencies of the pipeline are under less permissive licenses - in particular, RFdiffusion and BindCraft use [Rosetta/PyRosetta](https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md) which is **only free for Non-Commercial use**.

See also the [license section on the home page](../index.md#license) for commercial-use restrictions and citation guidance.
