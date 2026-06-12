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

## License

The `nf-binder-design` pipeline code is licensed under the MIT License.

> Note that some software dependencies of the pipeline are under less permissive licenses - in particular, RFdiffusion and BindCraft use [Rosetta/PyRosetta](https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md) which is **only free for Non-Commercial use**.

See also the [license section on the home page](../index.md#license) for commercial-use restrictions and citation guidance.
