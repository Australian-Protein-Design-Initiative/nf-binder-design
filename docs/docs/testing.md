# Testing

The pipeline uses [nf-test](https://www.nf-test.com/) for unit and integration
tests. Tests live in the `tests/` directory and mirror the layout of the
modules and workflows they cover.

## Install nf-test

`nf-test` is a single self-contained binary. Install it once and place it on
your `PATH`:

```bash
curl -fsSL https://code.askimed.com/install/nf-test | bash
mv nf-test "$HOME/.local/bin/"

nf-test version
```

It requires a working Nextflow installation (see [Setup](setup.md)).

## Layout

```
nf-test.config              # project-wide nf-test configuration
tests/
  nextflow.config           # Nextflow config used when running tests
                            # (forces local executor, disables containers)
  modules/local/...*.nf.test
  workflows/...*.nf.test
```

The project-level `nf-test.config` sets `testsDir = "tests"` so only files
under `tests/` are picked up by `nf-test test`.

The test-time `tests/nextflow.config` overrides the SLURM executor and the
Apptainer/Docker containers from the main `nextflow.config` so simple tests
can run on a developer workstation with no GPU and no scheduler.

## Run the tests

From the repository root:

```bash
# Run every test under tests/
nf-test test

# Run a single test file
nf-test test tests/modules/local/common/unique_id.nf.test

# Run by tag (tags are declared inside the .nf.test files)
nf-test test --tag unique_id

# Generate a TAP report (useful in CI)
nf-test test --tap report.tap
```

Use `--profile <name>` to apply a Nextflow profile from `nextflow.config`,
for example `--profile local` to use the local executor block defined in the
main config.

See the [nf-test running tests docs](https://www.nf-test.com/docs/running-tests/)
for the full list of options.

## Writing a new test

Generate a skeleton for a process or workflow:

```bash
# For a process (modules/)
nf-test generate process modules/local/common/unique_id.nf

# For a workflow
nf-test generate workflow workflows/rfd.nf
```

This creates a `*.nf.test` file under `tests/` mirroring the source path.
Edit the generated file to define the `when` block (inputs/params) and the
`then` block (assertions on `process.out` / `workflow.out`).

A minimal example, lifted from `tests/modules/local/common/unique_id.nf.test`:

```groovy
nextflow_process {

    name "Test Process UNIQUE_ID"
    script "modules/local/common/unique_id.nf"
    process "UNIQUE_ID"

    tag "unique_id"

    test("Should generate a 7-character base58 unique_id.txt") {

        when {
            params {
                outdir = "${outputDir}"
            }
            process {
                """
                """
            }
        }

        then {
            assert process.success
            with(process.out.id_file) {
                assert size() == 1
                def id = path(get(0)).text.trim()
                assert id.length() == 7
                assert id ==~ /[1-9A-HJ-NP-Za-km-z]{7}/
            }
        }
    }
}
```

Refer to the upstream documentation for details on writing
[process](https://www.nf-test.com/docs/testcases/nextflow_process/),
[workflow](https://www.nf-test.com/docs/testcases/nextflow_workflow/) and
[pipeline](https://www.nf-test.com/docs/testcases/nextflow_pipeline/) tests.

## Tips

- Prefer lightweight, non-GPU processes for fast unit tests (eg `UNIQUE_ID`,
  config-generation processes, simple Python helpers).
- Keep heavy GPU-bound processes (RFdiffusion, BindCraft, Boltz, BoltzGen,
  AF2 initial guess) out of the default test suite, or guard them behind tags
  so they are only run on machines with the required hardware and
  containers.
- Tag every test (`tag "<feature>"`) so suites can be selected with
  `nf-test test --tag <feature>`.
