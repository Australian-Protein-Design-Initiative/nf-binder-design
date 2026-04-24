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
                            # (includes ../nextflow.config, forces local executor)
  modules/local/...*.nf.test
  workflows/...*.nf.test
```

The project-level `nf-test.config` sets `testsDir = "tests"` so only files
under `tests/` are picked up by `nf-test test`.

The test-time `tests/nextflow.config` includes the repository root
`nextflow.config` via `includeConfig '../nextflow.config'`, so Apptainer stays
enabled and each process uses the `container` image declared on that process in
its module file (there is no single global test image). A follow-up
`process { executor = 'local' }` block overrides the default SLURM executor so
tests run on the current machine. Processes that do not declare a `container`
still execute on the host interpreter when Apptainer is enabled.

## Apptainer cache and image filenames

When Nextflow (or `apptainer pull docker://…`) pulls an image, the file on disk
is **not** named like the `docker://` URI. It is stored under
`NXF_APPTAINER_CACHEDIR` / `APPTAINER_CACHEDIR` (often `~/.apptainer/cache`) with
slashes turned into hyphens and the tag appended, for example:

| Module `container` directive | Typical cached filename |
|------------------------------|-------------------------|
| `ghcr.io/australian-protein-design-initiative/containers/rfdiffusion:pytorch2407` | `ghcr.io-australian-protein-design-initiative-containers-rfdiffusion-pytorch2407.img` |

Use that path for direct `apptainer exec` or debugging, e.g.:

```bash
APPTAINER_CACHE="${NXF_APPTAINER_CACHEDIR:-$HOME/.apptainer/cache}"
RFD="${APPTAINER_CACHE}/ghcr.io-australian-protein-design-initiative-containers-rfdiffusion-pytorch2407.img"
apptainer exec --nv "$RFD" python /app/RFdiffusion/scripts/run_inference.py --help
```

Nextflow still resolves the module `container` string to this file when the cache
is set; the table above is for **manual** `apptainer` use without `docker://`.

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
for example `--profile local` to merge the local executor block from the main
config (the test config already sets `executor = 'local'`; profiles are for
extra options such as queue tuning).

Containerised tests need Apptainer available and enough disk for image pulls;
see [Setup](setup.md) for setting `NXF_APPTAINER_CACHEDIR` and related variables.
The inherited `apptainer.runOptions` include `--nv`; GPU images need a host with
NVIDIA drivers unless you add a dedicated test profile that relaxes those
options for CPU-only CI. For on-disk image names and `apptainer exec` without
`docker://`, see [Apptainer cache and image filenames](#apptainer-cache-and-image-filenames) above.

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
- The `RFDIFFUSION` process tests (`tests/modules/local/rfd/rfdiffusion.nf.test`) need a
  visible GPU (`require_gpu = true` in the test params, matching the pipeline default),
  the rfdiffusion Apptainer image, and an extra minute or so per test. They are
  tagged `rfd` and `gpu`. Run only those with
  `nf-test test --tag rfd` (or `nf-test test tests/modules/local/rfd/rfdiffusion.nf.test`);
  to skip them, run
  `nf-test test --tag unique_id` (or any tag that excludes the RFD file). `RFDIFFUSION`
  picks up `inference.deterministic=true` from `tests/nextflow.config` via
  `process.withName: RFDIFFUSION { ext.args = '...' }` and a trailing append in
  `modules/local/rfd/rfdiffusion.nf`. A standalone reproducer lives under
  a full-`rfd` example: `examples/glp1-rfd-mini/` (runs `main.nf --method rfd` on the GLP-1 test structure).
- Keep other heavy GPU-bound processes (BindCraft, Boltz, BoltzGen, AF2 initial guess) out
  of the default test suite, or guard them behind tags, until you have a similar setup.
- Tag every test (`tag "<feature>"`) so suites can be selected with
  `nf-test test --tag <feature>`.
