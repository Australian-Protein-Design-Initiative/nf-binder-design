# Setup

## Prerequisites

### Install Nextflow

Follow the official [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html).

Nextflow requires:

- Java 11 or later
- Linux or macOS (Windows via WSL2)

The quick version is:

```bash
# Optional, if you don't have Java 17+
# Check your Java version
# java -version

# Install Java 17+, if you don't have it
# curl -s https://get.sdkman.io | bash
# sdk install java 17.0.10-tem

# Install Nextflow
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow $HOME/.local/bin/

# Ensure $HOME/.local/bin/ is on your PATH
# you should also add this line to your ~/.bashrc or ~/.zshrc if required
export PATH="$PATH:$HOME/.local/bin"

# Test it works
nextflow info
```

## Clone the repository using git

For development, it's often more convenient to clone the git repository directly:

```bash
git clone https://github.com/Australian-Protein-Design-Initiative/nf-binder-design
cd nf-binder-design

# See the help as a first test
nextflow run main.nf --help
```

## Alternative: pull using Nextflow

To pre-cache the `nf-binder-design` repository in `~/.nextflow/assets/`, run:

```bash
nextflow pull Australian-Protein-Design-Initiative/nf-binder-design
```

you can pull a specific version like:

```bash
nextflow pull -r 0.1.4 Australian-Protein-Design-Initiative/nf-binder-design
```

Test that you have successfully cached the repository by running:
```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design/main.nf --help
```

**Note:** The documentation examples currently assume you have `git clone`'d the repository, but can be adapted to use this alternative method. In this case, platform-specific config files will be in `~/.nextflow/assets/Australian-Protein-Design-Initiative/nf-binder-design/conf/platforms/` - in the future these will be turned into config profiles to simplify usage.

## Container Setup

These workflows use Apptainer containers for reproducibility. 
You'll typically need [Apptainer](https://apptainer.org/) installed - this is usually already available on modern HPC clusters.

Containers are automatically downloaded when running the pipeline - you will need ~60Gb+ of storage space for the containers.

Model weights / parameters are currently packaged inside the containers - this makes the containers large, but simplifies the workflow setup. In the future we will very likely provide more compact alternative containers without weights.

## Platform-Specific Configuration

### Local Workstation / single compute node

For local execution of a single compute node, use the `-profile local` flag when running workflows. See [`examples/pdl1-rfd/nextflow.dual-gpu.config`](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/blob/main/examples/pdl1-rfd/nextflow.dual-gpu.config) for a config example using dual-GPUs.

### HPC Clusters with SLURM

For running using the SLURM executor, use the `-profile slurm` flag (this is the default profile if not specified).

Site-specific configuration files are provided in `conf/platforms/`:

- `m3.config` - Monash M3 cluster
- `m3-bdi.config` - Monash M3 cluster with access to the `bdi` and `m3h` partitions
- `mlerp.config` - MLeRP cluster

These can be adapted to other HPC clusters - pull requests are welcome!

Use the `-c` flag to specify the path to a configuration file:

```bash
nextflow run main.nf -c conf/platforms/m3.config ...
```

#### HPC Cluster environment variables

Set up your Apptainer cache directory to store downloaded containers:

```bash
# Add to your ~/.bashrc or set before running
export APPTAINER_CACHEDIR=${HOME}/.apptainer/cache
export NXF_APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR}
```

Changing this location is particularly important on HPC clusters where you may have limited space in your home directory.

Ensure your temporary directory is set to a location with sufficient space, eg:

```bash
export TMPDIR=/scratch2/myproject/tmp
export NXF_TEMP=$TMPDIR
mkdir -p $TMPDIR
```

## Next Steps

See the [Usage](usage.md) page for examples and detailed workflow documentation.
