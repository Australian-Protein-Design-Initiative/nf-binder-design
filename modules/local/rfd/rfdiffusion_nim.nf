// Container set via conf/platforms/aws_batch_nims.config (rfdiffusion-nim-updated -
// same image NVIDIA publishes, rebuilt only to clear its baked-in ENTRYPOINT).
//
// The NIM image's default behavior is to start its own HTTP server and sit there
// forever - that's what ENTRYPOINT normally does here, and it would otherwise
// swallow whatever command Nextflow tries to run instead of letting it execute.
// With the entrypoint cleared, this script has to do that startup itself: launch
// the NIM's own server in the background, wait for it to report ready, call its
// API once, save the result, then exit - so the task still behaves like every
// other step in this pipeline (one container, runs, produces output, exits).
process RFDIFFUSION_NIM {
    publishDir "${params.outdir}/rfd/rfdiffusion", pattern: 'pdbs/*.pdb', mode: 'copy'

    input:
    path input_pdb
    val contigs
    val hotspot_res
    val design_index
    val unique_id

    output:
    path 'pdbs/*.pdb', emit: pdbs

    script:
    // Passed in via environment (see below) rather than interpolated directly
    // into the shell script, to avoid three-way quoting problems between
    // Nextflow/Groovy, bash, and the Python helper this script shells out to.
    """
    set -euo pipefail
    mkdir -p pdbs

    export RFD_INPUT_PDB="${input_pdb}"
    export RFD_CONTIGS="${contigs}"
    export RFD_HOTSPOT_RES="${hotspot_res}"
    export RFD_OUTPUT_PDB="pdbs/design_ppi_${unique_id}_${design_index}.pdb"

    # NIM downloads model weights from NVIDIA's NGC API on first startup and
    # needs this to authenticate. Read from the local shell environment (not
    # a --param) so it never gets written to params.json in the output bucket.
    export NGC_API_KEY="${System.getenv('NGC_API_KEY') ?: ''}"

    # Start the NIM's own server (normally done automatically by its entrypoint,
    # which this image build has cleared).
    /opt/nim/start_server.sh &

    for i in \$(seq 1 60); do
        if curl -sf http://localhost:8000/v1/health/ready > /dev/null 2>&1; then
            echo "NIM server ready after \${i} check(s)"
            break
        fi
        echo "Waiting for NIM server to be ready (\${i}/60)..."
        sleep 5
    done
    curl -sf http://localhost:8000/v1/health/ready

    rfdiffusion_nim_call.py
    """
}
