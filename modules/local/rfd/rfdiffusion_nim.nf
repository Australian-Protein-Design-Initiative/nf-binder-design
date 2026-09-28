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

    # The NIM server we start below is a real standing service - it never exits
    # on its own. Without cleanup, the container hangs forever after our own
    # work is done (confirmed: job completes the API call, then sits alive
    # until AWS Batch's job timeout kills it hours later). Two things are
    # needed, not just one: (1) the server's own worker processes fork into a
    # separate process group, so a plain process-group kill on exit doesn't
    # reach them - kill by matching the command line instead; (2) more
    # importantly, the server inherits our script's stdout/stderr, so even
    # after our script exits, Nextflow's own log-capturing wrapper can't see
    # end-of-stream until every process holding that pipe closes it -
    # including orphaned NIM workers regardless of process group. Redirecting
    # the server's output to its own file avoids that entirely.
    SELF_PID=\$\$
    trap '
        code=\$?
        for pid in \$(pgrep -g "\$SELF_PID" 2>/dev/null); do
            if [ "\$pid" != "\$SELF_PID" ]; then
                kill -9 "\$pid" 2>/dev/null || true
            fi
        done
        pkill -9 -f start_server 2>/dev/null || true
        exit "\$code"
    ' EXIT

    export RFD_INPUT_PDB="${input_pdb}"
    export RFD_CONTIGS="${contigs}"
    export RFD_HOTSPOT_RES="${hotspot_res}"
    export RFD_OUTPUT_PDB="pdbs/design_ppi_${unique_id}_${design_index}.pdb"

    # NIM downloads model weights from NVIDIA's NGC API on first startup and
    # needs this to authenticate. Read from the local shell environment (not
    # a --param) so it never gets written to params.json in the output bucket.
    export NGC_API_KEY="${System.getenv('NGC_API_KEY') ?: ''}"

    # Start the NIM's own server (normally done automatically by its entrypoint,
    # which this image build has cleared). Output redirected to its own file,
    # not inherited - see comment above for why.
    /opt/nim/start_server.sh > nim_server.log 2>&1 &

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
