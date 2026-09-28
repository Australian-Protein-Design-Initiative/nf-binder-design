// Container set via conf/platforms/aws_batch_nims.config (proteinmpnn-nim-updated).
// Same server-start/wait/call/exit pattern as rfdiffusion_nim.nf - see that
// file's header comment for why.
//
// IMPORTANT - output shape differs from the baseline DL_BINDER_DESIGN_PROTEINMPNN:
// the NIM's /predict endpoint returns designed sequences only (a multi-FASTA),
// not a fully-built PDB with the new sequence's side chains threaded onto the
// backbone the way the baseline module's output is. Emits `fasta`, not `pdbs`.
// AF2_INITIAL_GUESS as it exists today expects the latter - wiring this into the
// rest of the smoke test needs that gap resolved first (a separate threading
// step, or reworking what AF2_INITIAL_GUESS consumes). Not solved here.
process DL_BINDER_DESIGN_PROTEINMPNN_NIM {
    publishDir "${params.outdir}/rfd/proteinmpnn_nim", pattern: 'fasta/*.fasta', mode: 'copy'

    input:
    path backbone_pdb
    val design_chain
    val sampling_temp
    val design_index

    output:
    path 'fasta/*.fasta', emit: fasta

    script:
    """
    set -euo pipefail
    mkdir -p fasta

    # The NIM server we start below is a real standing service - it never exits
    # on its own. Without cleanup, the container hangs forever after our own
    # work is done until AWS Batch's job timeout kills it. Two things needed,
    # not just one (see rfdiffusion_nim.nf for the full story of how this was
    # confirmed): (1) the server's own worker processes fork into a separate
    # process group, so a plain process-group kill doesn't reach them - match
    # by command line too; (2) the server inherits our script's stdout/stderr,
    # so even after our script exits, Nextflow's log-capturing wrapper can't
    # see end-of-stream until every process holding that pipe closes it -
    # redirecting the server's output to its own file avoids that entirely.
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

    export PMPNN_BACKBONE_PDB="${backbone_pdb}"
    export PMPNN_DESIGN_CHAIN="${design_chain}"
    export PMPNN_SAMPLING_TEMP="${sampling_temp}"
    export PMPNN_OUTPUT_FASTA="fasta/${backbone_pdb.baseName}_${design_index}.fasta"

    # NIM downloads model weights from NVIDIA's NGC API on first startup and
    # needs this to authenticate. Read from the local shell environment (not
    # a --param) so it never gets written to params.json in the output bucket.
    export NGC_API_KEY="${System.getenv('NGC_API_KEY') ?: ''}"

    # Output redirected to its own file, not inherited - see comment above.
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

    proteinmpnn_nim_call.py
    """
}
