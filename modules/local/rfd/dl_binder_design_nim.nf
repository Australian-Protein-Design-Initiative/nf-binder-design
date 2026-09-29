// Uses proteinmpnn-nim-updated (ENTRYPOINT cleared, see aws_batch_nims.config).
// Emits sequences only (FASTA), not a threaded PDB - see THREAD_AND_RELAX for
// the step that builds a real structure from this output.
process DL_BINDER_DESIGN_PROTEINMPNN_NIM {
    publishDir "${params.outdir}/rfd/proteinmpnn_nim", pattern: 'fasta/*.fasta', mode: 'copy'

    input:
    path backbone_pdb
    val design_chain
    val sampling_temp
    val design_index

    output:
    path 'fasta/*.fasta', emit: fasta
    path backbone_pdb, emit: backbone

    script:
    """
    set -euo pipefail
    mkdir -p fasta

    # Kill the NIM server on exit so the container doesn't hang.
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
    export NGC_API_KEY="${System.getenv('NGC_API_KEY') ?: ''}"

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
