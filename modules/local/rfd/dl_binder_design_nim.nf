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

    export PMPNN_BACKBONE_PDB="${backbone_pdb}"
    export PMPNN_DESIGN_CHAIN="${design_chain}"
    export PMPNN_SAMPLING_TEMP="${sampling_temp}"
    export PMPNN_OUTPUT_FASTA="fasta/${backbone_pdb.baseName}_${design_index}.fasta"

    # NIM downloads model weights from NVIDIA's NGC API on first startup and
    # needs this to authenticate. Read from the local shell environment (not
    # a --param) so it never gets written to params.json in the output bucket.
    export NGC_API_KEY="${System.getenv('NGC_API_KEY') ?: ''}"

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

    proteinmpnn_nim_call.py
    """
}
