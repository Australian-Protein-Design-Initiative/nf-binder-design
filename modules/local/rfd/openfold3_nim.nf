// Uses openfold3-nim (ENTRYPOINT cleared, see aws_batch_nims.config).
// Replaces the baseline AF2_INITIAL_GUESS scoring step. Note this is not an
// "initial guess" fold: OpenFold3 predicts from sequence alone and cannot be
// seeded with the design's coordinates, so the input structure is only used to
// read the chain sequences off.
process OPENFOLD3_NIM {
    publishDir "${params.outdir}/rfd/openfold3_nim", pattern: 'pdbs/*.pdb', mode: 'copy'
    publishDir "${params.outdir}/rfd/openfold3_nim", pattern: 'scores/*.tsv', mode: 'copy'

    input:
    path design_pdb
    val diffusion_samples

    output:
    path 'pdbs/*.pdb', emit: pdbs
    path 'scores/*.tsv', emit: scores
    tuple path('pdbs/*.pdb'), path('scores/*.tsv'), emit: pdbs_with_scores

    script:
    """
    set -euo pipefail
    mkdir -p pdbs scores

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

    export OF3_INPUT_PDB="${design_pdb}"
    export OF3_OUTPUT_DIR="pdbs"
    export OF3_OUTPUT_PREFIX="${design_pdb.baseName}"
    export OF3_OUTPUT_SCORES="scores/${design_pdb.baseName}.of3_scores.tsv"
    export OF3_DIFFUSION_SAMPLES="${diffusion_samples}"
    export NGC_API_KEY="${System.getenv('NGC_API_KEY') ?: ''}"

    /opt/nim/start_server.sh > nim_server.log 2>&1 &

    # Longer budget than the RFdiffusion/ProteinMPNN NIMs - OpenFold3 pulls
    # ~15GB of model parameters into its cache on a cold start.
    for i in \$(seq 1 120); do
        if curl -sf http://localhost:8000/v1/health/ready > /dev/null 2>&1; then
            echo "NIM server ready after \${i} check(s)"
            break
        fi
        echo "Waiting for NIM server to be ready (\${i}/120)..."
        sleep 10
    done

    if ! curl -sf http://localhost:8000/v1/health/ready > /dev/null 2>&1; then
        echo "NIM server never became ready. Last 100 lines of nim_server.log:" >&2
        tail -n 100 nim_server.log >&2
        exit 1
    fi

    openfold3_nim_call.py
    """
}
