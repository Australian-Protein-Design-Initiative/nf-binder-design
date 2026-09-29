// Uses rfdiffusion-nim-updated (ENTRYPOINT cleared, see aws_batch_nims.config).
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
    """
    set -euo pipefail
    mkdir -p pdbs

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

    export RFD_INPUT_PDB="${input_pdb}"
    export RFD_CONTIGS="${contigs}"
    export RFD_HOTSPOT_RES="${hotspot_res}"
    export RFD_OUTPUT_PDB="pdbs/design_ppi_${unique_id}_${design_index}.pdb"
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

    rfdiffusion_nim_call.py
    """
}
