process BOLTZGEN_MERGE {
    tag "merge"

    container 'ghcr.io/australian-protein-design-initiative/containers/boltzgen:0.2.0'

    publishDir path: "${params.outdir}/boltzgen", pattern: '**', mode: 'copy'

    input:
    path batch_dirs

    output:
    path 'merged', type: 'dir', emit: merged_dir
    path 'gpu_stats.csv', optional: true, emit: gpu_stats

    script:
    """
    set -euo pipefail

    # Start GPU monitoring in background if enabled
    if [[ "${params.enable_gpu_stats}" == "true" ]]; then
        PARENT_DIR=\$(basename \$(dirname \$(pwd)))
        CURRENT_DIR=\$(basename \$(pwd))
        TASK_HASH="\${PARENT_DIR}/\${CURRENT_DIR}"
        TASK_HASH="\${TASK_HASH:0:9}"
        ${baseDir}/bin/monitor_gpu.py \
            --process-name "BOLTZGEN_MERGE" \
            --task-hash "\${TASK_HASH}" \
            --task-index "${task.index}" \
            --interval ${params.gpu_stats_interval} \
            --output gpu_stats.csv &
        GPU_MONITOR_PID=\${!:-}
        if [[ -n "\${GPU_MONITOR_PID}" ]]; then
            trap "kill \${GPU_MONITOR_PID} 2>/dev/null || true" EXIT
        fi
    fi

    # With _many_ batch_ directories, this will fail due to ARG_MAX limits with shell glob expansion
    #boltzgen merge batch_* \
    #    --output merged/ \
    #    ${task.ext.args ?: ''}

    # We need to do some shell-gymnastics to avoid ARG_MAX limits with shell glob expansion
    # Find batch directories, including symlinks to directories (Nextflow creates symlinks for inputs)
    # Use -L to follow symlinks, or check -d which works for both directories and symlinks to directories
    batch_dirs=()
    for item in \$(find . -maxdepth 1 -name 'batch_*' | sort); do
        if [ -d "\$item" ]; then
            batch_dirs+=("\$item")
        fi
    done
    
    if [ \${#batch_dirs[@]} -eq 0 ]; then
        echo "ERROR: No batch directories found."
        exit 1
    fi
    
    boltzgen merge "\${batch_dirs[@]}" \
        --output merged/ \
        ${task.ext.args ?: ''}

    if grep -q "No designs found to merge." .command.log; then
        echo "ERROR: No designs found to merge."
        exit 1
    fi
    """
}
