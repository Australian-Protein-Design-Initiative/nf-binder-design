process AF2_INITIAL_GUESS {
    container 'ghcr.io/australian-protein-design-initiative/containers/af2_initial_guess:nv-cuda12'

    publishDir "${params.outdir}/af2_initial_guess", pattern: 'pdbs/*.pdb', mode: 'copy'
    publishDir "${params.outdir}/af2_initial_guess", pattern: 'scores/*.cs', mode: 'copy'

    input:
    path 'input/*'

    output:
    tuple path('pdbs/*.pdb'), path('af2ig_scores.tsv'), emit: pdbs_with_scores
    path 'pdbs/*.pdb', emit: pdbs
    path 'scores/*.cs', emit: scores
    path 'gpu_stats.csv', optional: true, emit: gpu_stats

    script:
    """
    mkdir -p scores/

    # Find least-used GPU (by active processes and VRAM) and set CUDA_VISIBLE_DEVICES
    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    # Start GPU monitoring in background if enabled
    if [[ "${params.enable_gpu_stats}" == "true" ]]; then
        PARENT_DIR=\$(basename \$(dirname \$(pwd)))
        CURRENT_DIR=\$(basename \$(pwd))
        TASK_HASH="\${PARENT_DIR}/\${CURRENT_DIR}"
        TASK_HASH="\${TASK_HASH:0:9}"
        ${baseDir}/bin/monitor_gpu.py \
            --process-name "AF2_INITIAL_GUESS" \
            --task-hash "\${TASK_HASH}" \
            --task-index "${task.index}" \
            --interval ${params.gpu_stats_interval} \
            --output gpu_stats.csv &
        GPU_MONITOR_PID=\${!:-}
        if [[ -n "\${GPU_MONITOR_PID}" ]]; then
            trap "kill \${GPU_MONITOR_PID} 2>/dev/null || true" EXIT
        fi
    fi

    # Get first input PDB filename without extension
    PREFIX=\$(ls input/*.pdb | head -n1 | xargs basename | sed 's/\\.pdb\$//')

    python /app/dl_binder_design/af2_initial_guess/predict.py \
        -pdbdir input/ \
        -outpdbdir pdbs/ \
        -recycle ${params.af2ig_recycle} \
        -scorefilename scores/\${PREFIX}.scores.cs

    # Combine scores into a single TSV file
    python ${projectDir}/bin/af2_combine_scores.py scores/ -o af2ig_scores.tsv
    """
}
