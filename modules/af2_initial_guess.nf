process AF2_INITIAL_GUESS {
    container 'ghcr.io/australian-protein-design-initiative/containers/af2_initial_guess:nv-cuda12'

    publishDir "${params.outdir}/af2_initial_guess", pattern: 'pdbs/*.pdb', mode: 'copy'
    publishDir "${params.outdir}/af2_initial_guess", pattern: 'scores/*.cs', mode: 'copy'

    input:
    path 'input/*'

    output:
    path 'pdbs/*.pdb', emit: pdbs
    path 'scores/*.cs', emit: scores
    path 'af2ig_scores.tsv', emit: scores_tsv

    script:
    """
    mkdir -p scores/

    # Find least-used GPU (by active processes and VRAM) and set CUDA_VISIBLE_DEVICES
    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
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
