process AF2_INITIAL_GUESS {
    container 'ghcr.io/australian-protein-design-initiative/containers/af2_initial_guess:nv-cuda12'

    publishDir "${params.outdir}/rfd/af2_initial_guess", pattern: 'pdbs/*.pdb', mode: 'copy'
    publishDir "${params.outdir}/rfd/af2_initial_guess", pattern: 'scores/*.cs', mode: 'copy'

    input:
    path 'input/*'

    output:
    tuple path('pdbs/*.pdb'), path('af2ig_scores.tsv'), emit: pdbs_with_scores
    path 'pdbs/*.pdb', emit: pdbs
    path 'scores/*.cs', emit: scores

    script:
    """
    mkdir -p scores/

    # Claim a GPU for this task's lifetime (bin/gpu_lock.sh).
    if [[ -n "${params.gpu_devices}" ]]; then
        source ${projectDir}/bin/gpu_lock.sh
        nfbd_acquire_gpu "${params.gpu_devices}" "${params.gpu_lock_dir ?: workDir.toString() + '/.gpu_locks'}" ${task.ext.gpu_slots ?: params.gpu_slots_per_device} ${params.gpu_lock_timeout} || exit 1
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
