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

    # Claim a GPU for this task's lifetime, then record which card we got
    # (bin/gpu_lock.sh). The claim is required, and fails the task if it cannot
    # be made. The recording is diagnostic, and must never fail the task -- the
    # `|| true` also suspends `set -e` for the whole function body, so nothing
    # inside it can abort the script either.
    if [[ -n "${params.gpu_devices}" ]]; then
        source ${projectDir}/bin/gpu_lock.sh
        nfbd_acquire_gpu "${params.gpu_devices}" "${params.gpu_lock_dir ?: workDir.toString() + '/.gpu_locks'}" ${task.ext.gpu_slots ?: params.gpu_slots_per_device} ${params.gpu_lock_timeout} || exit 1
    else
        source ${projectDir}/bin/gpu_lock.sh || true
    fi
    nfbd_record_gpu_trace "${params.gpu_trace_dir ?: workDir.toString() + '/.gpu_trace'}" "${task.process}" || true

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
