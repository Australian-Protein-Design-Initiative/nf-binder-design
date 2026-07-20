// NEW, decoupled RF3 fold process for fold.nf - independent of the
// binder-design-coupled ROSETTAFOLD3 (modules/local/rfd3/rosettafold3.nf),
// which is hard-wired to a 2-body target+binder complex plus a binder
// postprocess step. This process folds whatever generic component list
// GENERATE_RF3_FOLD_INPUT built (monomer in Phase 1; multimer is Phase 2),
// with no postprocessing beyond what `rf3 fold` itself emits.
process RF3_FOLD {
    tag "${meta.id}"

    container 'oras://ghcr.io/australian-protein-design-initiative/containers/rc-foundry:0.2.0-weights'

    // Recursive glob publish - see modules/fold/af2/alphafold2.nf's comment
    // for why: publishDir pattern/saveAs only see the top-level output
    // *item* (a directory), never files nested inside it, so declaring the
    // output as "output/**" makes each nested file its own publish item.
    // rf3 fold writes to output/<name>/ where name == meta.id (see
    // generate_rf3_fold_input.nf --name '${meta.id}'), so just strip the
    // task-local "output/" prefix to land at <outdir>/rf3/<meta.id>/... .
    // (Replacing it with "${meta.id}/" instead would double-nest to
    // <outdir>/rf3/<meta.id>/<meta.id>/... .)
    publishDir(
        path: "${params.outdir}/rf3",
        mode: 'copy',
        saveAs: { filename -> filename.toString().replaceFirst(/^output\//, '') }
    )

    input:
    tuple val(meta), path(fasta), path(a3m), path(rf3_input_json)

    output:
    tuple val(meta), path('output/**'), emit: predictions
    tuple val(meta), path('output/**/*_summary_confidences.json'), emit: confidence_json

    script:
    """
    set -euo pipefail

    if [[ ${params.require_gpu} == "true" ]]; then
        if ! command -v nvidia-smi >/dev/null 2>&1; then
            echo "nvidia-smi not found / no NVIDIA driver detected! Failing fast rather than going slow (since --require_gpu=true; set --require_gpu false to bypass)"
            exit 1
        fi

        if [[ \$(nvidia-smi -L) =~ "No devices found" ]]; then
            echo "No GPU detected! Failing fast rather than going slow (since --require_gpu=true)"
            exit 1
        fi

        nvidia-smi
    fi

    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${projectDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    mkdir -p output

    rf3 fold \\
        inputs=${rf3_input_json} \\
        out_dir=output \\
        ckpt_path=${params.rf3_ckpt_path} \\
        num_steps=${params.rf3_num_steps} \\
        n_recycles=${params.rf3_n_recycles} \\
        diffusion_batch_size=${params.rf3_diffusion_batch_size} \\
        early_stopping_plddt_threshold=${params.rf3_early_stopping_plddt_threshold} \\
        annotate_b_factor_with_plddt=true \\
        one_model_per_file=true \\
        ${task.ext.args ?: ''}
    """
}
