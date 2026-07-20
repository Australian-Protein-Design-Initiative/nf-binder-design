// Protenix (AF3-style) folding process for fold.nf's 4th --methods engine
// (see plans/fold-nf-multi-method-folding.md Phase 3). Weights are baked into
// the container under /models/protenix/checkpoint (confirmed present -
// protenix_base_default_v1.0.0.pt et al. - and /models/protenix/common, so
// `protenix pred` needs no download at predict time, same as rf3's
// rc-foundry image).
process PROTENIX_FOLD {
    tag "${meta.id}"

    container 'oras://ghcr.io/australian-protein-design-initiative/containers/protenix:v2.0.0-weights'

    // Recursive glob publish - see modules/fold/rf3/rf3_fold.nf's comment for
    // why (publishDir pattern/saveAs only see the top-level output *item*,
    // never files nested inside it). protenix pred writes to
    // out_dir/<name>/seed_<seed>/predictions/... where name == meta.id (see
    // generate_protenix_input.nf's --name '${meta.id}'; confirmed against
    // runner/inference.py's infer_predict(), which uses each JSON job's
    // top-level "name" field as both dump dataset dir and sample_name), so
    // strip only the task-local "output/" prefix to land at
    // <outdir>/protenix/<meta.id>/... . (Replacing it with "${meta.id}/"
    // instead would double-nest to <outdir>/protenix/<meta.id>/<meta.id>/...,
    // the same trap rf3_fold.nf's comment warns about.)
    publishDir(
        path: "${params.outdir}/protenix",
        mode: 'copy',
        saveAs: { filename -> filename.toString().replaceFirst(/^output\//, '') }
    )

    input:
    tuple val(meta), path(fasta), path(a3m), path(protenix_input_json)

    output:
    tuple val(meta), path('output/**'), emit: predictions
    tuple val(meta), path('output/**/*_summary_confidence_sample_*.json'), emit: confidence_json

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

    protenix pred \\
        --input ${protenix_input_json} \\
        --out_dir output \\
        --seeds ${params.protenix_seeds} \\
        --cycle ${params.protenix_cycle} \\
        --step ${params.protenix_step} \\
        --sample ${params.protenix_sample} \\
        --model_name ${params.protenix_model_name} \\
        --use_msa ${params.protenix_use_msa} \\
        ${task.ext.args ?: ''}
    """
}
