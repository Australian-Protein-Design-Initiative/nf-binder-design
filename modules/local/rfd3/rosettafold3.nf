process ROSETTAFOLD3 {
    container 'ghcr.io/australian-protein-design-initiative/containers/rc-foundry:0.1.11-weights'
    // container "rosettacommons/foundry:0.1.9-weights"

    publishDir path: "${params.outdir}/rfd3/rosettafold3/output", pattern: 'output/*', mode: 'copy', saveAs: { it.replaceFirst('output/', '') }

    input:
    tuple val(meta), path(structure_cif)

    output:
    path 'output/*', emit: results
    path 'output/*/*_summary_confidences.json', emit: confidence_json
    tuple val(meta), path('output/*/*_model.cif'), emit: refolded_cif

    // TODO: support (and encourage) batching for lower process startup costs, inputs='folder/of/cifs'
    // TODO: support MSA input for target chain(s) to improve prediction quality
    // TODO: add params.rf3_early_stopping_plddt_threshold=0.5 as default for filtering low-confidence predictions
    // TODO: support JSON config mode with template_selection to template target chain(s)
    // TODO: support n_recycles, diffusion_batch_size, num_steps, seed, and consider if we can use a faster
    //       default num_steps=50 as a good compromise between speed and quality

    script:
    """
    set -euo pipefail

    if [[ ${params.require_gpu} == "true" ]]; then
       if [[ \$(nvidia-smi -L) =~ "No devices found" ]]; then
           echo "No GPU detected! Failing fast rather than going slow (since --require_gpu=true)"
            exit 1
        fi

        nvidia-smi
    fi

    # Find least-used GPU and set CUDA_VISIBLE_DEVICES
    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    mkdir -p output

    rf3 fold \
        inputs=${structure_cif} \
        out_dir=output \
        ckpt_path=${params.rf3_ckpt_path} \
        ${task.ext.args ?: ''}
    """
}
