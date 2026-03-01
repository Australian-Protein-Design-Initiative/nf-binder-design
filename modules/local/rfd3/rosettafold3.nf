process ROSETTAFOLD3 {
    container 'ghcr.io/australian-protein-design-initiative/containers/rc-foundry:0.1.11-weights'
    // container "rosettacommons/foundry:0.1.9-weights"

    publishDir path: "${params.outdir}/rfd3/rosettafold3/output", pattern: 'output/*', mode: 'copy', saveAs: { it.replaceFirst('output/', '') }

    input:
    tuple val(meta), path(structure_cif)
    val(uid)

    output:
    path 'output/*', emit: results
    path 'output/*/*_summary_confidences.json', emit: confidence_json
    tuple val(meta), path('output/*/*_model.cif'), emit: refolded_cif
    path 'rfd3_*_rf3_*.tsv', emit: scores

    // TODO: support MSA input for target chain(s) to improve prediction quality
    //       - this is important, since in my limited testing RosettaFold3 isn't performing very well
    //         on target protein prediction without an MSA (de novo binders seem to be predicted much better)
    //
    // TODO: support (and encourage) batching for lower process startup costs, inputs='folder/of/cifs'
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

    full_conf=\$(find output -maxdepth 2 -name '*_confidences.json' ! -name '*summary*' -print | head -1)
    summary_conf=\$(find output -maxdepth 2 -name '*_summary_confidences.json' -print | head -1)
    cif=\$(find output -maxdepth 2 -name '*_model.cif' -print | head -1)
    if [[ -n "\$full_conf" && -n "\$summary_conf" && -n "\$cif" ]]; then
        outdir=\$(dirname "\$cif")
        (cd "\$outdir" && python ${projectDir}/bin/ipsae.py \\
            --format rf3 \\
            --update-summary "\$(basename "\$summary_conf")" \\
            --binder-chain A --target-chain B \\
            "\$(basename "\$full_conf")" "\$(basename "\$cif")" 10 10)
    fi

    if [[ -n "\$summary_conf" ]]; then
        suffix=\$(basename "\$summary_conf" _summary_confidences.json | sed -n 's/.*\\(cif_b[0-9]*_d[0-9]*\\)/\\1/p')
        [[ -z "\$suffix" ]] && suffix=rf3
        python ${projectDir}/bin/rfd3/extract_rfd3_scores.py rosettafold3 "\$summary_conf" -o rfd3_${uid}_batch0_rf3_\${suffix}.tsv
    fi
    """
}
