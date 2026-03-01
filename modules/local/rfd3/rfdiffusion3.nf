process RFDIFFUSION3 {
    container 'ghcr.io/australian-protein-design-initiative/containers/rc-foundry:0.1.11-weights'
    // container "rosettacommons/foundry:0.1.9-weights"

    publishDir path: "${params.outdir}/rfd3/rfdiffusion3", pattern: 'output/*.cif.gz', mode: 'copy'
    publishDir path: "${params.outdir}/rfd3/rfdiffusion3", pattern: 'output/*.json', mode: 'copy'

    input:
    path config_json
    path input_pdb
    val design_name
    val batch_size
    val n_batches
    val design_startnum
    val step_scale
    val gamma_0
    val extra_args

    output:
    path 'output/*.cif.gz', emit: cifs
    path 'output/*.json', emit: json_metrics
    path 'rfd3_*_batch*_backbone.tsv', emit: scores

    script:
    // Always set global_prefix so all batches share the {design_name}_ prefix.
    // rfd3 appends the config key after global_prefix, so we include a trailing
    // separator to keep the batch index visually distinct.
    def global_prefix = "global_prefix=${design_name}_b${design_startnum}_"
    def extra_args_str = extra_args ?: ''
    def uid = design_name.replaceFirst(/^rfd3_/, '')
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

    # Rewrite input paths in config to use staged basenames (copy to avoid modifying staged symlink)
    ${projectDir}/bin/rfd3/stage_rfd3_config.py stage ${config_json} -o ${design_name}.json

    # Run RFDiffusion3
    rfd3 design \
        out_dir=output \
        inputs=${design_name}.json \
        diffusion_batch_size=${batch_size} \
        n_batches=${n_batches} \
        inference_sampler.step_scale=${step_scale} \
        inference_sampler.gamma_0=${gamma_0} \
        ${global_prefix} \
        ${extra_args_str} \
        ${task.ext.args ?: ''}

    for j in output/*.json; do
      [[ -f "\$j" ]] || continue
      batch=\$(basename "\$j" .json | sed -n 's/.*batch\\([0-9]*\\).*/\1/p; s/.*_b\\([0-9]*\\)_.*/\1/p' | head -1)
      [[ -z "\$batch" ]] && batch=0
      python ${projectDir}/bin/rfd3/extract_rfd3_scores.py rfdiffusion3 "\$j" -o rfd3_${uid}_batch\${batch}_backbone.tsv
    done
    """
}
