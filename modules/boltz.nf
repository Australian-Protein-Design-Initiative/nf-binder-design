process BOLTZ {
    tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/boltz:v2.1.1'
    publishDir "${params.outdir}/boltz_pulldown", mode: 'copy'

    input:
    tuple val(meta), path(yaml_file), path(target_msa), path(binder_msa)
    path(templates)
    val gpu_device

    output:
    path ("boltz_results_${meta.id}"), emit: results
    tuple val(meta), path("boltz_results_${meta.id}/predictions/${meta.id}/confidence_${meta.id}_model_0.json"), emit: confidence_json

    script:
    def use_msa_server_flag = params.use_msa_server ? "--use_msa_server" : ""
    def args = task.ext.args ?: ''
    """
    if [[ ${gpu_device} != "all" ]]; then
        export CUDA_VISIBLE_DEVICES=${gpu_device}
    fi

    # Boltz model weights are stored in our container
    export BOLTZ_CACHE=/app/boltz/cache

    # Create various tmp/cache directories that are expected to be in \$HOME by default
    export NUMBA_CACHE_DIR="\$(pwd)/.numba_cache"
    mkdir -p \$NUMBA_CACHE_DIR
    export XDG_CONFIG_HOME="\$(pwd)/.config"
    mkdir -p \$XDG_CONFIG_HOME
    export TRITON_CACHE_DIR="\$(pwd)/.triton_cache"
    mkdir -p \$TRITON_CACHE_DIR

    # Prevent Python from using ~/.local/lib/ packages mounted inside the container
    export PYTHONNOUSERSITE=1

    # We could autodetect if we have a GPU, but lets leave this up to task.ext.args
    # instead of using --accelerator \${ACCELERATOR} \
    # if nvidia-smi >/dev/null 2>&1; then
    #    ACCELERATOR=gpu
    # else
    #     ACCELERATOR=cpu
    # fi

    boltz predict \
        ${args} \
        ${use_msa_server_flag} \
        --preprocessing-threads ${task.cpus} \
        --num_workers ${task.cpus} \
        ${yaml_file}

        #--out_dir "results"
        #--output_format pdb
    """
    /*
    $ boltz predict --help

    boltz-latest.img boltz predict --help
    Usage: boltz predict [OPTIONS] DATA

    Run predictions with Boltz.

    Options:
    --out_dir PATH                  The path where to save the predictions.
    --cache PATH                    The directory where to download the data and
                                    model. Default is ~/.boltz, or $BOLTZ_CACHE
                                    if set.
    --checkpoint PATH               An optional checkpoint, will use the
                                    provided Boltz-1 model by default.
    --devices INTEGER               The number of devices to use for prediction.
                                    Default is 1.
    --accelerator [gpu|cpu|tpu]     The accelerator to use for prediction.
                                    Default is gpu.
    --recycling_steps INTEGER       The number of recycling steps to use for
                                    prediction. Default is 3.
    --sampling_steps INTEGER        The number of sampling steps to use for
                                    prediction. Default is 200.
    --diffusion_samples INTEGER     The number of diffusion samples to use for
                                    prediction. Default is 1.
    --max_parallel_samples INTEGER  The maximum number of samples to predict in
                                    parallel. Default is None.
    --step_scale FLOAT              The step size is related to the temperature
                                    at which the diffusion process samples the
                                    distribution. The lower the higher the
                                    diversity among samples (recommended between
                                    1 and 2). Default is 1.638 for Boltz-1 and
                                    1.5 for Boltz-2. If not provided, the
                                    default step size will be used.
    --write_full_pae                Whether to dump the pae into a npz file.
                                    Default is True.
    --write_full_pde                Whether to dump the pde into a npz file.
                                    Default is False.
    --output_format [pdb|mmcif]     The output format to use for the
                                    predictions. Default is mmcif.
    --num_workers INTEGER           The number of dataloader workers to use for
                                    prediction. Default is 2.
    --override                      Whether to override existing found
                                    predictions. Default is False.
    --seed INTEGER                  Seed to use for random number generator.
                                    Default is None (no seeding).
    --use_msa_server                Whether to use the MMSeqs2 server for MSA
                                    generation. Default is False.
    --msa_server_url TEXT           MSA server url. Used only if
                                    --use_msa_server is set.
    --msa_pairing_strategy TEXT     Pairing strategy to use. Used only if
                                    --use_msa_server is set. Options are
                                    'greedy' and 'complete'
    --use_potentials                Whether to not use potentials for steering.
                                    Default is False.
    --model [boltz1|boltz2]         The model to use for prediction. Default is
                                    boltz2.
    --method TEXT                   The method to use for prediction. Default is
                                    None.
    --preprocessing-threads INTEGER
                                    The number of threads to use for
                                    preprocessing. Default is 1.
    --affinity_mw_correction        Whether to add the Molecular Weight
                                    correction to the affinity value head.
    --sampling_steps_affinity INTEGER
                                    The number of sampling steps to use for
                                    affinity prediction. Default is 200.
    --diffusion_samples_affinity INTEGER
                                    The number of diffusion samples to use for
                                    affinity prediction. Default is 5.
    --affinity_checkpoint PATH      An optional checkpoint, will use the
                                    provided Boltz-1 model by default.
    --max_msa_seqs INTEGER          The maximum number of MSA sequences to use
                                    for prediction. Default is 8192.
    --subsample_msa                 Whether to subsample the MSA. Default is
                                    True.
    --num_subsampled_msa INTEGER    The number of MSA sequences to subsample.
                                    Default is 1024.
    --no_trifast                    Whether to not use trifast kernels for
                                    triangular updates. Default False
    --help                          Show this message and exit.
*/

}
