process RFD4_GENERATE {
    tag "batch_${batch_index}"

    container 'ghcr.io/australian-protein-design-initiative/containers/rfd4-proteina:0.0.3-aaccdd2_nv-cuda13-arm64'

    // saveAs strips the task-local 'output/' prefix so designs land directly in designs/
    publishDir path: "${params.outdir}/rfd4_proteina/designs", pattern: 'output/*', mode: 'copy', saveAs: { fn -> file(fn).name }
    publishDir path: "${params.outdir}/rfd4_proteina/designs", pattern: 'rfd4_designs_batch*.tsv', mode: 'copy'

    input:
    path config_spec
    path input_structures
    val design_name
    val batch_index
    val n_replicates
    val batch_size
    val seed
    val extra_args

    output:
    path 'output/*', emit: cifs
    path "rfd4_designs_batch${batch_index}.tsv", emit: design_table
    path 'rundir/inference_metrics/per_example_metrics.csv', emit: per_example_metrics, optional: true

    script:
    def nsteps_arg = params.rfd4_nsteps ? "generation.args.nsteps=${params.rfd4_nsteps}" : ''
    def recycle_arg = params.rfd4_recycle_strategy ? "recycle_strategy=${params.rfd4_recycle_strategy}" : ''
    def experiment = params.rfd4_inference_experiment ? params.rfd4_inference_experiment.toString() : ''
    def experiment_arg = experiment ? "inference_experiment=${experiment}" : ''
    // '~metrics' drops the metrics group so no refolding oracle is constructed, which otherwise
    // needs an AlphaFold 3 / OpenFold3 environment that is not in this image. Hydra errors if the
    // group was never composed, and the quickstart preset does not compose it, so it is omitted
    // there. Refolding is done downstream with Boltz-2 instead (--full_refold_with boltz).
    def metrics_arg = params.rfd4_metrics_tags
        ? "metrics_tags='[${params.rfd4_metrics_tags}]'"
        : (experiment == 'quickstart' ? '' : "'~metrics'")
    def compile_arg = rfd4IsFalse(params.rfd4_compile_lengths) ? 'compile_lengths=null' : ''
    def partial_arg = rfd4IsFalse(params.rfd4_partial_trajectory_fraction)
        ? ''
        : "generation.args.partial_trajectory_fraction_steps=${params.rfd4_partial_trajectory_fraction}"
    def inverse_folding_arg = rfd4IsFalse(params.rfd4_inverse_folding) ? '' : 'generation=inverse_folding'
    def extra_args_str = extra_args ?: ''
    def spec_name = "${design_name}_b${batch_index}.json"
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

    # Caches that default to \$HOME, which is not writable in every container runtime
    export TORCHINDUCTOR_CACHE_DIR="\$(pwd)/.inductor_cache"
    export TRITON_CACHE_DIR="\$(pwd)/.triton_cache"
    export CUEQ_TRITON_CACHE_DIR="\$(pwd)/.cuequivariance_triton_cache"
    export XDG_CONFIG_HOME="\$(pwd)/.config"
    mkdir -p "\$TORCHINDUCTOR_CACHE_DIR" "\$TRITON_CACHE_DIR" "\$CUEQ_TRITON_CACHE_DIR" "\$XDG_CONFIG_HOME"

    # Weights and Biases is a hard dependency of rfproteina but has nothing to log to here
    export WANDB_MODE=disabled

    # Rewrite input paths in the spec to the staged basenames
    ${projectDir}/bin/rfd4_proteina/rfd4_config.py stage ${config_spec} -o ${spec_name}

    rfd4-generate \\
        inputs=${spec_name} \\
        out_dir=rundir \\
        num_replicates=${n_replicates} \\
        batch_size=${batch_size} \\
        num_workers=${params.rfd4_num_workers} \\
        seed=${seed} \\
        hydra.run.dir=./hydra_outputs \\
        ${metrics_arg} \\
        ${experiment_arg} \\
        ${nsteps_arg} \\
        ${recycle_arg} \\
        ${compile_arg} \\
        ${partial_arg} \\
        ${inverse_folding_arg} \\
        ${extra_args_str} \\
        ${task.ext.args ?: ''}

    ${projectDir}/bin/rfd4_proteina/collect_rfd4_designs.py rundir \\
        --design-name "${design_name}" \\
        --batch ${batch_index} \\
        --seed ${seed} \\
        -o output \\
        --tsv rfd4_designs_batch${batch_index}.tsv
    """
}

// Compared as a string: Nextflow gives a boolean for `--flag false` but the literal string when
// the same value comes from a -params-file or a config.
boolean rfd4IsFalse(def value) {
    return value == null || value.toString().toLowerCase() in ['false', '0', '']
}

/*
Usage:

https://rfd4-proteina.org/latest/reference/cli-options.html
https://github.com/RosettaCommons/RFD4-Proteina

rfd4-generate is a Hydra application: every option is `key=value`, never `--key value`.

Key options used above:

  inputs                  Annotated CIF/BCIF, a JSON/YAML condition spec, or a directory of them
  out_dir                 Output root; designs land in <out_dir>/inference_samples/<example_id>/
  num_replicates          Designs generated per input example
  batch_size              Designs processed in parallel on the GPU
  seed                    Random seed; null (the default) means a fresh seed each run
  inference_experiment    Experiment preset composing generation, search and metrics.
                          `quickstart` skips metrics and the torch.compile warm-up.
  generation.args.nsteps  Integration steps in the sampling schedule (default 60)
  recycle_strategy        Trunk recycle schedule, e.g. thinking_medium, high/midpeak
  compile_lengths         Token buckets padded to for the compiled trunk; null runs it eagerly
  generation.args.partial_trajectory_fraction_steps
                          Partial diffusion: run only the final fraction of nsteps
  generation=inverse_folding
                          Hold the input backbone fixed and resample sequence/side chains
*/
