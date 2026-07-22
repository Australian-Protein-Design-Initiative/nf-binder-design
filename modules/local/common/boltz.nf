process BOLTZ {
    tag "${meta.id}${meta.fold_batch ? " batch${meta.fold_batch}" : ''}"
    container 'ghcr.io/australian-protein-design-initiative/containers/boltz:v2.2.1-2'
    // Recursive /** so each nested file is its own publish item. A directory
    // output is a single item - saveAs never sees *_model_N.cif inside it, so
    // the fold/predictions gather would publish nothing (same pattern as
    // modules/fold/af2/alphafold2.nf and modules/fold/rf3/rf3_fold.nf).
    publishDir(
        path: "${params.outdir}/${step_name}",
        mode: 'copy',
        saveAs: { filename ->
            meta.fold_namespaced ? "batch_${meta.fold_batch}/${filename}" : filename
        }
    )
    // Second publishDir (fold.nf only): gather per-sample structures into the
    // shared flat <outdir>/fold/predictions/ dir with a boltz_ prefix. Gated on
    // step_name (fold.nf passes 'fold/boltz'; boltz_pulldown passes
    // 'boltz_pulldown').
    publishDir(
        path: "${params.outdir}/fold/predictions",
        mode: 'copy',
        saveAs: { filename ->
            if (!step_name.toString().startsWith('fold/')) { return null }
            def bn = filename.toString().replaceFirst(/^.*\//, '')
            if (!(bn ==~ /.*_model_\d+\.(cif|pdb)/)) { return null }
            return meta.fold_namespaced ? "boltz_batch${meta.fold_batch}_${bn}" : "boltz_${bn}"
        }
    )

    input:
    tuple val(meta), path(yaml_file), path(target_msa), path(binder_msa)
    path templates
    val step_name
    // pdb | cif | mmcif — callers hardcode this (fold.nf -> cif, boltz_pulldown -> pdb).
    // File suffix is .cif for cif/mmcif; Boltz CLI uses --output_format mmcif for those.
    val output_format

    output:
    path ("boltz_results_${yaml_file.baseName}/**"), emit: results
    tuple val(meta), path("boltz_results_${yaml_file.baseName}/predictions/${yaml_file.baseName}/${yaml_file.baseName}_model_0.${ output_format.toString().toLowerCase() in ['cif', 'mmcif'] ? 'cif' : 'pdb' }"), emit: pdb
    // All per-sample structures (model_0..N-1); used by Engens / channel consumers.
    tuple val(meta), path("boltz_results_${yaml_file.baseName}/predictions/${yaml_file.baseName}/${yaml_file.baseName}_model_*.${ output_format.toString().toLowerCase() in ['cif', 'mmcif'] ? 'cif' : 'pdb' }"), emit: structure_all
    tuple val(meta), path("boltz_results_${yaml_file.baseName}/predictions/${yaml_file.baseName}/confidence_${yaml_file.baseName}_model_0.json"), emit: confidence_json
    // All per-sample confidence JSONs (model_0..N-1 when --diffusion_samples/-n_predictions > 1);
    // collapses to just model_0 for the default single-sample case, so this is
    // a non-breaking superset of confidence_json above.
    tuple val(meta), path("boltz_results_${yaml_file.baseName}/predictions/${yaml_file.baseName}/confidence_${yaml_file.baseName}_model_*.json"), emit: confidence_json_all
    tuple val(meta), path("boltz_results_${yaml_file.baseName}/predictions/${yaml_file.baseName}/*_ipsae.tsv"), emit: ipsae_tsv
    tuple val(meta), path("boltz_results_${yaml_file.baseName}/predictions/${yaml_file.baseName}/*_ipsae_byres.tsv"), emit: ipsae_byres_tsv

    script:
    def use_msa_server_flag = params.use_msa_server ? '--use_msa_server' : ''
    def args = task.ext.args ?: ''
    // fold.nf's --boltz_* knobs are resolved here in the module (reading params /
    // meta) rather than via a withName ext.args closure: ext.args set on the
    // subworkflow-qualified `BOLTZ_FOLD:BOLTZ` selector is silently dropped once
    // the m3 profile layers in its broad `withName: BOLTZ` selector (a Nextflow
    // profile-merge quirk). Per-job sample count and seed come from meta when
    // BOLTZ_FOLD fans --n_predictions across --boltz_batch_size jobs; params
    // remain unset for boltz_pulldown callers (no flag added - a no-op).
    def diffusion_samples = meta.fold_batch_size ?: params.boltz_batch_size
    def seed = meta.boltz_seed != null ? meta.boltz_seed : params.boltz_seed
    def fmt = output_format.toString().toLowerCase()
    def out_fmt = (fmt in ['cif', 'mmcif']) ? 'cif' : 'pdb'
    // Boltz CLI value for mmCIF is 'mmcif' (files are still written as .cif).
    def out_fmt_cli = out_fmt == 'cif' ? 'mmcif' : 'pdb'
    def boltz_flags = [
        params.boltz_recycling ? "--recycling_steps ${params.boltz_recycling}" : '',
        diffusion_samples ? "--diffusion_samples ${diffusion_samples}" : '',
        params.boltz_sampling_steps ? "--sampling_steps ${params.boltz_sampling_steps}" : '',
        seed ? "--seed ${seed}" : '',
    ].findAll { it }.join(' ')
    """
    # Find least-used GPU (by active processes and VRAM) and set CUDA_VISIBLE_DEVICES
    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
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

    BOLTZ_PREDICT_LOG=.boltz_predict_console.log
    rm -f "\$BOLTZ_PREDICT_LOG"
    set +e
    boltz predict \
        ${args} \
        ${boltz_flags} \
        ${use_msa_server_flag} \
        --preprocessing-threads ${task.cpus} \
        --num_workers ${task.cpus} \
        --output_format ${out_fmt_cli} \
        ${yaml_file} 2>&1 | tee "\$BOLTZ_PREDICT_LOG"
    boltz_rc=\${PIPESTATUS[0]}
    set -e
    if grep -qF 'ran out of memory, skipping batch' "\$BOLTZ_PREDICT_LOG"; then
        echo 'BOLTZ: Boltz logged GPU OOM (batch skipped); failing.' >&2
        exit 1
    fi
    if [[ "\$boltz_rc" -ne 0 ]]; then
        exit "\$boltz_rc"
    fi

    # Run ipSAE once per diffusion sample. Globbing pae_*.npz + structure into a
    # single ipsae.py call only works when there is exactly one model (ipsae.py
    # takes one pae + one structure positionally); with --diffusion_samples > 1
    # (e.g. via --n_predictions) there are N of each, so loop and pair each pae
    # npz with its matching structure and confidence JSON by model index. For the
    # single-sample default this loops exactly once - identical to before.
    for pae in boltz_results_*/predictions/*/pae_*.npz; do
        d=\$(dirname "\$pae")
        model=\$(basename "\$pae" .npz)      # pae_<name>_model_<i>
        model=\${model#pae_}                 # <name>_model_<i>
        ${projectDir}/bin/ipsae.py \\
            --update-summary "\$d/confidence_\${model}.json" \\
            --format boltz \\
            "\$pae" "\$d/\${model}.${out_fmt}" 10 10
    done
    """
}
