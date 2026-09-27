// AlphaFold3 inference for fold.nf's af3 --methods engine. MSAs come from the
// pipeline's shared MSA stage via the input JSON, so AF3's own data pipeline
// (and its ~630 GB of databases) is never run.
//
// The weights are not in the container (the AF3 terms forbid redistribution).
// --af3_model_dir is staged as a path input named af3_models, which makes the
// container engines auto-mount it wherever it lives on the host, and AF3 is
// pointed at it with --model_dir. AF3 reads af3.bin.zst directly.
process ALPHAFOLD3 {
    tag "${meta.id}${meta.fold_batch ? " batch${meta.fold_batch}" : ''}${meta.msa_depth_tag ? " msa${meta.msa_depth_tag}" : ''}"

    container 'ghcr.io/australian-protein-design-initiative/containers/alphafold3:3.0.4'

    // Recursive glob publish - see modules/fold/rf3/rf3_fold.nf. AF3 writes to
    // output/<sanitised name>/..., which is <meta.id>/ for ordinary ids.
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/af3",
        mode: 'copy',
        saveAs: { filename ->
            def rel = filename.toString().replaceFirst(/^output\//, '')
            def name = FoldNaming.af3Name(meta.id)
            if (meta.fold_namespaced && rel.startsWith("${name}/")) {
                def tail = rel.substring(name.length() + 1)
                def msa_bit = meta.msa_depth_tag ? "_msa_${meta.msa_depth_tag}" : ''
                return "${name}/batch_${meta.fold_batch}${msa_bit}/${tail}"
            }
            return rel
        }
    )
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/predictions",
        mode: 'copy',
        saveAs: { filename ->
            def bn = filename.toString().replaceFirst(/^.*\//, '')
            if (!(bn ==~ /.*_seed-\d+_sample-\d+_model\.cif/)) { return null }
            return "${FoldNaming.flatPrefix('af3', meta)}${bn}"
        }
    )
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/msa_ids",
        mode: 'copy',
        pattern: '*_ids.txt'
    )

    input:
    tuple val(meta), path(fasta), path(a3ms), path(af3_input_json)
    path af3_models, stageAs: 'af3_models'

    output:
    tuple val(meta), path('output/**'), emit: predictions
    path("*_ids.txt"), emit: msa_ids, optional: true
    tuple val(meta), path('output/**/*_summary_confidences.json'), emit: confidence_json

    script:
    def n_samples = meta.fold_batch_size ?: (params.af3_batch_size ?: 5)
    def do_subsample = meta.msa_max_seq != null
    def write_msa_ids = meta.msa_depth_tag != null
    def batch_bit = meta.fold_namespaced ? "batch${meta.fold_batch}_" : ''
    def msa_ids_file = write_msa_ids \
        ? "af3_${batch_bit}msa${meta.msa_depth_tag}_${meta.id}_ids.txt" \
        : ''
    def jax_cache_arg = params.af3_jax_cache_dir ? "--jax_compilation_cache_dir=${params.af3_jax_cache_dir}" : ''
    """
    set -euo pipefail

    # Optional CF-random-style MSA subsample (monomer only; the JSON references
    # chain_A_unpaired.a3m by basename, so replace the staged symlink in place)
    if [[ "${do_subsample}" == "true" ]]; then
        python3 ${projectDir}/bin/fold/subsample_a3m.py \\
            --a3m chain_A_unpaired.a3m \\
            --max-seq ${meta.msa_max_seq} \\
            --max-extra-seq ${meta.msa_max_extra_seq} \\
            --seed ${meta.msa_subsample_seed} \\
            -o subsampled.a3m \\
            --ids-output "${msa_ids_file}"
        rm -f chain_A_unpaired.a3m
        mv subsampled.a3m chain_A_unpaired.a3m
    elif [[ "${write_msa_ids}" == "true" ]]; then
        python3 ${projectDir}/bin/fold/subsample_a3m.py \\
            --a3m chain_A_unpaired.a3m \\
            --ids-only \\
            --ids-output "${msa_ids_file}"
    fi

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
    # (bin/gpu_lock.sh). See modules/fold/protenix/protenix_fold.nf.
    if [[ -n "${params.gpu_devices}" ]]; then
        source ${projectDir}/bin/gpu_lock.sh
        nfbd_acquire_gpu "${params.gpu_devices}" "${params.gpu_lock_dir ?: workDir.toString() + '/.gpu_locks'}" ${task.ext.gpu_slots ?: params.gpu_slots_per_device} ${params.gpu_lock_timeout} || exit 1
    else
        source ${projectDir}/bin/gpu_lock.sh || true
    fi
    nfbd_record_gpu_trace "${params.gpu_trace_dir ?: workDir.toString() + '/.gpu_trace'}" "${task.process}" || true

    # AF3 refuses to run on compute capability 7.x (V100, T4) unless XLA's custom
    # kernel fusion is disabled and flash attention uses the xla implementation;
    # otherwise it produces garbage structures. nvidia-smi ignores
    # CUDA_VISIBLE_DEVICES, so query the claimed card explicitly.
    flash_attention="${params.af3_flash_attention}"
    gpu_idx="\${CUDA_VISIBLE_DEVICES:-0}"
    gpu_idx="\${gpu_idx%%,*}"
    compute_cap=""
    if [[ "\${gpu_idx}" =~ ^[0-9]+\$ ]] && command -v nvidia-smi >/dev/null 2>&1; then
        compute_cap=\$(nvidia-smi -i "\${gpu_idx}" --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -n 1 || true)
    fi
    if [[ -n "\${compute_cap}" ]] && awk -v c="\${compute_cap}" 'BEGIN { exit !(c < 8.0) }'; then
        export XLA_FLAGS="\${XLA_FLAGS:-} --xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
        if [[ "\${flash_attention}" == "auto" ]]; then
            flash_attention=xla
        fi
    fi
    if [[ "\${flash_attention}" == "auto" ]]; then
        flash_attention=triton
    fi
    echo "GPU compute capability: \${compute_cap:-unknown}; flash attention: \${flash_attention}" >&2

    mkdir -p output

    python /app/alphafold/run_alphafold.py \\
        --json_path=${af3_input_json} \\
        --model_dir=af3_models \\
        --output_dir=output \\
        --force_output_dir \\
        --run_data_pipeline=false \\
        --num_diffusion_samples=${n_samples} \\
        --num_recycles=${params.af3_num_recycles} \\
        --flash_attention_implementation=\${flash_attention} \\
        ${jax_cache_arg} \\
        ${task.ext.args ?: ''}
    """
}
