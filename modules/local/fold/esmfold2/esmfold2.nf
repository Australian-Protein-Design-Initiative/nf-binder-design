// ESMFold2 inference for fold.nf's esmfold2 and esmfold2_fast --methods engines
// (meta.esmfold2_tool, set by ESMFOLD2_FOLD along with the model). The esm package
// ships no CLI, so this drives the documented Python API through
// bin/fold/run_esmfold2.py (see that script for the output layout).
//
// Weights are baked into the images under /weights/huggingface: biohub/ESMFold2
// (~25 GB incl. the ESMC-6B encoder) in _full_weights, biohub/ESMFold2-Fast in
// _fast_weights. --esmfold2_weights_dir points at an external HF cache instead.
process ESMFOLD2 {
    tag "${meta.id}${meta.fold_batch ? " batch${meta.fold_batch}" : ''}${meta.msa_depth_tag ? " msa${meta.msa_depth_tag}" : ''}"

    container "ghcr.io/australian-protein-design-initiative/containers/esmfold2:3.4.1.post1_nv-cuda13_${meta.esmfold2_tool == 'esmfold2_fast' ? 'fast' : 'full'}_weights"

    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/${meta.esmfold2_tool}",
        mode: 'copy',
        saveAs: { filename ->
            def rel = filename.toString().replaceFirst(/^output\//, '')
            if (meta.fold_namespaced) {
                def msa_bit = meta.msa_depth_tag ? "_msa_${meta.msa_depth_tag}" : ''
                return "${meta.id}/batch_${meta.fold_batch}${msa_bit}/${rel}"
            }
            return "${meta.id}/${rel}"
        }
    )
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/predictions",
        mode: 'copy',
        saveAs: { filename ->
            def bn = filename.toString().replaceFirst(/^.*\//, '')
            if (!(bn ==~ /.*_seed_\d+_sample_\d+_model\.cif/)) { return null }
            return "${FoldNaming.flatPrefix(meta.esmfold2_tool, meta)}${bn}"
        }
    )
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/msa_ids",
        mode: 'copy',
        pattern: '*_ids.txt'
    )

    input:
    tuple val(meta), path(fasta), path(a3ms)

    output:
    tuple val(meta), path('output/**'), emit: predictions
    path("*_ids.txt"), emit: msa_ids, optional: true
    tuple val(meta), path('output/*_summary_confidences.json'), emit: confidence_json

    script:
    def n_samples = meta.fold_batch_size ?: (params.esmfold2_batch_size ?: 5)
    def do_subsample = meta.msa_max_seq != null
    def write_msa_ids = meta.msa_depth_tag != null
    def batch_bit = meta.fold_namespaced ? "batch${meta.fold_batch}_" : ''
    def msa_ids_file = write_msa_ids \
        ? "${meta.esmfold2_tool}_${batch_bit}msa${meta.msa_depth_tag}_${meta.id}_ids.txt" \
        : ''
    def files = (a3ms instanceof List) ? a3ms : [a3ms]
    // chain_A/chain_B/... from GENERATE_ESMFOLD2_INPUT_COMPLEX, or the single
    // monomer a3m; sorted so chain order matches the FASTA record order.
    def a3m_names = files.collect { it.name }.sort()
    def a3m_arg = (meta.esmfold2_single_sequence || !a3m_names) \
        ? '--single-sequence' \
        : "--a3m ${a3m_names.collect { "'${it}'" }.join(' ')}"
    def opt_args = [
        params.esmfold2_num_loops ? "--num-loops ${params.esmfold2_num_loops}" : '',
        params.esmfold2_num_sampling_steps ? "--num-sampling-steps ${params.esmfold2_num_sampling_steps}" : '',
        params.esmfold2_msa_max_depth ? "--msa-max-depth ${params.esmfold2_msa_max_depth}" : '',
    ].findAll { it }.join(' ')
    def subsample_target = a3m_names ? a3m_names[0] : ''
    """
    set -euo pipefail

    if [[ -n "${params.esmfold2_weights_dir ?: ''}" ]]; then
        export HF_HOME="${params.esmfold2_weights_dir ?: ''}"
        export HF_HUB_OFFLINE=\${HF_HUB_OFFLINE:-1}
    fi

    # Optional CF-random-style MSA subsample (monomer only). The staged a3m is a
    # symlink into the MSA stage's work dir, so copy it before rewriting.
    if [[ -n "${subsample_target}" && ( "${do_subsample}" == "true" || "${write_msa_ids}" == "true" ) ]]; then
        cp -L "${subsample_target}" subsample_src.a3m
        if [[ "${do_subsample}" == "true" ]]; then
            python3 ${projectDir}/bin/fold/subsample_a3m.py \\
                --a3m subsample_src.a3m \\
                --max-seq ${meta.msa_max_seq} \\
                --max-extra-seq ${meta.msa_max_extra_seq} \\
                --seed ${meta.msa_subsample_seed} \\
                -o subsampled.a3m \\
                --ids-output "${msa_ids_file}"
            rm -f "${subsample_target}"
            mv subsampled.a3m "${subsample_target}"
        else
            python3 ${projectDir}/bin/fold/subsample_a3m.py \\
                --a3m subsample_src.a3m \\
                --ids-only \\
                --ids-output "${msa_ids_file}"
        fi
        rm -f subsample_src.a3m
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
    # (bin/gpu_lock.sh). See modules/local/fold/protenix/protenix_fold.nf.
    if [[ -n "${params.gpu_devices}" ]]; then
        source ${projectDir}/bin/gpu_lock.sh
        nfbd_acquire_gpu "${params.gpu_devices}" "${params.gpu_lock_dir ?: workDir.toString() + '/.gpu_locks'}" ${task.ext.gpu_slots ?: params.gpu_slots_per_device} ${params.gpu_lock_timeout} || exit 1
    else
        source ${projectDir}/bin/gpu_lock.sh || true
    fi
    nfbd_record_gpu_trace "${params.gpu_trace_dir ?: workDir.toString() + '/.gpu_trace'}" "${task.process}" || true

    python3 ${projectDir}/bin/fold/run_esmfold2.py \\
        --fasta ${fasta} \\
        --name '${meta.id}' \\
        --output-dir output \\
        --model '${meta.esmfold2_model}' \\
        --seed ${meta.esmfold2_seed} \\
        --num-diffusion-samples ${n_samples} \\
        --kernel-backend '${params.esmfold2_kernel_backend}' \\
        ${a3m_arg} \\
        ${opt_args} \\
        ${task.ext.args ?: ''}

    if ! compgen -G "output/*_model.cif" >/dev/null; then
        echo "ESMFold2 wrote no structures for ${meta.id} - see the log above" >&2
        exit 1
    fi
    """
}
