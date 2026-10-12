// Chai-1 inference for fold.nf's chai1 --methods engine, driven through
// bin/fold/run_chai1.py (the chai-lab CLI does not save the PAE). The weights
// are baked into the image under $CHAI_DOWNLOADS_DIR (/models/chai). MSAs come
// from the shared MSA stage and templates from --templates; Chai-1's MSA and
// template servers are never used.
process CHAI1 {
    tag "${meta.id}${meta.fold_batch ? " batch${meta.fold_batch}" : ''}${meta.msa_depth_tag ? " msa${meta.msa_depth_tag}" : ''}"

    container 'ghcr.io/australian-protein-design-initiative/containers/chai-1:2026-06-30_66c38d1_nv-cuda12_weights'

    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/chai1",
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
            return "${FoldNaming.flatPrefix('chai1', meta)}${bn}"
        }
    )
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/msa_ids",
        mode: 'copy',
        pattern: '*_ids.txt'
    )

    input:
    tuple val(meta), path(fasta), path(a3ms)
    path templates

    output:
    tuple val(meta), path('output/**'), emit: predictions
    path("*_ids.txt"), emit: msa_ids, optional: true
    tuple val(meta), path('output/*_summary_confidences.json'), emit: confidence_json

    script:
    def n_samples = meta.fold_batch_size ?: (params.chai1_batch_size ?: 5)
    def do_subsample = meta.msa_max_seq != null
    def write_msa_ids = meta.msa_depth_tag != null
    def batch_bit = meta.fold_namespaced ? "batch${meta.fold_batch}_" : ''
    def msa_ids_file = write_msa_ids ? "chai1_${batch_bit}msa${meta.msa_depth_tag}_${meta.id}_ids.txt" : ''
    // Monomer: the one shared a3m. Multimer: the Protenix bundle's
    // *.protenix_unpaired.a3m files (original headers, so taxids survive), whose
    // bundle order is chain order.
    def files = (a3ms instanceof List) ? a3ms : [a3ms]
    def a3m_names = (meta.n_chains ?: 1) > 1 \
        ? files.findAll { it.name.endsWith('.protenix_unpaired.a3m') }.collect { it.name } \
        : files.collect { it.name }
    def template_chains = (meta.template_chains ?: []) as List
    def templates_arg = params.templates \
        ? "--templates-dir ${templates}" + (template_chains ? " --template-chains ${template_chains.join(' ')}" : '') \
        : ''
    def opt_args = [
        meta.chai1_seed != null ? "--seed ${meta.chai1_seed}" : '',
        params.chai1_num_trunk_recycles ? "--num-trunk-recycles ${params.chai1_num_trunk_recycles}" : '',
        params.chai1_num_diffn_timesteps ? "--num-diffn-timesteps ${params.chai1_num_diffn_timesteps}" : '',
    ].findAll { it }.join(' ')
    def subsample_target = a3m_names ? a3m_names[0] : ''
    """
    set -euo pipefail

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

    python3 ${projectDir}/bin/fold/run_chai1.py \\
        --fasta ${fasta} \\
        --name '${meta.id}' \\
        --output-dir output \\
        --num-diffusion-samples ${n_samples} \\
        ${a3m_names ? "--a3m ${a3m_names.collect { "'${it}'" }.join(' ')}" : ''} \\
        ${templates_arg} \\
        ${opt_args} \\
        ${task.ext.args ?: ''}

    if ! compgen -G "output/*_model.cif" >/dev/null; then
        echo "Chai-1 wrote no structures for ${meta.id} - see the log above" >&2
        exit 1
    fi
    """
}
