// Protenix (AF3-style) folding process for fold.nf's 4th --methods engine
// (see plans/fold-nf-multi-method-folding.md Phase 3). Weights are baked into
// the container under /models/protenix/checkpoint (confirmed present -
// protenix_base_default_v1.0.0.pt et al. - and /models/protenix/common, so
// `protenix pred` needs no download at predict time, same as rf3's
// rc-foundry image).
process PROTENIX_FOLD {
    tag "${meta.id}${meta.fold_batch ? " batch${meta.fold_batch}" : ''}${meta.msa_depth_tag ? " msa${meta.msa_depth_tag}" : ''}"

    container 'oras://ghcr.io/australian-protein-design-initiative/containers/protenix:v2.0.0-weights'

    // Recursive glob publish - see modules/fold/rf3/rf3_fold.nf's comment for
    // why (publishDir pattern/saveAs only see the top-level output *item*,
    // never files nested inside it). protenix pred writes to
    // out_dir/<name>/seed_<seed>/predictions/... where name == meta.id (see
    // generate_protenix_input.nf's --name '${meta.id}'; confirmed against
    // runner/inference.py's infer_predict(), which uses each JSON job's
    // top-level "name" field as both dump dataset dir and sample_name), so
    // strip only the task-local "output/" prefix to land at
    // <outdir>/fold/protenix/<meta.id>/... . When fold.nf splits --n_predictions
    // across jobs, nest under batch_N/ so sample indices do not collide (seed
    // dirs alone are not enough when seeds are unset and collide).
    publishDir(
        path: "${params.outdir}/fold/protenix",
        mode: 'copy',
        saveAs: { filename ->
            def rel = filename.toString().replaceFirst(/^output\//, '')
            if (meta.fold_namespaced && rel.startsWith("${meta.id}/")) {
                def tail = rel.substring(meta.id.toString().length() + 1)
                def msa_bit = meta.msa_depth_tag ? "_msa_${meta.msa_depth_tag}" : ''
                return "${meta.id}/batch_${meta.fold_batch}${msa_bit}/${tail}"
            }
            return rel
        }
    )
    // Second publishDir: gather the per-sample mmCIF predictions into the shared
    // flat <outdir>/fold/predictions/ dir with a protenix_ prefix.
    publishDir(
        path: "${params.outdir}/fold/predictions",
        mode: 'copy',
        saveAs: { filename ->
            def bn = filename.toString().replaceFirst(/^.*\//, '')
            if (!(bn ==~ /.*_sample_\d+\.cif/)) { return null }
            def msa_bit = meta.msa_depth_tag ? "_msa${meta.msa_depth_tag}" : ''
            return meta.fold_namespaced \
                ? "protenix_batch${meta.fold_batch}${msa_bit}_${bn}" \
                : "protenix${msa_bit}_${bn}"
        }
    )
    // Sequence IDs used in each MSA depth job (when --msa_subsample is on).
    publishDir(
        path: "${params.outdir}/fold/msa_ids",
        mode: 'copy',
        pattern: '*_ids.txt'
    )

    input:
    tuple val(meta), path(fasta), path(a3m), path(protenix_input_json)

    output:
    tuple val(meta), path('output/**'), emit: predictions
    path("*_ids.txt"), emit: msa_ids, optional: true
    tuple val(meta), path('output/**/*_summary_confidence_sample_*.json'), emit: confidence_json

    script:
    // Per-job sample count from fold fan-out meta, else --protenix_batch_size /
    // default 5 when neither n_predictions nor batch size set a meta override.
    def sample = meta.fold_batch_size ?: (params.protenix_batch_size ?: 5)
    // Seed from meta when fold fans a pinned --protenix_seeds across batches;
    // else params; unset -> omit (protenix default, -resume-stable).
    def seeds = meta.protenix_seeds != null ? meta.protenix_seeds : params.protenix_seeds
    def seeds_arg = seeds ? "--seeds ${seeds}" : ''
    def do_subsample = meta.msa_max_seq != null
    def write_msa_ids = meta.msa_depth_tag != null
    def batch_bit = meta.fold_namespaced ? "batch${meta.fold_batch}_" : ''
    def msa_ids_file = write_msa_ids \
        ? "protenix_${batch_bit}msa${meta.msa_depth_tag}_${meta.id}_ids.txt" \
        : ''
    """
    set -euo pipefail

    # Optional CF-random-style MSA subsample (JSON uses a3m basename)
    if [[ "${do_subsample}" == "true" ]]; then
        python3 ${projectDir}/bin/subsample_a3m.py \
            --a3m "${a3m}" \
            --max-seq ${meta.msa_max_seq} \
            --max-extra-seq ${meta.msa_max_extra_seq} \
            --seed ${meta.msa_subsample_seed} \
            -o subsampled.a3m \
            --ids-output "${msa_ids_file}"
        rm -f "${a3m}"
        mv subsampled.a3m "${a3m}"
    elif [[ "${write_msa_ids}" == "true" ]]; then
        python3 ${projectDir}/bin/subsample_a3m.py \
            --a3m "${a3m}" \
            --ids-only \
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

    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${projectDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    mkdir -p output

    protenix pred \\
        --input ${protenix_input_json} \\
        --out_dir output \\
        ${seeds_arg} \\
        --cycle ${params.protenix_cycle} \\
        --step ${params.protenix_step} \\
        --sample ${sample} \\
        --model_name ${params.protenix_model_name} \\
        --use_msa ${params.protenix_use_msa} \\
        ${task.ext.args ?: ''}
    """
}
