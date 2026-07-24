// NEW, decoupled RF3 fold process for fold.nf - independent of the
// binder-design-coupled ROSETTAFOLD3 (modules/local/rfd3/rosettafold3.nf),
// which is hard-wired to a 2-body target+binder complex plus a binder
// postprocess step. This process folds whatever generic component list
// GENERATE_RF3_FOLD_INPUT built (monomer in Phase 1; multimer is Phase 2),
// with no postprocessing beyond what `rf3 fold` itself emits.
process RF3_FOLD {
    tag "${meta.id}${meta.fold_batch ? " batch${meta.fold_batch}" : ''}${meta.msa_depth_tag ? " msa${meta.msa_depth_tag}" : ''}"

    container 'oras://ghcr.io/australian-protein-design-initiative/containers/rc-foundry:0.2.0-weights'

    // Recursive glob publish - see modules/fold/af2/alphafold2.nf's comment
    // for why: publishDir pattern/saveAs only see the top-level output
    // *item* (a directory), never files nested inside it, so declaring the
    // output as "output/**" makes each nested file its own publish item.
    // rf3 fold writes to output/<name>/ where name == meta.id (see
    // generate_rf3_fold_input.nf --name '${meta.id}'), so just strip the
    // task-local "output/" prefix to land at <outdir>/fold/rf3/<meta.id>/... .
    // When fold.nf splits --n_predictions across jobs, nest under batch_N/ so
    // sample indices (which restart per job) do not clobber each other.
    publishDir(
        path: "${params.outdir}/fold/rf3",
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
    // Second publishDir: gather the per-sample mmCIF models into the shared flat
    // <outdir>/fold/predictions/ dir with an rf3_ prefix. Skip the top-level
    // merged best-model copy (no sample-N in its name).
    publishDir(
        path: "${params.outdir}/fold/predictions",
        mode: 'copy',
        saveAs: { filename ->
            def bn = filename.toString().replaceFirst(/^.*\//, '')
            if (!(bn ==~ /.*_sample-\d+_model\.cif/)) { return null }
            return "${FoldNaming.flatPrefix('rf3', meta)}${bn}"
        }
    )
    // Sequence IDs used in each MSA depth job (when --msa_subsample is on).
    publishDir(
        path: "${params.outdir}/fold/msa_ids",
        mode: 'copy',
        pattern: '*_ids.txt'
    )

    input:
    tuple val(meta), path(fasta), path(a3m), path(rf3_input_json)

    output:
    tuple val(meta), path('output/**'), emit: predictions
    tuple val(meta), path('output/**/*_summary_confidences.json'), emit: confidence_json
    path("*_ids.txt"), emit: msa_ids, optional: true

    script:
    // Per-job sample count from BOLTZ-style fan-out meta (fold.nf), else the
    // fold --rf3_batch_size / legacy default of 5 when neither n_predictions
    // nor batch size drove a meta override.
    def diffusion_batch_size = meta.fold_batch_size ?: (params.rf3_batch_size ?: 5)
    // Seed from meta when fold fans a pinned --rf3_seed across batches
    // (seed, seed+1, ...); else params; unset -> omit (RF3 draws its own,
    // task hash stays stable for -resume).
    def seed = meta.rf3_seed != null ? meta.rf3_seed : params.rf3_seed
    def seed_arg = seed ? "seed=${seed}" : ''
    def do_subsample = meta.msa_max_seq != null
    def write_msa_ids = meta.msa_depth_tag != null
    def batch_bit = meta.fold_namespaced ? "batch${meta.fold_batch}_" : ''
    def msa_ids_file = write_msa_ids \
        ? "rf3_${batch_bit}msa${meta.msa_depth_tag}_${meta.id}_ids.txt" \
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

    rf3 fold \\
        inputs=${rf3_input_json} \\
        out_dir=output \\
        ckpt_path=${params.rf3_ckpt_path} \\
        num_steps=${params.rf3_num_steps} \\
        n_recycles=${params.rf3_n_recycles} \\
        diffusion_batch_size=${diffusion_batch_size} \\
        ${seed_arg} \\
        early_stopping_plddt_threshold=${params.rf3_early_stopping_plddt_threshold} \\
        annotate_b_factor_with_plddt=true \\
        one_model_per_file=true \\
        ${task.ext.args ?: ''}
    """
}
