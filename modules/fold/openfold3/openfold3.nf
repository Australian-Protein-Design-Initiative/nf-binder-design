// OpenFold3 inference for fold.nf's openfold3 --methods engine. The weights are
// baked into the image under $OPENFOLD_CACHE (/models/openfold3), which
// run_openfold searches by default. MSAs come from the pipeline's shared MSA
// stage via the query JSON (--use-msa-server false), and templates are off.
process OPENFOLD3 {
    tag "${meta.id}${meta.fold_batch ? " batch${meta.fold_batch}" : ''}${meta.msa_depth_tag ? " msa${meta.msa_depth_tag}" : ''}"

    container 'ghcr.io/australian-protein-design-initiative/containers/openfold3:0.5.0_nv-cuda12_weights'

    // Recursive glob publish - see modules/fold/rf3/rf3_fold.nf. OpenFold3 writes
    // to output/<query name>/seed_<S>/..., which is <meta.id>/ for ordinary ids.
    // Its run-level files (experiment_config.json, summary.txt, ...) sit at the
    // output root, so they go under the query dir to stop tasks overwriting them.
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/openfold3",
        mode: 'copy',
        saveAs: { filename ->
            def rel = filename.toString().replaceFirst(/^output\//, '')
            def name = FoldNaming.openfold3Name(meta.id)
            def tail = rel.startsWith("${name}/") ? rel.substring(name.length() + 1) : rel
            if (meta.fold_namespaced) {
                def msa_bit = meta.msa_depth_tag ? "_msa_${meta.msa_depth_tag}" : ''
                return "${name}/batch_${meta.fold_batch}${msa_bit}/${tail}"
            }
            return "${name}/${tail}"
        }
    )
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/predictions",
        mode: 'copy',
        saveAs: { filename ->
            def bn = filename.toString().replaceFirst(/^.*\//, '')
            if (!(bn ==~ /.*_seed_\d+_sample_\d+_model\.cif/)) { return null }
            return "${FoldNaming.flatPrefix('openfold3', meta)}${bn}"
        }
    )
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/msa_ids",
        mode: 'copy',
        pattern: '*_ids.txt'
    )

    input:
    tuple val(meta), path(fasta), path(msa_dirs), path(query_json)

    output:
    tuple val(meta), path('output/**'), emit: predictions
    path("*_ids.txt"), emit: msa_ids, optional: true
    tuple val(meta), path('output/**/*_confidences_aggregated.json'), emit: confidence_json

    script:
    def n_samples = meta.fold_batch_size ?: (params.openfold3_batch_size ?: 5)
    def do_subsample = meta.msa_max_seq != null
    def write_msa_ids = meta.msa_depth_tag != null
    def batch_bit = meta.fold_namespaced ? "batch${meta.fold_batch}_" : ''
    def msa_ids_file = write_msa_ids \
        ? "openfold3_${batch_bit}msa${meta.msa_depth_tag}_${meta.id}_ids.txt" \
        : ''
    """
    set -euo pipefail

    # The image activates its pixi env (the only python3) in the Docker
    # entrypoint, which Nextflow bypasses. Its activation also points the Triton /
    # torch extension caches inside the read-only image; those are moved below.
    if [[ -f /opt/activate.sh ]]; then
        set +u
        source /opt/activate.sh
        set -u
    fi

    # Optional CF-random-style MSA subsample (monomer only). The staged MSA dir
    # is a symlink into the generator's work dir, so copy it before rewriting.
    if [[ "${do_subsample}" == "true" || "${write_msa_ids}" == "true" ]]; then
        msa_dir=\$(ls -d msa_* | head -n 1)
        cp -rL "\${msa_dir}" "\${msa_dir}.copy"
        rm -f "\${msa_dir}"
        mv "\${msa_dir}.copy" "\${msa_dir}"
        if [[ "${do_subsample}" == "true" ]]; then
            python3 ${projectDir}/bin/fold/subsample_a3m.py \\
                --a3m "\${msa_dir}/colabfold_main.a3m" \\
                --max-seq ${meta.msa_max_seq} \\
                --max-extra-seq ${meta.msa_max_extra_seq} \\
                --seed ${meta.msa_subsample_seed} \\
                -o subsampled.a3m \\
                --ids-output "${msa_ids_file}"
            mv subsampled.a3m "\${msa_dir}/colabfold_main.a3m"
        else
            python3 ${projectDir}/bin/fold/subsample_a3m.py \\
                --a3m "\${msa_dir}/colabfold_main.a3m" \\
                --ids-only \\
                --ids-output "${msa_ids_file}"
        fi
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

    kernel_cache="${params.openfold3_kernel_cache_dir ?: ''}"
    kernel_cache="\${kernel_cache:-\${PWD}/.kernel_cache}"
    export TRITON_CACHE_DIR="\${kernel_cache}/triton"
    export TORCH_EXTENSIONS_DIR="\${kernel_cache}/torch_extensions"
    mkdir -p "\${TRITON_CACHE_DIR}" "\${TORCH_EXTENSIONS_DIR}"

    # OpenFold3 takes seeds from the runner YAML, not the query JSON; the default
    # is a fixed [42], so each batch gets its own seed from OPENFOLD3_FOLD. Its
    # default of 10 DataLoader workers ignores the task's CPU allocation.
    printf 'experiment_settings:\\n  seeds: [%s]\\ndata_module_args:\\n  num_workers: %s\\n' \\
        "${meta.openfold3_seed}" "${task.cpus}" >runner.yml

    run_openfold predict \\
        --query-json ${query_json} \\
        --runner-yaml runner.yml \\
        --output-dir output \\
        --use-msa-server false \\
        --use-templates false \\
        --num-diffusion-samples ${n_samples} \\
        ${task.ext.args ?: ''}

    # A failed query is logged and skipped rather than failing the run.
    if ! compgen -G "output/*/seed_*/*_model.cif" >/dev/null; then
        echo "OpenFold3 wrote no structures for ${meta.id} - see the log above" >&2
        exit 1
    fi
    """
}
