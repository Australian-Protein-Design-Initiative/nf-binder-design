process BOLTZ_COMPARE_BINDER_MONOMER {
    // Batched tasks are tagged by their first design plus a count, so the trace stays
    // readable and a task is still identifiable when a whole batch fails.
    tag { (metas instanceof List && metas.size() > 1) ? "${metas.first().id} +${metas.size() - 1}" : "${(metas instanceof List ? metas.first() : metas).id}" }
    container 'ghcr.io/australian-protein-design-initiative/containers/boltz:v2.2.1-2'
    publishDir "${params.outdir}/boltz_refold/predict/binder_monomer", pattern: 'boltz_results_*', mode: 'copy'
    publishDir "${params.outdir}/boltz_refold/rmsd/aligned_rmsd_monomer_vs_af2ig", pattern: 'aligned_rmsd_monomer_vs_af2ig/*.pdb', mode: 'copy', saveAs: { file(it).name }
    publishDir "${params.outdir}/boltz_refold/rmsd/aligned_rmsd_monomer_vs_complex", pattern: 'aligned_rmsd_monomer_vs_complex/*.pdb', mode: 'copy', saveAs: { file(it).name }

    input:
    // Numbered staging for the same reason as BOLTZ_COMPARE_COMPLEX: a Nextflow collection
    // will not stage two files sharing a basename, and RFD3 hands every design's RF3 output
    // over as `model.cif`. Order is preserved, which is what pairs each staged file with its
    // meta below. Design extensions come alongside, since the wildcard discards them; the
    // Boltz complex outputs are always .pdb.
    tuple val(metas), path(af2ig_pdbs, stageAs: 'design_inputs/input*'), path(boltz_complex_pdbs, stageAs: 'complex_inputs/complex*.pdb'), val(design_exts)
    val binder_chain

    output:
    path ("boltz_results_*"), emit: results
    tuple val(metas), path('per_design'), emit: per_design_bundle
    path ("aligned_rmsd_monomer_vs_af2ig/*.pdb"), emit: aligned_pdbs_monomer_vs_af2ig, optional: true
    path ("aligned_rmsd_monomer_vs_complex/*.pdb"), emit: aligned_pdbs_monomer_vs_complex, optional: true

    script:
    def metaList = metas instanceof List ? metas : [metas]
    def designList = af2ig_pdbs instanceof List ? af2ig_pdbs : [af2ig_pdbs]
    def complexList = boltz_complex_pdbs instanceof List ? boltz_complex_pdbs : [boltz_complex_pdbs]

    // One line per design, consumed by the shell loops below. Built here rather than as a
    // generated block of shell so the script length does not grow with the batch size.
    def extList = design_exts instanceof List ? design_exts : [design_exts]
    def manifest = [metaList, designList, complexList, extList]
        .transpose()
        .collect { m, d, c, e -> "${m.id}\t${d}\t${c}\t${e}" }
        .join('\n')

    def output_transformed_flag = params.output_rmsd_aligned ? "--output-transformed aligned_rmsd_monomer_vs_af2ig/" : ''
    def output_transformed_flag_complex = params.output_rmsd_aligned ? "--output-transformed aligned_rmsd_monomer_vs_complex/" : ''
    def args = task.ext.args ?: ''

    """
    set -euo pipefail

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

    cat > .batch_manifest.tsv <<'NFBD_MANIFEST_EOF'
${manifest}
NFBD_MANIFEST_EOF

    # Step 1: Create one Boltz YAML per design for binder monomer only
    # (no target specified = monomer mode), all into a single input directory
    # Give every staged design its id and real extension back; tools downstream dispatch
    # on the extension and rmsd4all needs a per-design unique basename.
    mkdir -p boltz_inputs staged
    while IFS=\$'\\t' read -r design_id staged_design complex_pdb design_ext; do
        [[ -n "\$design_id" ]] || continue
        ln -s "\$(readlink -f "\$staged_design")" "staged/\${design_id}.\${design_ext}"
    done < .batch_manifest.tsv

    while IFS=\$'\\t' read -r design_id staged_design complex_pdb design_ext; do
        [[ -n "\$design_id" ]] || continue
        design_pdb="staged/\${design_id}.\${design_ext}"
        /usr/bin/python3 ${projectDir}/bin/create_boltz_yaml.py \\
            --binder_id "\${design_id}_binder_monomer" \\
            --binder_from_pdb "\$design_pdb" \\
            --binder_chains '${binder_chain}' \\
            --binder_msa empty \\
            --output_yaml "boltz_inputs/\${design_id}_monomer.yml"
    done < .batch_manifest.tsv

    # Step 2: Run Boltz prediction once for the whole batch; fail if log shows GPU OOM batch
    # skip (boltz may still exit 0). Passing the directory rather than a single YAML is what
    # amortises model load and process startup over every design in the batch.
    BOLTZ_PREDICT_LOG=.boltz_predict_console.log
    rm -f "\$BOLTZ_PREDICT_LOG"
    set +e
    boltz predict \\
        ${args} \\
        --preprocessing-threads ${task.cpus} \\
        --num_workers ${task.cpus} \\
        --output_format pdb \\
        --out_dir boltz_batch \\
        boltz_inputs 2>&1 | tee "\$BOLTZ_PREDICT_LOG"
    boltz_rc=\${PIPESTATUS[0]}
    set -e
    if grep -qF 'ran out of memory, skipping batch' "\$BOLTZ_PREDICT_LOG"; then
        echo 'BOLTZ_COMPARE_BINDER_MONOMER: Boltz logged GPU OOM (batch skipped); failing.' >&2
        exit 1
    fi
    if [[ "\$boltz_rc" -ne 0 ]]; then
        exit "\$boltz_rc"
    fi

    # Boltz writes a batch to boltz_batch/boltz_results_boltz_inputs/predictions/<design_id>_monomer/.
    # Reshape that into the per-design boltz_results_<design_id>_monomer/ layout this module has
    # always emitted, so publishDir patterns and every downstream path are unchanged whatever
    # the batch size is.
    mkdir -p per_design
    while IFS=\$'\\t' read -r design_id staged_design complex_pdb design_ext; do
        [[ -n "\$design_id" ]] || continue
        design_pdb="staged/\${design_id}.\${design_ext}"

        src="boltz_batch/boltz_results_boltz_inputs/predictions/\${design_id}_monomer"
        if [[ ! -d "\$src" ]]; then
            echo "BOLTZ_COMPARE_BINDER_MONOMER: Boltz produced no prediction for '\${design_id}'." >&2
            exit 1
        fi
        results_dir="boltz_results_\${design_id}_monomer"
        mkdir -p "\${results_dir}/predictions"
        mv "\$src" "\${results_dir}/predictions/\${design_id}_monomer"

        # msa/, processed/ and lightning_logs/ describe the whole batch rather than one design,
        # and this module has always published them inside boltz_results_<design_id>_monomer.
        # Hard-linked, so there are no duplicated bytes and at the default batch size of 1 the
        # published tree is exactly what it was before batching.
        for aux in msa processed lightning_logs; do
            aux_src="boltz_batch/boltz_results_boltz_inputs/\${aux}"
            [[ -d "\$aux_src" ]] && cp -al "\$aux_src" "\${results_dir}/\${aux}"
        done

        MONOMER_PDB="\${results_dir}/predictions/\${design_id}_monomer/\${design_id}_monomer_model_0.pdb"

        # Step 3: Run RMSD calculations
        # Run RMSD: monomer vs AF2IG binder (chain A)
        rm -rf fixed mobile
        mkdir -p fixed mobile

        ln -s "\$(readlink -f "\$design_pdb")" "fixed/\$(basename "\$design_pdb")"
        ln -s "\$(readlink -f "\$MONOMER_PDB")" "mobile/\$(basename "\$MONOMER_PDB")"

        # Monomer Boltz output uses chain A (create_boltz_yaml binder-only mode). AF2IG/RFD3 binder is ${binder_chain}.
        /usr/bin/python3 ${projectDir}/bin/rmsd4all.py \\
            --tm-score \\
            --superimpose-chains ${binder_chain} \\
            --mobile-superimpose-chains A \\
            --score-chains ${binder_chain} \\
            --mobile-score-chains A \\
            ${output_transformed_flag} \\
            fixed/ mobile/ > "rmsd_monomer_vs_af2ig_\${design_id}.tsv"

        # Run RMSD: monomer vs Boltz complex binder (Boltz complex always labels binder chain A; see boltz_compare_complex.nf)
        rm -rf fixed mobile
        mkdir -p fixed mobile

        ln -s "\$(readlink -f "\$complex_pdb")" "fixed/\$(basename "\$complex_pdb")"
        ln -s "\$(readlink -f "\$MONOMER_PDB")" "mobile/\$(basename "\$MONOMER_PDB")"

        /usr/bin/python3 ${projectDir}/bin/rmsd4all.py \\
            --tm-score \\
            --superimpose-chains A \\
            --mobile-superimpose-chains A \\
            --score-chains A \\
            --mobile-score-chains A \\
            ${output_transformed_flag_complex} \\
            fixed/ mobile/ > "rmsd_monomer_vs_complex_\${design_id}.tsv"

        # Step 4: Parse confidence JSON
        /usr/bin/python3 ${projectDir}/bin/parse_boltz_confidence.py \\
            --json "\${results_dir}/predictions/\${design_id}_monomer/confidence_\${design_id}_monomer_model_0.json" \\
            --id "\${design_id}_monomer" \\
            --binder "\${design_id}_binder_monomer" > "confidence_monomer_\${design_id}.tsv"

        # Per-design bundle. The workflow re-associates these with their meta by design id,
        # which is what keeps batched outputs correctly paired (see BOLTZ_REFOLD_CORE).
        # Hard links rather than copies: same filesystem, no duplicated bytes, and unlike
        # symlinks they stay valid wherever Nextflow stages the bundle.
        d="per_design/\${design_id}"
        mkdir -p "\$d"
        ln "\$MONOMER_PDB" "\$d/\${design_id}_monomer.pdb"
        ln "confidence_monomer_\${design_id}.tsv" "\$d/confidence.tsv"
        ln "rmsd_monomer_vs_af2ig_\${design_id}.tsv" "\$d/rmsd_monomer_vs_af2ig.tsv"
        ln "rmsd_monomer_vs_complex_\${design_id}.tsv" "\$d/rmsd_monomer_vs_complex.tsv"
    done < .batch_manifest.tsv

    rm -rf fixed mobile boltz_batch staged
    """
}
