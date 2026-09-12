process BOLTZ_COMPARE_COMPLEX {
    // Batched tasks are tagged by their first design plus a count, so the trace stays
    // readable and a task is still identifiable when a whole batch fails.
    tag { (metas instanceof List && metas.size() > 1) ? "${metas.first().id} +${metas.size() - 1}" : "${(metas instanceof List ? metas.first() : metas).id}" }
    container 'ghcr.io/australian-protein-design-initiative/containers/boltz:v2.2.1-2'
    publishDir "${params.outdir}/boltz_refold/predict/complex", pattern: 'boltz_results_*', mode: 'copy'
    publishDir "${params.outdir}/boltz_refold/rmsd/aligned_rmsd_target_aligned_binder", pattern: 'aligned_rmsd_target_aligned_binder/*.pdb', mode: 'copy', saveAs: { file(it).name }
    publishDir "${params.outdir}/boltz_refold/rmsd/aligned_rmsd_complex", pattern: 'aligned_rmsd_complex/*.pdb', mode: 'copy', saveAs: { file(it).name }

    input:
    // Staged as design_inputs/input1, input2, ... rather than under their own names:
    // a Nextflow collection will not stage two files with the same basename, and RFD3
    // hands every design's RF3 output over as `model.cif`. The numbered wildcard keeps
    // input order, which is what pairs each staged file with its meta below. The real
    // extensions come alongside in pdb_exts, since the wildcard discards them.
    tuple val(metas), path(pdbs, stageAs: 'design_inputs/input*'), val(pdb_exts)
    val binder_chain
    val target_chain
    val create_target_msa
    val use_msa_server
    path target_msa
    path binder_msa
    path templates
    path refold_target_fasta

    output:
    path ("boltz_results_*"), emit: results
    tuple val(metas), path('per_design'), emit: per_design_bundle
    path ("aligned_rmsd_target_aligned_binder/*.pdb"), emit: aligned_pdbs_target_aligned_binder, optional: true
    path ("aligned_rmsd_complex/*.pdb"), emit: aligned_pdbs_complex, optional: true

    script:
    def metaList = metas instanceof List ? metas : [metas]
    def pdbList = pdbs instanceof List ? pdbs : [pdbs]

    // One line per design, consumed by the shell loops below. Built here rather
    // than as a generated block of shell so the script length does not grow with
    // the batch size.
    def extList = pdb_exts instanceof List ? pdb_exts : [pdb_exts]
    def manifest = [metaList, pdbList, extList]
        .transpose()
        .collect { m, p, e -> "${m.id}\t${p}\t${e}" }
        .join('\n')

    def use_msa_server_flag = use_msa_server ? '--use_msa_server' : ''
    def templates_flag = templates ? "--templates '${templates}'" : ''

    def output_transformed_flag_target = params.output_rmsd_aligned ? "--output-transformed aligned_rmsd_target_aligned_binder/" : ''
    def output_transformed_flag_complex = params.output_rmsd_aligned ? "--output-transformed aligned_rmsd_complex/" : ''
    def args = task.ext.args ?: ''

    // Determine MSA flags - swap these because we're swapping target/binder roles
    def target_msa_flag = '--binder_msa empty'
    // Binder MSA (always empty)

    def binder_msa_flag = ''
    // Target MSA (may be created or from server)
    if (create_target_msa && use_msa_server) {
        binder_msa_flag = ''
    }
    else if (create_target_msa && !use_msa_server) {
        binder_msa_flag = ''
    }
    else {
        binder_msa_flag = '--target_msa empty'
    }

    // Target sequence for YAML "binder" slot (Boltz output chain B); see create_boltz_yaml.py protein ids A/B.
    // When the target comes from the design structure it differs per design, so it is read from the
    // loop variable rather than interpolated here.
    def target_source_args = ''
    if (refold_target_fasta.name != 'empty') {
        target_source_args = "--binder_from_fasta '${refold_target_fasta}' --binder_chains '${target_chain}'"
    }
    else {
        target_source_args = "--binder_from_pdb \"\$pdb_path\" --binder_chains '${target_chain}'"
    }

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

    # Create one Boltz YAML per design, all into a single input directory: binder sequence
    # goes to YAML target slot (Boltz chain A), target sequence to YAML binder slot (Boltz chain B)
    # Give every staged input its design id and real extension back. Tools downstream
    # dispatch on the extension, and rmsd4all needs a basename unique per design so its
    # `structure1` column does not collapse to one key across the batch.
    mkdir -p boltz_inputs staged
    while IFS=\$'\\t' read -r design_id staged_path pdb_ext; do
        [[ -n "\$design_id" ]] || continue
        ln -s "\$(readlink -f "\$staged_path")" "staged/\${design_id}.\${pdb_ext}"
    done < .batch_manifest.tsv

    while IFS=\$'\\t' read -r design_id staged_path pdb_ext; do
        [[ -n "\$design_id" ]] || continue
        pdb_path="staged/\${design_id}.\${pdb_ext}"
        /usr/bin/python3 ${projectDir}/bin/create_boltz_yaml.py \\
            --target_id "\${design_id}_binder" \\
            --binder_id "\${design_id}_target" \\
            --target_from_pdb "\$pdb_path" \\
            --target_chains '${binder_chain}' \\
            ${target_source_args} \\
            ${binder_msa_flag} \\
            ${target_msa_flag} \\
            --output_yaml "boltz_inputs/\${design_id}.yml" \\
            ${use_msa_server_flag} \\
            ${templates_flag}
    done < .batch_manifest.tsv

    # Run Boltz prediction once for the whole batch; fail if log shows GPU OOM batch skip
    # (boltz may still exit 0). Passing the directory rather than a single YAML is what
    # amortises model load and process startup over every design in the batch.
    BOLTZ_PREDICT_LOG=.boltz_predict_console.log
    rm -f "\$BOLTZ_PREDICT_LOG"
    set +e
    boltz predict \\
        ${args} \\
        --preprocessing-threads ${task.cpus} \\
        --num_workers ${task.cpus} \\
        --output_format pdb \\
        ${use_msa_server_flag} \\
        --out_dir boltz_batch \\
        boltz_inputs 2>&1 | tee "\$BOLTZ_PREDICT_LOG"
    boltz_rc=\${PIPESTATUS[0]}
    set -e
    if grep -qF 'ran out of memory, skipping batch' "\$BOLTZ_PREDICT_LOG"; then
        echo 'BOLTZ_COMPARE_COMPLEX: Boltz logged GPU OOM (batch skipped); failing.' >&2
        exit 1
    fi
    if [[ "\$boltz_rc" -ne 0 ]]; then
        exit "\$boltz_rc"
    fi

    # Boltz writes a batch to boltz_batch/boltz_results_boltz_inputs/predictions/<design_id>/.
    # Reshape that into the per-design boltz_results_<design_id>/predictions/<design_id>/ layout
    # this module has always emitted, so publishDir patterns and every downstream path are
    # unchanged whatever the batch size is.
    mkdir -p per_design
    while IFS=\$'\\t' read -r design_id staged_path pdb_ext; do
        [[ -n "\$design_id" ]] || continue
        pdb_path="staged/\${design_id}.\${pdb_ext}"

        src="boltz_batch/boltz_results_boltz_inputs/predictions/\${design_id}"
        if [[ ! -d "\$src" ]]; then
            echo "BOLTZ_COMPARE_COMPLEX: Boltz produced no prediction for '\${design_id}'." >&2
            exit 1
        fi
        mkdir -p "boltz_results_\${design_id}/predictions"
        mv "\$src" "boltz_results_\${design_id}/predictions/\${design_id}"

        # Boltz also writes msa/, processed/ and lightning_logs/ alongside predictions/, and
        # this module has always published them inside boltz_results_<design_id>. They describe
        # the whole batch rather than one design, so they are hard-linked into each design's
        # results directory: no duplicated bytes, and at the default batch size of 1 the
        # published tree is exactly what it was before batching.
        for aux in msa processed lightning_logs; do
            aux_src="boltz_batch/boltz_results_boltz_inputs/\${aux}"
            [[ -d "\$aux_src" ]] && cp -al "\$aux_src" "boltz_results_\${design_id}/\${aux}"
        done

        pred_dir="boltz_results_\${design_id}/predictions/\${design_id}"
        model_pdb="\${pred_dir}/\${design_id}_model_0.pdb"

        # Run RMSD calculations against a directory pair holding just this design.
        rm -rf fixed mobile
        mkdir -p fixed mobile

        # Stage the fixed structure under a design-id-aware basename so the rmsd4all `structure1`
        # column is unique per-design (RFD3 stages every per-design RF3 output as `model.cif`,
        # which collapses all rows to the same key and breaks downstream merges by structure1).
        ln -s "\$(readlink -f "\$pdb_path")" "fixed/\${design_id}.\${pdb_ext}"
        ln -s "\$(readlink -f "\$model_pdb")" "mobile/\$(basename "\$model_pdb")"

        # Boltz complex PDB chain IDs are always A,B from create_boltz_yaml.py: YAML "target" slot is
        # chain A (here: input binder sequence), YAML "binder" slot is chain B (input target sequence).
        # Input design structure keeps original chain IDs (${target_chain}=target, ${binder_chain}=binder).
        # Mobile must use A=binder and B=target, not the input chain letters.
        # Target-aligned binder RMSD: superimpose on target, score binder
        /usr/bin/python3 ${projectDir}/bin/rmsd4all.py \\
            --tm-score \\
            --superimpose-chains ${target_chain} \\
            --mobile-superimpose-chains B \\
            --score-chains ${binder_chain} \\
            --mobile-score-chains A \\
            ${output_transformed_flag_target} \\
            fixed/ mobile/ > "rmsd_target_aligned_binder_\${design_id}.tsv"

        # Complex RMSD: Boltz mobile is A=binder, B=target. Fixed side keeps input chain IDs; when those
        # match Boltz (e.g. rfd AF2ig A=binder B=target), order aligns. Otherwise sequence-based alignment
        # still pairs homologous chains; target-aligned row above is the robust pose metric.
        /usr/bin/python3 ${projectDir}/bin/rmsd4all.py \\
            --tm-score \\
            --superimpose-chains ${binder_chain},${target_chain} \\
            --mobile-superimpose-chains A,B \\
            --score-chains ${binder_chain},${target_chain} \\
            --mobile-score-chains A,B \\
            ${output_transformed_flag_complex} \\
            fixed/ mobile/ > "rmsd_complex_\${design_id}.tsv"

        /usr/bin/python3 ${projectDir}/bin/ipsae.py \\
            --update-summary "\${pred_dir}/confidence_\${design_id}_model_0.json" \\
            --binder-chain A \\
            --target-chain B \\
            --format boltz \\
            "\${pred_dir}"/pae*.npz \\
            "\${pred_dir}"/*.pdb \\
            10 10

        # Parse confidence JSON
        /usr/bin/python3 ${projectDir}/bin/parse_boltz_confidence.py \\
            --json "\${pred_dir}/confidence_\${design_id}_model_0.json" \\
            --id "\${design_id}" \\
            --target "\${design_id}_target" \\
            --binder "\${design_id}_binder" > "confidence_\${design_id}.tsv"

        # Per-design bundle. The workflow re-associates these with their meta by design id,
        # which is what keeps batched outputs correctly paired (see BOLTZ_REFOLD_CORE).
        # Hard links rather than copies: same filesystem, no duplicated bytes, and unlike
        # symlinks they stay valid wherever Nextflow stages the bundle.
        d="per_design/\${design_id}"
        mkdir -p "\$d"
        ln "\$model_pdb" "\$d/\${design_id}_complex.pdb"
        ln "rmsd_target_aligned_binder_\${design_id}.tsv" "\$d/rmsd_target_aligned_binder.tsv"
        ln "rmsd_complex_\${design_id}.tsv" "\$d/rmsd_complex.tsv"
        ln "confidence_\${design_id}.tsv" "\$d/confidence.tsv"
        for f in "\${pred_dir}"/*_ipsae.tsv; do
            [[ -e "\$f" ]] && ln "\$f" "\$d/ipsae.tsv"
        done
        for f in "\${pred_dir}"/*_ipsae_byres.tsv; do
            [[ -e "\$f" ]] && ln "\$f" "\$d/ipsae_byres.tsv"
        done
    done < .batch_manifest.tsv

    rm -rf fixed mobile boltz_batch staged
    """
}
