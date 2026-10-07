// Match --templates structure files to the chains being folded (one task per
// run). Engine input generators look templates up in fold_templates/index.json
// by chain sequence; see bin/fold/match_templates.py.
process MATCH_FOLD_TEMPLATES {
    tag "${templates instanceof List ? templates.size() : 1} template file(s)"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}",
        mode: 'copy',
        saveAs: { filename -> filename == 'fold_templates' ? 'templates' : null },
    )

    input:
    path templates, stageAs: 'templates_in/*'
    path queries, stageAs: 'queries/query*.fasta'

    output:
    path 'fold_templates', type: 'dir', emit: templates

    script:
    """
    python ${projectDir}/bin/fold/match_templates.py \\
        --templates templates_in/* \\
        --queries queries/* \\
        --min-identity ${params.template_min_identity} \\
        --min-coverage ${params.template_min_coverage} \\
        --max-per-chain ${params.template_max_per_chain} \\
        -o fold_templates
    """
}
