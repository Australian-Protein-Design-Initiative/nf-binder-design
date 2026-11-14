// Pivot the refold's long tables into one row per design, so a design workflow can
// merge them into combined_scores.tsv on the design id. Confidence comes in per
// (target, binder, engine), RMSD per (design, engine, model) and interface scores per
// refolded complex; all three leave here as <engine>-prefixed columns.
//
// The per-sample table comes in as well, to choose which sample the structural columns
// describe: the best-ranked one, so the RMSD and interface values in a row all belong
// to the same predicted coordinates.
process REFOLD_SCORES_WIDE {
    tag "refold_scores"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.7'

    publishDir path: "${params.outdir}/${params.refold_publish_dir ?: 'refold'}", mode: 'copy'

    input:
    path summary_tsv
    path rmsd_tsvs, stageAs: 'rmsd/*'
    path interface_tsvs, stageAs: 'interface/*'
    path sample_scores_tsv
    val metrics

    output:
    path 'refold_scores.tsv', emit: scores

    script:
    def long_args = (rmsd_tsvs instanceof List ? rmsd_tsvs : [rmsd_tsvs])
        .findAll { it.name != 'empty' }
        .collect { "--long '${it}'" }
        .join(' ')
    def interface_args = (interface_tsvs instanceof List ? interface_tsvs : [interface_tsvs])
        .findAll { it.name != 'empty' }
        .collect { "--bindcraft '${it}'" }
        .join(' ')
    """
    set -euo pipefail

    python ${projectDir}/bin/refold/summary_to_wide.py '${summary_tsv}' \\
        --metrics '${metrics}' \\
        ${long_args} \\
        ${interface_args} \\
        --samples '${sample_scores_tsv}' \\
        --rank-metric '${params.consensus_metric}' \\
        --structure-stat '${params.refold_structure_stat}' \\
        -o refold_scores.tsv
    """
}
