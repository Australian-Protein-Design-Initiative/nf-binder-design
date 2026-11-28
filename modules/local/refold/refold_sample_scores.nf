// One row per refolded sample, as opposed to REFOLD_SCORES_WIDE's one row per design.
//
// The wide table reduces each engine's samples to a single value per metric, which is
// what makes it rankable and also what discards the spread. This one keeps every
// sample's own confidence, RMSD and interface scores, so any summary statistic over
// the samples can be computed from it afterwards.
process REFOLD_SAMPLE_SCORES {
    tag "refold_sample_scores"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.8'

    publishDir path: "${params.outdir}/${params.refold_publish_dir ?: 'refold'}", mode: 'copy'

    input:
    path scores_tsv
    path rmsd_tsvs, stageAs: 'rmsd/*'
    path interface_tsvs, stageAs: 'interface/*'

    output:
    path 'refold_sample_scores.tsv', emit: scores

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

    python ${projectDir}/bin/refold/sample_scores.py '${scores_tsv}' \\
        ${long_args} \\
        ${interface_args} \\
        -o refold_sample_scores.tsv
    """
}
