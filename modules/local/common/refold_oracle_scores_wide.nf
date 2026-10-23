// Pivot the long per-(target, binder, tool) pulldown summary into one row per design,
// so a design workflow can merge it into combined_scores.tsv on the design id.
process REFOLD_ORACLE_SCORES_WIDE {
    tag "refold_oracle_scores"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    publishDir path: "${params.outdir}/${params.fold_publish_dir ?: 'refold_oracles'}", mode: 'copy'

    input:
    path summary_tsv
    val metrics

    output:
    path 'refold_oracle_scores.tsv', emit: scores

    script:
    """
    set -euo pipefail

    python ${projectDir}/bin/refold/oracle_summary_to_wide.py '${summary_tsv}' \\
        --metrics '${metrics}' \\
        -o refold_oracle_scores.tsv
    """
}
