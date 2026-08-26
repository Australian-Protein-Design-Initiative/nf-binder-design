process FOLD_PULLDOWN_REPORTING {
    publishDir "${params.outdir}/${params.fold_publish_dir ?: 'fold_pulldown'}", mode: 'copy'

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    path('fold_pulldown_reporting.qmd')
    path('fold_pulldown_scores.tsv')
    path('fold_pulldown_summary.tsv')

    output:
    path('fold_pulldown_report.html')

    script:
    """
    set -euo pipefail
    quarto render fold_pulldown_reporting.qmd --execute-dir \${PWD} --output - >fold_pulldown_report.html
    """
}
