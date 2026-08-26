process FOLD_PULLDOWN_MERGE_SCORES {
    tag "fold_pulldown_scores"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    publishDir "${params.outdir}/${params.fold_publish_dir ?: 'fold_pulldown'}", mode: 'copy'

    input:
    path scores_tsv
    path pairs_tsv

    output:
    path('fold_pulldown_scores.tsv'), emit: scores
    path('fold_pulldown_summary.tsv'), emit: summary

    script:
    """
    python3 ${projectDir}/bin/fold_pulldown_summarise.py \
        --scores '${scores_tsv}' \
        --pairs '${pairs_tsv}' \
        --scores-out fold_pulldown_scores.tsv \
        --summary-out fold_pulldown_summary.tsv
    """
}
