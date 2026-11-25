#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_CLUSTER_TABLE process

Rewrites foldseek's representative/member cluster TSV as one row per design, keyed on
`id`, so it can be merged into a combined_scores.tsv.

Separate from FOLDSEEK_CLUSTER because the foldseek image carries no python3.

Input:
    path cluster_tsv  - foldseek's *_cluster.tsv

Output:
    cluster_table     - one row per design: rep, cluster id, size, is_rep
*/

process FOLDSEEK_CLUSTER_TABLE {
    tag "foldseek_cluster_table"
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.7'

    publishDir path: "${params.outdir}/foldseek_cluster", mode: 'copy', pattern: 'cluster_table.tsv'

    input:
    path cluster_tsv

    output:
    path 'cluster_table.tsv', emit: cluster_table

    script:
    """
    python3 ${projectDir}/bin/foldseek/cluster_table.py ${cluster_tsv} -o cluster_table.tsv
    """
}
