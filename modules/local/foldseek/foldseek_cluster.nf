#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_CLUSTER process

All-vs-all structural clustering of a run's own designs, answering "how many
structurally distinct binders did this run produce?". This is a different job from
FOLDSEEK_SEARCH, which looks for database homologs -- no database is involved here.

`easy-cluster` clusters binder chains; `easy-multimercluster` clusters whole complexes
and so requires every chain to match, which makes it much stricter for a set of designs
that share one target.

foldseek release 10 has no `easy-interfacecluster`; add it as a third mode once a
release carrying it is available.

Input:
    path structures   - Binder chain PDBs, or whole complexes with --foldseek_cluster_multimer

Output:
    cluster_tsv       - foldseek's two-column representative/member table
    rep_seq           - representative sequences

FOLDSEEK_CLUSTER_TABLE turns cluster_tsv into the per-design table. It is a separate
process because this image has no python3.
*/

process FOLDSEEK_CLUSTER {
    tag "foldseek_cluster"
    container 'quay.io/biocontainers/foldseek:10.941cd33--h5021889_1'

    publishDir path: "${params.outdir}/foldseek_cluster", mode: 'copy',
        pattern: '{clusters.tsv,cluster_rep_seq.fasta,cluster_report.tsv}'

    input:
    path structures, stageAs: 'structures/*'

    output:
    path 'clusters.tsv', emit: cluster_tsv
    path 'cluster_rep_seq.fasta', emit: rep_seq, optional: true
    path 'cluster_report.tsv', emit: cluster_report, optional: true

    script:
    def multimer = params.foldseek_cluster_multimer as boolean
    def subcommand = multimer ? 'easy-multimercluster' : 'easy-cluster'
    def mode_args = multimer
        ? "--multimer-tm-threshold ${params.foldseek_cluster_multimer_tm} " +
          "--chain-tm-threshold ${params.foldseek_cluster_chain_tm} " +
          "--interface-lddt-threshold ${params.foldseek_cluster_interface_lddt}"
        : "--alignment-type ${params.foldseek_cluster_alignment_type} " +
          "--cov-mode ${params.foldseek_cluster_cov_mode} " +
          "-c ${params.foldseek_cluster_coverage} " +
          "--min-seq-id ${params.foldseek_cluster_min_seq_id} " +
          "--tmscore-threshold ${params.foldseek_cluster_tmscore}"
    def extra_args = params.foldseek_cluster_extra_args ?: ''
    """
    set -euo pipefail

    foldseek ${subcommand} \\
        structures/ \\
        cluster \\
        tmp \\
        ${mode_args} \\
        --threads ${task.cpus} \\
        ${extra_args} \\
        ${task.ext.args ?: ''}

    # foldseek names its outputs <prefix>_cluster.tsv / _rep_seq.fasta, and
    # easy-multimercluster adds <prefix>_cluster_report.
    mv cluster_cluster.tsv clusters.tsv
    # easy-multimercluster's report has no extension; the rep_seq fasta is already named right
    if [[ -f cluster_cluster_report ]]; then mv cluster_cluster_report cluster_report.tsv; fi
    """
}
