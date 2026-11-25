#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_CLUSTER_DESIGNS subworkflow

Structural clustering of a run's own designs, to measure how many distinct binders it
produced. Unrelated to FOLDSEEK_SEARCH, which finds database homologs -- the two share
only the binder-chain extraction step, and have separate flag namespaces.

Inputs:
    ch_designs   Channel of design complex structures (one path per design)

Parameters (from params scope):
    params.foldseek_cluster_max        Designs to cluster [default: false, meaning all]
    params.foldseek_cluster_multimer   Cluster whole complexes rather than binder chains
    params.foldseek_cluster_*          Thresholds, see conf/refold_fold_stack.config

Outputs:
    cluster_tsv      foldseek's representative/member table
    cluster_table    One row per design, joinable onto combined_scores.tsv on `id`
    rep_seq          Representative sequences
*/

// Aliased because the calling workflow may already have invoked the unaliased process for
// the database search, and a process instance can only be invoked once per run.
include { FOLDSEEK_PREPARE_QUERIES as FOLDSEEK_CLUSTER_PREPARE_QUERIES } from '../../modules/local/foldseek/foldseek_prepare_queries'
include { FOLDSEEK_CLUSTER } from '../../modules/local/foldseek/foldseek_cluster'
include { FOLDSEEK_CLUSTER_TABLE } from '../../modules/local/foldseek/foldseek_cluster_table'

workflow FOLDSEEK_CLUSTER_DESIGNS {

    take:
    ch_designs   // channel of design structure files

    main:

    def cluster_max = params.foldseek_cluster_max
    def ch_capped = cluster_max ? ch_designs.take(cluster_max as int) : ch_designs

    if (params.foldseek_cluster_multimer) {
        // easy-multimercluster needs every chain, so the complexes go in as they are
        FOLDSEEK_CLUSTER(ch_capped.collect())
    }
    else {
        FOLDSEEK_CLUSTER_PREPARE_QUERIES(ch_capped.collect(), 0)
        FOLDSEEK_CLUSTER(FOLDSEEK_CLUSTER_PREPARE_QUERIES.out.design_chains.collect())
    }

    FOLDSEEK_CLUSTER_TABLE(FOLDSEEK_CLUSTER.out.cluster_tsv)

    emit:
    cluster_tsv = FOLDSEEK_CLUSTER.out.cluster_tsv
    cluster_table = FOLDSEEK_CLUSTER_TABLE.out.cluster_table
    rep_seq = FOLDSEEK_CLUSTER.out.rep_seq
}
