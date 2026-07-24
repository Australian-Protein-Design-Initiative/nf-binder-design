// Merge the per-tool fold score TSVs (af2/rf3/protenix already canonical; boltz
// native) into the master <outdir>/fold/fold_scores.tsv. bin/merge_fold_scores.py
// auto-detects each input's schema, normalizes to the canonical columns
// (equivalent scores share a name, plddt on 0-1, asymmetric per-chain-pair
// scores dropped), and concatenates one row per generated structure. CPU-only,
// local executor.
process FOLD_MERGE_SCORES {
    tag "fold_scores"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    publishDir "${params.outdir}/fold", mode: 'copy'

    input:
    path(tsvs)

    output:
    path('fold_scores.tsv'), emit: tsv

    script:
    def inputs = tsvs.collect { "--input '${it}'" }.join(' ')
    """
    python3 ${projectDir}/bin/merge_fold_scores.py ${inputs} -o fold_scores.tsv
    """
}
