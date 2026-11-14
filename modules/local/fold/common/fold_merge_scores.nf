// Merge the per-tool fold score TSVs (af2/rf3/protenix already canonical; boltz
// native) into the master <outdir>/fold/fold_scores.tsv. bin/fold/merge_fold_scores.py
// auto-detects each input's schema, normalizes to the canonical columns
// (equivalent scores share a name, plddt on 0-1, asymmetric per-chain-pair
// scores dropped), and concatenates one row per generated structure. CPU-only,
// local executor.
process FOLD_MERGE_SCORES {
    tag "fold_scores"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.7'

    // publish_label distinguishes two fold runs that share an output directory: the
    // refold path folds the complexes and then the binder monomers, and without a label
    // the second run's fold_scores.tsv would overwrite the first's.
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}",
        mode: 'copy',
        saveAs: { filename -> publish_label ? "${publish_label}_${filename}" : filename },
        enabled: (params.method != 'fold_pulldown' && params.fold_publish_dir != 'fold_pulldown')
    )

    input:
    path(tsvs)
    val publish_label

    output:
    path('fold_scores.tsv'), emit: tsv

    script:
    def inputs = tsvs.collect { "--input '${it}'" }.join(' ')
    """
    python3 ${projectDir}/bin/fold/merge_fold_scores.py ${inputs} -o fold_scores.tsv
    """
}
