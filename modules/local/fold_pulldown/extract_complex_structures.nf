// Read binder sequences (and optionally the target) out of supplied target-binder
// complexes, for --method fold_pulldown --complex_structures.
//
// One task for the whole set: the target chain resolution compares against every
// --targets record, and deriving the targets from the complexes needs them all in
// hand to deduplicate by sequence.
process EXTRACT_COMPLEX_STRUCTURES {
    tag "complex_structures"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.8'

    publishDir path: "${params.outdir}/${params.fold_publish_dir ?: 'fold_pulldown'}/input", mode: 'copy'

    input:
    path complexes, stageAs: 'complexes/*'
    val binder_chains
    val target_chains
    // Empty placeholder means "no --targets given": the targets are then read back out
    // of the complexes instead, and --target_chains auto has nothing to match against.
    path targets

    output:
    path 'binders.fasta', emit: binders
    path 'chains.tsv', emit: chains
    path 'targets.fasta', emit: targets, optional: true

    script:
    def target_chains_arg = target_chains ? "--target-chains '${target_chains}'" : ''
    def has_targets = targets.name != 'no_pulldown_targets.fasta'
    def targets_arg = has_targets ? "--targets '${targets}'" : ''
    def targets_out = has_targets ? '' : '--targets-out targets.fasta'
    """
    set -euo pipefail

    python ${projectDir}/bin/pulldown/complex_structures.py complexes/ \\
        --binder-chains '${binder_chains}' \\
        ${target_chains_arg} \\
        ${targets_arg} \\
        ${targets_out} \\
        --binders-out binders.fasta \\
        --chains-out chains.tsv \\
        --verbose
    """
}
