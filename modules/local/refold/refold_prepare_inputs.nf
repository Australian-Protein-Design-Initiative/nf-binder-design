// Extract the target and per-design binder sequences from a set of design complexes,
// as the targets.fasta / binders.fasta pair FOLD_PULLDOWN_CORE takes. One task for
// the whole set, so the target is read once and every design is checked against it.
process REFOLD_PREPARE_INPUTS {
    tag "refold_inputs"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.7'

    publishDir path: "${params.outdir}/${params.refold_publish_dir ?: 'refold'}/input", mode: 'copy'

    input:
    path designs, stageAs: 'designs/*'
    val binder_chain
    val target_chain
    val target_id

    output:
    path 'targets.fasta', emit: targets
    path 'binders.fasta', emit: binders

    script:
    """
    set -euo pipefail

    python ${projectDir}/bin/refold/prepare_inputs.py designs/ \\
        --binder-chain '${binder_chain}' \\
        --target-chain '${target_chain}' \\
        --target-id '${target_id}' \\
        --targets-out targets.fasta \\
        --binders-out binders.fasta
    """
}
