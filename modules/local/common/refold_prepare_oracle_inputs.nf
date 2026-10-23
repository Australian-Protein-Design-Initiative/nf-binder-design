// Extract the target and per-design binder sequences from a set of design complexes,
// as the targets.fasta / binders.fasta pair FOLD_PULLDOWN_CORE takes.
process REFOLD_PREPARE_ORACLE_INPUTS {
    tag "refold_oracle_inputs"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    publishDir path: "${params.outdir}/${params.fold_publish_dir ?: 'refold_oracles'}/input", mode: 'copy'

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

    python ${projectDir}/bin/refold/prepare_oracle_inputs.py designs/ \\
        --binder-chain '${binder_chain}' \\
        --target-chain '${target_chain}' \\
        --target-id '${target_id}' \\
        --targets-out targets.fasta \\
        --binders-out binders.fasta
    """
}
