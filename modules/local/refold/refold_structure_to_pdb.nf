// mmCIF -> PDB for the scoring tools that read PDB only (BindCraft's interface
// metrics go through PyRosetta's pose_from_pdb). Most engines now write mmCIF, so
// this sits between the refold and those scores. The output is named for the design,
// the engine, the fold batch and the sample, which is what the score tables key on.
// Both the batch and the sample are needed: an engine handed a list of seeds returns
// several samples from one batch, while an engine run once per seed writes every
// batch's sample under the same label.
process REFOLD_STRUCTURE_TO_PDB {
    tag "${meta.id} ${meta.fold_tool} b${meta.fold_batch} ${meta.fold_sample}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.7'

    input:
    tuple val(meta), path(structure, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}__${meta.fold_tool}__b${meta.fold_batch}__${meta.fold_sample}.pdb"), emit: pdb

    script:
    """
    set -euo pipefail

    python ${projectDir}/bin/structure_to_pdb.py '${structure}' \\
        -o '${meta.id}__${meta.fold_tool}__b${meta.fold_batch}__${meta.fold_sample}.pdb'
    """
}
