// RF3/RFD3 output is fixed: chain A = target, chain B = binder
process RFD3_RMSD {
    tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    input:
    tuple val(meta), path(design_cif), path(refolded_cif)

    output:
    tuple val(meta), path("rmsd_target_aligned_binder_${meta.id}.tsv"), emit: rmsd_target_aligned_binder
    tuple val(meta), path("rmsd_complex_${meta.id}.tsv"), emit: rmsd_complex
    tuple val(meta), path("rmsd_binder_aligned_binder_${meta.id}.tsv"), emit: rmsd_binder_aligned_binder
    tuple val(meta), path("rmsd_target_aligned_target_${meta.id}.tsv"), emit: rmsd_target_aligned_target

    script:
    """
    set -euo pipefail

    mkdir -p fixed/ mobile/
    ln -s "\$(readlink -f ${design_cif})" "fixed/\$(basename ${design_cif})"
    ln -s "\$(readlink -f ${refolded_cif})" "mobile/\$(basename ${refolded_cif})"

    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains A \\
        --mobile-superimpose-chains A \\
        --score-chains B \\
        --mobile-score-chains B \\
        fixed/ mobile/ > rmsd_target_aligned_binder_${meta.id}.tsv

    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains A,B \\
        --mobile-superimpose-chains A,B \\
        --score-chains A,B \\
        --mobile-score-chains A,B \\
        fixed/ mobile/ > rmsd_complex_${meta.id}.tsv

    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains B \\
        --mobile-superimpose-chains B \\
        --score-chains B \\
        --mobile-score-chains B \\
        fixed/ mobile/ > rmsd_binder_aligned_binder_${meta.id}.tsv

    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains A \\
        --mobile-superimpose-chains A \\
        --score-chains A \\
        --mobile-score-chains A \\
        fixed/ mobile/ > rmsd_target_aligned_target_${meta.id}.tsv
    """
}
