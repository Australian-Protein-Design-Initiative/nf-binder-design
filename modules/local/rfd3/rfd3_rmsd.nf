// Refolded RF3 CIF chains (default A = target, B = binder with our RF3 JSON layout)
process RFD3_RMSD {
    tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(design_cif), path(refolded_cif), val(target_chain), val(binder_chain)

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
        --superimpose-chains ${target_chain} \\
        --mobile-superimpose-chains ${target_chain} \\
        --score-chains ${binder_chain} \\
        --mobile-score-chains ${binder_chain} \\
        fixed/ mobile/ > rmsd_target_aligned_binder_${meta.id}.tsv

    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${target_chain},${binder_chain} \\
        --mobile-superimpose-chains ${target_chain},${binder_chain} \\
        --score-chains ${target_chain},${binder_chain} \\
        --mobile-score-chains ${target_chain},${binder_chain} \\
        fixed/ mobile/ > rmsd_complex_${meta.id}.tsv

    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${binder_chain} \\
        --mobile-superimpose-chains ${binder_chain} \\
        --score-chains ${binder_chain} \\
        --mobile-score-chains ${binder_chain} \\
        fixed/ mobile/ > rmsd_binder_aligned_binder_${meta.id}.tsv

    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${target_chain} \\
        --mobile-superimpose-chains ${target_chain} \\
        --score-chains ${target_chain} \\
        --mobile-score-chains ${target_chain} \\
        fixed/ mobile/ > rmsd_target_aligned_target_${meta.id}.tsv
    """
}
