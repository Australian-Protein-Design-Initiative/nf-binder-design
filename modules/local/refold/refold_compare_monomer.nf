// Compare a binder folded on its own with the same binder in the design and in the
// refolded complex.
//
// The pair of questions BOLTZ_COMPARE_BINDER_MONOMER asked, with no engine assumed:
// does the binder hold its designed fold without the target (monomer vs design), and
// does binding change it (monomer vs refolded complex). A binder that only adopts its
// fold against the target is a different proposition from one that is stable alone.
//
// The monomer and the complex are paired by fold batch, so each task compares one
// monomer sample with the complex sample from the same seed. Pairing every monomer
// sample with every complex sample would give a grid rather than one number per
// sample, and across seeds the samples have no meaningful order to pair on.
process REFOLD_COMPARE_MONOMER {
    tag "${meta.id} ${meta.fold_tool} b${meta.fold_batch} ${meta.fold_sample}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.7'

    input:
    tuple val(meta),
          path(design, stageAs: 'design/*'),
          path(complex_model, stageAs: 'complex/*'),
          path(monomers, stageAs: 'monomers/*')
    val binder_chain        // chain id in the design
    val pred_binder_chain   // binder chain id in the refolded complex
    val monomer_chain       // chain id in the folded monomer

    output:
    tuple val(meta), path('refold_monomer_rmsd.tsv'), emit: rmsd

    script:
    """
    set -euo pipefail

    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${binder_chain} \\
        --mobile-superimpose-chains ${monomer_chain} \\
        --score-chains ${binder_chain} \\
        --mobile-score-chains ${monomer_chain} \\
        design/ monomers/ > rmsd_monomer_vs_design.tsv

    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${pred_binder_chain} \\
        --mobile-superimpose-chains ${monomer_chain} \\
        --score-chains ${pred_binder_chain} \\
        --mobile-score-chains ${monomer_chain} \\
        complex/ monomers/ > rmsd_monomer_vs_complex.tsv

    python ${projectDir}/bin/refold/rmsd_rows.py \\
        --id '${meta.id}' \\
        --tool '${meta.fold_tool}' \\
        --batch '${meta.fold_batch}' \\
        --model '${meta.fold_sample}' \\
        --table 'monomer_vs_design=rmsd_monomer_vs_design.tsv' \\
        --table 'monomer_vs_complex=rmsd_monomer_vs_complex.tsv' \\
        -o refold_monomer_rmsd.tsv
    """
}
