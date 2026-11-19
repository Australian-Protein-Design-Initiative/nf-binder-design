// Compare a refolded complex with the design it was folded from.
//
// Two comparisons, the pair BOLTZ_COMPARE_COMPLEX used to compute inline:
//   - superimpose on the target, score the binder. This is the pose metric: did the
//     binder come back in the same place against the target.
//   - superimpose and score the whole complex.
//
// Engine-agnostic: it reads structures, so it does not care which predictor wrote
// them, and it needs no GPU, so it runs in the utils container rather than in an
// engine image. The refold path pairs the design with every model the engine
// produced for it, so a run with --n_predictions > 1 gets a row per model. The refold
// hands over one sample per task, and the rows are keyed by (design, engine, batch,
// model) to match the confidence rows the fold stack writes for the same samples.
process REFOLD_COMPARE_COMPLEX {
    tag "${meta.id} ${meta.fold_tool} b${meta.fold_batch} ${meta.fold_sample}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.8'

    publishDir(
        path: "${params.outdir}/${params.refold_publish_dir ?: 'refold'}/rmsd/aligned",
        pattern: 'aligned_*/*.pdb',
        mode: 'copy',
        saveAs: { file(it).name }
    )

    input:
    // The design is staged under its own name (unique per design, and rmsd4all's
    // structure1 column keys on it); the models keep theirs, which carry the pair id
    // and sample number.
    tuple val(meta), path(design, stageAs: 'design/*'), path(models, stageAs: 'models/*')
    val binder_chain        // chain id in the design
    val target_chain        // chain id in the design
    val pred_binder_chain   // chain id in the refolded complex
    val pred_target_chain   // chain id in the refolded complex

    output:
    tuple val(meta), path('refold_rmsd.tsv'), emit: rmsd
    path 'aligned_*/*.pdb', emit: aligned, optional: true

    script:
    def aligned_target = params.output_rmsd_aligned ? '--output-transformed aligned_target_aligned_binder/' : ''
    def aligned_complex = params.output_rmsd_aligned ? '--output-transformed aligned_complex/' : ''
    // Only the target-aligned comparison: DockQ is defined on a receptor and a ligand,
    // and the whole-complex comparison puts both chains on both sides.
    def dockq = params.refold_dockq ? '--dockq' : ''
    """
    set -euo pipefail

    # Target-aligned binder RMSD: superimpose on the target, score the binder.
    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${target_chain} \\
        --mobile-superimpose-chains ${pred_target_chain} \\
        --score-chains ${binder_chain} \\
        --mobile-score-chains ${pred_binder_chain} \\
        ${dockq} \\
        ${aligned_target} \\
        design/ models/ > rmsd_target_aligned_binder.tsv

    # Whole-complex RMSD. Chain order is given explicitly on both sides rather than
    # left to sequence-based matching, which can pair the binder with the target when
    # the two are similar in length.
    python ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${binder_chain},${target_chain} \\
        --mobile-superimpose-chains ${pred_binder_chain},${pred_target_chain} \\
        --score-chains ${binder_chain},${target_chain} \\
        --mobile-score-chains ${pred_binder_chain},${pred_target_chain} \\
        ${aligned_complex} \\
        design/ models/ > rmsd_complex.tsv

    python ${projectDir}/bin/refold/rmsd_rows.py \\
        --id '${meta.id}' \\
        --tool '${meta.fold_tool}' \\
        --batch '${meta.fold_batch}' \\
        --model '${meta.fold_sample}' \\
        --table 'target_aligned_binder=rmsd_target_aligned_binder.tsv' \\
        --table 'complex=rmsd_complex.tsv' \\
        -o refold_rmsd.tsv
    """
}
