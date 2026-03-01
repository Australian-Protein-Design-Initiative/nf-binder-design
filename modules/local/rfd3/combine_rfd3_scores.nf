process COMBINE_RFD3_SCORES {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    publishDir path: "${params.outdir}/rfd3", pattern: 'combined_scores.tsv', mode: 'copy'

    input:
    path rf3_scores, stageAs: 'rf3_scores/*'
    path rfd3_scores, stageAs: 'rfd3_scores/*'
    tuple path(rmsd_target_aligned_binder_tsv), path(rmsd_complex_tsv), path(rmsd_binder_aligned_binder_tsv), path(rmsd_target_aligned_target_tsv)

    output:
    path 'combined_scores.tsv', emit: combined_scores

    script:
    """
    set -euo pipefail

    awk 'NR==1 || FNR>1' rf3_scores/*.tsv > rf3_all.tsv
    awk 'NR==1 || FNR>1' rfd3_scores/*.tsv > rfd3_all.tsv

    python ${projectDir}/bin/merge_scores.py \\
      rf3_all.tsv rfd3_all.tsv \\
      --keys backbone_id \\
      --sort-by pair_pae_min \\
      --first-column id,filename,pair_pae_min,ranking_score,iptm,plddt \\
      -o combined_scores.tsv

    if [[ -s ${rmsd_target_aligned_binder_tsv} && -s ${rmsd_complex_tsv} && -s ${rmsd_binder_aligned_binder_tsv} && -s ${rmsd_target_aligned_target_tsv} ]]; then
      csvtk -t cut -b -f structure1,rmsd_all ${rmsd_target_aligned_binder_tsv} > tmp_rmsd_target_aligned_binder.tsv
      python ${projectDir}/bin/merge_scores.py \\
        combined_scores.tsv tmp_rmsd_target_aligned_binder.tsv \\
        --keys filename,structure1 \\
        --column-prefix refold_rmsd_target_aligned_binder_ \\
        --drop-columns 'refold_rmsd_target_aligned_binder_structure.*' \\
        -o step1.tsv

      csvtk -t cut -b -f structure1,rmsd_all ${rmsd_complex_tsv} > tmp_rmsd_complex.tsv
      python ${projectDir}/bin/merge_scores.py \\
        step1.tsv tmp_rmsd_complex.tsv \\
        --keys filename,structure1 \\
        --column-prefix refold_rmsd_complex_ \\
        --drop-columns 'refold_rmsd_complex_structure.*' \\
        -o step2.tsv

      csvtk -t cut -b -f structure1,rmsd_all ${rmsd_binder_aligned_binder_tsv} > tmp_rmsd_binder_aligned_binder.tsv
      python ${projectDir}/bin/merge_scores.py \\
        step2.tsv tmp_rmsd_binder_aligned_binder.tsv \\
        --keys filename,structure1 \\
        --column-prefix refold_rmsd_binder_aligned_binder_ \\
        --drop-columns 'refold_rmsd_binder_aligned_binder_structure.*' \\
        -o step3.tsv

      csvtk -t cut -b -f structure1,rmsd_all ${rmsd_target_aligned_target_tsv} > tmp_rmsd_target_aligned_target.tsv
      python ${projectDir}/bin/merge_scores.py \\
        step3.tsv tmp_rmsd_target_aligned_target.tsv \\
        --keys filename,structure1 \\
        --column-prefix refold_rmsd_target_aligned_target_ \\
        --drop-columns 'refold_rmsd_target_aligned_target_structure.*' \\
        -o combined_scores.tsv
    fi
    """
}
