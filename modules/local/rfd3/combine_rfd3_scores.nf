process COMBINE_RFD3_SCORES {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    publishDir path: "${params.outdir}/rfd3", pattern: 'combined_scores.tsv', mode: 'copy'

    input:
    path rf3_scores, stageAs: 'rf3_scores/*'
    path rfd3_scores, stageAs: 'rfd3_scores/*'

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
    """
}
