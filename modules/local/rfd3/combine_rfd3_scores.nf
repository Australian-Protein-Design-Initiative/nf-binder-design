process COMBINE_RFD3_SCORES {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.8'

    publishDir path: "${params.outdir}/rfd3", pattern: 'combined_scores.tsv', mode: 'copy'
    publishDir path: "${params.outdir}/rfd3", pattern: 'binders.fasta', mode: 'copy'

    input:
    tuple path(rf3_scores),
          path(rfd3_scores),
          path(rmsd_target_aligned_binder_tsv, stageAs: 'rmsd_target_aligned_binder.tsv'),
          path(rmsd_complex_tsv, stageAs: 'rmsd_complex.tsv'),
          path(rmsd_binder_aligned_binder_tsv, stageAs: 'rmsd_binder_aligned_binder.tsv'),
          path(rmsd_target_aligned_target_tsv, stageAs: 'rmsd_target_aligned_target.tsv'),
          path(refold_scores, stageAs: 'refold_scores.tsv'),
          val(binder_seq_chain),
          path(mpnn_cifs, stageAs: 'cifs/*')

    output:
    path 'combined_scores.tsv', emit: combined_scores
    path 'binders.fasta', emit: binders_fasta

    script:
    """
    set -euo pipefail

    python ${projectDir}/bin/merge_scores.py \\
      ${rf3_scores} ${rfd3_scores} \\
      --keys backbone_id \\
      --sort-by pair_pae_min \\
      --first-column id,filename,pair_pae_min,ranking_score,iptm,plddt \\
      -o combined_scores.tsv

    python ${projectDir}/bin/pdb_to_fasta.py cifs/ --tsv --chains "${binder_seq_chain}" -o sequences.tsv

    seq_rows=\$(awk 'NR>1 {n++} END {print n+0}' sequences.tsv 2>/dev/null || echo 0)
    if [[ -s sequences.tsv && "\${seq_rows}" -ge 1 ]]; then
      csvtk replace -t -f filename -p '\\.cif\$' -r '_rf3.cif' sequences.tsv -o sequences_for_merge.tsv
      mv sequences_for_merge.tsv sequences.tsv
      python ${projectDir}/bin/merge_scores.py \\
        combined_scores.tsv sequences.tsv \\
        --keys filename \\
        --first-column id,filename,sequence,length,chain,pair_pae_min,ranking_score,iptm,plddt \\
        -o with_sequence.tsv
      mv with_sequence.tsv combined_scores.tsv
    fi

    if [[ -s ${rmsd_target_aligned_binder_tsv} && -s ${rmsd_complex_tsv} && -s ${rmsd_binder_aligned_binder_tsv} && -s ${rmsd_target_aligned_target_tsv} ]]; then
      python ${projectDir}/bin/rfd3/merge_rmsd_columns.py \\
        combined_scores.tsv \\
        "${rmsd_target_aligned_binder_tsv}=rf3_rmsd_target_aligned_binder_" \\
        "${rmsd_complex_tsv}=rf3_rmsd_complex_" \\
        "${rmsd_binder_aligned_binder_tsv}=rf3_rmsd_binder_aligned_binder_" \\
        "${rmsd_target_aligned_target_tsv}=rf3_rmsd_target_aligned_target_" \\
        --keys filename,structure1 \\
        --strip-suffix '(\\.pdb|\\.cif)\$' \\
        --merge-key-replace-basename '_rf3\\.cif\$' '.cif' \\
        --value-columns rmsd_all \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    # Refold scores: already one row per design, with engine-prefixed columns for
    # confidence, RMSD against the design and interface metrics (REFOLD_SCORES_WIDE).
    # This replaces the four separate Boltz tables that BOLTZ_REFOLD_CORE produced, and
    # carries every engine in --refold_with rather than Boltz alone.
    #
    # Keyed on merge_id: this table's `filename` is <design>_rf3.cif, while the refold
    # names each row by the design id it started from.
    if [[ -s ${refold_scores} ]]; then
      csvtk mutate -t -f filename -n merge_id -p '^(.+)\$' combined_scores.tsv \\
        | csvtk replace -t -f merge_id -p '^.*/' -r '' \\
        | csvtk replace -t -f merge_id -p '_rf3\\.cif' -r '' \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv

      python ${projectDir}/bin/merge_scores.py \\
        combined_scores.tsv ${refold_scores} \\
        --keys merge_id,id \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv

      csvtk -t cut -f -merge_id combined_scores.tsv -o combined_scores.tmp.tsv \\
        && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s combined_scores.tsv ]] && grep -q 'sequence' combined_scores.tsv; then
      python ${projectDir}/bin/rfd3/tsv_to_binders_fasta.py combined_scores.tsv \\
        --scores '${params.refold_fasta_score_columns}' \\
        -o binders.fasta
    else
      touch binders.fasta
    fi
    """
}
