process COMBINE_RFD3_SCORES {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    publishDir path: "${params.outdir}/rfd3", pattern: 'combined_scores.tsv', mode: 'copy'
    publishDir path: "${params.outdir}/rfd3", pattern: 'binders.fasta', mode: 'copy'

    input:
    tuple path(rf3_scores),
          path(rfd3_scores),
          path(rmsd_target_aligned_binder_tsv, stageAs: 'rmsd_target_aligned_binder.tsv'),
          path(rmsd_complex_tsv, stageAs: 'rmsd_complex.tsv'),
          path(rmsd_binder_aligned_binder_tsv, stageAs: 'rmsd_binder_aligned_binder.tsv'),
          path(rmsd_target_aligned_target_tsv, stageAs: 'rmsd_target_aligned_target.tsv'),
          path(boltz_scores_complex, stageAs: 'boltz_scores_complex.tsv'),
          path(boltz_scores_monomer, stageAs: 'boltz_scores_monomer.tsv'),
          path(boltz_rmsd_target_aligned_binder, stageAs: 'boltz_rmsd_target_aligned_binder.tsv'),
          path(boltz_rmsd_monomer_vs_complex, stageAs: 'boltz_rmsd_monomer_vs_complex.tsv'),
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

    if [[ -s ${boltz_scores_complex} ]] || [[ -s ${boltz_scores_monomer} ]] \\
        || [[ -s ${boltz_rmsd_target_aligned_binder} ]] || [[ -s ${boltz_rmsd_monomer_vs_complex} ]]; then
      csvtk mutate -t -f filename -n merge_id -p '^(.+)\$' combined_scores.tsv \\
        | csvtk replace -t -f merge_id -p '^.*/' -r '' \\
        | csvtk replace -t -f merge_id -p '_rf3\\.cif' -r '' \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s ${boltz_scores_complex} ]]; then
      cp ${boltz_scores_complex} tmp_boltz_complex.tsv
      python ${projectDir}/bin/merge_scores.py \\
        combined_scores.tsv tmp_boltz_complex.tsv \\
        --keys merge_id,id \\
        --column-prefix boltz_complex_ \\
        --drop-columns 'boltz_complex_target,boltz_complex_binder,boltz_complex_state,boltz_complex_ptm,boltz_complex_iptm,boltz_complex_ligand_iptm,boltz_complex_protein_iptm,boltz_complex_has_clash' \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s ${boltz_scores_monomer} ]]; then
      csvtk replace -t -f id -p "_monomer\$" -r "" ${boltz_scores_monomer} > tmp_boltz_monomer.tsv
      python ${projectDir}/bin/merge_scores.py \\
        combined_scores.tsv tmp_boltz_monomer.tsv \\
        --keys merge_id,id \\
        --column-prefix boltz_monomer_ \\
        --drop-columns 'boltz_monomer_target,boltz_monomer_binder,boltz_monomer_state,boltz_monomer_ptm,boltz_monomer_ligand_iptm,boltz_monomer_protein_iptm,boltz_monomer_has_clash' \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s ${boltz_rmsd_target_aligned_binder} ]]; then
      python ${projectDir}/bin/rfd3/merge_rmsd_columns.py \\
        combined_scores.tsv \\
        "${boltz_rmsd_target_aligned_binder}=boltz_target_aligned_binder_" \\
        --keys merge_id,structure1 \\
        --strip-suffix '(\\.pdb|\\.cif)\$' \\
        --value-columns rmsd_pruned,rmsd_all \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s ${boltz_rmsd_monomer_vs_complex} ]]; then
      python ${projectDir}/bin/rfd3/merge_rmsd_columns.py \\
        combined_scores.tsv \\
        "${boltz_rmsd_monomer_vs_complex}=boltz_monomer_vs_complex_" \\
        --keys merge_id,structure1 \\
        --strip-suffix '(_model_0\\.pdb|_monomer)\$' \\
        --value-columns rmsd_pruned,rmsd_all \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s combined_scores.tsv ]] && grep -q 'merge_id' combined_scores.tsv; then
      csvtk -t cut -f -merge_id combined_scores.tsv -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s combined_scores.tsv ]] && grep -q 'sequence' combined_scores.tsv; then
      python ${projectDir}/bin/rfd3/tsv_to_binders_fasta.py combined_scores.tsv -o binders.fasta
    else
      touch binders.fasta
    fi
    """
}
