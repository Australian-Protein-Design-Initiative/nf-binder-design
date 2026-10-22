process COMBINE_RFD4_SCORES {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    publishDir path: "${params.outdir}/rfd4_proteina", pattern: 'combined_scores.tsv', mode: 'copy'
    publishDir path: "${params.outdir}/rfd4_proteina", pattern: 'binders.fasta', mode: 'copy'

    input:
    tuple path(design_table, stageAs: 'rfd4_designs.tsv'),
          path(boltz_scores_complex, stageAs: 'boltz_scores_complex.tsv'),
          path(boltz_scores_monomer, stageAs: 'boltz_scores_monomer.tsv'),
          path(boltz_rmsd_target_aligned_binder, stageAs: 'boltz_rmsd_target_aligned_binder.tsv'),
          path(boltz_rmsd_monomer_vs_complex, stageAs: 'boltz_rmsd_monomer_vs_complex.tsv'),
          val(binder_chain),
          path(design_cifs, stageAs: 'cifs/*')

    output:
    path 'combined_scores.tsv', emit: combined_scores
    path 'binders.fasta', emit: binders_fasta

    script:
    """
    set -euo pipefail

    cp ${design_table} combined_scores.tsv

    # RFD4-Proteina co-designs the sequence, so it is read straight off the generated complex.
    python ${projectDir}/bin/pdb_to_fasta.py cifs/ --tsv --chains "${binder_chain}" -o sequences.tsv

    seq_rows=\$(awk 'NR>1 {n++} END {print n+0}' sequences.tsv 2>/dev/null || echo 0)
    if [[ -s sequences.tsv && "\${seq_rows}" -ge 1 ]]; then
      python ${projectDir}/bin/merge_scores.py \\
        combined_scores.tsv sequences.tsv \\
        --keys filename \\
        --strip-suffix '(\\.cif\\.zst|\\.cif\\.gz|\\.cif)\$' \\
        --sort-by id \\
        --first-column id,filename,sequence,length \\
        -o with_sequence.tsv
      mv with_sequence.tsv combined_scores.tsv
    fi

    if [[ -s ${boltz_scores_complex} ]]; then
      python ${projectDir}/bin/merge_scores.py \\
        combined_scores.tsv ${boltz_scores_complex} \\
        --keys id \\
        --column-prefix boltz_complex_ \\
        --sort-by boltz_complex_pair_pae_min \\
        --first-column id,filename,sequence,length,boltz_complex_pair_pae_min \\
        --drop-columns 'boltz_complex_target,boltz_complex_binder,boltz_complex_state,boltz_complex_ptm,boltz_complex_ligand_iptm,boltz_complex_protein_iptm,boltz_complex_has_clash' \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s ${boltz_scores_monomer} ]]; then
      csvtk replace -t -f id -p "_monomer\$" -r "" ${boltz_scores_monomer} >tmp_boltz_monomer.tsv
      python ${projectDir}/bin/merge_scores.py \\
        combined_scores.tsv tmp_boltz_monomer.tsv \\
        --keys id \\
        --column-prefix boltz_monomer_ \\
        --sort-by id \\
        --first-column id,filename,sequence,length \\
        --drop-columns 'boltz_monomer_target,boltz_monomer_binder,boltz_monomer_state,boltz_monomer_ptm,boltz_monomer_ligand_iptm,boltz_monomer_protein_iptm,boltz_monomer_has_clash' \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s ${boltz_rmsd_target_aligned_binder} ]]; then
      python ${projectDir}/bin/rfd3/merge_rmsd_columns.py \\
        combined_scores.tsv \\
        "${boltz_rmsd_target_aligned_binder}=boltz_target_aligned_binder_" \\
        --keys id,structure1 \\
        --strip-suffix '(\\.pdb|\\.cif)\$' \\
        --value-columns rmsd_pruned,rmsd_all \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s ${boltz_rmsd_monomer_vs_complex} ]]; then
      python ${projectDir}/bin/rfd3/merge_rmsd_columns.py \\
        combined_scores.tsv \\
        "${boltz_rmsd_monomer_vs_complex}=boltz_monomer_vs_complex_" \\
        --keys id,structure1 \\
        --strip-suffix '(_model_0\\.pdb|_monomer)\$' \\
        --value-columns rmsd_pruned,rmsd_all \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
    fi

    if [[ -s combined_scores.tsv ]] && grep -q 'sequence' combined_scores.tsv; then
      python ${projectDir}/bin/rfd3/tsv_to_binders_fasta.py combined_scores.tsv \\
        --scores boltz_complex_pair_pae_min,boltz_complex_plddt \\
        -o binders.fasta
    else
      touch binders.fasta
    fi
    """
}
