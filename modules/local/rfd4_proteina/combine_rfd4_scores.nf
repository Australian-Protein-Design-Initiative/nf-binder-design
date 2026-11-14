process COMBINE_RFD4_SCORES {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.7'

    publishDir path: "${params.outdir}/rfd4_proteina", pattern: 'combined_scores.tsv', mode: 'copy'
    publishDir path: "${params.outdir}/rfd4_proteina", pattern: 'binders.fasta', mode: 'copy'

    input:
    tuple path(design_table, stageAs: 'rfd4_designs.tsv'),
          path(refold_scores, stageAs: 'refold_scores.tsv'),
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

    # Refold scores: already one row per design, with engine-prefixed columns for
    # confidence, RMSD against the design and interface metrics (REFOLD_SCORES_WIDE).
    if [[ -s ${refold_scores} ]]; then
      python ${projectDir}/bin/merge_scores.py \\
        combined_scores.tsv ${refold_scores} \\
        --keys id \\
        --sort-by id \\
        --first-column id,filename,sequence,length \\
        -o combined_scores.tmp.tsv && mv combined_scores.tmp.tsv combined_scores.tsv
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
