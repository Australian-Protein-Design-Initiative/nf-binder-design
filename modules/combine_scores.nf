process COMBINE_SCORES {
    container 'ghcr.io/australian-protein-design-initiative/containers/mdanalysis:2.8.0'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path 'scores/*'
    path 'extra_scores.tsv'
    path 'pdbs/*'

    output:
    path 'af2_initial_guess_scores.tsv', emit: af2_initial_guess_scores
    path 'extra_scores.tsv', emit: extra_scores
    path 'shape_scores.tsv', emit: shape_scores
    path 'combined_scores.tsv', emit: combined_scores
    path 'binders.fasta', emit: binders_fasta

    script:
    """
    # Run the af2 score aggregation script
    python ${projectDir}/bin/af2_combine_scores.py scores --output af2_initial_guess_scores.tsv

    # Run the shape score calculation script (Rg, Dmax, asphericity, Stokes Radius, chain, length, sequence)
    pushd pdbs
      python ${projectDir}/bin/calculate_shape_scores.py --chain A *.pdb >../shape_scores.tsv
    popd

    # Merge both tables
    python ${projectDir}/bin/merge_scores.py \
      af2_initial_guess_scores.tsv \
      shape_scores.tsv \
      extra_scores.tsv \
        >combined_scores.tsv

    # Output FASTA sequences of binders, with scores in the header
    python ${projectDir}/bin/pdb_to_fasta.py \
        --scores-table combined_scores.tsv \
        --scores pae_interaction,rg,length \
        --chain A \
        pdbs/*.pdb \
      >binders.fasta
    """
}
