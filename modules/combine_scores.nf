process COMBINE_SCORES {
  container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

  publishDir "${params.outdir}", pattern: 'combined_scores.tsv', mode: 'copy'
  publishDir "${params.outdir}", pattern: 'binders.fasta', mode: 'copy'

  input:
  path 'af2ig_scores/*'
  path 'extra_scores.tsv'
  // TODO: boltz_scores must be optional, detected by dummy files ?
  path 'boltz_scores_complex.tsv'
  path 'rmsd_monomer_vs_complex.tsv'
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
    python ${projectDir}/bin/af2_combine_scores.py af2ig_scores --output af2_initial_guess_scores.tsv

    # Run the shape score calculation script (Rg, Dmax, asphericity, Stokes Radius, chain, length, sequence)
    pushd pdbs
      python ${projectDir}/bin/calculate_shape_scores.py --chain A *.pdb >../shape_scores.tsv
    popd

    # Merge both tables
    python ${projectDir}/bin/merge_scores.py \
      af2_initial_guess_scores.tsv \
      shape_scores.tsv \
      extra_scores.tsv \
      --keys description,filename \
      --strip-suffix '(\\.pdb|_model_0\\.pdb)\$' \
        >prerefold_combined_scores.tsv

    # Only run this if the boltz scores and rmsd scores are non-empty files
    if [[ -s boltz_scores_complex.tsv && -s rmsd_monomer_vs_complex.tsv ]]; then
      # Merge boltz refolding scores with combined scores
      python ${projectDir}/bin/merge_scores.py \
        prerefold_combined_scores.tsv boltz_scores_complex.tsv \
        --column-prefix boltz_ \
        --first-column filename,description,pae_interaction,plddt_binder,boltz_confidence_score,boltz_iptm \
        --keys id,description \
        >refold_combined_scores.tsv

      python ${projectDir}/bin/merge_scores.py \
        refold_combined_scores.tsv rmsd_monomer_vs_complex.tsv \
        --column-prefix boltz_monomer_vs_complex_ \
        --first-column filename,description,pae_interaction,plddt_binder,boltz_confidence_score,boltz_iptm,boltz_monomer_vs_complex_rmsd_all \
        --keys structure1,description \
        --strip-suffix '(_model_0\\.pdb|_monomer)\$' \
        >combined_scores.tsv
    else
      cp prerefold_combined_scores.tsv combined_scores.tsv && \
      rm prerefold_combined_scores.tsv
    fi
  

    # Output FASTA sequences of binders, with scores in the header
    python ${projectDir}/bin/pdb_to_fasta.py \
        --scores-table combined_scores.tsv \
        --scores pae_interaction,plddt_binder,rg,length \
        --chain A \
        pdbs/*.pdb \
      >binders.fasta
    """
}
