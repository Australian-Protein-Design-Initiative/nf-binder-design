process COMBINE_SCORES {
  container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.8'

  publishDir "${params.outdir}/rfd", pattern: 'combined_scores.tsv', mode: 'copy'
  publishDir "${params.outdir}/rfd", pattern: 'binders.fasta', mode: 'copy'

  input:
  path 'af2ig_scores/*'
  path 'extra_scores.tsv'
  // REFOLD's wide table, or an empty placeholder when nothing was refolded. One file
  // rather than the three Boltz-specific tables this took before: the refold reduces
  // every engine to one row per design, with engine-prefixed columns.
  path 'refold_scores.tsv'
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

    # Refold scores, when --refold_with named any engines. Already one row per design,
    # keyed by the design id, which is `description` on this side.
    if [[ -s refold_scores.tsv ]]; then
      python ${projectDir}/bin/merge_scores.py \
        prerefold_combined_scores.tsv refold_scores.tsv \
        --keys description,id \
        --first-column filename,description,pae_interaction,plddt_binder \
        >combined_scores.tsv
    else
      cp prerefold_combined_scores.tsv combined_scores.tsv && \
      rm prerefold_combined_scores.tsv
    fi
  

    # Output FASTA sequences of binders, with scores in the header
    python ${projectDir}/bin/pdb_to_fasta.py \
        --scores-table combined_scores.tsv \
        --scores pae_interaction,plddt_binder,rg,length \
        --chains A \
        pdbs/*.pdb \
      >binders.fasta
    """
}
