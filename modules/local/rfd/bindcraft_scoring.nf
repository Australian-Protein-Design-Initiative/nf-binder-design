process BINDCRAFT_SCORING {
  container 'ghcr.io/australian-protein-design-initiative/containers/bindcraft:0366085-nv-cuda120'

  publishDir path: { "${params.outdir}/${extra_scores_publish_dir}" }, pattern: '*.tsv', mode: 'copy'

  input:
  path pdb_file
  val binder_chain
  val advanced_settings_preset
  val extra_scores_publish_dir

  output:
  path '*.tsv', emit: scores

  script:
  """
    /opt/conda/envs/BindCraft/bin/python ${projectDir}/bin/bindcraft/bindcraft_scoring.py \
      --format tsv \
      --output ${pdb_file.simpleName}.tsv \
      --binder-chain ${binder_chain} \
      --dalphaball-path /app/BindCraft/functions/DAlphaBall.gcc \
      --dssp-path /app/BindCraft/functions/dssp \
      ${pdb_file}

    # --omit-aas is left at its default, which is already CX. Note that a \$-less
    # reference to a param here would be interpolated by Groovy even inside a shell
    # comment, and warn about an undefined parameter on every run.
    """
}
