process BINDCRAFT_SCORING {
    container 'ghcr.io/australian-protein-design-initiative/containers/bindcraft:nv-cuda12'

    publishDir path: "${params.outdir}/af2_initial_guess/extra_scores/", pattern: '*.tsv', mode: 'copy'

    input:
    path pdb_file
    val binder_chain
    val advanced_settings_preset

    output:
    path '*.tsv', emit: scores

    script:
    def advanced_settings_file = advanced_settings_preset ? "/app/BindCraft/settings_advanced/${advanced_settings_preset}.json" : '/app/BindCraft/settings_advanced/default_4stage_multimer.json'

    """
    /opt/conda/envs/BindCraft/bin/python ${projectDir}/bin/bindcraft_scoring.py \
      --format tsv \
      --output ${pdb_file.simpleName}.tsv \
      --binder-chain ${binder_chain} \
      --dalphaball-path /app/BindCraft/functions/DAlphaBall.gcc \
      --dssp-path /app/BindCraft/functions/dssp \
      ${pdb_file}

    # We don't use this, just the default (which is already CX)
    # --omit-aas ${params.pmpnn_omit_aas}
    """
}
