process GENERATE_RF3_INPUT_JSON {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(structure_cif), path(target_msa), path(target_templates), val(target_chain), val(binder_chain)

    output:
    tuple val(meta), path(structure_cif), path(target_msa), path(target_templates), path('rf3.json'), emit: with_json

    script:
    def use_msa = (target_msa.name != 'empty_target_msa' && target_msa.name != 'boltz_will_make_target_msa') ? 'true' : 'false'
    def target_msa_arg = use_msa == 'true' ? "--target-msa ${target_msa}" : ''
    """
    set -euo pipefail
    python ${projectDir}/bin/rfd3/make_rf3_input_spec.py json \\
        --structure-cif ${structure_cif} \\
        --target-chain ${target_chain} \\
        --binder-chain ${binder_chain} \\
        --binder-sequence-chain ${binder_chain} \\
        --pdb-to-fasta ${projectDir}/bin/pdb_to_fasta.py \\
        --name "${meta.id}_rf3_config" \\
        ${target_msa_arg} \\
        --template-structure ${target_templates} \\
        --template-selection ${target_chain} \\
        --basename-paths \\
        -o rf3.json
    """
}
