process GENERATE_RF3_INPUT_JSON {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(metas), path(structure_cifs, stageAs: 'cifs/*'), path(target_msa), path(target_templates), val(target_chain), val(binder_chain), val(binder_sequence_chain)

    output:
    tuple val(metas), path(structure_cifs), path(target_msa), path(target_templates), path('rf3.json'), emit: with_json

    script:
    def metaList = metas instanceof List ? metas : [metas]
    def cifList = structure_cifs instanceof List ? structure_cifs : [structure_cifs]
    def use_msa = (target_msa.name != 'empty_target_msa' && target_msa.name != 'boltz_will_make_target_msa') ? 'true' : 'false'
    def target_msa_arg = use_msa == 'true' ? "--target-msa ${target_msa}" : ''
    def structArgs = cifList.collect { "'" + it.toString().replace("'", "'\\''") + "'" }.join(' ')
    def nameArgs = metaList.collect { "'" + "${it.id}_rf3".replace("'", "'\\''") + "'" }.join(' ')
    """
    set -euo pipefail
    python ${projectDir}/bin/rfd3/make_rf3_input_spec.py json-batch \\
        --structure-cifs ${structArgs} \\
        --names ${nameArgs} \\
        --target-chain ${target_chain} \\
        --binder-chain ${binder_chain} \\
        --binder-sequence-chain ${binder_sequence_chain} \\
        --pdb-to-fasta ${projectDir}/bin/pdb_to_fasta.py \\
        ${target_msa_arg} \\
        --template-structure ${target_templates} \\
        --template-selection ${target_chain} \\
        --basename-paths \\
        -o rf3.json
    """
}
