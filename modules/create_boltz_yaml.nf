// Helper function to sanitize strings for filenames
def sanitize(name) {
    return name.replaceAll(/[^a-zA-Z0-9_.-]/, '_')
}

// Complex mode: both target and binder
process CREATE_BOLTZ_YAML {
    tag "${target_meta.id}_and_${binder_meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    input:
    tuple val(target_meta), path(target_msa), val(binder_meta), path(binder_msa)
    path(templates)

    output:
    tuple val(meta), path(yaml), path(target_msa), path(binder_msa)

    script:
    def id = "${sanitize(target_meta.id)}_and_${sanitize(binder_meta.id)}"
    def use_msa_server_flag = params.use_msa_server ? '--use_msa_server' : ''
    def templates_flag = params.templates ? "--templates '${templates}'" : ''
    def target_chain_flag = params.target_chains ? "--target_chain ${params.target_chains}" : ''
    def binder_chain_flag = params.binder_chains ? "--binder_chain ${params.binder_chains}" : ''

    /* we use 'flag files' empty_target_msa and boltz_will_make_target_msa
       to indicate wheather we want no MSA (empty), or allow Boltz to generate one
       via a remote mmseqs2 server (boltz_will_make_target_msa).
       If another (.a3m) file is provided, we use that as the MSA - this will generally
       come from a mmseqs2 search task in the pipeline (eg MMSEQS_COLABFOLDSEARCH).

       An alternative approach here might be to have mutliple variants of
       CREATE_BOLTZ_YAML (CREATE_BOLTZ_YAML_NO_MSA, CREATE_BOLTZ_YAML_MSA,
       CREATE_BOLTZ_YAML_BINDER_MSA_ONLY, CREATE_BOLTZ_YAML_TARGET_MSA_ONLY) and a
       conditional in the main pipeline that selects the appropriate variant.
    */
    def target_msa_flag = ''
    if (target_msa.name == 'empty_target_msa') {
        target_msa_flag = '--target_msa empty'
    } else if (target_msa.name == 'boltz_will_make_target_msa') {
        target_msa_flag = ''
    } else {
        target_msa_flag = "--target_msa '${target_msa}'"
    }

    def binder_msa_flag = ''
    if (binder_msa.name == 'empty_binder_msa') {
        binder_msa_flag = '--binder_msa empty'
    } else if (binder_msa.name == 'boltz_will_make_binder_msa') {
        binder_msa_flag = ''
    } else {
        binder_msa_flag = "--binder_msa '${binder_msa}'"
    }

    meta = [id: id, type: 'complex', target: sanitize(target_meta.id), binder: sanitize(binder_meta.id)]
    yaml = "${id}.yml"
    """
    ${projectDir}/bin/create_boltz_yaml.py \
        --target_id '${target_meta.id}' \
        --target_seq '${target_meta.seq}' \
        --binder_id '${binder_meta.id}' \
        --binder_seq '${binder_meta.seq}' \
        ${target_msa_flag} \
        ${binder_msa_flag} \
        ${target_chain_flag} \
        ${binder_chain_flag} \
        --output_yaml '${yaml}' \
        ${use_msa_server_flag} \
        ${templates_flag}
    """
}

// Monomer mode: single protein (target or binder)
process CREATE_BOLTZ_YAML_MONOMER {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    input:
    tuple val(meta), path(msa), val(type)
    path(templates)

    output:
    tuple val(meta), path(yaml), path(msa)

    script:
    def id = "${sanitize(meta.id)}"
    def use_msa_server_flag = params.use_msa_server ? '--use_msa_server' : ''
    def templates_flag = params.templates ? "--templates '${templates}'" : ''
    def yaml_name = "${id}.yml"

    // Handle MSA flags based on type and MSA file name
    def target_msa_flag = ''
    def binder_msa_flag = ''
    if (type == 'target') {
        if (msa.name == 'empty_target_msa') {
            target_msa_flag = '--target_msa empty'
        } else if (msa.name == 'boltz_will_make_target_msa') {
            target_msa_flag = ''
        } else {
            target_msa_flag = "--target_msa '${msa}'"
        }
    } else if (type == 'binder') {
        if (msa.name == 'empty_binder_msa') {
            binder_msa_flag = '--binder_msa empty'
        } else if (msa.name == 'boltz_will_make_binder_msa') {
            binder_msa_flag = ''
        } else {
            binder_msa_flag = "--binder_msa '${msa}'"
        }
    }

    // Handle chain flags based on type
    def target_chain_flag = ''
    def binder_chain_flag = ''
    if (type == 'target') {
        target_chain_flag = params.target_chains ? "--target_chain '${params.target_chains.split(',')[0]}'" : ''
    } else if (type == 'binder') {
        binder_chain_flag = params.binder_chains ? "--binder_chain '${params.binder_chains.split(',')[0]}'" : ''
    }

    // Handle protein arguments based on type
    def target_args = ''
    def binder_args = ''
    if (type == 'target') {
        target_args = "--target_id '${meta.id}' --target_seq '${meta.seq}'"
    } else if (type == 'binder') {
        binder_args = "--binder_id '${meta.id}' --binder_seq '${meta.seq}'"
    }

    // Update meta with type information
    meta = [id: id, type: type,
            target : (type == 'target' ? sanitize(meta.id) : ''),
            binder : (type == 'binder' ? sanitize(meta.id) : '')]
    yaml = yaml_name
    """
    ${projectDir}/bin/create_boltz_yaml.py \
        ${target_args} \
        ${binder_args} \
        ${target_msa_flag} \
        ${binder_msa_flag} \
        ${target_chain_flag} \
        ${binder_chain_flag} \
        --output_yaml '${yaml}' \
        ${use_msa_server_flag} \
        ${templates_flag}
    """
}
