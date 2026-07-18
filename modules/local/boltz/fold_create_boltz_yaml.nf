// Generic single-target Boltz YAML for fold.nf. Adapted from the (previously
// unwired) CREATE_BOLTZ_YAML_MONOMER in modules/local/common/create_boltz_yaml.nf,
// but takes the FASTA file directly (via bin/create_boltz_yaml.py's
// --binder_from_fasta) rather than a meta.seq string, since fold.nf's meta
// contract (see fold.nf §1) doesn't carry a sequence field.
//
// Phase 1 scope: monomer only, so this always emits a single "id: [A]"
// protein entry (bin/create_boltz_yaml.py's binder-only mode). Multimer
// (multiple protein entries, one id per chain) is Phase 2 work.

process FOLD_CREATE_BOLTZ_YAML {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3m)

    output:
    tuple val(meta), path(yaml), path(a3m), emit: yaml

    script:
    yaml = "${meta.id}.yml"
    def use_msa_server_flag = params.use_msa_server ? '--use_msa_server' : ''
    """
    ${projectDir}/bin/create_boltz_yaml.py \
        --binder_id '${meta.id}' \
        --binder_from_fasta '${fasta}' \
        --binder_msa '${a3m}' \
        --output_yaml '${yaml}' \
        ${use_msa_server_flag}
    """
}
