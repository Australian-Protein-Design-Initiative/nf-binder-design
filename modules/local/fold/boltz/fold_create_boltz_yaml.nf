// Monomer Boltz-2 YAML for fold.nf: a single "id: [A]" protein entry with the
// shared a3m as its msa:. Uses the same generator as the multimer path so chain
// templates (--templates) are handled in one place.

process FOLD_CREATE_BOLTZ_YAML {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3m)
    path templates

    output:
    tuple val(meta), path(yaml), path(a3m), emit: yaml

    script:
    yaml = "${meta.id}.yml"
    def template_chains = (meta.template_chains ?: []) as List
    // Boltz force holds a templated chain within --boltz_template_threshold of the
    // template. Off by default for fold (guide, don't over-bias an unknown), on for
    // fold_pulldown (chain structures are known; only the pose is predicted).
    def force = params.boltz_template_force != null ? params.boltz_template_force : (params.method == 'fold_pulldown')
    def templates_flag = params.templates ? [
        "--templates '${templates}'",
        template_chains ? "--template_chains ${template_chains.join(' ')}" : '',
        force ? "--template_force --template_threshold ${params.boltz_template_threshold}" : '',
    ].findAll { it }.join(' ') : ''
    """
    ${projectDir}/bin/fold/make_boltz_complex_yaml.py \
        --fasta ${fasta} \
        --msa ${a3m} \
        --output_yaml ${yaml} \
        ${templates_flag}
    """
}
