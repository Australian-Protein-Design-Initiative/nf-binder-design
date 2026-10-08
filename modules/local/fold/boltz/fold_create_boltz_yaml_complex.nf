// Multimer Boltz-2 YAML for fold.nf: one `protein:` entry per chain, each with
// its own `msa:` pointing at that chain's key,sequence CSV (rendered by
// bin/fold/msa_taxonomy.py --tool boltz; Boltz pairs rows across chains sharing a
// taxid key). With --use_msa_server the msa: field is omitted and Boltz fetches
// + pairs its own MSA (the staged CSVs are then unused). CSVs arrive as a
// chain-ordered list from FOLD_MSA.
process FOLD_CREATE_BOLTZ_YAML_COMPLEX {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(csvs)
    path templates

    output:
    tuple val(meta), path(yaml), path(csvs), emit: yaml

    script:
    yaml = "${meta.id}.yml"
    def files = (csvs instanceof List) ? csvs : [csvs]
    def msa_arg = files.collect { it.name }.join(' ')
    def use_msa_server_flag = params.use_msa_server ? '--use_msa_server' : ''
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
    // meta.query_only_chains (fold_pulldown contract): chain letters that must
    // stay query-only (msa: empty) even under --use_msa_server, so the server
    // does not fetch an MSA for a chain --create_binder_msa false disabled it for.
    def query_only_chains = (meta.query_only_chains ?: []) as List
    def query_only_flag = query_only_chains ? "--query_only_chains ${query_only_chains.join(' ')}" : ''
    """
    ${projectDir}/bin/fold/make_boltz_complex_yaml.py \
        --fasta ${fasta} \
        --msa ${msa_arg} \
        --output_yaml ${yaml} \
        ${use_msa_server_flag} \
        ${templates_flag} \
        ${query_only_flag}
    """
}
