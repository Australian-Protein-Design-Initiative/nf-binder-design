// Multimer Boltz-2 YAML for fold.nf: one `protein:` entry per chain, each with
// its own `msa:` pointing at that chain's key,sequence CSV (rendered by
// bin/msa_taxonomy.py --tool boltz; Boltz pairs rows across chains sharing a
// taxid key). With --use_msa_server the msa: field is omitted and Boltz fetches
// + pairs its own MSA (the staged CSVs are then unused). CSVs arrive as a
// chain-ordered list from FOLD_MSA.
process FOLD_CREATE_BOLTZ_YAML_COMPLEX {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(csvs)

    output:
    tuple val(meta), path(yaml), path(csvs), emit: yaml

    script:
    yaml = "${meta.id}.yml"
    def files = (csvs instanceof List) ? csvs : [csvs]
    def msa_arg = files.collect { it.name }.join(' ')
    def use_msa_server_flag = params.use_msa_server ? '--use_msa_server' : ''
    """
    ${projectDir}/bin/make_boltz_complex_yaml.py \
        --fasta ${fasta} \
        --msa ${msa_arg} \
        --output_yaml ${yaml} \
        ${use_msa_server_flag}
    """
}
