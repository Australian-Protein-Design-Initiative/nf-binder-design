// Resolve the input structure path(s) a condition spec refers to.
//
// Runs as a process rather than on the driver because a spec may be YAML, and PyYAML is
// available in the container but not necessarily on the Nextflow head node.
process PARSE_RFD4_CONFIG {
    container 'ghcr.io/australian-protein-design-initiative/containers/rfd4-proteina:0.0.3_nv-cuda13-arm64'

    input:
    path config_spec
    val config_dir

    output:
    stdout

    script:
    """
    ${projectDir}/bin/rfd4_proteina/rfd4_config.py parse-inputs ${config_spec} --config-dir "${config_dir}"
    """
}
