// RFD3 output: chain A = target, chain B = binder. Filter uses binder chain B for metrics (e.g. Rg).
process FILTER_DESIGNS {
    tag "filter_${cif.name}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    publishDir "${params.outdir}/rfd3/filter/scores", mode: 'copy', pattern: '*.scores.tsv'
    publishDir "${params.outdir}/rfd3/filter/filtered", mode: 'copy', pattern: 'accepted/*'
    publishDir "${params.outdir}/rfd3/filter/filtered", mode: 'copy', pattern: 'rejected/*'

    input:
    path cif
    val filters
    val binder_chains
    val step

    output:
    path "accepted/${cif.name}", emit: accepted, optional: true
    path "rejected/${cif.name}", emit: rejected, optional: true
    path "${cif.baseName}.scores.tsv", emit: scores

    script:
    def filter_args_str = filters ? filters.split(';').collect { "--filter \"${it}\"" }.join(' ') : ''
    """
    set -euo pipefail

    input_file="${cif}"
    input_name=\$(basename "${cif}")
    filter_input="\${input_name}"

    if [[ "\${input_name}" == *.gz ]]; then
        gunzip -c "\${input_file}" >"\${input_name%.gz}"
        filter_input="\${input_name%.gz}"
    fi

    filter_designs.py \\
        "\${filter_input}" \\
        ${filter_args_str} \\
        --binder-chains ${binder_chains} \\
        --collect-in . \\
        --output ${cif.baseName}.scores.tsv

    if [[ "\${input_name}" == *.gz ]]; then
        if [[ -f "accepted/\${filter_input}" ]]; then
            gzip -f "accepted/\${filter_input}"
        fi
        if [[ -f "rejected/\${filter_input}" ]]; then
            gzip -f "rejected/\${filter_input}"
        fi
    fi
    """
}
