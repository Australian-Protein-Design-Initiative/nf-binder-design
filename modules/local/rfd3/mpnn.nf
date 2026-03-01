process MPNN {
    container 'ghcr.io/australian-protein-design-initiative/containers/rc-foundry:0.1.11-weights'
    // container "rosettacommons/foundry:0.1.9-weights"

    publishDir path: "${params.outdir}/rfd3/mpnn", pattern: 'output/*.cif', mode: 'copy'
    publishDir path: "${params.outdir}/rfd3/mpnn", pattern: 'output/*.fa', mode: 'copy'

    input:
    path structure_file
    val mpnn_args

    output:
    path 'output/*.cif', emit: cifs
    path 'output/*.fa', emit: fastas

    script:
    """
    set -euo pipefail

    mkdir -p output

    mpnn \
        --structure_path "${structure_file}" \
        --out_directory output \
        ${mpnn_args} \
        ${task.ext.args ?: ''}
    """
}
