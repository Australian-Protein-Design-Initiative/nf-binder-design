process RMSD4ALL {
    tag "${out_name}"
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.4'
    publishDir "${params.outdir}/boltz_pulldown/rmsd", mode: 'copy'

    input:
    path 'query_dir/*'
    val vs_dir
    val chains
    val out_name

    output:
    path (out_name), emit: rmsd_tsv

    script:
    def chains_flag = chains ? "--vs-chains ${chains}" : ''
    def out = out_name
    def mode = 'rmsd_pruned'
    //def mode = 'rmsd_all'
    """
    set -euo pipefail

    # Run RMSD calculation: query vs provided directory
    ${projectDir}/bin/rmsd4all.py query_dir/ --matrix ${mode} --vs "${vs_dir}" ${chains_flag} > "${out}"
    """
}
