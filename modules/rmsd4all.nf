process RMSD4ALL {
    tag "${out_name}"
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.4'
    publishDir "${params.outdir}/${step_name}/rmsd", mode: 'copy'

    input:
    path 'fixed/*'
    path 'mobile/*'
    val superimpose_chains
    val score_chains
    val out_name
    val skip_tm
    val step_name

    output:
    path (out_name), emit: rmsd_tsv

    script:
    def tm_flag = skip_tm ? '' : '--tm-score'
    def superimpose_flag = superimpose_chains ? "--superimpose-chains ${superimpose_chains} --mobile-superimpose-chains ${superimpose_chains}" : ''
    def score_flag = score_chains ? "--score-chains ${score_chains} --mobile-score-chains ${score_chains}" : ''
    """
    set -euo pipefail

    # Run RMSD calculation: query vs provided directory
    ${projectDir}/bin/rmsd4all.py ${tm_flag} ${superimpose_flag} ${score_flag} fixed/ mobile/ > "${out_name}"
    """
}
