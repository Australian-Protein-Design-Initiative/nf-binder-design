process BC2_REPORTING {
    tag 'report'

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    publishDir "${params.outdir}/bc2", mode: 'copy'

    input:
    path summary_csv
    path ranked_csv
    path refolded_csv
    path trajectories_csv

    output:
    path 'bc2_report.html'

    script:
    def qmd_file = "${projectDir}/assets/bc2_reporting.qmd"
    """
    export XDG_CACHE_HOME="./.cache"
    export XDG_DATA_HOME="./.local/share"
    export JUPYTER_RUNTIME_DIR="./.jupyter"
    export XDG_RUNTIME_DIR="/tmp"

    cp ${qmd_file} .

    quarto render bc2_reporting.qmd --execute-dir \${PWD} --output - >bc2_report.html
    """
}
