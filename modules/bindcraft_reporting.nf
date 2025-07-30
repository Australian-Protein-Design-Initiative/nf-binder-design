nextflow.enable.dsl = 2

process BINDCRAFT_REPORTING {
    publishDir "${params.outdir}/bindcraft", mode: 'copy'

    container 'ghcr.io/australian-protein-design-initiative/containers/mdanalysis:2.9.0'

    input:
    path('failure_csv.csv')
    path('final_design_stats.csv')
    path('mpnn_design_stats.csv')
    path('trajectory_stats.csv')

    output:
    path('bindcraft_report.html')

    script:
    def qmd_file = "${projectDir}/assets/bindcraft_reporting.qmd"
    """
    export XDG_CACHE_HOME="./.cache"
    export XDG_DATA_HOME="./.local/share"
    export JUPYTER_RUNTIME_DIR="./.jupyter"
    export XDG_RUNTIME_DIR="/tmp"
    quarto render ${qmd_file} --execute-dir \${PWD} --output - >bindcraft_report.html
    """
}
