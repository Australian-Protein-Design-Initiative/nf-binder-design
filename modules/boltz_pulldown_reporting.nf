nextflow.enable.dsl=2

process BOLTZ_PULLDOWN_REPORTING {    
    publishDir "${params.outdir}/boltz_pulldown", mode: 'copy'

    container "ghcr.io/australian-protein-design-initiative/containers/mdanalysis:2.9.0"

    input:
    path('boltz_pulldown_report.qmd')
    path('boltz_pulldown.tsv')

    output:
    path("boltz_pulldown_report.html")

    script:
    def qmd_file = "${projectDir}/assets/boltz_pulldown_reporting.qmd"
    """
    export XDG_CACHE_HOME="./.cache"
    export XDG_DATA_HOME="./.local/share"
    quarto render ${qmd_file} --execute-dir \${PWD} --output - >boltz_pulldown_report.html
    """
}
