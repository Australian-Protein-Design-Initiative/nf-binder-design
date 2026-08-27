nextflow.enable.dsl = 2

process BOLTZ_PULLDOWN_REPORTING {
    publishDir "${params.outdir}/boltz_pulldown", mode: 'copy'

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    path('boltz_pulldown_report.qmd')
    path('boltz_pulldown.tsv')

    output:
    path('boltz_pulldown_report.html')

    script:
    def qmd_file = "${projectDir}/assets/boltz_pulldown_reporting.qmd"
    """
    export XDG_CACHE_HOME="./.cache"
    export XDG_DATA_HOME="./.local/share"
    export JUPYTER_RUNTIME_DIR="./.jupyter"
    export XDG_RUNTIME_DIR="/tmp"
    mkdir -p "\${XDG_CACHE_HOME}" "\${XDG_DATA_HOME}/quarto/logs" "\${JUPYTER_RUNTIME_DIR}"

    cp ${qmd_file} .

    quarto render boltz_pulldown_reporting.qmd --execute-dir \${PWD} --output - >boltz_pulldown_report.html
    """
}
