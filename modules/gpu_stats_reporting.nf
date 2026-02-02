nextflow.enable.dsl = 2

process GPU_STATS_REPORTING {
    publishDir "${params.outdir}/logs", mode: 'copy'

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    input:
    path('gpu_stats.csv')

    output:
    path('gpu_stats_report.html')

    script:
    def qmd_file = "${projectDir}/assets/gpu_stats.qmd"
    """
    export XDG_CACHE_HOME="./.cache"
    export XDG_DATA_HOME="./.local/share"

    cp ${qmd_file} .

    quarto render gpu_stats.qmd --execute-dir \${PWD} --output - >gpu_stats_report.html
    """
}
