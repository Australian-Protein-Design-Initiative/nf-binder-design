nextflow.enable.dsl = 2

process ENGENS {
    tag "${meta.id}"

    container 'docker://ghcr.io/australian-protein-design-initiative/containers/engens:2024-11-22_b097882_py311'

    publishDir(
        path: "${params.outdir}/engens/${meta.id}",
        mode: 'copy'
    )

    input:
    path qmd
    tuple val(meta), path(structures)

    output:
    tuple val(meta), path('clusters.html'), emit: report
    tuple val(meta), path('clustering'), emit: clustering, optional: true

    script:
    def clustering = params.engens_clustering ?: 'gmm'
    def dimred = params.engens_dimred ?: 'umap'
    def min_structures = params.engens_min_structures ?: 3
    def max_clusters = params.engens_max_clusters ?: 10
    def gmm_ic = params.engens_gmm_ic ?: 'aic'
    def seed = (params.engens_seed != null && !(params.engens_seed instanceof Boolean)) \
        ? params.engens_seed : ''
    """
    set -euo pipefail

    export XDG_CACHE_HOME="./.cache"
    export XDG_DATA_HOME="./.local/share"
    export MPLCONFIGDIR="./.matplotlib"
    mkdir -p "\${XDG_CACHE_HOME}" "\${XDG_DATA_HOME}" "\${MPLCONFIGDIR}"

    mkdir -p structures
    for f in ${structures}; do
        ln -s "\$(realpath "\$f")" "structures/\$(basename "\$f")"
    done

    export ENGENS_INPUT_DIR="\${PWD}/structures"
    export ENGENS_TARGET_ID="${meta.id}"
    export ENGENS_DIMRED="${dimred}"
    export ENGENS_CLUSTERING="${clustering}"
    export ENGENS_MIN_STRUCTURES="${min_structures}"
    export ENGENS_MAX_CLUSTERS="${max_clusters}"
    export ENGENS_GMM_IC="${gmm_ic}"
    export ENGENS_SEED="${seed}"
    export ENGENS_OUTDIR="\${PWD}"

    # Apptainer --cleanenv drops image ENV (MAMBA activate, LD_LIBRARY_PATH)
    # and keeps the host PATH, so conda/python are not visible to Quarto.
    export PATH="/opt/conda/bin:/usr/local/bin:\${PATH}"
    export QUARTO_PYTHON="/opt/conda/bin/python3"
    export LD_LIBRARY_PATH="/opt/conda/lib\${LD_LIBRARY_PATH:+:\${LD_LIBRARY_PATH}}"
    export HOME="\${PWD}"
    export IPYTHONDIR="\${PWD}/.ipython"

    quarto render ${qmd} \\
        --execute-dir "\${PWD}" \\
        --output clusters.html
    """
}
