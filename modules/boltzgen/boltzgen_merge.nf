process BOLTZGEN_MERGE {
    tag "merge"

    container 'ghcr.io/australian-protein-design-initiative/containers/boltzgen:0.2.0'

    publishDir path: "${params.outdir}/boltzgen", pattern: '**', mode: 'copy'

    input:
    path batch_dirs

    output:
    path 'merged', type: 'dir', emit: merged_dir

    script:
    """
    set -euo pipefail

    # With _many_ batch_ directories, this will fail due to ARG_MAX limits with shell glob expansion
    #boltzgen merge batch_* \
    #    --output merged/ \
    #    ${task.ext.args ?: ''}

    # We need to do some shell-gymnastics to avoid ARG_MAX limits with shell glob expansion
    # Find batch directories, including symlinks to directories (Nextflow creates symlinks for inputs)
    # Use -L to follow symlinks, or check -d which works for both directories and symlinks to directories
    batch_dirs=()
    for item in \$(find . -maxdepth 1 -name 'batch_*' | sort); do
        if [ -d "\$item" ]; then
            batch_dirs+=("\$item")
        fi
    done
    
    if [ \${#batch_dirs[@]} -eq 0 ]; then
        echo "ERROR: No batch directories found."
        exit 1
    fi
    
    boltzgen merge "\${batch_dirs[@]}" \
        --output merged/ \
        ${task.ext.args ?: ''}

    if grep -q "No designs found to merge." .command.log; then
        echo "ERROR: No designs found to merge."
        exit 1
    fi
    """
}
