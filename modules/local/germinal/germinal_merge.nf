process GERMINAL_MERGE {
    tag 'merge'

    publishDir "${params.outdir}/germinal", mode: 'copy', overwrite: true

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    input:
    path 'batch_dirs/*'

    output:
    path 'all_trajectories.csv', optional: true, emit: all_trajectories
    path 'accepted_designs.csv', optional: true, emit: accepted_designs
    path 'failure_counts.csv', optional: true, emit: failure_counts
    path 'config', type: 'dir', optional: true, emit: config
    path 'accepted', type: 'dir', optional: true, emit: accepted
    path 'trajectories', type: 'dir', optional: true, emit: trajectories
    path 'redesign_candidates', type: 'dir', optional: true, emit: redesign_candidates

    script:
    """
    set -euo pipefail

    batch_dirs=()
    for item in \$(find batch_dirs -mindepth 1 -maxdepth 1 | sort -V); do
        if [ -d "\$item" ]; then
            batch_dirs+=("\$item")
        fi
    done

    if [ \${#batch_dirs[@]} -eq 0 ]; then
        echo "ERROR: No batch directories found."
        exit 1
    fi

    ${baseDir}/bin/germinal/germinal_merge.py "\${batch_dirs[@]}" -o .
    """
}
