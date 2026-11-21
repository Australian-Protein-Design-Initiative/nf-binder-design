process BC2_CAMPAIGN_OUTPUT {
    tag 'campaign_output'

    container 'ghcr.io/australian-protein-design-initiative/containers/bindcraft2:v1.0.1'

    // `campaign_output` rewrites each campaign's ranking in place. Under the
    // default symlink staging that writes back into the BINDCRAFT2 work
    // directories, mutating upstream task outputs; copy so the gather cannot
    // corrupt what a -resume would reuse.
    stageInMode 'copy'

    publishDir path: "${params.outdir}/bc2", pattern: '*.csv', mode: 'copy'
    publishDir path: "${params.outdir}/bc2/ranked", pattern: 'ranked/*.cif', mode: 'copy',
        saveAs: { filename -> java.nio.file.Paths.get(filename).getFileName().toString() }

    input:
    path batch_dirs, stageAs: 'campaigns/*'

    output:
    path 'summary.csv',           emit: summary
    path 'ranked.csv',            optional: true, emit: ranked
    path 'refolded.csv',          optional: true, emit: refolded
    path 'trajectories.csv',      optional: true, emit: trajectories
    path 'campaigns/*/summary.csv', optional: true, emit: per_campaign_summaries
    path 'ranked/*.cif',          optional: true, emit: ranked_cifs

    script:
    """
    set -euo pipefail

    # BC2 writes caches under \$HOME by default; keep them in the task directory.
    export XDG_CACHE_HOME="\$(pwd)/.cache"
    mkdir -p "\$XDG_CACHE_HOME"

    # campaign_output treats each immediate subdirectory of the parent as one
    # campaign, rewrites each ranking, and writes a combined summary at the
    # parent. Passing the parent rather than a glob of the batch directories
    # also avoids ARG_MAX on runs with hundreds of batches.
    bindcraft campaign_output campaigns/ ${task.ext.args ?: ''}

    if [[ -f campaigns/summary.csv ]]; then
        cp campaigns/summary.csv summary.csv
    else
        # A single batch is treated as one campaign, so no combined file is
        # written at the parent; use that campaign's own summary.
        found=\$(find campaigns -mindepth 2 -maxdepth 2 -name summary.csv | head -n 1)
        if [[ -z "\${found}" ]]; then
            echo "ERROR: bindcraft campaign_output produced no summary.csv" >&2
            exit 1
        fi
        cp "\${found}" summary.csv
    fi

    # campaign_output does not concatenate the per-stage tables across
    # campaigns, so do it here while every batch is already staged.
    python3 ${projectDir}/bin/bindcraft2/merge_bc2_tables.py campaigns \\
        --ranked ranked.csv \\
        --refolded refolded.csv \\
        --trajectories trajectories.csv

    # Flatten the accepted structures into one directory. Done here rather than
    # in BINDCRAFT2 because campaign_output re-ranks, so only now is 3_Ranked
    # authoritative. Design names embed a hash of the recipe, which includes the
    # per-batch campaign_seed, so they do not collide across batches; prefix on
    # the off chance that two do, rather than silently overwrite.
    mkdir -p ranked
    find campaigns -mindepth 3 -maxdepth 3 -path '*/3_Ranked/*.cif' -print0 |
        while IFS= read -r -d '' cif; do
            campaign=\$(basename "\$(dirname "\$(dirname "\${cif}")")")
            name=\$(basename "\${cif}")
            if [[ -e "ranked/\${name}" ]]; then
                name="\${campaign}_\${name}"
            fi
            cp "\${cif}" "ranked/\${name}"
        done
    """
}
