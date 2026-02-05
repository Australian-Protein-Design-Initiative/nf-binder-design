process AF2IG_SCORE_FILTER {
    tag "af2ig_score_filter"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    publishDir "${params.outdir}/af2_initial_guess/filtered", mode: 'copy'

    input:
    path scores_tsv
    path pdb_files
    val filters

    output:
    path "accepted/*.pdb", emit: accepted, optional: true
    path "rejected/*.pdb", emit: rejected, optional: true
    path "af2ig_filtered_scores.tsv", emit: scores

    script:
    def filter_args_str = filters ? filters.split(';').collect { "--filter \"${it.trim()}\"" }.join(' ') : ''
    """
    # Create a list of PDB files for filter_designs.py
    echo "Creating PDB list..."
    find . -name "*.pdb" > pdb_list.txt

    # Run filter_designs.py
    filter_designs.py \\
        --pdb-list pdb_list.txt \\
        --scores ${scores_tsv}:description \\
        ${filter_args_str} \\
        --binder-chains A \\
        --collect-in . \\
        --pass-col-name af2ig_pass_filter \\
        --output af2ig_filtered_scores.tsv

    # Clean up
    rm pdb_list.txt
    """
}
