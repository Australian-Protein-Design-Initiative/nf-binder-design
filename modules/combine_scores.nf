process COMBINE_SCORES {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path 'scores/*'
    
    output:
    path "combined_scores.tsv", emit: combined
    
    script:
    """
    # Run the combination script
    python /scratch/projects/apdi/nf-binder-design/bin/af2_combine_scores.py scores --output combined_scores.tsv
    """
} 