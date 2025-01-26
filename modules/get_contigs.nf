process GET_CONTIGS {
    input:
    tuple path(input_pdb), val(chain_id)
    
    output:
    tuple path(input_pdb), env(contigs), emit: contigs
    
    script:
    def target_contigs_arg = params.target_contigs != 'auto' ? "--target_contigs '${params.target_contigs}'" : ''
    """
    # Get chain ranges and construct contigs string
    contigs=\$(python ${projectDir}/bin/get_contigs.py ${input_pdb} ${chain_id} ${target_contigs_arg})
    """
} 