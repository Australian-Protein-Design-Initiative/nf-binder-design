process RENUMBER_RESIDUES {
    input:
    tuple path(input_pdb), val(binder_chains)
    
    output:
    path "${input_pdb.baseName}_renum.pdb", emit: renumbered_pdb
    
    script:
    """
    # Get chain ranges and construct contigs string
    python ${projectDir}/bin/renumber_chains.py ${input_pdb} --binder-chains ${binder_chains} >${input_pdb.baseName}_renum.pdb
    """
}