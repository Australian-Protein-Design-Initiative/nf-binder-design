process PDB_TO_FASTA {
    //tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.4'

    input:
    path(pdb_file)
    val(chains)  //  string like: A,B

    output:
    path "${fasta_file}"

    script:
    clean_chains = chains ? chains.replaceAll(",", "") : ""
    chain_suffix = clean_chains ? "_${clean_chains}" : ""
    fasta_file = "${pdb_file.baseName}${chain_suffix}.fasta"
    """
    python ${projectDir}/bin/pdb_to_fasta.py ${pdb_file} --chains ${chains} >${fasta_file}
    """
}