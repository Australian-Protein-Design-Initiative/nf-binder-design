#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_PREPARE_QUERIES process

Extracts the design chain from protein structure files.
Accepts PDB, mmCIF, and gzipped variants. Always outputs PDB format.
Useful for extracting just the binder chain from design complexes before
running FoldSeek structural search.

Input:
    path input_files      - One or more structure files (.pdb, .cif, .pdb.gz, .cif.gz)
    val n_binder_chains   - Number of trailing chains that make up the binder, for
                            writers that place every target chain before it. 0 falls
                            back to taking the single chain with fewest atoms, which
                            is what a one-target/one-binder complex needs.

Output:
    path 'design_chains/*.pdb' - Individual PDB files containing only the design chain
*/

process FOLDSEEK_PREPARE_QUERIES {
    tag "prepare_foldseek"
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.7'

    input:
    path input_files
    val n_binder_chains

    output:
    path 'design_chains/*.pdb', emit: design_chains

    script:
    def chain_arg = (n_binder_chains as Integer) > 0 ? "--last-chains ${n_binder_chains}" : ''
    """
    python3 ${projectDir}/bin/foldseek/extract_design_chain.py design_chains ${chain_arg} ${input_files}
    """
}
