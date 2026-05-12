#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_PREPARE_QUERIES process

Extracts the design (shorter) chain from protein structure files.
Accepts PDB, mmCIF, and gzipped variants. Always outputs PDB format.
Useful for extracting just the binder chain from design complexes before
running FoldSeek structural search.

Input:
    path input_files - One or more structure files (.pdb, .cif, .pdb.gz, .cif.gz)

Output:
    path 'design_chains/*.pdb' - Individual PDB files containing only the design chain
*/

process FOLDSEEK_PREPARE_QUERIES {
    tag "prepare_foldseek"
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    input:
    path input_files

    output:
    path 'design_chains/*.pdb', emit: design_chains

    script:
    """
    python3 ${projectDir}/bin/foldseek/extract_design_chain.py design_chains ${input_files}
    """
}
