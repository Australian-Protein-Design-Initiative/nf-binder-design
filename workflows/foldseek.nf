#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FoldSeek standalone workflow

Runs FoldSeek structural similarity search against a database.
Can be used standalone via `--method foldseek` or the FOLDSEEK_SEARCH
subworkflow can be imported into other workflows.

Usage:
    nextflow run main.nf --method foldseek --input_pdbs 'designs/*.pdb'
    nextflow run main.nf --method foldseek --input_pdbs 'designs/*.pdb' --foldseek_database PDB
    nextflow run main.nf --method foldseek --input_pdbs 'designs/*.pdb' --foldseek_use_webserver true
    nextflow run main.nf --method foldseek --input_pdbs 'designs/*.pdb' --foldseek_databases_path /data/foldseek --foldseek_database CATH50
*/

// Workflow-level parameters
params.input_pdbs = false
params.foldseek_database = 'CATH50'
params.foldseek_databases_path = false
params.foldseek_use_webserver = false
params.foldseek_mode = '3diaa'
params.foldseek_maxaccept = 1
params.foldseek_gzip_output = false
params.foldseek_include_html_output = false
params.foldseek_cath_names_path = false

include { FOLDSEEK_SEARCH } from '../subworkflows/local/foldseek_search'

workflow FOLDSEEK {

    main:

    // Show help if no input specified
    if (params.input_pdbs == false) {
        log.info(
            """
        ==================================================================
        FOLDSEEK STRUCTURAL SEARCH
        ==================================================================

        Required arguments:
            --input_pdbs                 Glob pattern or path to input structure files
                                         (PDB, mmCIF; plain or gzipped)

        Database options:
            --foldseek_database          Database to search [default: CATH50]
                                         Local search (default):
                                           - CATH50 (default) — CATH domain DB, combined AF2+PDB at 50% seq.id.
                                           - PDB              — Protein Data Bank
                                           - Alphafold/UniProt — Full AF2 database (~700GB)
                                           - Alphafold/UniProt50 — AF2 clustered at 50% seq.id.
                                           - Alphafold/UniProt50-minimal — AF2 50% clusters (reps only)
                                           - Alphafold/Proteome — AF2 proteomes
                                           - Alphafold/Swiss-Prot — AF2 Swiss-Prot
                                           - ESMAtlas30       — ESM metagenomic atlas at 30% seq.id.
                                           - BFMD             — Big Fantastic Multimer Database
                                           - BFVD             — Big Fantastic Virus Database
                                         Remote search (--foldseek_use_webserver true):
                                           - pdb100           — PDB (complex-aware, includes multimers)
                                           - afdb50           — AlphaFold/UniProt50 (~37M structures)
                                           - afdb-swissprot   — AlphaFold/Swiss-Prot
                                           - afdb-proteome    — AlphaFold/Proteomes
                                           - cath50           — CATH50 domain database
                                           - bfmd             — BFMD multimers
                                           - BFVD             — BFVD viruses
                                           - mgnify_esm30     — MGnify-ESM30 metagenomic atlas
            --foldseek_databases_path    Path to local databases directory.
                                         Each database in a subdirectory: {path}/{db_name}/
                                         If unset, database is auto-downloaded.
            --foldseek_use_webserver     Use FoldSeek web API instead of local search [default: false]

        Search options:
            --foldseek_mode              Search mode [default: 3diaa]
            --foldseek_maxaccept         Max accepted alignments per query [default: 1]
                                         Use 0 or unset for unlimited hits

        Output options:
            --foldseek_gzip_output       Gzip TSV output [default: false]
            --foldseek_include_html_output  Include HTML report [default: false]
            --foldseek_cath_names_path    Path to local cath-names.txt (auto-downloaded if unset)
                                          CATH annotation is automatic when database name starts with "CATH"

        Examples:
            # Local search with auto-downloaded CATH50
            nextflow run main.nf --method foldseek --input_pdbs 'designs/*.pdb'

            # Local search with user-provided database
            nextflow run main.nf --method foldseek --input_pdbs 'designs/*.pdb' \\
                --foldseek_databases_path /data/foldseek --foldseek_database CATH50

            # Remote search against PDB
            nextflow run main.nf --method foldseek --input_pdbs 'designs/*.pdb' \\
                --foldseek_use_webserver true --foldseek_database pdb100

        """.stripIndent()
        )
        exit(1)
    }

    // Collect input PDBs
    ch_pdbs = Channel.fromPath(params.input_pdbs)
                     .map { pdb -> pdb }

    // Count inputs for metadata
    ch_count = ch_pdbs.count()
    ch_count.subscribe { n -> log.info "[FoldSeek] Found ${n} input PDB files" }

    // Collect all PDBs into a single list and create metadata
    ch_grouped = ch_pdbs.collect()
                        .map { pdbs -> [[id: 'foldseek_query', count: pdbs.size()], pdbs] }

    // Split into meta and pdb channels for the subworkflow
    ch_meta = ch_grouped.map { meta, pdbs -> meta }
    ch_pdb_list = ch_grouped.map { meta, pdbs -> pdbs }

    // Actually, the subworkflow expects separate channels:
    // ch_pdbs: individual PDB files, ch_meta: metadata
    // The subworkflow internally combines them
    ch_input_pdbs = Channel.fromPath(params.input_pdbs)
    ch_input_meta = Channel.of([id: 'foldseek_query'])

    FOLDSEEK_SEARCH(ch_input_pdbs, ch_input_meta)

    // Output (published via publishDir in processes)
}