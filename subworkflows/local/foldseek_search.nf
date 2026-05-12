#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_SEARCH subworkflow

Structural similarity search using FoldSeek against a database.
Automatically annotates results with CATH hierarchy descriptions when
using a CATH database (database name starts with "CATH").

Inputs:
    ch_pdbs    Channel of PDB files (individual paths)
    ch_meta    Channel of metadata map (val)

Parameters (from params scope):
    params.foldseek_use_webserver       Boolean  Use remote API [default: false]
    params.foldseek_databases_path      Path     Local databases base dir
    params.foldseek_database            String   Database name [default: 'CATH50']
    params.foldseek_cath_names_path     Path     Local path to pre-downloaded cath-names.txt

Outputs:
    tsv              TSV file with header + alignment results (or .tsv.gz if gzip enabled)
    tsv_annotated    Annotated TSV with CATH columns (empty channel if not a CATH database)
    html             Pretty HTML alignment report
    all_output       All output files from the search
    db_name          Name of the database searched
*/

include { FOLDSEEK_DOWNLOAD_DB } from '../../modules/local/foldseek/foldseek_download_db'
include { FOLDSEEK_DOWNLOAD_CATH_NAMES } from '../../modules/local/foldseek/foldseek_download_cath_names'
include { FOLDSEEK_LOCAL_SEARCH } from '../../modules/local/foldseek/foldseek_local_search'
include { FOLDSEEK_REMOTE_SEARCH } from '../../modules/local/foldseek/foldseek_remote_search'
include { FOLDSEEK_ANNOTATE_CATH } from '../../modules/local/foldseek/foldseek_annotate_cath'

workflow FOLDSEEK_SEARCH {

    take:
    ch_pdbs      // channel of PDB files (individual paths)
    ch_meta      // channel of metadata map (val)

    main:

    def use_remote = params.foldseek_use_webserver ?: false
    def db_name = params.foldseek_database ?: 'CATH50'
    def db_base_path = params.foldseek_databases_path ?: false
    def annotate_cath = db_name.toUpperCase().startsWith('CATH')
    def cath_names_path = params.foldseek_cath_names_path ?: false

    // Collect all PDBs into a single list, create input tuple with meta
    ch_input = ch_pdbs.collect().map { pdb_list -> tuple([id: 'foldseek_query'], pdb_list) }

    // --- CATH names resolution (only when annotation is enabled) ---
    ch_cath_names = Channel.empty()
    if (annotate_cath && !use_remote) {
        if (cath_names_path) {
            // Use pre-downloaded cath-names.txt
            def cn_file = new File("${cath_names_path}")
            if (cn_file.exists()) {
                ch_cath_names = Channel.fromPath(cath_names_path).first()
            } else {
                log.warn "[FoldSeek] --foldseek_cath_names_path set to ${cath_names_path} but file not found. Downloading instead."
                FOLDSEEK_DOWNLOAD_CATH_NAMES(Channel.of('trigger'))
                ch_cath_names = FOLDSEEK_DOWNLOAD_CATH_NAMES.out.names_file
            }
        } else if (db_base_path) {
            // Check if cath-names.txt exists alongside the databases
            def cn_file = new File("${db_base_path}/cath-names.txt")
            if (cn_file.exists()) {
                ch_cath_names = Channel.fromPath("${db_base_path}/cath-names.txt").first()
            } else {
                FOLDSEEK_DOWNLOAD_CATH_NAMES(Channel.of('trigger'))
                ch_cath_names = FOLDSEEK_DOWNLOAD_CATH_NAMES.out.names_file
            }
        } else {
            FOLDSEEK_DOWNLOAD_CATH_NAMES(Channel.of('trigger'))
            ch_cath_names = FOLDSEEK_DOWNLOAD_CATH_NAMES.out.names_file
        }
    }

    if (use_remote) {
        FOLDSEEK_REMOTE_SEARCH(
            ch_input,
            Channel.of(db_name),
            Channel.of(params.foldseek_mode ?: '3diaa')
        )

    } else {
        ch_db_dir = Channel.empty()
        ch_db_name = Channel.of(db_name)

        if (db_base_path) {
            // Check if database already exists at the specified path
            def db_path = new File("${db_base_path}/${db_name}")
            if (db_path.exists() && db_path.isDirectory()) {
                ch_db_dir = Channel.fromPath("${db_base_path}/${db_name}").first()
            } else {
                FOLDSEEK_DOWNLOAD_DB(Channel.of(db_name))
                ch_db_dir = FOLDSEEK_DOWNLOAD_DB.out.db
            }
        } else {
            FOLDSEEK_DOWNLOAD_DB(Channel.of(db_name))
            ch_db_dir = FOLDSEEK_DOWNLOAD_DB.out.db
        }

        FOLDSEEK_LOCAL_SEARCH(
            ch_input,
            ch_db_dir,
            ch_db_name
        )

        // --- CATH annotation (only for local search with CATH-like databases) ---
        if (annotate_cath) {
            // Combine the foldseek TSV (may be .tsv or .tsv.gz) with meta for annotation
            ch_foldseek_tsv = FOLDSEEK_LOCAL_SEARCH.out.tsv.mix(FOLDSEEK_LOCAL_SEARCH.out.tsv_gz)

            FOLDSEEK_ANNOTATE_CATH(
                ch_foldseek_tsv.map { f -> tuple([id: 'cath_annotate'], f) },
                ch_cath_names,
                ch_db_name
            )
        }
    }

    emit:
    tsv = use_remote ? FOLDSEEK_REMOTE_SEARCH.out.tsv_summary : FOLDSEEK_LOCAL_SEARCH.out.tsv.mix(FOLDSEEK_LOCAL_SEARCH.out.tsv_gz)
    tsv_annotated = annotate_cath ? FOLDSEEK_ANNOTATE_CATH.out.tsv_annotated.mix(FOLDSEEK_ANNOTATE_CATH.out.tsv_annotated_gz) : Channel.empty()
    html = use_remote ? Channel.empty() : FOLDSEEK_LOCAL_SEARCH.out.html
    all_output = use_remote ? FOLDSEEK_REMOTE_SEARCH.out.tsv_summary.collect() : FOLDSEEK_LOCAL_SEARCH.out.tsv.mix(FOLDSEEK_LOCAL_SEARCH.out.tsv_gz).mix(FOLDSEEK_LOCAL_SEARCH.out.html).collect()
    db_name_out = use_remote ? FOLDSEEK_REMOTE_SEARCH.out.db_name_out : FOLDSEEK_LOCAL_SEARCH.out.db_name_out
}
