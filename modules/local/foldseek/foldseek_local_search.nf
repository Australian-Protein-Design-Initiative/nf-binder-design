#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_LOCAL_SEARCH process

Runs a local FoldSeek structural search against a pre-built database.
Produces TSV with column headers and optionally HTML pretty-printed alignments.
*/

process FOLDSEEK_LOCAL_SEARCH {
    tag "${meta.id}"
    container 'quay.io/biocontainers/foldseek:10.941cd33--h5021889_1'
    publishDir "${params.outdir}/foldseek/${db_name}", mode: 'copy',
        pattern: '{foldseek_results.tsv,foldseek_results.tsv.gz,foldseek_results.html}'

    input:
    tuple val(meta), path(query_pdbs, stageAs: 'queries/*')
    path db_dir
    val db_name

    output:
    path "foldseek_results.tsv", emit: tsv, optional: true
    path "foldseek_results.tsv.gz", emit: tsv_gz, optional: true
    path "foldseek_results.html", emit: html, optional: true
    path "query_db*", emit: query_db
    path "result*", emit: result_db
    val db_name, emit: db_name_out

    script:
    def args = task.ext.args ?: ''
    def max_accept = params.foldseek_maxaccept != null ? params.foldseek_maxaccept : 1
    def gzip_output = params.foldseek_gzip_output ?: false
    def include_html = params.foldseek_include_html_output ?: false
    def format_columns = 'query,target,fident,alnlen,mismatch,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,qtmscore,ttmscore,alntmscore,rmsd,lddt,prob'
    """
    # Create query database from all staged structures (PDB, mmCIF, with optional .gz)
    foldseek createdb queries/* query_db

    # Find the database prefix inside the directory.
    # Names like CATH50 land directly under ${db_dir}, but names containing a
    # "/" (e.g. Alphafold/UniProt50) are nested in a subdir, and may be nested
    # one level deeper again if this is a published/cached DB (see saveAs
    # note in foldseek_download_db.nf). Try the deterministic path first,
    # then fall back to a recursive search for the base .dbtype file.
    DB_PREFIX="${db_dir}/${db_name}"
    if [ ! -f "\${DB_PREFIX}.dbtype" ]; then
        DB_PREFIX=""
        for f in \$(find ${db_dir} -name '*.dbtype'); do
            base=\$(basename "\$f" .dbtype)
            if [[ ! "\$base" =~ _ca\$ && ! "\$base" =~ _ss\$ && ! "\$base" =~ _h\$ && ! "\$base" =~ _clu\$ ]]; then
                DB_PREFIX="\$(dirname "\$f")/\$base"
                break
            fi
        done
    fi

    if [ -z "\$DB_PREFIX" ]; then
        echo "ERROR: Could not find FoldSeek database prefix in ${db_dir}"
        ls -la ${db_dir}/
        exit 1
    fi

    echo "Using database: \$DB_PREFIX"

    # Run exhaustive search with exact TM-scores
    # --max-accept limits accepted alignments per query (default: 1 = top hit)
    foldseek search \\
        query_db \\
        \$DB_PREFIX \\
        result \\
        ./tmp_search \\
        --exhaustive-search \\
        -s 9.0 \\
        --tmalign-fast 0 \\
        -a \\
        --max-seqs 500 \\
        --max-accept ${max_accept} \\
        ${args}

    # TSV with column headers (--format-mode 4) and full structural metrics
    foldseek convertalis \\
        query_db \\
        \$DB_PREFIX \\
        result \\
        foldseek_results.tsv \\
        --format-mode 4 \\
        --format-output ${format_columns} \\
        --exact-tmscore 1

    ${gzip_output ? '''
    # Gzip TSV output
    gzip foldseek_results.tsv
    ''' : ''}

    ${include_html ? '''
    # Pretty HTML alignments (--format-mode 3)
    foldseek convertalis \\
        query_db \\
        \\$DB_PREFIX \\
        result \\
        foldseek_results.html \\
        --format-mode 3 \\
        --exact-tmscore 1
    ''' : ''}
    """
}
