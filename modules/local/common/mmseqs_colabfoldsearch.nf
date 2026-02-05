// based on: https://github.com/nf-core/proteinfold/blob/master/modules/local/mmseqs_colabfoldsearch.nf

process MMSEQS_COLABFOLDSEARCH {
    tag "$meta.id"

    container "quay.io/nf-core/proteinfold_colabfold:1.1.1"

    input:
    tuple val(meta), path(fasta)
    path colabfold_envdb
    path uniref30

    output:
    tuple val(meta), path("**.a3m"), emit: a3m

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p db
    ln -r -s $uniref30/uniref30_* ./db
    ln -r -s $colabfold_envdb/colabfold_envdb* ./db

    /localcolabfold/colabfold-conda/bin/colabfold_search \\
        $args \\
        --threads $task.cpus ${fasta} \\
        ./db \\
        "result/"
    """
    /*
    usage: colabfold_search [-h] [--prefilter-mode {0,1,2}] [-s S] [--db1 DB1] [--db2 DB2] [--db3 DB3] [--db4 DB4] [--use-env {0,1}] [--use-env-pairing {0,1}]
                        [--use-templates {0,1}] [--filter {0,1}] [--mmseqs MMSEQS] [--expand-eval EXPAND_EVAL] [--align-eval ALIGN_EVAL] [--diff DIFF]
                        [--qsc QSC] [--max-accept MAX_ACCEPT] [--pairing_strategy PAIRING_STRATEGY] [--db-load-mode DB_LOAD_MODE] [--unpack {0,1}]
                        [--threads THREADS]                                                                                                                        
                        query dbbase base                                                                                                                          
                                        
    positional arguments:                                                                                                                                              
        query                 fasta files with the queries.
        dbbase                The path to the database and indices you downloaded and created with setup_databases.sh
        base                  Directory for the results (and intermediate files)

    options:
        -h, --help            show this help message and exit
        --prefilter-mode {0,1,2}
                                Prefiltering algorithm to use: 0: k-mer (high-mem), 1: ungapped (high-cpu), 2: exhaustive (no prefilter, very slow). See wiki for more
                                details: https://github.com/sokrypton/ColabFold/wiki#colabfold_search (default: 0)
        -s S                  MMseqs2 sensitivity. Lowering this will result in a much faster search but possibly sparser MSAs. By default, the k-mer threshold is
                                directly set to the same one of the server, which corresponds to a sensitivity of ~8. (default: None)
        --db1 DB1             UniRef database (default: uniref30_2302_db)
        --db2 DB2             Templates database (default: .)
        --db3 DB3             Environmental database (default: colabfold_envdb_202108_db)
        --db4 DB4             Environmental pairing database (default: spire_ctg10_2401_db)
        --use-env {0,1}       Use --db3 (default: 1)
        --use-env-pairing {0,1}
                                Use --db4 (default: 0)
        --use-templates {0,1}
                                Use --db2 (default: 0)
        --filter {0,1}        Filter the MSA by pre-defined align_eval, qsc, max_accept (default: 1)
        --mmseqs MMSEQS       Location of the mmseqs binary. (default: mmseqs)
        --expand-eval EXPAND_EVAL
                                e-val threshold for 'expandaln'. (default: inf)
        --align-eval ALIGN_EVAL
                                e-val threshold for 'align'. (default: 10)
        --diff DIFF           filterresult - Keep at least this many seqs in each MSA block. (default: 3000)
        --qsc QSC             filterresult - reduce diversity of output MSAs using min score thresh. (default: -20.0)
        --max-accept MAX_ACCEPT
                                align - Maximum accepted alignments before alignment calculation for a query is stopped. (default: 1000000)
        --pairing_strategy PAIRING_STRATEGY
                                pairaln - Pairing strategy. (default: 0)
        --db-load-mode DB_LOAD_MODE
                                Database preload mode 0: auto, 1: fread, 2: mmap, 3: mmap+touch (default: 0)
        --unpack {0,1}        Unpack results to loose files or keep MMseqs2 databases. (default: 1)
        --threads THREADS     Number of threads to use. (default: 64)
  */
}