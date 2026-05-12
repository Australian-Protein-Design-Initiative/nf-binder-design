#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_REMOTE_SEARCH process

Submits PDB queries to the FoldSeek web API (https://search.foldseek.com)
and retrieves JSON results.

Available remote databases: pdb100, afdb50, afdb-swissprot, afdb-proteome,
                             cath50, bfmd, BFVD, mgnify_esm30
*/

process FOLDSEEK_REMOTE_SEARCH {
    tag "${meta.id}"
    container 'quay.io/biocontainers/foldseek:10.941cd33--h5021889_1'
    publishDir "${params.outdir}/foldseek/${db_name}", mode: 'copy'

    input:
    tuple val(meta), path(query_pdbs, stageAs: 'queries/*')
    val db_name
    val mode

    output:
    path "foldseek_remote_results/*.json", emit: json_results
    path "foldseek_remote_results.tsv", emit: tsv_summary
    val db_name, emit: db_name_out

    script:
    def args = task.ext.args ?: ''
    def search_mode = mode ?: '3diaa'
    """
    mkdir -p foldseek_remote_results

    # Submit all PDBs and collect ticket IDs
    for pdb in queries/*.pdb; do
        base=\$(basename "\$pdb" .pdb)
        echo "Submitting \$base to FoldSeek web API..."

        RESPONSE=\$(curl -s -X POST "https://search.foldseek.com/api/ticket" \\
            -F "q=@\$pdb" \\
            -F "database[]=${db_name}" \\
            -F "mode=${search_mode}" \\
            ${args})

        TICKET=\$(echo "\$RESPONSE" | grep -oP '"id":"\\K[^"]+' 2>/dev/null || true)
        STATUS=\$(echo "\$RESPONSE" | grep -oP '"status":"\\K[^"]+' 2>/dev/null || echo "ERROR")

        if [ -z "\$TICKET" ] || [ "\$STATUS" = "ERROR" ]; then
            echo "ERROR: Failed to submit \$base"
            echo "\$RESPONSE"
            exit 1
        fi

        echo "\$base \$TICKET" >> tickets.txt
        echo "  Submitted: ticket=\$TICKET status=\$STATUS"
    done

    # Poll until all jobs complete (max 30 min)
    MAX_POLLS=180
    SLEEP_INTERVAL=10
    POLL_COUNT=0

    while [ \$POLL_COUNT -lt \$MAX_POLLS ]; do
        ALL_DONE=true
        while IFS=' ' read -r design ticket; do
            STATUS=\$(curl -s "https://search.foldseek.com/api/ticket/\$ticket" \\
                | grep -oP '"status":"\\K[^"]+' 2>/dev/null || echo "PENDING")

            if [ "\$STATUS" != "COMPLETE" ]; then
                ALL_DONE=false
            fi
        done < tickets.txt

        if [ "\$ALL_DONE" = true ]; then
            echo "All jobs complete!"
            break
        fi

        echo "Poll \$((POLL_COUNT+1))/\$MAX_POLLS - some jobs still running..."
        sleep \$SLEEP_INTERVAL
        POLL_COUNT=\$((POLL_COUNT+1))
    done

    if [ \$POLL_COUNT -ge \$MAX_POLLS ]; then
        echo "ERROR: Timeout waiting for FoldSeek results"
        exit 1
    fi

    # Download results
    while IFS=' ' read -r design ticket; do
        curl -s "https://search.foldseek.com/api/result/\$ticket/0" \\
            > "foldseek_remote_results/\${design}.json"
        echo "Downloaded results for \$design"
    done < tickets.txt

    # Convert JSON to TSV using jq
    echo "query_design\\tdatabase\\ttarget\\tseqId\\talnLength\\tprob\\teval\\tscore\\tqLen\\tdbLen" > foldseek_remote_results.tsv
    for jf in foldseek_remote_results/*.json; do
        design=\$(basename "\$jf" .json)
        jq -r --arg d "\$design" '
            .results[]?.alignments[0][]? |
            [$d, .database // "?", .target // "", .seqId // 0, .alnLength // 0,
             .prob // 0, .eval // 0, .score // 0, .qLen // 0, .dbLen // 0] |
            @tsv
        ' "\$jf" >> foldseek_remote_results.tsv 2>/dev/null || true
    done

    # Fallback: if jq is not available, just ship the JSON files
    if [ ! -s foldseek_remote_results.tsv ]; then
        echo "WARNING: Could not parse JSON results (jq may not be available)"
        echo "query_design\\tdatabase\\ttarget" > foldseek_remote_results.tsv
    fi
    """
}
