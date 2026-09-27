// Monomer AlphaFold3 input JSON for fold.nf: one protein chain with the shared
// a3m as its unpaired MSA (paired MSA is query-only). Runs after batch fan-out
// because the per-batch seed lives in the JSON (AF3 has no --seed CLI flag).
process GENERATE_AF3_INPUT {
    tag "${meta.id}${meta.fold_batch ? " batch${meta.fold_batch}" : ''}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3m)

    output:
    tuple val(meta), path(fasta), path('chain_*_{unpaired,paired}.a3m'), path('af3_input.json'), emit: with_json

    script:
    """
    python ${projectDir}/bin/fold/make_af3_input.py \\
        --fasta ${fasta} \\
        --name '${meta.id}' \\
        --a3m ${a3m} \\
        --seed ${meta.af3_seed} \\
        -o af3_input.json
    """
}
