// Assemble an AF2 multimer precomputed-MSA tree from per-chain ColabFold a3ms
// (FOLD_MSA's multimer per-chain search, --msa_method mmseqs2_colabfold), so
// --methods af2_mono works on a ColabFold-searched multimer input. AF2's native
// multimer pairing ('af2', not af2_mono) still requires --msa_method jackhmmer_af2
// (FoldValidation enforces this) - this process only ever feeds af2_mono, whose
// monomer-weights chain-break trick needs one MSA per chain, not a paired one.
//
// Layout written: <base_id>/msas/<CHAIN>/colabfold.a3m for each chain, then
// features.pkl built the same way as the native route (af2_multimer_features_from_msas.py,
// which now recognises colabfold.a3m as an MSA source - see A2 in the fold review).
process AF2_STAGE_COLABFOLD_MULTIMER_MSAS {
    tag "${meta.id}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/ghcr.io-australian-protein-design-initiative-containers-alphafold2-2.3.2-custom.img'

    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/af2/msas",
        mode: 'copy',
        saveAs: { filename ->
            def rel = filename.toString()
            return rel.startsWith("${meta.id}/") ? rel : null
        }
    )

    input:
    tuple val(meta), path(fasta), val(chain_ids), path(chain_a3ms)

    output:
    tuple val(meta), path(fasta), path("${meta.id}"), emit: msas

    script:
    def a3m_files = (chain_a3ms instanceof List) ? chain_a3ms : [chain_a3ms]
    def stage_cmds = chain_ids.withIndex().collect { cid, i ->
        "mkdir -p \"${meta.id}/msas/${cid}\"\n    cp \"${a3m_files[i]}\" \"${meta.id}/msas/${cid}/colabfold.a3m\""
    }.join('\n    ')
    """
    set -euo pipefail
    ${stage_cmds}

    python ${projectDir}/bin/fold/af2_multimer_features_from_msas.py \\
        --fasta "${fasta}" \\
        --msas-dir "${meta.id}"
    """
}
