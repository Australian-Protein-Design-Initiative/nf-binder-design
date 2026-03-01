process EXTRACT_RFD3_BACKBONE_SCORES {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    input:
    tuple path(json_file), val(uid)

    output:
    path '*.tsv', emit: scores

    script:
    def batchMatcher = (json_file.name =~ /batch(\d+)|_b(\d+)_/)
    def batch = batchMatcher ? (batchMatcher[0][1] ?: batchMatcher[0][2]) : '0'
    def out_tsv = "rfd3_${uid}_batch${batch}_backbone.tsv"
    """
    set -euo pipefail
    python ${projectDir}/bin/rfd3/extract_rfd3_scores.py rfdiffusion3 ${json_file} -o ${out_tsv}
    """
}

process EXTRACT_RF3_SCORES {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'

    input:
    tuple path(confidence_json), val(uid)

    output:
    path '*.tsv', emit: scores

    script:
    def batchMatcher = (confidence_json.name =~ /batch(\d+)|_b(\d+)_/)
    def batch = batchMatcher ? (batchMatcher[0][1] ?: batchMatcher[0][2]) : '0'
    def suffixMatcher = (confidence_json.name =~ /(cif_b\d+_d\d+)/)
    def suffix = suffixMatcher ? suffixMatcher[0][1] : 'rf3'
    def out_tsv = "rfd3_${uid}_batch${batch}_rf3_${suffix}.tsv"
    """
    set -euo pipefail
    python ${projectDir}/bin/rfd3/extract_rfd3_scores.py rosettafold3 ${confidence_json} -o ${out_tsv}
    """
}
