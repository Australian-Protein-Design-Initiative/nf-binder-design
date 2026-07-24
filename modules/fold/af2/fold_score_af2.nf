// Score one AF2 run for fold.nf: bin/score_af2_run.py runs ipSAE (bin/ipsae.py)
// on each structure the run publishes to fold/predictions/ and emits normalized
// TSV rows (via bin/parse_fold_confidence.py). One process per AF2 run; rows for
// all kept models go to stdout, collected into af2_fold_scores.tsv by the
// ALPHAFOLD2 subworkflow. CPU-only (ipSAE is numpy); runs on the local executor.
process FOLD_SCORE_AF2 {
    tag "${meta.id} run${meta.af2_run}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(run_dir), val(pred_prefix)

    output:
    stdout

    script:
    def keep = meta.af2_keep_models ?: 'all'
    def no_relax = (params.af2_no_relax instanceof Boolean ? params.af2_no_relax
        : params.af2_no_relax.toString().toBoolean()) ? '--no-relax' : ''
    """
    python3 ${projectDir}/bin/score_af2_run.py \
        --run-dir "${run_dir}" \
        --id "${meta.id}" \
        --pred-prefix "${pred_prefix}" \
        --keep-models ${keep} \
        ${no_relax}
    """
}
