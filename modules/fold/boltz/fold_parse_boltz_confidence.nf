// fold.nf-local confidence parser: reuses bin/parse_boltz_confidence.py (the
// same script boltz_pulldown.nf's PARSE_BOLTZ_CONFIDENCE_JSON uses) without
// the target/binder metadata columns that script's binder-design callers add
// - fold.nf's meta has no target/binder split (it's a single fold target).
//
// Runs once per diffusion sample (Boltz model_0..N-1) so every one of the
// --n_predictions structures gets a row in boltz_fold_scores.tsv, not just the
// top-ranked model_0. ipSAE is deliberately not merged here: fold.nf is
// monomer-only (Phase 1) and ipSAE is an inter-chain interface metric, so it
// carries no signal for a single-chain fold.
process FOLD_PARSE_BOLTZ_CONFIDENCE {
    tag "${meta.id}${meta.fold_batch ? "_batch${meta.fold_batch}" : ''}_model_${model}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), val(model), path(json_file)

    output:
    stdout

    script:
    // When fold.nf splits --n_predictions across jobs, model indices restart
    // per batch - qualify the id so boltz_fold_scores.tsv rows stay unique.
    def row_id = meta.fold_namespaced ? "${meta.id}_batch${meta.fold_batch}" : meta.id
    """
    python3 ${projectDir}/bin/parse_boltz_confidence.py \
        --json "${json_file}" \
        --id "${row_id}" \
        --model "${model}"
    """
}
