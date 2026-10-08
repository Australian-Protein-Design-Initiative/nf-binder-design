// fold.nf-local confidence parser: reuses bin/parse_boltz_confidence.py (the
// same script the Boltz refold comparison modules use) without
// the target/binder metadata columns that script's binder-design callers add
// - fold.nf's meta has no target/binder split (it's a single fold target).
//
// Runs once per diffusion sample (Boltz model_0..N-1) so every one of the
// --n_predictions structures gets a row in boltz_fold_scores.tsv, not just the
// top-ranked model_0. ipSAE is merged from BOLTZ.out.ipsae_tsv (ipsae.py runs
// unconditionally in boltz.nf, once per sample); for a monomer fold there are
// no chain pairs, so the merge is a no-op (Type==min row absent).
process FOLD_PARSE_BOLTZ_CONFIDENCE {
    tag "${meta.id}${meta.fold_batch ? "_batch${meta.fold_batch}" : ''}_model_${model}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), val(model), val(original_file), val(predictions_file), path(json_file), path(ipsae_tsv)

    output:
    stdout

    script:
    // id stays plain meta.id (matching every other engine) - batch/msa_depth
    // columns (not a mangled id) are what keep rows unique across
    // --n_predictions batches / --msa_subsample depth jobs; fold_pulldown's
    // summariser looks ids up in pairs.tsv, so a "${id}_batchN" id broke that.
    """
    python3 ${projectDir}/bin/parse_boltz_confidence.py \
        --json "${json_file}" \
        --id "${meta.id}" \
        --model "${model}" \
        --batch "${meta.fold_batch ?: ''}" \
        --msa-depth "${meta.msa_depth_tag ?: ''}" \
        --original-file "${original_file}" \
        --predictions-file "${predictions_file}" \
        --merge-ipsae "${ipsae_tsv}"
    """
}
