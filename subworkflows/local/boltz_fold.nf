/*
BOLTZ_FOLD: generic single-target Boltz-2 folding for fold.nf.

Phase 1 scope: monomer only. Multimer (one `protein:` entry per FASTA record,
each chain's own `msa:`) is Phase 2 work (see plans/
fold-nf-multi-method-folding.md §5).
*/

include { FOLD_CREATE_BOLTZ_YAML } from '../../modules/fold/boltz/fold_create_boltz_yaml'
include { BOLTZ } from '../../modules/local/common/boltz'
include { FOLD_PARSE_BOLTZ_CONFIDENCE } from '../../modules/fold/boltz/fold_parse_boltz_confidence'

workflow BOLTZ_FOLD {
    take:
    ch_for_boltz // tuple(meta, fasta, a3m)

    main:
    FOLD_CREATE_BOLTZ_YAML(ch_for_boltz)

    ch_templates = params.templates ? file(params.templates) : file("${projectDir}/assets/dummy_files/empty_templates")

    // BOLTZ's process signature is shared with boltz_pulldown.nf's
    // target+binder complex mode, so it always expects two MSA paths to
    // stage (their names never appear in the script body - they only need to
    // be physically present alongside the YAML, since the YAML's `msa:`
    // field is what boltz predict actually reads). fold.nf has one MSA, so
    // the unused "target" slot gets the same empty-placeholder file
    // boltz_pulldown.nf uses for its own target-msa-less branches.
    ch_boltz_input = FOLD_CREATE_BOLTZ_YAML.out.yaml.map { meta, yaml, msa ->
        [meta, yaml, file("${projectDir}/assets/dummy_files/empty_target_msa"), msa]
    }

    // step_name is BOLTZ's publishDir subdir under params.outdir; use 'boltz'
    // so predictions land at <outdir>/boltz/ alongside boltz_fold_scores.tsv
    // (matching af2's <outdir>/af2/ and rf3's <outdir>/rf3/ layout).
    BOLTZ(ch_boltz_input, ch_templates, 'boltz')

    FOLD_PARSE_BOLTZ_CONFIDENCE(BOLTZ.out.confidence_json.join(BOLTZ.out.ipsae_tsv))

    ch_tsv = FOLD_PARSE_BOLTZ_CONFIDENCE.out.collectFile(
        name: 'boltz_fold_scores.tsv',
        storeDir: "${params.outdir}/boltz",
        keepHeader: true,
        skip: 1,
    )

    emit:
    predictions = BOLTZ.out.pdb
    confidence_json = BOLTZ.out.confidence_json
    tsv = ch_tsv
}
