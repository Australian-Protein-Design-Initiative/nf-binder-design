/*
REFOLD_ORACLES: score binder designs by co-folding them against their target with
any set of folding oracles.

PROTOTYPE. Written for --method rfd4_proteina, but deliberately free of anything
rfd4-specific, so the other design workflows can adopt it. It takes designs as
tuple(meta, structure) and gives back one row per design.

Where BOLTZ_REFOLD_CORE is tied to Boltz-2 (its RMSD, ipSAE and BindCraft metrics
come out of Boltz-specific processes), this reuses FOLD_PULLDOWN_CORE, so the
oracle set is whatever FOLD_PREDICT supports: af2, af2_mono, boltz, rf3, protenix,
af3, openfold3, chai1, esmfold2, esmfold2_fast. Each design's binder sequence is
co-folded against the target as a one-target pulldown, and the per-(binder, tool)
summary is pivoted into tool-prefixed columns.

The trade-off against BOLTZ_REFOLD_CORE: this gives iptm, ipSAE and the cross-tool
consensus z-score for any oracle, but not the RMSD-to-design or binder-monomer
metrics, which need the design structure rather than just its sequence. The two
are complementary, and a run can use both (--full_refold_with boltz for the Boltz
metrics, --refold_oracles for breadth across oracles).

The fold-family parameters govern it (--n_predictions, --msa_method,
--create_target_msa, --use_msa_server, and each engine's own flags), not the
--full_refold_* parameters, which stay with the Boltz path.

Parameters are NOT declared here: Nextflow does not evaluate top-level `params.x = ...`
statements in an included script, so a calling workflow must declare them itself. The block
this subworkflow needs is in workflows/rfd4_proteina.nf under "Oracle refold scoring"; copy
it when adopting this subworkflow elsewhere.
*/

include { FOLD_PULLDOWN_CORE } from './fold_pulldown_core'
include { REFOLD_PREPARE_ORACLE_INPUTS } from '../../modules/local/common/refold_prepare_oracle_inputs'
include { REFOLD_ORACLE_SCORES_WIDE } from '../../modules/local/common/refold_oracle_scores_wide'

workflow REFOLD_ORACLES {

    take:
    ch_designs      // tuple(meta, structure): designs to score
    binder_chain    // String: binder chain id in the design
    target_chain    // String: target chain id in the design
    refold_max      // Integer or false: cap on how many designs to score
    oracles         // List<String>: FOLD_PREDICT method names
    msa_method      // 'jackhmmer_af2' | 'mmseqs2_colabfold'

    main:

    ch_selected = refold_max
        ? ch_designs.take(refold_max as int)
        : ch_designs

    // One task for the whole set: the target is read once and every binder sequence
    // comes out of the same pass, which also makes the "do the designs agree on the
    // target" check possible at all.
    REFOLD_PREPARE_ORACLE_INPUTS(
        ch_selected.map { meta, structure -> structure }.collect(),
        binder_chain,
        target_chain,
        params.refold_oracle_target_id,
    )

    // Re-read the FASTAs as (meta, single-record file) pairs, the shape
    // FOLD_PULLDOWN_CORE takes. This mirrors how --method fold_pulldown reads
    // --targets and --binders.
    ch_targets = REFOLD_PREPARE_ORACLE_INPUTS.out.targets
        .splitFasta(record: [id: true, seqString: true])
        .map { record -> [id: record.id, seq: record.seqString] }
        .merge(REFOLD_PREPARE_ORACLE_INPUTS.out.targets.splitFasta(file: true))

    ch_binders = REFOLD_PREPARE_ORACLE_INPUTS.out.binders
        .splitFasta(record: [id: true, seqString: true])
        .map { record -> [id: record.id, seq: record.seqString] }
        .merge(REFOLD_PREPARE_ORACLE_INPUTS.out.binders.splitFasta(file: true))

    FOLD_PULLDOWN_CORE(
        ch_targets,
        ch_binders,
        oracles,
        msa_method,
        // Query FASTAs for --templates matching; binders only with --binder_templates,
        // as in --method fold_pulldown.
        params.binder_templates
            ? REFOLD_PREPARE_ORACLE_INPUTS.out.targets.mix(REFOLD_PREPARE_ORACLE_INPUTS.out.binders)
            : REFOLD_PREPARE_ORACLE_INPUTS.out.targets,
    )

    REFOLD_ORACLE_SCORES_WIDE(
        FOLD_PULLDOWN_CORE.out.summary,
        params.refold_oracle_metrics,
    )

    emit:
    scores = REFOLD_ORACLE_SCORES_WIDE.out.scores
    pulldown_scores = FOLD_PULLDOWN_CORE.out.scores
    pulldown_summary = FOLD_PULLDOWN_CORE.out.summary
    predictions = FOLD_PULLDOWN_CORE.out.predictions
}
