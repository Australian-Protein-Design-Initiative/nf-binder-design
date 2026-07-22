/*
ENGENS_CLUSTER: conformational clustering of predicted structures from fold.nf.

Takes a channel of [meta, structure_files] (already uniquely named), stages
them if needed, and renders assets/engens/engens-analysis.qmd to
results/engens/<id>/clusters.html with representative PDBs under
clustering/<featurizer>/<gmm|km|hdbscan>/conformations/
(e.g. residue_mindist, backbone_torsions-residue_mindist).
*/

include { ENGENS_STAGE_STRUCTURE } from '../../modules/local/engens/engens_stage_structure'
include { ENGENS } from '../../modules/local/engens/engens_cluster'

// Filter a prediction emit (often a nested directory dump) down to the
// structure files that fold.nf also publishes into fold/predictions/.
def engensStructurePaths(files) {
    def list = (files instanceof List) ? files.flatten() : [files]
    def no_relax = params.af2_no_relax instanceof Boolean \
        ? params.af2_no_relax \
        : params.af2_no_relax.toString().toBoolean()
    def keep_best = params.af2_keep_models == 'best'
    return list.findAll { f ->
        def bn = f.getName().toString()
        if (!(bn ==~ /.*\.(cif|pdb)$/) || (bn ==~ /.*empty\.(cif|pdb)$/)) {
            return false
        }
        // AF2: match the flat fold/predictions/ publish selection.
        if (bn ==~ /relaxed_model_.*\.cif/) {
            return !no_relax
        }
        if (bn ==~ /unrelaxed_model_.*\.cif/) {
            return no_relax && !keep_best
        }
        if (bn == 'ranked_0.cif') {
            return no_relax && keep_best
        }
        // Boltz: <id>_model_N.cif|pdb
        if (bn ==~ /.*_model_\d+\.(cif|pdb)/ && !(bn ==~ /.*_sample.*/)) {
            return true
        }
        // RF3: *_sample-N_model.cif
        if (bn ==~ /.*_sample-\d+_model\.cif/) {
            return true
        }
        // Protenix: *_sample_N.cif
        if (bn ==~ /.*_sample_\d+\.cif/) {
            return true
        }
        return false
    }
}

// Collision-free basename matching fold/predictions/ publish naming
// (includes _msa<depth> when --msa_subsample fans out depth jobs).
def engensUniqueName(meta, path) {
    def bn = path.getName().toString()
    def msa_bit = meta.msa_depth_tag ? "_msa${meta.msa_depth_tag}" : ''
    if (bn ==~ /relaxed_model_.*\.cif/ || bn ==~ /unrelaxed_model_.*\.cif/ || bn == 'ranked_0.cif') {
        def run = meta.af2_run != null ? "run${meta.af2_run}" : ''
        def run_bit = run ? "_${run}" : ''
        return "af2_${meta.id}${run_bit}${msa_bit}_${bn}"
    }
    if (bn ==~ /.*_model_\d+\.(cif|pdb)/ && !(bn ==~ /.*_sample.*/)) {
        return meta.fold_namespaced \
            ? "boltz_batch${meta.fold_batch}${msa_bit}_${bn}" \
            : "boltz${msa_bit}_${bn}"
    }
    if (bn ==~ /.*_sample-\d+_model\.cif/) {
        return meta.fold_namespaced \
            ? "rf3_batch${meta.fold_batch}${msa_bit}_${bn}" \
            : "rf3${msa_bit}_${bn}"
    }
    if (bn ==~ /.*_sample_\d+\.cif/) {
        return meta.fold_namespaced \
            ? "protenix_batch${meta.fold_batch}${msa_bit}_${bn}" \
            : "protenix${msa_bit}_${bn}"
    }
    def batch = meta.fold_batch != null ? "_batch${meta.fold_batch}" : ''
    return "${meta.id}${batch}${msa_bit}_${bn}"
}

workflow ENGENS_CLUSTER {
    take:
    ch_predictions // tuple(meta, path|list) - mixed method prediction emits

    main:
    // Flatten -> filter to structures -> rename uniquely -> group by target id.
    ch_named = ch_predictions
        .flatMap { meta, files ->
            engensStructurePaths(files).collect { f ->
                [meta.id.toString(), engensUniqueName(meta, f), f]
            }
        }

    ENGENS_STAGE_STRUCTURE(ch_named)

    ch_grouped = ENGENS_STAGE_STRUCTURE.out.staged
        .groupTuple()
        .map { id, paths -> [[id: id], paths] }

    ch_qmd = Channel.value(file("${projectDir}/assets/engens/engens-analysis.qmd"))
    ENGENS(ch_qmd, ch_grouped)

    emit:
    report = ENGENS.out.report
    clustering = ENGENS.out.clustering
}
