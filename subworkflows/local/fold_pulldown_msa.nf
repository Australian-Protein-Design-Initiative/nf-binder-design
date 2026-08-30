/*
FOLD_PULLDOWN_MSA: build per-sequence MSAs once, then fan out target x binder
pairs into the channel shapes FOLD_PREDICT expects.

No cross-chain MSA pairing: binders are treated as having no useful homologs.
MSA cost is O(N_targets + N_binders), not O(N x M).
*/

include { ALPHAFOLD2_JACKHMMER_MSA as JACKHMMER_TARGET } from '../../modules/fold/af2/alphafold2_jackhmmer_msa'
include { ALPHAFOLD2_JACKHMMER_MSA as JACKHMMER_BINDER } from '../../modules/fold/af2/alphafold2_jackhmmer_msa'
include { AF2_MSAS_TO_A3M as AF2_MSAS_TO_A3M_TARGET } from '../../modules/fold/af2/af2_msas_to_a3m'
include { AF2_MSAS_TO_A3M as AF2_MSAS_TO_A3M_BINDER } from '../../modules/fold/af2/af2_msas_to_a3m'
include { MMSEQS_COLABFOLDSEARCH as MMSEQS_TARGET } from '../../modules/local/common/mmseqs_colabfoldsearch'
include { MMSEQS_COLABFOLDSEARCH as MMSEQS_BINDER } from '../../modules/local/common/mmseqs_colabfoldsearch'
include { ANNOTATE_MSA as ANNOTATE_MSA_TARGET } from '../../modules/fold/common/annotate_msa'
include { ANNOTATE_MSA as ANNOTATE_MSA_BINDER } from '../../modules/fold/common/annotate_msa'
include { SINGLE_SEQ_A3M as SINGLE_SEQ_A3M_TARGET } from '../../modules/fold/common/single_seq_a3m'
include { SINGLE_SEQ_A3M as SINGLE_SEQ_A3M_BINDER } from '../../modules/fold/common/single_seq_a3m'
include { FOLD_ASSEMBLE_AF2_MULTIMER_MSAS } from '../../modules/fold/af2/fold_assemble_af2_multimer_msas'

def sanitize(name) {
    return name.toString().replaceAll(/[^a-zA-Z0-9_.-]/, "_")
}

def pickPrimaryA3m(meta, a3m) {
    def files = (a3m instanceof List) ? a3m : [a3m]
    return files.find { it.name == "${meta.id}.a3m" } \
        ?: files.find { it.toString().contains('result') } \
        ?: files[0]
}

workflow FOLD_PULLDOWN_MSA {
    take:
    ch_targets // tuple(meta, fasta)  meta: [id, seq]
    ch_binders // tuple(meta, fasta)
    methods    // List<String>
    msa_method // 'jackhmmer_af2' | 'mmseqs2_colabfold'

    main:
    // af2_mono consumes exactly the same assembled AF2 input as af2 - it differs only in
    // how features.pkl is built inside the task - so it must switch this on too. Gating on
    // 'af2' alone made `--methods af2_mono` (without af2) skip FOLD_ASSEMBLE_AF, leaving
    // ALPHAFOLD2_MONO with an empty channel: zero predictions, exit 0, no warning.
    def need_af2 = ('af2' in methods) || ('af2_mono' in methods)
    def need_annotate = ('boltz' in methods) || ('rf3' in methods) || ('protenix' in methods)
    def pub = "${params.fold_publish_dir ?: 'fold'}/msa"
    def empty_msa = file("${projectDir}/assets/dummy_files/empty")

    // ---------------- target a3m (+ optional AF2 msas dir) ----------------
    ch_target_a3m = Channel.empty()
    ch_target_af2_msas = Channel.empty() // tuple(target_id, msas_dir)

    if (params.create_target_msa && msa_method == 'jackhmmer_af2') {
        ch_tin = ch_targets.map { meta, fasta ->
            [meta + [n_chains: 1, af2_force_monomer_msa: true], fasta]
        }
        JACKHMMER_TARGET(ch_tin)
        ch_target_af2_msas = JACKHMMER_TARGET.out.msa
            .map { meta, _fasta, msas -> [meta.id.toString(), msas] }
        AF2_MSAS_TO_A3M_TARGET(JACKHMMER_TARGET.out.msa)
        ch_target_a3m = AF2_MSAS_TO_A3M_TARGET.out.a3m
    }
    else if (params.create_target_msa && msa_method == 'mmseqs2_colabfold') {
        def envdb = params.use_remote_server ? file("${projectDir}/assets/dummy_files/empty") : file(params.colabfold_envdb)
        def uniref30_db = params.use_remote_server ? file("${projectDir}/assets/dummy_files/empty") : file(params.uniref30)
        ch_tin = ch_targets.map { meta, fasta -> [meta + [n_chains: 1], fasta] }
        MMSEQS_TARGET(ch_tin, params.use_remote_server, envdb, uniref30_db, "${pub}/mmseqs2_colabfold")
        ch_target_a3m = ch_tin.join(MMSEQS_TARGET.out.a3m).map { meta, fasta, a3m ->
            [meta, fasta, pickPrimaryA3m(meta, a3m)]
        }
    }
    else {
        SINGLE_SEQ_A3M_TARGET(ch_targets.map { meta, fasta -> [meta + [n_chains: 1], fasta] })
        ch_target_a3m = SINGLE_SEQ_A3M_TARGET.out.a3m
    }

    // ---------------- binder a3m ----------------
    ch_binder_a3m = Channel.empty()

    if (params.create_binder_msa && msa_method == 'jackhmmer_af2') {
        ch_bin = ch_binders.map { meta, fasta ->
            [meta + [n_chains: 1, af2_force_monomer_msa: true], fasta]
        }
        JACKHMMER_BINDER(ch_bin)
        AF2_MSAS_TO_A3M_BINDER(JACKHMMER_BINDER.out.msa)
        ch_binder_a3m = AF2_MSAS_TO_A3M_BINDER.out.a3m
    }
    else if (params.create_binder_msa && msa_method == 'mmseqs2_colabfold') {
        def envdb2 = params.use_remote_server ? file("${projectDir}/assets/dummy_files/empty") : file(params.colabfold_envdb)
        def uniref30_db2 = params.use_remote_server ? file("${projectDir}/assets/dummy_files/empty") : file(params.uniref30)
        ch_bin = ch_binders.map { meta, fasta -> [meta + [n_chains: 1], fasta] }
        MMSEQS_BINDER(ch_bin, params.use_remote_server, envdb2, uniref30_db2, "${pub}/mmseqs2_colabfold")
        ch_binder_a3m = ch_bin.join(MMSEQS_BINDER.out.a3m).map { meta, fasta, a3m ->
            [meta, fasta, pickPrimaryA3m(meta, a3m)]
        }
    }
    else {
        SINGLE_SEQ_A3M_BINDER(ch_binders.map { meta, fasta -> [meta + [n_chains: 1], fasta] })
        ch_binder_a3m = SINGLE_SEQ_A3M_BINDER.out.a3m
    }

    // ---------------- pair meta + FASTA (always) ----------------
    ch_pair_meta = ch_targets.combine(ch_binders)
        .map { tmeta, _tfasta, bmeta, _bfasta ->
            def tid = sanitize(tmeta.id)
            def bid = sanitize(bmeta.id)
            def pair_id = "${tid}_and_${bid}"
            def pair_meta = [
                id: pair_id,
                target: tid,
                binder: bid,
                n_chains: 2,
                target_seq: tmeta.seq,
                binder_seq: bmeta.seq,
            ]
            [pair_meta, tmeta.id.toString(), bmeta.id.toString()]
        }

    ch_pair_fasta = ch_pair_meta
        .map { pmeta, tid, bid ->
            def content = ">${pmeta.target}\n${pmeta.target_seq}\n>${pmeta.binder}\n${pmeta.binder_seq}\n"
            [pmeta.id, content]
        }
        .collectFile { id, content ->
            ["${id}.fasta", content]
        }
        .map { f -> [f.baseName, f] }

    ch_pairs_base = ch_pair_meta
        .map { pmeta, tid, bid -> [pmeta.id, pmeta, tid, bid] }
        .join(ch_pair_fasta)
        .map { _id, pmeta, tid, bid, fasta -> [pmeta, fasta, tid, bid] }

    ch_pairs_tsv = ch_pairs_base.map { pmeta, _fasta, _tid, _bid ->
        "${pmeta.id}\t${pmeta.target}\t${pmeta.binder}\n"
    }

    // ---------------- annotate + engine channels ----------------
    ch_for_boltz = Channel.empty()
    ch_for_rf3 = Channel.empty()
    ch_for_protenix = Channel.empty()

    if (need_annotate) {
        ANNOTATE_MSA_TARGET(
            ch_target_a3m.map { meta, fasta, a3m ->
                [meta + [chain_id: 'A', chain_index: 0], fasta, a3m]
            }
        )
        ANNOTATE_MSA_BINDER(
            ch_binder_a3m.map { meta, fasta, a3m ->
                [meta + [chain_id: 'B', chain_index: 1], fasta, a3m]
            }
        )

        ch_t_rend = ANNOTATE_MSA_TARGET.out.rendered
            .map { meta, rf3, pp, pu, bc -> [meta.id.toString(), rf3, pp, pu, bc] }
        ch_b_rend = ANNOTATE_MSA_BINDER.out.rendered
            .map { meta, rf3, pp, pu, bc -> [meta.id.toString(), rf3, pp, pu, bc] }

        // combine(by: 0), NOT join(): join() matches each key exactly once and drops
        // the rest, but this is inherently one-to-many - N_targets x N_binders pairs
        // share only N_targets rendered target MSAs and N_binders binder MSAs. With
        // join(), 27 pairs collapsed to 3 (one per target) and then to 1 (one per
        // binder), so the pulldown silently folded a SINGLE complex regardless of
        // input size. combine(by: 0) emits every pair whose key matches.
        ch_engine_pairs = ch_pairs_base
            .map { pmeta, fasta, tid, bid -> [tid, pmeta, fasta, bid] }
            .combine(ch_t_rend, by: 0)
            .map { tid, pmeta, fasta, bid, trf3, tpp, tpu, tbc ->
                [bid, pmeta, fasta, trf3, tpp, tpu, tbc]
            }
            .combine(ch_b_rend, by: 0)
            .map { bid, pmeta, fasta, trf3, tpp, tpu, tbc, brf3, bpp, bpu, bbc ->
                [pmeta, fasta, trf3, tpp, tpu, tbc, brf3, bpp, bpu, bbc]
            }

        ch_for_boltz = ch_engine_pairs.map { pmeta, fasta, _trf3, _tpp, _tpu, tbc, _brf3, _bpp, _bpu, bbc ->
            [pmeta, fasta, [tbc, bbc]]
        }
        ch_for_rf3 = ch_engine_pairs.map { pmeta, fasta, trf3, _tpp, _tpu, _tbc, brf3, _bpp, _bpu, _bbc ->
            [pmeta, fasta, [trf3, brf3]]
        }
        ch_for_protenix = ch_engine_pairs.map { pmeta, fasta, _trf3, tpp, tpu, _tbc, _brf3, bpp, bpu, _bbc ->
            [pmeta, fasta, [tpp, bpp, tpu, bpu]]
        }
    }

    // ---------------- AF2 assemble ----------------
    ch_af2_in = Channel.empty()
    if (need_af2) {
        // Always stage both optional inputs; unused slot is the empty stub.
        // Prefer jackhmmer msas dir when present; else ColabFold/mmseqs2 a3m.
        ch_target_a3m_by_id = ch_target_a3m.map { meta, _fasta, a3m ->
            [meta.id.toString(), a3m]
        }

        if (params.create_target_msa && msa_method == 'jackhmmer_af2') {
            // combine(by: 0) not join() - see ch_engine_pairs above; every binder
            // paired with a given target must reuse that target's MSA.
            ch_assemble_in = ch_pairs_base
                .map { pmeta, fasta, tid, _bid -> [tid, pmeta, fasta] }
                .combine(ch_target_af2_msas, by: 0)
                .map { _tid, pmeta, fasta, msas ->
                    [pmeta, fasta, msas, empty_msa]
                }
        }
        else if (params.create_target_msa && msa_method == 'mmseqs2_colabfold') {
            ch_assemble_in = ch_pairs_base
                .map { pmeta, fasta, tid, _bid -> [tid, pmeta, fasta] }
                .combine(ch_target_a3m_by_id, by: 0)
                .map { _tid, pmeta, fasta, a3m ->
                    [pmeta, fasta, empty_msa, a3m]
                }
        }
        else {
            ch_assemble_in = ch_pairs_base.map { pmeta, fasta, _tid, _bid ->
                [pmeta, fasta, empty_msa, empty_msa]
            }
        }
        FOLD_ASSEMBLE_AF2_MULTIMER_MSAS(ch_assemble_in)
        def a3m_stub = file("${projectDir}/assets/dummy_files/empty")
        ch_af2_in = FOLD_ASSEMBLE_AF2_MULTIMER_MSAS.out.msas
            .map { meta, fasta, msas -> [meta, fasta, msas, a3m_stub] }
    }

    emit:
    af2_in = ch_af2_in
    for_boltz = ch_for_boltz
    for_rf3 = ch_for_rf3
    for_protenix = ch_for_protenix
    pairs_rows = ch_pairs_tsv
}
