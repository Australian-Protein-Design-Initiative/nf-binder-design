/*
FOLD_MSA centralises MSA generation for fold.nf and adapts each engine's
required format from a single --msa_method (see plans/
fold-nf-multi-method-folding.md §3). It only runs the stages the requested
--methods actually need.

Monomer path (meta.n_chains == 1): unchanged from Phase 1 - one unpaired a3m
shared by Boltz/RF3/Protenix, plus AF2's native msas dir.

Multimer path (meta.n_chains > 1, plans/fold-nf-multimer-paired-msa.md §4):
  - Boltz/RF3/Protenix: split the complex into per-chain FASTAs, run the
    per-chain MSA search, then ANNOTATE_MSA renders each engine's native paired
    format via bin/fold/msa_taxonomy.py. The per-chain rendered files are grouped
    back per complex (chain order) into per-tool bundles.
  - AF2: fed the WHOLE complex to its native multimer MSA pipeline (jackhmmer +
    internal species pairing against the 2021 uniprot DB); no bespoke pairing.
    That run already searches every chain against the same DBs, so when AF2 is
    selected under jackhmmer_af2 the per-chain a3ms are taken from its
    msas/<chain>/ dirs rather than searching each chain a second time.

The monomer path is kept byte-for-byte identical so -resume caches unchanged;
the multimer processes sit on separate (aliased) invocations and are inert on
monomer-only runs (their input channels are empty).
*/

include { ALPHAFOLD2_JACKHMMER_MSA } from '../../modules/local/fold/af2/alphafold2_jackhmmer_msa'
include { ALPHAFOLD2_JACKHMMER_MSA as JACKHMMER_MSA_PERCHAIN } from '../../modules/local/fold/af2/alphafold2_jackhmmer_msa'
include { ALPHAFOLD2_JACKHMMER_MSA as JACKHMMER_MSA_COMPLEX } from '../../modules/local/fold/af2/alphafold2_jackhmmer_msa'
include { MMSEQS_COLABFOLDSEARCH } from '../../modules/local/common/mmseqs_colabfoldsearch'
include { MMSEQS_COLABFOLDSEARCH as MMSEQS_COLABFOLDSEARCH_PERCHAIN } from '../../modules/local/common/mmseqs_colabfoldsearch'
include { COLABFOLD_A3M_TO_AF2_MSAS } from '../../modules/local/fold/af2/colabfold_a3m_to_af2_msas'
include { AF2_MSAS_TO_A3M } from '../../modules/local/fold/af2/af2_msas_to_a3m'
include { AF2_MSAS_TO_A3M as AF2_MSAS_TO_A3M_PERCHAIN } from '../../modules/local/fold/af2/af2_msas_to_a3m'
include { SPLIT_COMPLEX_FASTA } from '../../modules/local/fold/common/split_complex_fasta'
include { ANNOTATE_MSA } from '../../modules/local/fold/common/annotate_msa'
include { AF2_STAGE_COLABFOLD_MULTIMER_MSAS } from '../../modules/local/fold/af2/af2_stage_colabfold_multimer_msas'

workflow FOLD_MSA {
    take:
    ch_input   // tuple(meta, fasta)
    methods    // List<String>, subset of ['af2', 'af2_mono', 'boltz', 'rf3', 'protenix', 'af3', 'openfold3', 'esmfold2', 'esmfold2_fast']
    msa_method // 'jackhmmer_af2' | 'mmseqs2_colabfold'

    main:
    // af2_mono consumes the same per-chain AF2 msas dir as af2; only the in-task
    // features.pkl assembly differs. See FOLD_PULLDOWN_MSA for the same gate.
    def need_af2_msas = ('af2' in methods) || ('af2_mono' in methods)
    // a3m needed for Boltz/RF3/Protenix/AF3/OpenFold3/ESMFold2, and for AF2 when
    // --msa_subsample is on (shallow jobs rebuild features.pkl from a subsampled a3m).
    // ESMFold2 wants none when --esmfold2_single_sequence folds it from sequence alone.
    def esmfold2_no_msa = ('esmfold2' in methods) && params.esmfold2_single_sequence
    def need_esmfold2_msa = ('esmfold2' in methods) && !params.esmfold2_single_sequence
    def need_a3m = ('boltz' in methods) || ('rf3' in methods) || ('protenix' in methods) || ('af3' in methods) \
        || ('openfold3' in methods) || need_esmfold2_msa \
        || (need_af2_msas && MsaSubsample.isEnabled(params.msa_subsample))
    // Boltz/RF3/Protenix/AF3/OpenFold3/ESMFold2 need per-chain paired MSAs on the multimer path.
    def need_paired = ('boltz' in methods) || ('rf3' in methods) || ('protenix' in methods) || ('af3' in methods) \
        || ('openfold3' in methods) || need_esmfold2_msa
    // af2_mono (monomer-weights chain-break) on a ColabFold-searched multimer input
    // needs one plain a3m per chain too, same split as need_paired, but AF2's own
    // multimer pairing ('af2') under mmseqs2_colabfold is rejected up front by
    // FoldValidation - only af2_mono reaches this branch.
    def need_af2_colabfold_multi = need_af2_msas && msa_method == 'mmseqs2_colabfold'
    def need_chain_split = need_paired || need_af2_colabfold_multi
    def af2_complex_search = need_af2_msas && msa_method == 'jackhmmer_af2'

    ch_mono = ch_input.filter { meta, fasta -> (meta.n_chains ?: 1) == 1 }
    ch_multi = ch_input.filter { meta, fasta -> (meta.n_chains ?: 1) > 1 }

    ch_af2_msas_mono = Channel.empty()
    ch_a3m_mono = Channel.empty()

    // ==================== MONOMER (unchanged Phase-1 path) ====================
    if (msa_method == 'jackhmmer_af2') {
        // Both gates are false only when nothing consumes an MSA at all
        // (single-sequence ESMFold2 as the sole engine), so skip the search rather
        // than run it for nobody. The msa_method branch itself must stay, or an
        // otherwise-valid method falls through to the unknown-method error below.
        if (need_af2_msas || need_a3m) {
            ALPHAFOLD2_JACKHMMER_MSA(ch_mono)
            ch_af2_msas_mono = ALPHAFOLD2_JACKHMMER_MSA.out.msa // tuple(meta, fasta, msas_dir)

            if (need_a3m) {
                AF2_MSAS_TO_A3M(ch_af2_msas_mono)
                ch_a3m_mono = AF2_MSAS_TO_A3M.out.a3m // tuple(meta, fasta, a3m)
            }
        }
    }
    else if (msa_method == 'mmseqs2_colabfold') {
        // See the historical comment block below for why the DBs are passed as
        // dummy files in remote-server mode.
        def envdb = params.use_remote_server ? file("${projectDir}/assets/dummy_files/empty") : file(params.colabfold_envdb)
        def uniref30_db = params.use_remote_server ? file("${projectDir}/assets/dummy_files/empty") : file(params.uniref30)
        MMSEQS_COLABFOLDSEARCH(
            ch_mono,
            params.use_remote_server,
            envdb,
            uniref30_db,
            "${params.fold_publish_dir ?: 'fold'}/msa/mmseqs2_colabfold",
        )
        ch_a3m_mono = ch_mono.join(MMSEQS_COLABFOLDSEARCH.out.a3m).map { meta, fasta, a3m ->
            def files = (a3m instanceof List) ? a3m : [a3m]
            def primary = files.find { it.name == "${meta.id}.a3m" } ?: files.find { it.toString().contains('result') } ?: files[0]
            [meta, fasta, primary]
        } // tuple(meta, fasta, a3m)

        if (need_af2_msas) {
            COLABFOLD_A3M_TO_AF2_MSAS(ch_a3m_mono)
            ch_af2_msas_mono = COLABFOLD_A3M_TO_AF2_MSAS.out.msas // tuple(meta, fasta, msas_dir)
        }
    }
    else {
        error("FOLD_MSA: unknown msa_method '${msa_method}'")
    }

    // ==================== MULTIMER (paired-MSA path) ====================
    ch_rf3_multi = Channel.empty()
    ch_protenix_multi = Channel.empty()
    ch_boltz_multi = Channel.empty()
    ch_af2_msas_multi = Channel.empty()

    // AF2 multimer uses its own native multimer MSA pipeline on the whole
    // complex (jackhmmer + internal pairing). fold.nf guarantees af2 multimer
    // only runs under --msa_method jackhmmer_af2 against a uniprot/-bearing DB.
    if (af2_complex_search) {
        JACKHMMER_MSA_COMPLEX(ch_multi)
        ch_af2_msas_multi = JACKHMMER_MSA_COMPLEX.out.msa
    }

    if (need_chain_split) {
        // 1. Split each complex into per-chain FASTAs, one search unit each.
        SPLIT_COMPLEX_FASTA(ch_multi)
        ch_chain = SPLIT_COMPLEX_FASTA.out.chains.flatMap { meta, files ->
            def fs = ((files instanceof List) ? files : [files]).sort { it.name }
            fs.withIndex().collect { f, i ->
                def chain_letter = ((('A' as char) as int) + i) as char
                def chain_meta = meta + [
                    id: f.baseName,
                    base_id: meta.id,
                    chain_index: i,
                    chain_id: "${chain_letter}",
                    total_chains: meta.n_chains,
                    n_chains: 1,
                ]
                [chain_meta, f]
            }
        }

        // 2. Per-chain MSA search (same route as the monomer path, one query each).
        ch_chain_a3m = Channel.empty()
        if (af2_complex_search) {
            // Reuse the complex run's per-chain MSAs. AF2 writes msas/<chain>/ only
            // for the first chain with each sequence, so repeated chains (homomers)
            // point at that chain's dir (chain_id_map.json maps chain -> sequence).
            ch_chain_from_complex = ch_chain.map { cm, f -> [cm.base_id, cm, f] }
                .combine(JACKHMMER_MSA_COMPLEX.out.msa.map { m, _fa, dir -> [m.id, dir] }, by: 0)
                .map { _id, cm, f, dir ->
                    def chain_map = new groovy.json.JsonSlurper().parseText(dir.resolve('msas/chain_id_map.json').text)
                    def seq = chain_map[cm.chain_id]?.sequence
                    def msa_chain = chain_map.keySet().sort().find { chain_map[it].sequence == seq } ?: cm.chain_id
                    [cm + [af2_msa_chain: msa_chain], f, dir]
                }
            AF2_MSAS_TO_A3M_PERCHAIN(ch_chain_from_complex)
            ch_chain_a3m = AF2_MSAS_TO_A3M_PERCHAIN.out.a3m.map { cm, f, a3m ->
                def m = cm.findAll { k, _v -> k != 'af2_msa_chain' }
                [m, f, a3m]
            } // tuple(chain_meta, chain_fasta, a3m)
        }
        else if (msa_method == 'jackhmmer_af2') {
            JACKHMMER_MSA_PERCHAIN(ch_chain)
            AF2_MSAS_TO_A3M_PERCHAIN(JACKHMMER_MSA_PERCHAIN.out.msa)
            ch_chain_a3m = AF2_MSAS_TO_A3M_PERCHAIN.out.a3m // tuple(chain_meta, chain_fasta, a3m)
        }
        else if (msa_method == 'mmseqs2_colabfold') {
            def envdb2 = params.use_remote_server ? file("${projectDir}/assets/dummy_files/empty") : file(params.colabfold_envdb)
            def uniref30_db2 = params.use_remote_server ? file("${projectDir}/assets/dummy_files/empty") : file(params.uniref30)
            MMSEQS_COLABFOLDSEARCH_PERCHAIN(
                ch_chain,
                params.use_remote_server,
                envdb2,
                uniref30_db2,
                "${params.fold_publish_dir ?: 'fold'}/msa/mmseqs2_colabfold",
            )
            ch_chain_a3m = ch_chain.join(MMSEQS_COLABFOLDSEARCH_PERCHAIN.out.a3m).map { meta, fasta, a3m ->
                def files = (a3m instanceof List) ? a3m : [a3m]
                def primary = files.find { it.name == "${meta.id}.a3m" } ?: files.find { it.toString().contains('result') } ?: files[0]
                [meta, fasta, primary]
            }
        }

        if (need_paired) {
            // 3. Render each chain into every engine's native paired format.
            ANNOTATE_MSA(ch_chain_a3m)

            // 4. Group per complex (chain order) into per-tool bundles. A plain
            //    groupTuple(by: 0) buffers until the WHOLE channel completes (it
            //    has no way to know a group is done), which serialises every
            //    multimer engine behind the slowest chain of the slowest complex.
            //    groupKey(base_id, total_chains) tells it the exact group size, so
            //    each complex's bundle emits as soon as its own chains land.
            //    Reorder by chain_index since group order is arrival order, not
            //    chain order.
            // The grouping key (element 0) is a GroupKey; base_id is carried alongside
            // it (plain String) so the downstream .join() below matches on a plain
            // value rather than depending on GroupKey's equality semantics.
            ch_grouped = ANNOTATE_MSA.out.rendered
                .map { cm, rf3, pp, pu, bc -> [groupKey(cm.base_id, cm.total_chains), cm.base_id, cm.chain_index, rf3, pp, pu, bc] }
                .groupTuple()
                .map { _gkey, base_ids, idxs, rf3s, pps, pus, bcs ->
                    def base_id = base_ids[0]
                    def order = (0..<idxs.size()).toList().sort { idxs[it] }
                    [
                        base_id,
                        order.collect { rf3s[it] },
                        order.collect { pps[it] },
                        order.collect { pus[it] },
                        order.collect { bcs[it] },
                    ]
                }

            // 5. Rejoin the complex fasta and split into per-tool channels.
            ch_bundle = ch_multi.map { meta, fasta -> [meta.id, meta, fasta] }
                .join(ch_grouped)
                .map { id, meta, fasta, rf3o, ppo, puo, bco -> [meta, fasta, rf3o, ppo, puo, bco] }

            ch_rf3_multi = ch_bundle.map { meta, fasta, rf3o, ppo, puo, bco -> [meta, fasta, rf3o] }
            // Protenix takes a combined paired+unpaired list; the COMPLEX generator
            // splits it by filename suffix.
            ch_protenix_multi = ch_bundle.map { meta, fasta, rf3o, ppo, puo, bco -> [meta, fasta, ppo + puo] }
            ch_boltz_multi = ch_bundle.map { meta, fasta, rf3o, ppo, puo, bco -> [meta, fasta, bco] }
        }

        // 6. af2_mono under mmseqs2_colabfold: stage the SAME per-chain a3ms into
        // an AF2 multimer-format msas/<CHAIN>/ tree (one plain a3m per chain, no
        // pairing) so the chain-break trick has MSAs to read. AF2's own multimer
        // pairing ('af2') never reaches this branch under mmseqs2_colabfold -
        // FoldValidation rejects that combination up front.
        if (need_af2_colabfold_multi) {
            ch_chain_a3m_grouped = ch_chain_a3m
                .map { cm, _f, a3m -> [groupKey(cm.base_id, cm.total_chains), cm.base_id, cm.chain_index, cm.chain_id, a3m] }
                .groupTuple()
                .map { _gkey, base_ids, idxs, chain_ids, a3ms ->
                    def order = (0..<idxs.size()).toList().sort { idxs[it] }
                    [base_ids[0], order.collect { chain_ids[it] }, order.collect { a3ms[it] }]
                }

            ch_af2_colabfold_multi_input = ch_multi.map { meta, fasta -> [meta.id, meta, fasta] }
                .join(ch_chain_a3m_grouped)
                .map { _id, meta, fasta, chain_ids, a3ms -> [meta, fasta, chain_ids, a3ms] }

            AF2_STAGE_COLABFOLD_MULTIMER_MSAS(ch_af2_colabfold_multi_input)
            ch_af2_msas_multi = AF2_STAGE_COLABFOLD_MULTIMER_MSAS.out.msas
        }
    }


    emit:
    af2_msas = ch_af2_msas_mono.mix(ch_af2_msas_multi) // tuple(meta, fasta, msas_dir)
    a3m = ch_a3m_mono                                  // monomer only (multimer disables subsample)
    for_boltz = ch_a3m_mono.mix(ch_boltz_multi)        // monomer: (meta,fasta,a3m); multimer: (meta,fasta,[csv...])
    for_rf3 = ch_a3m_mono.mix(ch_rf3_multi)            // multimer: (meta,fasta,[rf3_a3m...])
    for_protenix = ch_a3m_mono.mix(ch_protenix_multi)  // multimer: (meta,fasta,[paired...+unpaired...])
    for_af3 = ch_a3m_mono.mix(ch_protenix_multi)       // same bundle; AF3 re-renders pairing from the unpaired a3m
    for_openfold3 = ch_a3m_mono.mix(ch_protenix_multi) // same bundle; OpenFold3 re-renders pairing likewise
    // Single-sequence mode skips the MSA stage entirely for this engine, so it
    // cannot draw on ch_a3m_mono (which is not even built when nothing else needs
    // an a3m) - it takes the raw input FASTA plus a placeholder instead.
    for_esmfold2 = esmfold2_no_msa \
        ? ch_input.map { meta, fasta -> [meta, fasta, file("${projectDir}/assets/dummy_files/empty")] } \
        : ch_a3m_mono.mix(ch_protenix_multi)
    // ESMFold2-Fast has no MSA encoder, so it never takes (or triggers) an MSA.
    for_esmfold2_fast = ch_input.map { meta, fasta -> [meta, fasta, file("${projectDir}/assets/dummy_files/empty")] }
}
