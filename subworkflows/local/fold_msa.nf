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
    format via bin/msa_taxonomy.py. The per-chain rendered files are grouped
    back per complex (chain order) into per-tool bundles.
  - AF2: fed the WHOLE complex to its native multimer MSA pipeline (jackhmmer +
    internal species pairing against the 2021 uniprot DB); no bespoke pairing.

The monomer path is kept byte-for-byte identical so -resume caches unchanged;
the multimer processes sit on separate (aliased) invocations and are inert on
monomer-only runs (their input channels are empty).
*/

include { ALPHAFOLD2_JACKHMMER_MSA } from '../../modules/fold/af2/alphafold2_jackhmmer_msa'
include { ALPHAFOLD2_JACKHMMER_MSA as JACKHMMER_MSA_PERCHAIN } from '../../modules/fold/af2/alphafold2_jackhmmer_msa'
include { ALPHAFOLD2_JACKHMMER_MSA as JACKHMMER_MSA_COMPLEX } from '../../modules/fold/af2/alphafold2_jackhmmer_msa'
include { MMSEQS_COLABFOLDSEARCH } from '../../modules/local/common/mmseqs_colabfoldsearch'
include { MMSEQS_COLABFOLDSEARCH as MMSEQS_COLABFOLDSEARCH_PERCHAIN } from '../../modules/local/common/mmseqs_colabfoldsearch'
include { COLABFOLD_A3M_TO_AF2_MSAS } from '../../modules/fold/af2/colabfold_a3m_to_af2_msas'
include { AF2_MSAS_TO_A3M } from '../../modules/fold/af2/af2_msas_to_a3m'
include { AF2_MSAS_TO_A3M as AF2_MSAS_TO_A3M_PERCHAIN } from '../../modules/fold/af2/af2_msas_to_a3m'
include { SPLIT_COMPLEX_FASTA } from '../../modules/fold/common/split_complex_fasta'
include { ANNOTATE_MSA } from '../../modules/fold/common/annotate_msa'

workflow FOLD_MSA {
    take:
    ch_input   // tuple(meta, fasta)
    methods    // List<String>, subset of ['af2', 'boltz', 'rf3', 'protenix']
    msa_method // 'jackhmmer_af2' | 'mmseqs2_colabfold'

    main:
    def need_af2_msas = 'af2' in methods
    // a3m needed for Boltz/RF3/Protenix, and for AF2 when --msa_subsample is on
    // (shallow jobs rebuild features.pkl from a subsampled a3m).
    def need_a3m = ('boltz' in methods) || ('rf3' in methods) || ('protenix' in methods) \
        || (need_af2_msas && MsaSubsample.isEnabled(params.msa_subsample))
    // Boltz/RF3/Protenix need per-chain paired MSAs on the multimer path.
    def need_paired = ('boltz' in methods) || ('rf3' in methods) || ('protenix' in methods)

    ch_mono = ch_input.filter { meta, fasta -> (meta.n_chains ?: 1) == 1 }
    ch_multi = ch_input.filter { meta, fasta -> (meta.n_chains ?: 1) > 1 }

    ch_af2_msas_mono = Channel.empty()
    ch_a3m_mono = Channel.empty()

    // ==================== MONOMER (unchanged Phase-1 path) ====================
    if (msa_method == 'jackhmmer_af2') {
        ALPHAFOLD2_JACKHMMER_MSA(ch_mono)
        ch_af2_msas_mono = ALPHAFOLD2_JACKHMMER_MSA.out.msa // tuple(meta, fasta, msas_dir)

        if (need_a3m) {
            AF2_MSAS_TO_A3M(ch_af2_msas_mono)
            ch_a3m_mono = AF2_MSAS_TO_A3M.out.a3m // tuple(meta, fasta, a3m)
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
            'fold/msa/mmseqs2_colabfold',
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

    if (need_paired) {
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
                    n_chains: 1,
                ]
                [chain_meta, f]
            }
        }

        // 2. Per-chain MSA search (same route as the monomer path, one query each).
        ch_chain_a3m = Channel.empty()
        if (msa_method == 'jackhmmer_af2') {
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
                'fold/msa/mmseqs2_colabfold',
            )
            ch_chain_a3m = ch_chain.join(MMSEQS_COLABFOLDSEARCH_PERCHAIN.out.a3m).map { meta, fasta, a3m ->
                def files = (a3m instanceof List) ? a3m : [a3m]
                def primary = files.find { it.name == "${meta.id}.a3m" } ?: files.find { it.toString().contains('result') } ?: files[0]
                [meta, fasta, primary]
            }
        }

        // 3. Render each chain into every engine's native paired format.
        ANNOTATE_MSA(ch_chain_a3m)

        // 4. Group per complex (chain order) into per-tool bundles. groupTuple
        //    buffers to channel completion; reorder by chain_index since group
        //    order is arrival order, not chain order.
        ch_grouped = ANNOTATE_MSA.out.rendered
            .map { cm, rf3, pp, pu, bc -> [cm.base_id, cm.chain_index, rf3, pp, pu, bc] }
            .groupTuple(by: 0)
            .map { base_id, idxs, rf3s, pps, pus, bcs ->
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

    // AF2 multimer uses its own native multimer MSA pipeline on the whole
    // complex (jackhmmer + internal pairing). fold.nf guarantees af2 multimer
    // only runs under --msa_method jackhmmer_af2 against a uniprot/-bearing DB.
    if (need_af2_msas && msa_method == 'jackhmmer_af2') {
        JACKHMMER_MSA_COMPLEX(ch_multi)
        ch_af2_msas_multi = JACKHMMER_MSA_COMPLEX.out.msa
    }

    emit:
    af2_msas = ch_af2_msas_mono.mix(ch_af2_msas_multi) // tuple(meta, fasta, msas_dir)
    a3m = ch_a3m_mono                                  // monomer only (multimer disables subsample)
    for_boltz = ch_a3m_mono.mix(ch_boltz_multi)        // monomer: (meta,fasta,a3m); multimer: (meta,fasta,[csv...])
    for_rf3 = ch_a3m_mono.mix(ch_rf3_multi)            // multimer: (meta,fasta,[rf3_a3m...])
    for_protenix = ch_a3m_mono.mix(ch_protenix_multi)  // multimer: (meta,fasta,[paired...+unpaired...])
}
