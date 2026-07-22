/*
FOLD_MSA centralises MSA generation for fold.nf and adapts each engine's
required format from a single --msa_method (see plans/
fold-nf-multi-method-folding.md §3). It only runs the stages the requested
--methods actually need.

Phase 1 scope: monomer only. Paired multimer MSAs (plan §4) are deliberately
not implemented here yet.
*/

include { ALPHAFOLD2_JACKHMMER_MSA } from '../../modules/fold/af2/alphafold2_jackhmmer_msa'
include { MMSEQS_COLABFOLDSEARCH } from '../../modules/local/common/mmseqs_colabfoldsearch'
include { COLABFOLD_A3M_TO_AF2_MSAS } from '../../modules/fold/af2/colabfold_a3m_to_af2_msas'
include { AF2_MSAS_TO_A3M } from '../../modules/fold/af2/af2_msas_to_a3m'

workflow FOLD_MSA {
    take:
    ch_input   // tuple(meta, fasta)
    methods    // List<String>, subset of ['af2', 'boltz', 'rf3', 'protenix']
    msa_method // 'jackhmmer_af2' | 'mmseqs2_colabfold'

    main:
    def need_af2_msas = 'af2' in methods
    def need_a3m = ('boltz' in methods) || ('rf3' in methods) || ('protenix' in methods)

    ch_af2_msas = Channel.empty()
    ch_a3m = Channel.empty()

    if (msa_method == 'jackhmmer_af2') {
        // Native for AF2 - always run if af2 is requested. If af2 is NOT
        // requested but boltz/rf3/protenix are, we still need to run it once
        // as the only route to an a3m under this --msa_method (there is no
        // standalone jackhmmer-only a3m stage).
        ALPHAFOLD2_JACKHMMER_MSA(ch_input)
        ch_af2_msas = ALPHAFOLD2_JACKHMMER_MSA.out.msa // tuple(meta, fasta, msas_dir)

        if (need_a3m) {
            AF2_MSAS_TO_A3M(ch_af2_msas)
            ch_a3m = AF2_MSAS_TO_A3M.out.a3m // tuple(meta, fasta, a3m)
        }
    }
    else if (msa_method == 'mmseqs2_colabfold') {
        // fold.nf validates --colabfold_envdb/--uniref30 are set before this
        // subworkflow runs whenever --use_remote_server is false (see
        // fold.nf's help/validation block) - no M3 default local ColabFold DB
        // path was found during this feature's investigation (searched
        // /mnt/datasets and /mnt/reference; only raw hh-suite-format uniref30
        // exists, not colabfold_search's expected mmseqs2 DB layout).
        // MMSEQS_COLABFOLDSEARCH declares the two DBs as `path` inputs, so they
        // must be real path values even in remote-server mode (where the script
        // ignores them). Pass the repo's empty dummy file when --use_remote_server
        // is set - mirrors workflows/rfd3.nf. Passing `false` (a Boolean) instead
        // fails with "Not a valid path value type: java.lang.Boolean".
        def envdb = params.use_remote_server ? file("${projectDir}/assets/dummy_files/empty") : file(params.colabfold_envdb)
        def uniref30_db = params.use_remote_server ? file("${projectDir}/assets/dummy_files/empty") : file(params.uniref30)
        MMSEQS_COLABFOLDSEARCH(
            ch_input,
            params.use_remote_server,
            envdb,
            uniref30_db,
            'fold/msa/mmseqs2_colabfold',
        )
        // MMSEQS_COLABFOLDSEARCH re-emits meta but drops fasta; rejoin it so
        // downstream converters/consumers keep the harmonized (meta, fasta,
        // a3m) contract. Its output glob (**.a3m) can emit MULTIPLE files -
        // notably the remote path writes uniref.a3m + bfd.mgnify30...a3m plus
        // the MERGED per-target ${meta.id}.a3m - so select the single merged
        // a3m (query-first, all hits) rather than passing the whole list (which
        // downstream stringifies into one broken space-joined path). Mirrors the
        // primary-a3m pick in workflows/rfd3.nf.
        ch_a3m = ch_input.join(MMSEQS_COLABFOLDSEARCH.out.a3m).map { meta, fasta, a3m ->
            def files = (a3m instanceof List) ? a3m : [a3m]
            def primary = files.find { it.name == "${meta.id}.a3m" } ?: files.find { it.toString().contains('result') } ?: files[0]
            [meta, fasta, primary]
        } // tuple(meta, fasta, a3m)

        if (need_af2_msas) {
            COLABFOLD_A3M_TO_AF2_MSAS(ch_a3m)
            ch_af2_msas = COLABFOLD_A3M_TO_AF2_MSAS.out.msas // tuple(meta, fasta, msas_dir)
        }
    }
    else {
        error("FOLD_MSA: unknown msa_method '${msa_method}'")
    }

    emit:
    af2_msas = ch_af2_msas  // tuple(meta, fasta, msas_dir) - empty channel if af2 wasn't requested
    for_boltz = ch_a3m      // tuple(meta, fasta, a3m) - empty channel if none of boltz/rf3/protenix requested
    for_rf3 = ch_a3m
    for_protenix = ch_a3m
}
