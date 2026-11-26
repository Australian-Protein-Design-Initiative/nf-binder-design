#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Fold pulldown: co-fold every binder against every target with one or more of
AF2 / AF2-monomer (af2_mono) / Boltz-2 / RF3 / Protenix / AlphaFold3 / OpenFold3,
then summarise interface scores.

Usage via main.nf:
  nextflow run main.nf --method fold_pulldown \
      --targets targets.fasta --binders binders.fasta \
      --methods boltz,rf3 --create_target_msa true
*/

// Every parameter this workflow accepts is declared in conf/fold_pulldown.config, which
// nextflow.config includes when --method selects this entrypoint. They cannot be
// declared here: Nextflow does not evaluate top-level `params.x = ...` in an
// included script, so the subworkflows and modules below would see null.

include { COMPARE_TO_REFERENCE } from '../subworkflows/local/compare_to_reference'
include { EXTRACT_COMPLEX_STRUCTURES } from '../modules/local/fold_pulldown/extract_complex_structures'
include { FOLD_PULLDOWN_CORE } from '../subworkflows/local/fold_pulldown_core'
include { FOLD_PULLDOWN_REPORTING } from '../modules/local/common/fold_pulldown_reporting'
include { foldIdsExtractFastaIds; foldIdsValidatePulldownIds } from '../modules/local/fold/common/fold_ids'
include { foldValidationParseMethods; foldValidationValidate } from '../modules/local/fold/common/fold_validation'
include { msaSubsampleIsEnabled } from '../modules/local/fold/common/msa_subsample'

workflow FOLD_PULLDOWN {

    main:

    def have_binders = params.binders || params.complex_structures
    // --targets is optional only with --complex_structures, which can read the
    // target back out of the complexes' own target chain.
    def have_targets = params.targets || params.complex_structures
    if (params.help || !have_targets || !have_binders) {
        log.info(
            """
        ==================================================================
        FOLD PULLDOWN PIPELINE
        ==================================================================

        Co-fold every binder against every target with one or more structure
        predictors (AF2, Boltz-2, RF3, Protenix, AlphaFold3, OpenFold3), then summarise interface
        scores (iptm, ipsae) across replicates and models.

        Required arguments:
            --targets             FASTA file of target sequences. Optional with
                                  --complex_structures, which can read the target back
                                  out of each complex's target chain.
            --binders             FASTA file of binder sequences (or --complex_structures)

        Supplied complexes (alternative to --binders):
            --complex_structures  Target-binder complexes (.pdb/.cif; dir or glob). The
                                  binder chain of each is folded against every target,
                                  and the complex is the reference the prediction is
                                  compared against: target-aligned binder RMSD, TM-score,
                                  DockQ and interface quality, in complex_rmsd.tsv.
            --binder_chains       Binder chain id(s) in those complexes, comma-separated.
                                  Several means several binder records per complex.
                                  Required with --complex_structures.
            --target_chains       Chain used as the superposition reference. A single
                                  chain id, or 'auto' to pick it per complex by matching
                                  against --targets. Default: the one non-binder chain.
                                  A multi-chain target is not supported: the pulldown
                                  folds two-chain pairs, so there is nothing to
                                  superimpose a second target chain onto.

        Optional arguments:
            --outdir              Output directory [default: ${params.outdir}]
            --methods             Comma-separated af2,af2_mono,boltz,rf3,protenix,af3,openfold3,chai1,esmfold2,esmfold2_fast [default: ${params.methods}]
                                   af2      = AlphaFold2 multimer (every pair here is a 2-chain complex,
                                              so af2 always uses AF2's native multimer weights/pipeline).
                                   af2_mono = AF2 MONOMER weights on a concatenated complex, chains
                                              separated only by a residue_index jump. Shares af2's
                                              MSAs; only features.pkl differs. Without an initial
                                              guess the monomer models often fail to dock at all and
                                              only ranking separates the good pose, so consider
                                              --af2_keep_models best. Not an independent engine:
                                              it shares weights lineage with af2.
            --af2_chain_break_offset  residue_index jump per chain break, must exceed AF2's
                                   relative-position clip of 32 [default: ${params.af2_chain_break_offset}]
            --msa_method          jackhmmer_af2|mmseqs2_colabfold [default: ${params.msa_method}]
            --create_target_msa   Build MSA for each target [default: ${params.create_target_msa}]
            --create_binder_msa   Build MSA for each binder [default: ${params.create_binder_msa}]
            --n_predictions       Structures per complex per method [default: ${params.n_predictions}]
            --use_msa_server      Boltz fetches its own MSA [default: ${params.use_msa_server}]
            --templates           Template structures (.pdb/.cif; dir or glob), matched to target chains
                                  by sequence alignment; not used by esmfold2 [default: ${params.templates}]
            --target_msa          External A3M for the target chain, used instead of
                                  searching for one. Its first record must be the target
                                  sequence. Mutually exclusive with --create_target_msa.
                                  [default: ${params.target_msa}]
            --binder_templates    Also match templates to binder chains [default: ${params.binder_templates}]
            --template_min_identity / --template_min_coverage / --template_min_aligned / --template_max_per_chain
                                  Matching thresholds [default: ${params.template_min_identity} / ${params.template_min_coverage} / ${params.template_min_aligned} / ${params.template_max_per_chain}]
            --boltz_template_force  Hold templated chains near the template (Boltz force) [default: true]
            --boltz_template_threshold  Boltz force threshold in Angstrom [default: ${params.boltz_template_threshold}]
            --skip_engens         Skip EnGens clustering [default: ${params.skip_engens}]

            Every method-specific flag from --method fold --help (--af2_*, --boltz_*,
            --rf3_*, --protenix_*, --af3_*, --openfold3_*, --chai1_*, --esmfold2_*) applies here too, e.g.
            --af2_keep_models (which of AF2's 5 models/run to keep toward
            --n_predictions) [default: ${params.af2_keep_models}]. 'all' keeps all 5
            models from AF2's single run, matching --n_predictions 5 at no extra cost;
            --method fold defaults to 'best'.

            Ranking (summary table):
            --consensus_metric    ipsae|iptm; metric averaged into consensus_z [default: ${params.consensus_metric}]
            --z_stat              max|mean; per-complex statistic over samples that is
                                   standardised [default: ${params.z_stat}]
            --z_scope             target|global; standardise within (target, tool) or over all
                                   targets together [default: ${params.z_scope}]
            --min_pool            flag z-score pools smaller than this as z_pool_small, since
                                   over k complexes |z| cannot exceed (k-1)/sqrt(k)
                                   [default: ${params.min_pool}]

            AF2 (--methods includes af2 or af2_mono) needs the 2021 DB snapshot with uniprot/:
            --af2_db_path         [default: ${params.af2_db_path}]
            --af2_uniref30_subpath / --af2_mgnify_subpath / --af2_uniprot_subpath /
            --af2_pdb_seqres_subpath / --af2_pdb70_subpath   Per-DB subpaths under
                                   --af2_db_path, as for --method fold (see --method fold --help)
            --af2_data_dir        Directory containing params/ (model weights); bundled in the
                                   container by default [default: ${params.af2_data_dir}]

            AlphaFold3 (--methods includes af3; weights are NOT bundled - see models/download_af3_weights.sh):
            --af3_model_dir       Directory holding exactly one af3.bin.zst / af3.bin
                                   [default: ${params.af3_model_dir}]
            --af3_batch_size / --af3_seeds / --af3_num_recycles / --af3_flash_attention /
            --af3_jax_cache_dir   As for --method fold (see --method fold --help)
            --af3_paired_msa / --af3_templates / --af3_run_data_pipeline / --af3_db_dir
                                  As for --method fold. -profile af3_germinal_parity sets
                                  the combination Germinal uses (needs --af3_db_dir).

            OpenFold3 (--methods includes openfold3; weights are bundled in the container):
            --openfold3_batch_size / --openfold3_seeds /
            --openfold3_kernel_cache_dir   As for --method fold (see --method fold --help)

            Chai-1 (--methods includes chai1; weights are bundled in the container):
            --chai1_batch_size / --chai1_seeds / --chai1_num_trunk_recycles /
            --chai1_num_diffn_timesteps    As for --method fold (see --method fold --help)

            ESMFold2 (--methods includes esmfold2 and/or esmfold2_fast; weights are bundled in the containers):
            --esmfold2_weights_dir / --esmfold2_batch_size /
            --esmfold2_seeds / --esmfold2_single_sequence / --esmfold2_num_loops /
            --esmfold2_num_sampling_steps / --esmfold2_msa_max_depth /
            --esmfold2_kernel_backend      As for --method fold (see --method fold --help)

            ColabFold MSA (--msa_method mmseqs2_colabfold):
            --use_remote_server   Query ColabFold API [default: ${params.use_remote_server}]
            --uniref30 / --colabfold_envdb   Local ColabFold DB paths

            --gpu_devices           GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
            --gpu_slots_per_device  Concurrent tasks allowed per GPU [default: ${params.gpu_slots_per_device}]
            --gpu_lock_timeout      Seconds a task waits for a free GPU [default: ${params.gpu_lock_timeout}]

        Example:
            nextflow run main.nf --method fold_pulldown \\
                --targets targets.fasta --binders binders.fasta \\
                --methods boltz,rf3,protenix --create_target_msa true \\
                --msa_method jackhmmer_af2 -profile slurm,m3

        """.stripIndent()
        )
        exit(1)
    }

    if (params.complex_structures && params.binders) {
        error(
            'fold_pulldown: --complex_structures and --binders are two sources for the same binder set; give one.'
        )
    }
    if (params.complex_structures && !params.binder_chains) {
        error(
            'fold_pulldown: --complex_structures needs --binder_chains to say which chain of each complex is the binder.'
        )
    }
    if (params.binder_chains && !params.complex_structures) {
        error('fold_pulldown: --binder_chains only applies with --complex_structures.')
    }
    if (params.target_chains && !params.complex_structures) {
        error('fold_pulldown: --target_chains only applies with --complex_structures.')
    }
    if (params.target_chains == 'auto' && !params.targets) {
        error(
            'fold_pulldown: --target_chains auto picks the target chain by matching it against the --targets records, so it needs --targets.'
        )
    }

    def methods = foldValidationParseMethods(params.methods)
    // AF2 in fold_pulldown assembles the target MSA itself (jackhmmer dir or
    // ColabFold/mmseqs2 a3m via FOLD_ASSEMBLE_AF2_MULTIMER_MSAS) and always
    // folds the binder chain query-only, so neither fold_validation's generic
    // "AF2 multimer requires jackhmmer_af2" error nor its ColabFold "will run
    // UNPAIRED" warning apply here - opts.pulldown: true skips both.
    def (errors, warnings) = foldValidationValidate(params, methods, [
        hasMultimer: true,
        af2DbPath: params.af2_db_path,
        pulldown: true,
    ])

    warnings.each { log.warn("fold_pulldown: ${it}") }
    if (errors) {
        error("fold_pulldown: ${errors.join('\nfold_pulldown: ')}")
    }

    if (msaSubsampleIsEnabled(params.msa_subsample)) {
        error("fold_pulldown: --msa_subsample is not supported (multimer pairs only).")
    }

    // With --complex_structures the binder FASTA is written by a task, so there is no
    // file to read ids out of here. EXTRACT_COMPLEX_STRUCTURES rejects duplicate binder
    // ids itself, which is the check this would have done.
    if (!params.complex_structures) {
        def target_ids = foldIdsExtractFastaIds(file(params.targets))
        def binder_ids = foldIdsExtractFastaIds(file(params.binders))
        def id_errors = foldIdsValidatePulldownIds(target_ids, binder_ids)
        if (id_errors) {
            error("fold_pulldown: ${id_errors.join('\nfold_pulldown: ')}")
        }
    }

    ch_complex_chains = Channel.empty()
    if (params.complex_structures) {
        EXTRACT_COMPLEX_STRUCTURES(
            Channel.fromPath(params.complex_structures, checkIfExists: true).collect(),
            params.binder_chains,
            params.target_chains ?: '',
            params.targets
                ? Channel.fromPath(params.targets, checkIfExists: true)
                : file("${projectDir}/assets/dummy_files/no_pulldown_targets.fasta"),
        )
        ch_binders_fasta = EXTRACT_COMPLEX_STRUCTURES.out.binders
        // Without --targets the targets are the complexes' own target chains, read back
        // out and deduplicated by sequence.
        ch_targets_fasta = params.targets
            ? Channel.fromPath(params.targets)
            : EXTRACT_COMPLEX_STRUCTURES.out.targets
        ch_complex_chains = EXTRACT_COMPLEX_STRUCTURES.out.chains
    }
    else {
        ch_binders_fasta = Channel.fromPath(params.binders)
        ch_targets_fasta = Channel.fromPath(params.targets)
    }

    ch_targets_meta = ch_targets_fasta
        .splitFasta(record: [id: true, seqString: true])
        .map { record -> [id: record.id, seq: record.seqString] }

    ch_binders_meta = ch_binders_fasta
        .splitFasta(record: [id: true, seqString: true])
        .map { record -> [id: record.id, seq: record.seqString] }

    ch_targets_fasta_paths = ch_targets_fasta.splitFasta(file: true)
    ch_binders_fasta_paths = ch_binders_fasta.splitFasta(file: true)

    ch_targets = ch_targets_meta.merge(ch_targets_fasta_paths)
    ch_binders = ch_binders_meta.merge(ch_binders_fasta_paths)

    // Binders are only templated with --binder_templates (see meta.template_chains).
    FOLD_PULLDOWN_CORE(
        ch_targets,
        ch_binders,
        methods,
        params.msa_method,
        params.binder_templates \
            ? ch_targets_fasta.mix(ch_binders_fasta)
            : ch_targets_fasta,
    )

    // --complex_structures: compare every predicted complex back to the complex its
    // binder came from. With several --targets each binder yields one comparison per
    // (target, engine, sample), all against that one reference.
    if (params.complex_structures) {
        // Keyed back to the input files by name rather than to the published copies:
        // publishDir is asynchronous, so the published path is not guaranteed to exist
        // when the comparison task is staged.
        ch_complex_files = Channel
            .fromPath(params.complex_structures, checkIfExists: true)
            .map { f -> tuple(f.name, f) }

        ch_references = ch_complex_chains
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.complex, row.binder_id, row.binder_chain, row.target_chain) }
            .combine(ch_complex_files, by: 0)
            .map { _name, binder_id, binder_chain, target_chain, structure ->
                tuple(binder_id, structure, binder_chain, target_chain)
            }

        COMPARE_TO_REFERENCE(
            FOLD_PULLDOWN_CORE.out.predictions,
            ch_references,
            params.fold_publish_dir ?: 'fold_pulldown',
            'complex_rmsd.tsv',
            true,
        )
    }

    FOLD_PULLDOWN_REPORTING(
        file("${projectDir}/assets/fold_pulldown_reporting.qmd"),
        FOLD_PULLDOWN_CORE.out.scores,
        FOLD_PULLDOWN_CORE.out.summary,
    )

    emit:
    scores = FOLD_PULLDOWN_CORE.out.scores
    summary = FOLD_PULLDOWN_CORE.out.summary
}
