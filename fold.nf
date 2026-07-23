#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Standalone multi-method structure folding dispatcher: predicts structures for
FASTA inputs with any combination of --methods af2,boltz,rf3,protenix, sharing
a single MSA-generation stage (FOLD_MSA) with a selectable --msa_method. This
deliberately overlaps workflows/boltz_pulldown.nf; it is a standalone
prototype for now and is expected to be folded into a renamed/expanded
boltz_pulldown inside main.nf later (not in this work). See
plans/fold-nf-multi-method-folding.md for the full design.

Multimer: a multi-record FASTA folds as a protein complex (one record = one
chain, in file order -> chain IDs A, B, C, ...; homo-oligomers = repeated
records). Each engine consumes a taxonomically-paired MSA in its own native
format, built by the shared bin/msa_taxonomy.py renderer - see
plans/fold-nf-multimer-paired-msa.md. Header-derived pairing needs the
--msa_method jackhmmer_af2 route (rich UniProt/UniRef headers); the ColabFold
route emits taxonomy-less headers, so ColabFold multimer should use
--use_msa_server (Boltz) instead. AF2 multimer needs the 2021 DB snapshot
(--af2_db_path .../alphafold_20211129, which has uniprot/ + pdb_seqres/).

Usage:
  nextflow run fold.nf --input 'input/*.fasta' --outdir results \
      --methods af2,boltz,rf3,protenix --msa_method jackhmmer_af2 -profile slurm,m3
*/

// Default parameters
params.help = false
params.input = false
params.outdir = 'results'
params.methods = 'af2'
params.msa_method = 'jackhmmer_af2'

// Total predicted structures per input, per method. Default 5 matches the
// typical out-of-the-box behaviour of these tools (AF2's 5 trained models; the
// diffusion engines' usual sample counts). Split across GPU jobs by each
// method's --*_batch_size (samples per invocation); when a method's batch size
// is unset, all N land in ONE job (amortize model load). AF2 has no
// --af2_batch_size (generation batch is fixed at 5; see --af2_keep_models).
params.n_predictions = 5

// Seeds are left UNSET by default: each engine draws its own fresh random seed
// when no seed is passed, so repeated runs already explore different stochastic
// trajectories without us injecting anything. Injecting a random int here
// instead would break -resume - every re-parse would draw a new value, change
// the task command, and bust the cache. Pin any of them (--protenix_seeds /
// --boltz_seed / --rf3_seed / --af2_random_seed) on the CLI for
// reproducibility.

// --- AF2 (monomer-verified route: jackhmmer MSA -> ALPHAFOLD2 predict) ---
params.af2_db_path = '/mnt/datasets/alphafold/alphafold_20240229'
params.af2_model_preset = 'monomer_ptm' // monomer|monomer_ptm|monomer_casp14|multimer
params.af2_db_preset = 'full_dbs'
params.af2_max_template_date = '2024-01-01'
params.af2_random_seed = false
params.af2_num_predictions_per_model = 1 // multimer only
// DB sub-paths (relative to --af2_db_path). Defaults match the monomer
// 20240229 layout so the monomer command is byte-identical. AF2 MULTIMER needs
// the 2021 snapshot (alphafold_20211129), whose HHblits DB is uniclust30 (not
// uniref30) - override --af2_uniref30_subpath accordingly, e.g.
// 'uniclust30/uniclust30_2018_08/uniclust30_2018_08'. uniprot/ + pdb_seqres/
// are multimer-only and present in the 2021 snapshot.
params.af2_uniref30_subpath = 'uniref30/UniRef30_2021_03'
params.af2_uniprot_subpath = 'uniprot/uniprot.fasta'
params.af2_pdb_seqres_subpath = 'pdb_seqres/pdb_seqres.txt'
// mgnify filename differs by snapshot (20240229: mgy_clusters_2022_05.fa;
// 2021 snapshot: mgy_clusters_2018_12.fa).
params.af2_mgnify_subpath = 'mgnify/mgy_clusters_2022_05.fa'
// AF2 model weights (params/) source. --data_dir only locates params/; the
// genetic DBs are passed as separate flags. Defaults to --af2_db_path (monomer
// path unchanged). AF2 MULTIMER needs multimer_v3 weights, which the 2021 DB
// snapshot lacks (it ships only v1 multimer params) - so for multimer point
// --af2_data_dir at a v3-params snapshot (e.g. alphafold_20240229) while the DB
// paths stay on the 2021 snapshot.
params.af2_data_dir = false
// AF2 has no in-invocation sampling knob (a monomer run always emits its 5
// trained models). --af2_keep_models chooses which of those 5 we retain, which
// also sets the effective kept-per-run used to meet --n_predictions:
//   all  = keep all 5 models/run  -> ceil(n_predictions / 5) runs
//   best = keep only the top-ranked -> n_predictions runs (5x more GPU work)
// There is no --af2_batch_size: the generation batch is fixed at 5.
// --models_to_relax follows keep_models (relax best or all) unless --af2_no_relax.
params.af2_keep_models = 'all' // all|best
params.af2_no_relax = false // skip Amber relax; publish unrelaxed into fold/predictions/

// Gather every engine's per-sample mmCIF prediction into a single flat
// <outdir>/fold/predictions/ dir with a tool prefix (af2_ / boltz_ / rf3_ /
// protenix_). BOLTZ's shared module gates that second publishDir on step_name
// starting with 'fold/' (fold.nf passes 'fold/boltz').
params.fold_predictions_dir = "${params.outdir}/fold/predictions"
// ColabFold MSA publish subdir name under workflow_publish_dir. fold.nf leaves
// the on-disk 'result/' name (published at fold/msa/mmseqs2_colabfold/result/);
// boltz_pulldown also uses 'result' via nextflow.config.
params.colabfold_msa_publish_name = 'result'

// --- Boltz ---
params.use_msa_server = false // Boltz's own MMseqs2 server, independent of --msa_method
params.templates = false
params.boltz_recycling = false // --recycling_steps override (default: boltz's own default, currently 3)
params.boltz_batch_size = false // samples per Boltz job (--diffusion_samples); unset + --n_predictions => one job of N
params.boltz_sampling_steps = false // --sampling_steps override (default: 200)
params.boltz_seed = false // Boltz --seed (unset -> Boltz draws its own random seed; pin for reproducibility)

// --- RF3 (same defaults as workflows/rfd3.nf) ---
params.rf3_ckpt_path = '/models/foundry/rf3_foundry_01_24_latest_remapped.ckpt'
params.rf3_num_steps = 50 // default for rf3 cli is 200, however 50 is faster with no difference in quality
params.rf3_n_recycles = 10
// fold.nf meaning: diffusion samples per RF3_FOLD job. Distinct from rfd3.nf's
// --rf3_batch_size (MPNN designs per binder-design ROSETTAFOLD3 task).
params.rf3_batch_size = false // samples per RF3 job (diffusion_batch_size); unset => --n_predictions in one job
params.rf3_early_stopping_plddt_threshold = 0.5 // exits early if mean pLDDT < 0.5 after the first recycle
params.rf3_seed = false // RF3 hydra `seed=` (unset -> RF3 draws its own random seed; pin for reproducibility)

// --- Protenix (defaults confirmed against `protenix pred --help` in
// protenix:v2.0.0-weights on 2026-07-17 - identical to that build's own CLI
// defaults) ---
params.protenix_seeds = false // Comma-separated --seeds (unset -> protenix uses its own default; pin for reproducibility)
params.protenix_cycle = 10 // pairformer cycles (protenix's own --cycle default)
params.protenix_step = 200 // diffusion steps (protenix's own --step default)
params.protenix_batch_size = false // samples per Protenix job (--sample); unset + --n_predictions => one job of N
params.protenix_model_name = 'protenix_base_default_v1.0.0' // checkpoint baked into the container at /models/protenix/checkpoint
params.protenix_use_msa = true // false forces protenix to ignore our a3m and predict MSA-free

// --- MSA subsample (CF-random-style shallow random MSA per predict task) ---
// false = off; true = CF-random default depths; or "1:2,4:8,..." custom list.
params.msa_subsample = false
params.msa_subsample_include_full = true // also keep one full-MSA job when subsampling

// --- EnGens (post-prediction conformational clustering; on by default) ---
params.skip_engens = false // skip EnGens clustering / report
params.engens_dimred = 'umap' // umap only for static predicted ensembles
params.engens_clustering = 'hdbscan' // hdbscan (default); gmm, km, or comma-separated combinations
params.engens_min_structures = 3 // skip clustering below this many usable structures
params.engens_max_clusters = 10 // upper bound for auto cluster-count search
params.engens_gmm_ic = 'aic' // aic|bic for GMM information-criterion selection
params.engens_seed = false // optional RNG seed for UMAP / clustering

// --- ColabFold MSA (--msa_method mmseqs2_colabfold) ---
params.use_remote_server = false // query the ColabFold MMseqs2 API instead of a local DB search
params.uniref30 = false
params.colabfold_envdb = false

// require_gpu / gpu_devices / gpu_allocation_detect_process_regex are already
// defined with pipeline-wide defaults in nextflow.config

include { FOLD_MSA } from './subworkflows/local/fold_msa'
include { ALPHAFOLD2 } from './subworkflows/local/alphafold2'
include { BOLTZ_FOLD } from './subworkflows/local/boltz_fold'
include { ROSETTAFOLD3_FOLD } from './subworkflows/local/rosettafold3_fold'
include { PROTENIX_FOLD } from './subworkflows/local/protenix_fold'
include { ENGENS_CLUSTER } from './subworkflows/local/engens'

def VALID_METHODS = ['af2', 'boltz', 'rf3', 'protenix']
def VALID_MSA_METHODS = ['jackhmmer_af2', 'mmseqs2_colabfold']

workflow {
    if (params.help || params.input == false) {
        log.info(
            """
        ==================================================================
        🧬 FOLD WORKFLOW (multi-method structure prediction) 🧬
        ==================================================================

        Predict structures for one or more FASTA files with any combination of
        AlphaFold2, Boltz-2, RosettaFold3 and Protenix, sharing one
        MSA-generation stage.

        Multimer: a multi-record FASTA folds as a protein complex (one record
        = one chain -> chain IDs A, B, C, ...; homo-oligomers = repeated
        records; up to 26 chains). Header-derived MSA pairing needs
        --msa_method jackhmmer_af2; ColabFold multimer should use
        --use_msa_server. AF2 multimer needs the 2021 DB snapshot
        (see --af2_db_path).

        Required arguments:
            --input                            Single FASTA file, glob, or directory of FASTA files.
                                                Each file is one prediction unit (multi-record = one
                                                complex, folded as chains A, B, C, ...).

        Optional arguments:
            --outdir                           Output directory [default: ${params.outdir}]
            --methods                          Comma-separated list of af2,boltz,rf3,protenix [default: ${params.methods}]
            --msa_method                       jackhmmer_af2|mmseqs2_colabfold [default: ${params.msa_method}]
            --n_predictions                    Total structures per input, per method. Split across jobs by
                                                --boltz_batch_size / --rf3_batch_size / --protenix_batch_size
                                                (AF2 uses --af2_keep_models; generation batch is fixed at 5).
                                                With no method batch size, all N are one job.
                                                [default: ${params.n_predictions}]

            AlphaFold2 (--methods includes af2):
            --af2_db_path                       AlphaFold2 database directory [default: ${params.af2_db_path}]
            --af2_model_preset                  monomer|monomer_ptm|monomer_casp14|multimer [default: ${params.af2_model_preset}]
            --af2_db_preset                     full_dbs|reduced_dbs [default: ${params.af2_db_preset}]
            --af2_max_template_date             Maximum template release date [default: ${params.af2_max_template_date}]
            --af2_random_seed                   Fix the data pipeline's random seed [default: unset -> AF2 draws its own random seed per run]
            --af2_keep_models                   Which of AF2's 5 models/run to keep toward --n_predictions
                                                (also which to relax, unless --af2_no_relax):
                                                'all'  = keep 5/run  -> ceil(N/5) runs;
                                                'best' = keep 1/run  -> N runs [default: ${params.af2_keep_models}]
            --af2_no_relax                      Skip Amber relaxation; publish unrelaxed models to
                                                fold/predictions/ [default: ${params.af2_no_relax}]

            Boltz-2 (--methods includes boltz):
            --use_msa_server                   Use Boltz's own MMseqs2 MSA server instead of the shared --msa_method a3m [default: ${params.use_msa_server}]
            --templates                        Templates directory with .cif files [default: ${params.templates}]
            --boltz_recycling                  Boltz --recycling_steps override [default: boltz's own default]
            --boltz_batch_size                 Samples per Boltz job (--diffusion_samples). Splits
                                                --n_predictions into ceil(N/batch) jobs; unset => one job of N.
            --boltz_sampling_steps              Boltz --sampling_steps override [default: boltz's own default]
            --boltz_seed                        Boltz --seed [default: unset -> Boltz draws its own random seed]

            RosettaFold3 (--methods includes rf3):
            --rf3_ckpt_path                     RF3 checkpoint path [default: ${params.rf3_ckpt_path}]
            --rf3_num_steps                     [default: ${params.rf3_num_steps}]
            --rf3_n_recycles                    [default: ${params.rf3_n_recycles}]
            --rf3_batch_size                    Samples per RF3 job (diffusion_batch_size). Splits
                                                --n_predictions into ceil(N/batch) jobs; unset => one job of N.
            --rf3_early_stopping_plddt_threshold [default: ${params.rf3_early_stopping_plddt_threshold}]
            --rf3_seed                          RF3 hydra seed= [default: unset -> RF3 draws its own random seed]

            Protenix (--methods includes protenix):
            --protenix_seeds                    Single seed (or comma-separated only when not batching) [default: unset -> protenix's own default]
            --protenix_cycle                    Pairformer cycles [default: ${params.protenix_cycle}]
            --protenix_step                     Diffusion steps [default: ${params.protenix_step}]
            --protenix_batch_size               Samples per Protenix job (--sample). Splits
                                                --n_predictions into ceil(N/batch) jobs; unset => one job of N.
            --protenix_model_name               Checkpoint name (baked into the container) [default: ${params.protenix_model_name}]
            --protenix_use_msa                  Feed our shared a3m to Protenix; false predicts MSA-free [default: ${params.protenix_use_msa}]

            MSA subsample (optional CF-random-style shallow MSA per predict task):
            --msa_subsample                     false (default), true (CF-random depths
                                                1:2,2:4,4:8,8:16,16:32,32:64,64:128), or a custom
                                                comma-separated max_seq:max_extra_seq list.
                                                Depths with max_seq >= MSA size are skipped.
                                                [default: ${params.msa_subsample}]
            --msa_subsample_include_full        When subsampling, also keep one full-MSA job
                                                (AF2 reuses templated features.pkl for that job)
                                                [default: ${params.msa_subsample_include_full}]

            ColabFold MSA (--msa_method mmseqs2_colabfold):
            --use_remote_server                 Query the ColabFold MMseqs2 API instead of a local DB search [default: ${params.use_remote_server}]
            --uniref30                          UniRef30 database path (local search only) [default: ${params.uniref30}]
            --colabfold_envdb                   ColabFold environment database path (local search only) [default: ${params.colabfold_envdb}]
                                                Note: no M3 default path exists for these local ColabFold DBs
                                                (see plans/fold-nf-multi-method-folding.md); use --use_remote_server
                                                true if you haven't set them up yourself.

            EnGens (runs by default after prediction; UMAP + HDBSCAN clustering):
            --skip_engens                       Skip EnGens clustering / clusters.html report [default: ${params.skip_engens}]
            --engens_clustering                 hdbscan (default), gmm, km, or comma-separated
                                                combinations [default: ${params.engens_clustering}]
            --engens_dimred                     Dimensionality reduction (umap) [default: ${params.engens_dimred}]
            --engens_min_structures             Minimum usable structures before clustering [default: ${params.engens_min_structures}]
            --engens_max_clusters               Upper bound for auto cluster-count search [default: ${params.engens_max_clusters}]
            --engens_gmm_ic                     GMM information criterion: aic|bic [default: ${params.engens_gmm_ic}]
            --engens_seed                       Optional RNG seed for UMAP / clustering [default: unset]
                                                To cluster an existing folder of .cif/.pdb only (no folding),
                                                use the standalone engens.nf workflow instead.

            --gpu_devices                        GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
            --gpu_allocation_detect_process_regex  Regex pattern to detect busy GPU processes [default: ${params.gpu_allocation_detect_process_regex}]

        Example:
            nextflow run fold.nf --input 'input/*.fasta' --outdir results \\
                --methods af2,boltz,rf3,protenix --msa_method jackhmmer_af2 -profile slurm,m3

        """.stripIndent()
        )
        exit(1)
    }

    def methods = params.methods.toString().split(',').collect { it.trim().toLowerCase() }
    methods.each { m ->
        if (!(m in VALID_METHODS)) {
            error("fold.nf: unknown --methods entry '${m}' (valid: ${VALID_METHODS.join(', ')})")
        }
    }
    if (!(params.msa_method in VALID_MSA_METHODS)) {
        error("fold.nf: unknown --msa_method '${params.msa_method}' (valid: ${VALID_MSA_METHODS.join(', ')})")
    }
    if ('af2' in methods) {
        if (!(params.af2_keep_models in ['all', 'best'])) {
            error("fold.nf: --af2_keep_models must be 'all' or 'best' (got '${params.af2_keep_models}')")
        }
        // 'all' keeps 5 models/run, so n_predictions that isn't a multiple of 5
        // rounds UP to the next multiple (an extra run's worth of models).
        if (params.n_predictions && params.af2_keep_models == 'all' && (params.n_predictions as int) % 5 != 0) {
            def n_runs = ((params.n_predictions as int) + 4).intdiv(5)
            log.warn(
                "fold.nf: AF2 --af2_keep_models=all emits 5 models per run; --n_predictions " +
                "${params.n_predictions} is not a multiple of 5, so AF2 will produce ${n_runs * 5} " +
                "models across ${n_runs} runs (nearest multiple of 5 >= n_predictions)."
            )
        }
    }
    ['boltz_batch_size', 'rf3_batch_size', 'protenix_batch_size'].each { pname ->
        def v = params[pname]
        // Default false is Boolean; a CLI 0 is Integer/String - do not use Groovy
        // truthiness (0 == false) or we would silently treat 0 as unset.
        if (v != null && !(v instanceof Boolean)) {
            def n = v as int
            if (n < 1) {
                error("fold.nf: --${pname} must be >= 1 (got '${v}')")
            }
        }
    }
    if (params.n_predictions && (params.n_predictions as int) < 1) {
        error("fold.nf: --n_predictions must be >= 1 (got '${params.n_predictions}')")
    }
    def engens_clustering = params.engens_clustering.toString().split(',').collect { it.trim().toLowerCase() }
    engens_clustering.each { c ->
        if (!(c in ['gmm', 'km', 'hdbscan'])) {
            error("fold.nf: unknown --engens_clustering entry '${c}' (valid: gmm, km, hdbscan)")
        }
    }
    if (!(params.engens_dimred.toString().toLowerCase() in ['umap'])) {
        error("fold.nf: --engens_dimred must be 'umap' (got '${params.engens_dimred}')")
    }
    if (!(params.engens_gmm_ic.toString().toLowerCase() in ['aic', 'bic'])) {
        error("fold.nf: --engens_gmm_ic must be 'aic' or 'bic' (got '${params.engens_gmm_ic}')")
    }
    if ((params.engens_min_structures as int) < 2) {
        error("fold.nf: --engens_min_structures must be >= 2 (got '${params.engens_min_structures}')")
    }
    if ((params.engens_max_clusters as int) < 2) {
        error("fold.nf: --engens_max_clusters must be >= 2 (got '${params.engens_max_clusters}')")
    }
    // Validate --msa_subsample (false | true | "1:2,4:8,...")
    if (MsaSubsample.isEnabled(params.msa_subsample)) {
        def raw = (params.msa_subsample instanceof Boolean \
            || params.msa_subsample.toString().trim().toLowerCase() in ['true', '1', 'yes']) \
            ? MsaSubsample.MSA_DEPTH_DEFAULTS \
            : params.msa_subsample.toString().trim()
        raw.split(',').each { part ->
            def pdepth = part.trim()
            if (!pdepth) {
                return
            }
            if (!(pdepth ==~ /^\d+:\d+$/)) {
                error(
                    "fold.nf: invalid --msa_subsample depth '${pdepth}' " +
                    "(expected max_seq:max_extra_seq, or true for CF-random defaults)"
                )
            }
            def bits = pdepth.split(':')
            if ((bits[0] as int) < 1) {
                error("fold.nf: --msa_subsample max_seq must be >= 1 (got '${pdepth}')")
            }
        }
    }
    if (params.msa_method == 'mmseqs2_colabfold' && !params.use_remote_server && !(params.uniref30 && params.colabfold_envdb)) {
        error(
            "fold.nf: --msa_method mmseqs2_colabfold needs --uniref30 and --colabfold_envdb " +
            "(local ColabFold DBs) or --use_remote_server true. No M3 default local ColabFold " +
            "DB path exists yet (see plans/fold-nf-multi-method-folding.md) - pass " +
            "--use_remote_server true to query the ColabFold API instead."
        )
    }

    // Accept a single file, a glob, or a directory of FASTA files (one FASTA
    // file = one prediction unit). Resolved
    // synchronously via the `file()` DSL function (which returns a List<Path>
    // for a glob pattern) rather than a lazy Channel.fromPath, so the
    // chain-count validation below can fail the whole run up front, before any
    // process is scheduled. NOTE: Channel.fromPath(...).toList().getVal()
    // looks equivalent but deadlocks here - the dataflow scheduler isn't
    // pumped yet at this point in workflow-body execution, so avoid it for
    // this kind of eager, pre-channel validation.
    def p = params.input
    def resolved = file(p).isDirectory() ? file("${p}/*.{fasta,fa,faa}") : file(p)
    List input_paths = (resolved instanceof List) ? resolved : [resolved]

    if (!input_paths) {
        error("fold.nf: no FASTA files found for --input '${p}'")
    }

    // One FASTA file = one complex. Each record = one chain (file order ->
    // chain IDs A, B, C, ...). Validate eagerly (file() already read the paths
    // synchronously above) so bad inputs fail the whole run before any process
    // is scheduled, rather than mis-folding.
    def MAX_CHAINS = 26 // chain IDs A..Z
    def chain_counts = [:]
    input_paths.each { f ->
        int n_chains = 0
        int empty_records = 0
        boolean in_record = false
        int seq_len = 0
        f.eachLine { line ->
            def l = line.trim()
            if (l.startsWith('>')) {
                if (in_record && seq_len == 0) {
                    empty_records += 1
                }
                n_chains += 1
                in_record = true
                seq_len = 0
            }
            else if (in_record) {
                seq_len += l.length()
            }
        }
        if (in_record && seq_len == 0) {
            empty_records += 1
        }

        if (n_chains == 0) {
            error("fold.nf: ${f} contains no FASTA records.")
        }
        if (n_chains > MAX_CHAINS) {
            error(
                "fold.nf: ${f} has ${n_chains} FASTA records but at most ${MAX_CHAINS} chains " +
                "(A-Z) are supported this round. Reduce the chain count."
            )
        }
        if (empty_records > 0) {
            error("fold.nf: ${f} has ${empty_records} FASTA record(s) with an empty sequence.")
        }
        chain_counts[f.toString()] = n_chains
    }

    def has_multimer = chain_counts.values().any { it > 1 }
    // --msa_subsample is an AF2/monomer-centric depth sweep that operates on a
    // single per-target a3m; multimer feeds per-chain MSA bundles instead, so
    // the two don't combine this round.
    if (has_multimer && MsaSubsample.isEnabled(params.msa_subsample)) {
        error("fold.nf: --msa_subsample is not supported for multimer inputs (monomer only).")
    }
    // ColabFold a3m headers carry no taxonomy, so cross-chain pairing can't be
    // derived from them offline (see plans/fold-nf-multimer-paired-msa.md §0).
    if (has_multimer && params.msa_method == 'mmseqs2_colabfold' && !params.use_msa_server) {
        log.warn(
            "fold.nf: multimer input with --msa_method mmseqs2_colabfold - ColabFold a3m " +
            "headers carry no taxonomy, so RF3/Protenix/Boltz will run UNPAIRED. Use " +
            "--msa_method jackhmmer_af2, or --use_msa_server true (Boltz fetches + pairs itself)."
        )
    }
    // AF2 multimer runs its own native jackhmmer + species-pairing pipeline; the
    // ColabFold bridge has no multimer feature path, so require jackhmmer_af2.
    if (has_multimer && 'af2' in methods && params.msa_method != 'jackhmmer_af2') {
        error(
            "fold.nf: AF2 multimer requires --msa_method jackhmmer_af2 (AF2's native " +
            "multimer MSA pipeline). Drop af2 from --methods for ColabFold multimer runs."
        )
    }
    // AF2 multimer needs the 2021 snapshot's uniprot/ all-seqs DB (the default
    // 20240229 snapshot is monomer-only). Fail fast before scheduling.
    if (has_multimer && 'af2' in methods) {
        def uniprot_dir = file("${params.af2_db_path}/uniprot")
        if (!uniprot_dir.exists()) {
            error(
                "fold.nf: AF2 multimer needs a uniprot/ DB under --af2_db_path, absent in " +
                "'${params.af2_db_path}'. Point --af2_db_path at the 2021 snapshot (e.g. " +
                "/mnt/datasets/alphafold/alphafold_20211129), which has uniprot/ + pdb_seqres/."
            )
        }
    }

    // meta.id is the FASTA stem, matching AF2's own per-target output
    // directory naming - the MSA -> predict wiring in
    // subworkflows/local/alphafold2.nf relies on that matching. meta.n_chains
    // drives every downstream multimer branch.
    ch_input = Channel.fromList(input_paths).map { f -> [[id: f.baseName, n_chains: chain_counts[f.toString()]], f] }

    FOLD_MSA(ch_input, methods, params.msa_method)

    // Collect prediction emits for optional EnGens clustering. Empty channels
    // from unused methods are fine to mix - they contribute nothing.
    ch_af2_pred = Channel.empty()
    ch_boltz_pred = Channel.empty()
    ch_rf3_pred = Channel.empty()
    ch_protenix_pred = Channel.empty()

    if ('af2' in methods) {
        // Hybrid: always stage msas_dir (templated features.pkl for full jobs);
        // join real a3m when --msa_subsample is on, else a dummy stub.
        def a3m_stub = file("${projectDir}/assets/dummy_files/empty")
        if (MsaSubsample.isEnabled(params.msa_subsample)) {
            ch_af2_in = FOLD_MSA.out.af2_msas
                .join(FOLD_MSA.out.a3m.map { meta, fasta, a3m -> [meta, a3m] })
                .map { meta, fasta, msas, a3m -> [meta, fasta, msas, a3m] }
        }
        else {
            ch_af2_in = FOLD_MSA.out.af2_msas.map { meta, fasta, msas -> [meta, fasta, msas, a3m_stub] }
        }
        ALPHAFOLD2(ch_af2_in)
        ch_af2_pred = ALPHAFOLD2.out.predictions
    }
    if ('boltz' in methods) {
        BOLTZ_FOLD(FOLD_MSA.out.for_boltz)
        ch_boltz_pred = BOLTZ_FOLD.out.predictions
    }
    if ('rf3' in methods) {
        ROSETTAFOLD3_FOLD(FOLD_MSA.out.for_rf3)
        ch_rf3_pred = ROSETTAFOLD3_FOLD.out.predictions
    }
    if ('protenix' in methods) {
        PROTENIX_FOLD(FOLD_MSA.out.for_protenix)
        ch_protenix_pred = PROTENIX_FOLD.out.predictions
    }

    if (!params.skip_engens) {
        ch_engens_in = ch_af2_pred
            .mix(ch_boltz_pred, ch_rf3_pred, ch_protenix_pred)
        ENGENS_CLUSTER(ch_engens_in)
    }

    ///////////////////////////////////////////////////////////////////////////
    workflow.onComplete = {
        def fold_params_json = [:]

        fold_params_json['params'] = ParamsHelper.paramsToMap(params)

        fold_params_json['workflow'] = [
            name: workflow.manifest.name,
            version: workflow.manifest.version,
            revision: workflow.revision ?: null,
            commit: workflow.commitId ?: null,
            runName: workflow.runName,
            start: workflow.start.format('yyyy-MM-dd HH:mm:ss'),
            complete: workflow.complete.format('yyyy-MM-dd HH:mm:ss'),
            duration: workflow.duration,
            success: workflow.success,
        ]

        def output_file = "${params.outdir}/fold/params.json"
        def json_string = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(fold_params_json))

        new File(output_file).parentFile.mkdirs()
        new File(output_file).text = json_string

        log.info("Parameters saved to: ${output_file}")
    }
}
