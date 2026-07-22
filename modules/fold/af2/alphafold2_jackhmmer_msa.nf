process ALPHAFOLD2_JACKHMMER_MSA {
    tag "${meta.id}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/alphafold_cuda12_upstream-c77e5d2_custom-57618c5.sif'

    // Split publish: raw msas/ (shared with a3m conversion) under
    // fold/msa/jackhmmer_af2/; AF2-only features.pkl under fold/af2/msas/.
    // Emit "${meta.id}/**" so each nested file is a top-level publish item
    // (pattern/saveAs cannot see inside a directory output item). The directory
    // itself stays on the msa channel for AF2 predict / AF2_MSAS_TO_A3M.
    publishDir(
        path: "${params.outdir}/fold/msa/jackhmmer_af2",
        mode: 'copy',
        saveAs: { filename ->
            def rel = filename.toString()
            return rel.startsWith("${meta.id}/msas/") ? rel : null
        }
    )
    publishDir(
        path: "${params.outdir}/fold/af2/msas",
        mode: 'copy',
        saveAs: { filename ->
            def rel = filename.toString()
            return rel == "${meta.id}/features.pkl" ? rel : null
        }
    )

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("${meta.id}"), emit: msa
    path "${meta.id}/**", emit: msa_files

    script:
    // All DB flags below are required on both the MSA and predict stages -
    // run_alphafold.py has no defaults for them and fails flag parsing if any
    // are omitted, even though only this (CPU) stage actually reads them.
    def d = params.af2_db_path
    def db_flags = [
        "--data_dir=${d}",
        "--uniref90_database_path=${d}/uniref90/uniref90.fasta",
        "--mgnify_database_path=${d}/mgnify/mgy_clusters_2022_05.fa",
        "--bfd_database_path=${d}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
        "--uniref30_database_path=${d}/uniref30/UniRef30_2021_03",
        "--pdb70_database_path=${d}/pdb70/pdb70",
        "--template_mmcif_dir=${d}/pdb_mmcif/mmcif_files",
        "--obsolete_pdbs_path=${d}/pdb_mmcif/obsolete.dat",
        "--max_template_date=${params.af2_max_template_date}",
        "--db_preset=${params.af2_db_preset}",
        "--model_preset=${params.af2_model_preset}",
    ].join(' ')
    """
    # AlphaFold names its per-target output directory after the FASTA stem, which
    # is exactly meta.id (see fold.nf's input channel), so the predict stage can
    # find this directory again with no extra bookkeeping.
    python /app/alphafold/run_alphafold.py \
        --fasta_paths=${fasta} \
        --output_dir=\$PWD \
        --generate_msas_only=true \
        --use_precomputed_msas=false \
        --use_gpu_relax=false \
        --models_to_relax=none \
        ${db_flags}
    """
}
