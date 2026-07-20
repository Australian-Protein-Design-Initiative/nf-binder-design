process ALPHAFOLD2_JACKHMMER_MSA {
    tag "${meta.id}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/alphafold_cuda12_upstream-c77e5d2_custom-57618c5.sif'

    // Publish the per-target directory (msas/ + features.pkl) as
    // <outdir>/af2/msas/${meta.id}/... . pattern/saveAs operate on the top-level
    // output *item* (here the "${meta.id}" directory and the re-emitted input
    // FASTA), NOT the files nested inside the directory - so "${meta.id}" matches
    // the directory item (published recursively) while excluding the stray
    // "${meta.id}.fasta" file item.
    publishDir(
        path: "${params.outdir}/af2/msas",
        mode: 'copy',
        pattern: "${meta.id}"
    )

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("${meta.id}"), emit: msa

    script:
    // All DB flags below are required on both the MSA and predict stages -
    // run_alphafold.py has no defaults for them and fails flag parsing if any
    // are omitted, even though only this (CPU) stage actually reads them.
    def d = params.alphafold2_db_path
    def db_flags = [
        "--data_dir=${d}",
        "--uniref90_database_path=${d}/uniref90/uniref90.fasta",
        "--mgnify_database_path=${d}/mgnify/mgy_clusters_2022_05.fa",
        "--bfd_database_path=${d}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
        "--uniref30_database_path=${d}/uniref30/UniRef30_2021_03",
        "--pdb70_database_path=${d}/pdb70/pdb70",
        "--template_mmcif_dir=${d}/pdb_mmcif/mmcif_files",
        "--obsolete_pdbs_path=${d}/pdb_mmcif/obsolete.dat",
        "--max_template_date=${params.alphafold2_max_template_date}",
        "--db_preset=${params.alphafold2_db_preset}",
        "--model_preset=${params.alphafold2_model_preset}",
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
