/*
AlphaFold3 input-JSON and run-mode resolution, shared by the monomer and
complex generators, the ALPHAFOLD3 process and ALPHAFOLD3_FOLD, so they cannot
drift apart.

Three independent params control how AF3 is fed and run:

  --af3_paired_msa          true (default)  one paired a3m per chain, referenced
                                            by pairedMsaPath
                            false           "pairedMsa": "" on every chain, i.e.
                                            no cross-chain pairing

  --af3_templates           inline (default) embed templates matched by --templates
                            none             "templates": []
                            search           omit the key so AF3's own data pipeline
                                             searches pdb_seqres/mmcif_files.
                                             Requires --af3_run_data_pipeline true.

  --af3_run_data_pipeline   false (default)  --run_data_pipeline=false, inference only
                            true             data pipeline runs; needs --af3_db_dir

The combination that reproduces how Germinal calls AF3 (no pairing, AF3's own
template search, data pipeline on) is packaged as `-profile af3_germinal_parity`
- see conf/af3_germinal_parity.config. It is a profile rather than a parameter
because nothing in the pipeline needs to know that the combination has a name.
*/

def af3InputTemplateModes() {
    return ['inline', 'none', 'search']
}

boolean af3InputPairedMsa(params) {
    return params.af3_paired_msa as boolean
}

String af3InputTemplatesMode(params) {
    return (params.af3_templates ?: 'inline').toString()
}

boolean af3InputRunDataPipeline(params) {
    return params.af3_run_data_pipeline as boolean
}

// The --paired-msa-mode / --templates-mode flags for make_af3_input.py.
String af3InputModeArgs(params) {
    def args = []
    if (!af3InputPairedMsa(params)) {
        args << '--paired-msa-mode empty'
    }
    def mode = af3InputTemplatesMode(params)
    if (mode != 'inline') {
        args << "--templates-mode ${mode}"
    }
    return args.join(' ')
}

boolean af3InputLiteralSeeds(params) {
    return params.af3_seeds && params.af3_seeds.toString().contains(',')
}

/*
Seeds for one AF3 job. AF3 has no CLI seed flag - modelSeeds lives in the
JSON - and every seed in that list runs --num_diffusion_samples structures.

A single base seed (the default) keeps the historical behaviour: batch i
gets base+i, one seed per job. A comma list is taken literally and ALL of
its seeds go into every job's JSON, which is how Germinal runs AF3
(3 seeds for its initial fold, 5 for its final one).
*/
List af3InputSeedsForBatch(params, int batch_index) {
    if (af3InputLiteralSeeds(params)) {
        return params.af3_seeds.toString().split(',').collect { s -> s.trim() as int }
    }
    def base = params.af3_seeds ? (params.af3_seeds.toString().trim() as int) : 1
    return [base + batch_index]
}
