/*
FOLD_TEMPLATES: match --templates structures to the query chains (fold.nf /
fold_pulldown.nf). Emits a value channel holding the matched-templates
directory, or the empty_templates placeholder when --templates is unset, so
engine input generators can always take it as a path input.
*/

include { MATCH_FOLD_TEMPLATES } from '../../modules/local/fold/common/match_fold_templates'

def templateFiles(spec) {
    def given = file(spec.toString())
    def files = (!(given instanceof List) && given.isDirectory()) \
        ? file("${spec}/*.{pdb,ent,cif,mmcif,pdb.gz,ent.gz,cif.gz,mmcif.gz}") \
        : given
    return (files instanceof List ? files : [files]).findAll { it.exists() }
}

workflow FOLD_TEMPLATES {
    take:
    ch_query_fastas // FASTA path(s) whose records may be templated

    main:
    if (params.templates) {
        def files = templateFiles(params.templates)
        if (!files) {
            error("--templates ${params.templates}: no .pdb / .cif files found")
        }
        MATCH_FOLD_TEMPLATES(files, ch_query_fastas.collect())
        ch_templates = MATCH_FOLD_TEMPLATES.out.templates
    }
    else {
        ch_templates = Channel.value(file("${projectDir}/assets/dummy_files/empty_templates"))
    }

    emit:
    templates = ch_templates
}
