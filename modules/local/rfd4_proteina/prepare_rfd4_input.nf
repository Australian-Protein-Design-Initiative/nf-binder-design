// RFD4-Proteina only reads .cif, .cif.gz and .bcif inputs, so a PDB is converted here.
// An input that is already CIF passes through untouched, so its own numbering is preserved.
process PREPARE_RFD4_INPUT {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    publishDir path: "${params.outdir}/rfd4_proteina/input", pattern: 'prepared/*', mode: 'copy'

    input:
    path structure

    output:
    path 'prepared/*', emit: cif

    script:
    """
    set -euo pipefail

    mkdir -p prepared
    name="\$(basename "${structure}")"

    case "\${name}" in
        *.cif|*.cif.gz|*.bcif)
            cp "${structure}" "prepared/\${name}"
            ;;
        *)
            # Strip .gz before the format suffix, so foo.pdb.gz becomes foo.cif rather than foo.pdb.cif
            stem="\${name%.gz}"
            stem="\${stem%.*}"
            ${projectDir}/bin/rfd4_proteina/pdb_to_rfd4_cif.py "${structure}" -o "prepared/\${stem}.cif"
            ;;
    esac
    """
}
