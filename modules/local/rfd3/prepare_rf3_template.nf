// Prepares the target structure for use as the single RF3 template:
// Resolve contigs from contigs_string (params mode) or from contigs_or_config file (config mode: JSON or plain).
// 0) Convert CIF to PDB if needed (gemmi)
// 1) Trim to contigs so the template matches what RFDiffusion3 sees (trim_to_contigs.py; supports RFD3 v3 contig syntax)
// 2) Rename all chains to target_chain (gemmi) so the template chain_id matches RFD3 output for the target polymer
process PREPARE_RF3_TEMPLATE {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple path(structure), path(contigs_or_config), val(contigs_string), val(target_chain)

    output:
    path('template_rf3.pdb'), emit: pdb
    publishDir path: "${params.outdir}/rfd3/rf3_template", pattern: 'template_rf3.pdb', mode: 'copy'

    script:
    def esc = contigs_string?.toString()?.replace("'", "'\"'\"'") ?: ''
    """
    set -euo pipefail

    # Resolve contigs: from value (params mode) or from file (config mode)
    if [[ -n '${esc}' ]]; then
        echo -n '${esc}' > contigs.txt
        contigs=\$(cat contigs.txt)
    elif [[ "${contigs_or_config}" == *.json ]]; then
        contigs=\$(python3 -c "import json; c=json.load(open('${contigs_or_config}')); v=list(c.values())[0] if c else {}; print(v.get('contig',''))")
    else
        contigs=\$(cat "${contigs_or_config}" | tr -d '\\n' | sed 's#^[[:space:]]*##;s#[[:space:]]*\$##')
    fi

    # Step 0: Convert CIF to PDB if needed (trim_to_contigs uses Biopython PDB parser)
    work="${structure}"
    if [[ "${structure}" == *.cif ]] || [[ "${structure}" == *.cif.gz ]]; then
        gemmi convert "${structure}" work.pdb
        work=work.pdb
    fi

    # Step 1: Trim to contigs (skip if no contigs so we keep full structure)
    if [[ -n "\${contigs}" ]]; then
        python3 ${projectDir}/bin/trim_to_contigs.py "\${work}" "\${contigs}" -o trimmed.pdb
    else
        cp "\${work}" trimmed.pdb
    fi

    # Step 2: Rename all chains to target_chain (matches RFD3 target polymer letter from contig order)
    chains=\$(gemmi residues -c -s -s -s trimmed.pdb 2>/dev/null | tail -n +2 | awk '{print \$1}' | sort -u || true)
    rename_args=""
    for ch in \$chains; do
        if [[ "\$ch" != "${target_chain}" ]]; then
            rename_args="\${rename_args} --rename-chain=\${ch}:${target_chain}"
        fi
    done
    if [[ -n "\${rename_args}" ]]; then
        gemmi convert trimmed.pdb template_rf3.pdb \${rename_args}
    else
        cp trimmed.pdb template_rf3.pdb
    fi
    """
}
