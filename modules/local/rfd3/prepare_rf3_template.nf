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

// Chain rename only for a user-supplied RF3 template (no contig trim, preserve PDB vs CIF).
process RENAME_RF3_TEMPLATE_CHAINS {
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple path(structure), val(target_chain)

    output:
    path('template_rf3.*'), emit: structure
    publishDir path: "${params.outdir}/rfd3/rf3_template", pattern: 'template_rf3.*', mode: 'copy'

    script:
    def name = structure.name
    def ext = (name.endsWith('.cif') || name.endsWith('.cif.gz')) ? 'cif' : 'pdb'
    def tc = target_chain.toString()
    """
    set -euo pipefail

    # Renaming every chain to the target letter merges multiple polymers into one chain_id and
    # breaks RF3 (duplicate res_id / ambiguous bonds). Single chain: rename if needed. Multiple:
    # keep only the chain that already matches the RFD3 target chain (gemmi MMDB select //CHAIN/).
    readarray -t chain_arr < <(gemmi residues -c -s -s -s '${structure}' 2>/dev/null | tail -n +2 | awk '{print \$1}' | sort -u)
    n=\${#chain_arr[@]}
    if [[ "\${n}" -eq 0 ]]; then
        echo "RENAME_RF3_TEMPLATE_CHAINS: no polymer chains found in '${structure}'" >&2
        exit 1
    fi
    if [[ "\${n}" -eq 1 ]]; then
        single="\${chain_arr[0]}"
        if [[ "\${single}" == "${tc}" ]]; then
            cp '${structure}' template_rf3.${ext}
        else
            gemmi convert --rename-chain="\${single}:${tc}" '${structure}' template_rf3.${ext}
        fi
    else
        has_tc=0
        for c in "\${chain_arr[@]}"; do
            if [[ "\${c}" == "${tc}" ]]; then has_tc=1; break; fi
        done
        if [[ "\${has_tc}" -eq 0 ]]; then
            echo "RENAME_RF3_TEMPLATE_CHAINS: multi-chain template [\${chain_arr[*]}] has no chain ${tc}. Put the RF3 target on chain ${tc}, or use a single-chain structure." >&2
            exit 1
        fi
        gemmi convert --select '//${tc}/' '${structure}' template_rf3.${ext}
    fi
    """
}
