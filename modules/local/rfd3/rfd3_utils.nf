// Utility functions for RFD3 workflows

import groovy.transform.Field

@Field
def AA_1TO3 = [
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
    'X': 'UNK',
]

/**
 * Resolve a parameter with dual naming (legacy pmpnn_* and new mpnn_*).
 * The modern (new) value takes precedence if explicitly set.
 */
def resolveParam(legacy, modern) {
    if (modern != null && modern != false) {
        return modern
    }
    return legacy
}

/**
 * Convert a one-letter amino acid omit string (eg "CX") to
 * the mpnn --omit JSON list format (eg '["CYS","UNK"]').
 */
def convertOmitAas(String omitStr) {
    if (!omitStr) {
        return null
    }
    def threeLetterList = omitStr.toUpperCase().collect { ch ->
        AA_1TO3.containsKey(ch) ? "\"${AA_1TO3[ch]}\"" : null
    }.findAll { it != null }

    return "[${threeLetterList.join(',')}]"
}

/**
 * Build the MPNN CLI argument string from resolved params.
 */
def buildMpnnArgs(params) {
    def model_type = resolveParam(null, params.mpnn_model_type)
    def legacy_weights = resolveParam(null, params.mpnn_legacy_weights)
    def designed_chains = resolveParam(null, params.mpnn_designed_chains)
    def batch_size = resolveParam(params.pmpnn_seqs_per_struct, params.mpnn_batch_size)
    def temperature = resolveParam(params.pmpnn_temperature, params.mpnn_temperature)
    def structure_noise = resolveParam(params.pmpnn_augment_eps, params.mpnn_structure_noise)
    def checkpoint_path = resolveParam(params.pmpnn_weights, params.mpnn_checkpoint_path)

    // Resolve omit: new mpnn_omit can be either 1-letter (auto-converted) or already 3-letter JSON list
    def omit_raw = resolveParam(params.pmpnn_omit_aas, params.mpnn_omit)

    def args = []
    args << "--model_type ${model_type}"
    args << "--is_legacy_weights ${legacy_weights ? 'True' : 'False'}"
    args << "--designed_chains ${designed_chains}"
    args << "--batch_size ${batch_size}"
    args << "--temperature ${temperature}"
    args << "--structure_noise ${structure_noise}"
    args << "--write_structures True"
    args << "--write_fasta True"

    if (omit_raw) {
        def omit_str = omit_raw.toString()
        // If it looks like a JSON list already, pass through; otherwise convert from 1-letter
        if (omit_str.startsWith('[')) {
            args << "--omit '${omit_str}'"
        } else {
            def converted = convertOmitAas(omit_str)
            if (converted) {
                args << "--omit '${converted}'"
            }
        }
    }

    args << "--checkpoint_path ${checkpoint_path ?: '/weights/proteinmpnn_v_48_020.pt'}"

    return args.join(' ')
}
