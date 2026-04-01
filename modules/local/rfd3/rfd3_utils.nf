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
 * Normalise contig string for RFD3: if it looks like v1 style (starts with [ and ends with ]),
 * translate to v3 style (comma-separated, e.g. A18-132,/0,65-120). Otherwise pass through.
 */
def normaliseContigToV3(String contig) {
    def s = contig?.trim() ?: ''
    if (s.length() >= 2 && s.startsWith('[') && s.endsWith(']')) {
        s = s[1..-2].trim()
        s = s.replace(' ', ',')
        s = s.replace('/', ',/')
        return s
    }
    return contig
}

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
 * First chain ID from resolved designed_chains (e.g. "A,B" -> "A").
 */
def mpnnDesignedChainsFirst(String resolvedDesignedChains) {
    def s = resolvedDesignedChains?.toString()?.trim()
    if (!s) {
        return 'B'
    }
    return s.split(',')[0].trim()
}

/**
 * RFD3 target/binder chain letters (polymer order: first polymer A, second B, …).
 * Two-polymer contigs only. With explicit --mpnn_designed_chains, passes first letter as --binder.
 */
List<String> resolveRfd3TargetBinderChains(projectDir, params) {
    // ProcessBuilder needs java.lang.String elements, not GString (avoids arraycopy type mismatch).
    def cmdList = new ArrayList<String>()
    cmdList.add('python3')
    cmdList.add("${projectDir}/bin/rfd3/stage_rfd3_config.py".toString())
    cmdList.add('infer-rfd3-chain-pair')
    if (params.rfd3_config) {
        cmdList.add('--rfd3-config')
        cmdList.add(file(params.rfd3_config).toString())
    } else if (params.contigs?.toString()?.trim()) {
        cmdList.add('--contig')
        cmdList.add(params.contigs.toString().trim())
    } else {
        throw new Exception('Cannot infer RFD3 target/binder chains without --rfd3_config or --contigs')
    }
    def mpnnRaw = params.mpnn_designed_chains?.toString()?.trim()
    if (mpnnRaw && !mpnnRaw.equalsIgnoreCase('auto')) {
        cmdList.add('--binder')
        cmdList.add(mpnnDesignedChainsFirst(mpnnRaw).toString())
    }
    def pb = new ProcessBuilder(cmdList)
    pb.redirectErrorStream(true)
    def proc = pb.start()
    proc.waitFor()
    def out = proc.inputStream.text.trim()
    if (proc.exitValue() != 0) {
        throw new Exception("stage_rfd3_config.py infer-rfd3-chain-pair failed: ${out}")
    }
    def line = out ? out.split(/\r?\n/)[0].trim() : ''
    if (!line || !line.contains(',')) {
        throw new Exception("stage_rfd3_config.py infer-rfd3-chain-pair produced invalid output: ${line}")
    }
    def parts = line.split(',', 2)
    return [parts[0].trim().toString(), parts[1].trim().toString()]
}

/**
 * True if the user set a non-default workflow param (not null, not Groovy false, not string "false").
 */
boolean mpnnWorkflowParamSet(def v) {
    if (v == null || v == false) {
        return false
    }
    def s = v.toString().trim()
    return s && !s.equalsIgnoreCase('false')
}

/**
 * Rewrite params.mpnn_weights_noise to canonical string tokens (005, 010, 020, 030) after Nextflow/CLI coercion.
 * Call at workflow entry before validation so logs and params.json dumps match tier names.
 */
def canonicalizeMpnnWeightsNoiseParam(params) {
    if (!mpnnWorkflowParamSet(params.mpnn_weights_noise)) {
        return
    }
    def n = normalizeMpnnWeightsNoiseToken(params.mpnn_weights_noise)
    if (n != null) {
        params.mpnn_weights_noise = n
    }
}

/**
 * Normalise --mpnn_weights_noise: accept 005/010/020/030 strings and integer 5/10/20/30 after CLI coercion.
 */
String normalizeMpnnWeightsNoiseToken(def v) {
    if (!mpnnWorkflowParamSet(v)) {
        return null
    }
    def s = v.toString().trim()
    def allowed = ['005', '010', '020', '030'] as Set
    if (allowed.contains(s)) {
        return s
    }
    if (s ==~ /^\d+$/) {
        def i = s.toInteger()
        if (i in [5, 10, 20, 30]) {
            return String.format('%03d', i)
        }
    }
    return null
}

/**
 * Validate --mpnn_preset and --mpnn_weights_noise before the workflow runs.
 */
def validateRfd3MpnnPresetParams(params) {
    def allowedPresets = ['vanilla', 'soluble', 'hyper'] as Set

    if (mpnnWorkflowParamSet(params.mpnn_preset)) {
        def p = params.mpnn_preset.toString().trim().toLowerCase()
        if (!allowedPresets.contains(p)) {
            throw new Exception("--mpnn_preset must be one of vanilla, soluble, hyper (got '${params.mpnn_preset}')")
        }
    }

    if (mpnnWorkflowParamSet(params.mpnn_weights_noise)) {
        def n = normalizeMpnnWeightsNoiseToken(params.mpnn_weights_noise)
        if (n == null) {
            throw new Exception("--mpnn_weights_noise must be one of 005, 010, 020, 030 (got '${params.mpnn_weights_noise}')")
        }
    }

    if (mpnnWorkflowParamSet(params.mpnn_preset)
            && params.mpnn_preset.toString().trim().equalsIgnoreCase('hyper')
            && mpnnWorkflowParamSet(params.mpnn_weights_noise)
            && normalizeMpnnWeightsNoiseToken(params.mpnn_weights_noise) == '005') {
        throw new Exception('--mpnn_weights_noise 005 is not valid with --mpnn_preset hyper (no v48_005 hyper checkpoint; use 010, 020, 030, or set --mpnn_checkpoint_path to a custom file)')
    }
}

/**
 * Resolve MPNN checkpoint path and model settings. Explicit checkpoint wins; otherwise preset + weights_noise.
 * Vanilla/soluble filenames use _002.pt for CLI noise 005 (0.05 Å training noise).
 */
Map resolveMpnnCheckpoint(params) {
    def explicit = resolveParam(params.pmpnn_weights, params.mpnn_checkpoint_path)
    def userModelType = resolveParam(null, params.mpnn_model_type)
    def userLegacy = resolveParam(null, params.mpnn_legacy_weights)

    if (explicit != null && explicit != false && explicit.toString().trim()) {
        return [
            checkpoint_path: explicit.toString().trim(),
            model_type    : userModelType,
            legacy_weights: userLegacy,
        ]
    }

    def presetSet = mpnnWorkflowParamSet(params.mpnn_preset)
    def noiseSet = mpnnWorkflowParamSet(params.mpnn_weights_noise)
    def noise = noiseSet ? normalizeMpnnWeightsNoiseToken(params.mpnn_weights_noise) : '020'

    String relPath
    if (!presetSet && !noiseSet) {
        relPath = 'proteinmpnn_v_48_020.pt'
    } else if (!presetSet && noiseSet) {
        def fileNoise = noise == '005' ? '002' : noise
        relPath = "proteinmpnn_v_48_${fileNoise}.pt"
    } else {
        def preset = params.mpnn_preset.toString().trim().toLowerCase()
        switch (preset) {
            case 'vanilla':
                def fileNoise = noise == '005' ? '002' : noise
                relPath = "proteinmpnn_v_48_${fileNoise}.pt"
                break
            case 'soluble':
                def fileNoiseS = noise == '005' ? '002' : noise
                relPath = "solublempnn_v_48_${fileNoiseS}.pt"
                break
            case 'hyper':
                relPath = "v48_${noise}_epoch300_hyper.pt"
                break
            default:
                throw new Exception("invalid mpnn_preset: ${preset}")
        }
    }

    return [
        checkpoint_path: "/models/foundry/${relPath}",
        model_type    : 'protein_mpnn',
        legacy_weights: true,
    ]
}

/**
 * Build the MPNN CLI argument string from resolved params.
 */
def buildMpnnArgs(params, String resolvedDesignedChains) {
    def resolved = resolveMpnnCheckpoint(params)
    def model_type = resolved.model_type
    def legacy_weights = resolved.legacy_weights
    def checkpoint_path = resolved.checkpoint_path

    def batch_size = resolveParam(params.pmpnn_seqs_per_struct, params.mpnn_batch_size)
    def temperature = resolveParam(params.pmpnn_temperature, params.mpnn_temperature)
    def structure_noise = resolveParam(params.pmpnn_augment_eps, params.mpnn_structure_noise)

    // Resolve omit: new mpnn_omit can be either 1-letter (auto-converted) or already 3-letter JSON list
    def omit_raw = resolveParam(params.pmpnn_omit_aas, params.mpnn_omit)

    def args = []
    args << "--model_type ${model_type}"
    args << "--is_legacy_weights ${legacy_weights ? 'True' : 'False'}"
    args << "--designed_chains ${resolvedDesignedChains}"
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

    args << "--checkpoint_path ${checkpoint_path}"

    return args.join(' ')
}
