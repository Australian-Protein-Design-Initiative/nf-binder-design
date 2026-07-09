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
 * The modern (new) value takes precedence if explicitly set; otherwise the
 * legacy value is used if explicitly set; otherwise fallback (the real
 * default) applies. Both legacy and modern params must default to `false`
 * (the "not set by user" sentinel) for this to work correctly.
 */
def resolveParam(legacy, modern, fallback = null) {
    if (modern != null && modern != false) {
        return modern
    }
    if (legacy != null && legacy != false) {
        return legacy
    }
    return fallback
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
 * Convert RFDiffusion v1 contig format to v3 (comma-separated, /0 as its own token).
 */
def convertV1ContigsToV3(String contigs) {
    def s = contigs?.trim() ?: ''
    s = s.replaceAll(/^\[+/, '')
    s = s.replaceAll(/\]+$/, '')
    s = s.trim()
    if (s.contains(',') && !s.contains(' ')) {
        return s
    }
    def v3Parts = []
    s.split(/\s+/).each { part ->
        if (part.contains('/0')) {
            // Limit -1: Java/Groovy regex split drops trailing empties; Python str.split does not.
            def segments = part.split('/0', -1)
            segments.eachWithIndex { seg, i ->
                if (seg) {
                    v3Parts << seg
                }
                if (i < segments.size() - 1) {
                    v3Parts << '/0'
                }
            }
        } else {
            v3Parts << part
        }
    }
    return v3Parts.join(',')
}

def normaliseContigForChainInfer(String contig) {
    def s = contig?.trim() ?: ''
    if (!s) {
        return ''
    }
    if (s.length() >= 2 && s.startsWith('[') && s.endsWith(']')) {
        return normaliseContigToV3(s)
    }
    if (s.contains(',') && !s.contains(' ')) {
        return s
    }
    return convertV1ContigsToV3(s)
}

def splitPolymersV3(String v3Contig) {
    def tokens = v3Contig.trim().split(',').collect { it.trim() }.findAll { it }
    def polymers = []
    def current = []
    tokens.each { t ->
        if (t == '/0') {
            polymers << current
            current = []
        } else {
            current << t
        }
    }
    polymers << current
    return polymers
}

def stripSlashSuffix(String part) {
    def idx = part.indexOf('/')
    return idx >= 0 ? part.substring(0, idx).trim() : part.trim()
}

boolean isBinderPolymer(List segments) {
    if (!segments) {
        return false
    }
    return segments.every { seg ->
        def s = stripSlashSuffix(seg as String)
        if (!s) {
            return false
        }
        if (s ==~ /^([A-Za-z]+)(-?\d+)-(-?\d+)$/) {
            return false
        }
        return s ==~ /^\d+-\d+$/
    }
}

String inferBinderChainFromContig(String contig) {
    def v3 = normaliseContigForChainInfer(contig)
    if (!v3) {
        throw new Exception('Empty contig after normalisation')
    }
    def polymers = splitPolymersV3(v3).findAll { it }
    def binderIndices = []
    polymers.eachWithIndex { p, i ->
        if (isBinderPolymer(p)) {
            binderIndices << i
        }
    }
    if (binderIndices.size() != 1) {
        throw new Exception(
            "Expected exactly one binder polymer (length-only segments); found ${binderIndices.size()} in contig '${contig}' (v3='${v3}')"
        )
    }
    return String.valueOf((char) ('A'.charAt(0) + binderIndices[0]))
}

List<String> resolveRfd3TargetBinderChainLetters(String contig, String explicitBinder = null) {
    def v3 = normaliseContigForChainInfer(contig)
    if (!v3) {
        throw new Exception('Empty contig after normalisation')
    }
    def polymers = splitPolymersV3(v3).findAll { it }
    if (polymers.size() != 2) {
        throw new Exception(
            "Expected exactly two contig polymers for RFD3 target/binder chain IDs; got ${polymers.size()} in contig '${contig}' (v3='${v3}')"
        )
    }
    int bi
    String binder
    if (explicitBinder) {
        binder = explicitBinder.trim().split(',')[0].trim().toUpperCase()
        if (binder.length() != 1 || !(binder in ['A', 'B'])) {
            throw new Exception("For a two-polymer contig, explicit binder chain must be A or B; got '${explicitBinder}'")
        }
        bi = ((int) binder.charAt(0)) - ((int) 'A'.charAt(0))
    } else {
        binder = inferBinderChainFromContig(contig)
        bi = ((int) binder.charAt(0)) - ((int) 'A'.charAt(0))
    }
    def ti = 1 - bi
    def target = String.valueOf((char) ('A'.charAt(0) + ti))
    return [target, binder]
}

String contigFromRfd3ConfigPath(String configPath) {
    def pathObj = file(configPath)
    if (!pathObj.exists()) {
        throw new Exception("RFD3 config not found: ${configPath}")
    }
    def name = pathObj.name.toLowerCase()
    if (name.endsWith('.json')) {
        def data = new groovy.json.JsonSlurper().parse(pathObj)
        if (!(data instanceof Map) || data.isEmpty()) {
            throw new Exception("RFD3 config must be a non-empty JSON object: ${configPath}")
        }
        def spec = data.values().iterator().next()
        if (!(spec instanceof Map)) {
            throw new Exception("First config entry must be an object: ${configPath}")
        }
        def c = spec.contig
        if (!(c instanceof String)) {
            throw new Exception("contig must be a string in ${configPath}")
        }
        return c.trim()
    }
    if (name.endsWith('.yaml') || name.endsWith('.yml')) {
        def yaml = Class.forName('org.yaml.snakeyaml.Yaml').newInstance()
        def data = yaml.load(pathObj.text)
        if (!(data instanceof Map) || data.isEmpty()) {
            throw new Exception("RFD3 config must be a non-empty mapping: ${configPath}")
        }
        def spec = data.values().iterator().next()
        if (!(spec instanceof Map)) {
            throw new Exception("First config entry must be a mapping: ${configPath}")
        }
        def c = spec.contig
        if (!(c instanceof String)) {
            throw new Exception("contig must be a string in ${configPath}")
        }
        return c.trim()
    }
    return pathObj.text.split(/\s+/).join(' ').trim()
}

/**
 * Resolved input structure path(s) from an RFD3 JSON config (driver-side; JSON only).
 */
List<String> extractRfd3InputPaths(String configPath, String configDir = null) {
    def pathObj = file(configPath)
    def baseDir = configDir ? file(configDir) : pathObj.parent
    if (!pathObj.exists()) {
        throw new Exception("Config file not found: ${configPath}")
    }
    if (!pathObj.name.toLowerCase().endsWith('.json')) {
        throw new Exception(
            "Driver-side input path parsing requires a JSON config (got ${configPath}). Use .json or provide --input_pdb in params mode."
        )
    }
    def config = new groovy.json.JsonSlurper().parse(pathObj)
    def seen = [] as LinkedHashSet
    def result = []
    config.each { key, spec ->
        if (spec instanceof Map && spec.containsKey('input')) {
            def raw = spec.input.toString()
            def resolvedFile = new File(raw).isAbsolute() ? file(raw) : file("${baseDir}/${raw}")
            def resolved = resolvedFile.toString()
            if (!seen.contains(resolved)) {
                seen << resolved
                result << resolved
            }
        }
    }
    return result
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
List<String> resolveRfd3TargetBinderChains(params) {
    String contig
    if (params.rfd3_config) {
        contig = contigFromRfd3ConfigPath(params.rfd3_config.toString())
    } else if (params.contigs?.toString()?.trim()) {
        contig = params.contigs.toString().trim()
    } else {
        throw new Exception('Cannot infer RFD3 target/binder chains without --rfd3_config or --contigs')
    }
    def mpnnRaw = params.mpnn_designed_chains?.toString()?.trim()
    def explicitBinder = (mpnnRaw && !mpnnRaw.equalsIgnoreCase('auto'))
        ? mpnnDesignedChainsFirst(mpnnRaw)
        : null
    return resolveRfd3TargetBinderChainLetters(contig, explicitBinder)
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
 * Canonical mpnn_weights_noise token (005, 010, 020, 030) after Nextflow/CLI coercion, or null if unset.
 * Does not mutate params; pass the return value through validation and MPNN resolution explicitly.
 */
def canonicalizeMpnnWeightsNoiseParam(params) {
    if (!mpnnWorkflowParamSet(params.mpnn_weights_noise)) {
        return null
    }
    return normalizeMpnnWeightsNoiseToken(params.mpnn_weights_noise)
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
def validateRfd3MpnnPresetParams(params, def mpnnWeightsNoise = null) {
    def allowedPresets = ['vanilla', 'soluble', 'hyper'] as Set

    if (mpnnWorkflowParamSet(params.mpnn_preset)) {
        def p = params.mpnn_preset.toString().trim().toLowerCase()
        if (!allowedPresets.contains(p)) {
            throw new Exception("--mpnn_preset must be one of vanilla, soluble, hyper (got '${params.mpnn_preset}')")
        }
    }

    if (mpnnWorkflowParamSet(params.mpnn_weights_noise)) {
        def n = mpnnWeightsNoise != null ? mpnnWeightsNoise : normalizeMpnnWeightsNoiseToken(params.mpnn_weights_noise)
        if (n == null) {
            throw new Exception("--mpnn_weights_noise must be one of 005, 010, 020, 030 (got '${params.mpnn_weights_noise}')")
        }
    }

    if (mpnnWorkflowParamSet(params.mpnn_preset)
            && params.mpnn_preset.toString().trim().equalsIgnoreCase('hyper')
            && mpnnWorkflowParamSet(params.mpnn_weights_noise)
            && (mpnnWeightsNoise ?: normalizeMpnnWeightsNoiseToken(params.mpnn_weights_noise)) == '005') {
        throw new Exception('--mpnn_weights_noise 005 is not valid with --mpnn_preset hyper (no v48_005 hyper checkpoint; use 010, 020, 030, or set --mpnn_checkpoint_path to a custom file)')
    }
}

/**
 * Resolve MPNN checkpoint path and model settings. Explicit checkpoint wins; otherwise preset + weights_noise.
 * Vanilla/soluble filenames use _002.pt for CLI noise 005 (0.05 Å training noise).
 */
Map resolveMpnnCheckpoint(params, def mpnnWeightsNoise = null) {
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
    def noise = noiseSet
        ? (mpnnWeightsNoise ?: normalizeMpnnWeightsNoiseToken(params.mpnn_weights_noise))
        : '020'

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
def buildMpnnArgs(params, String resolvedDesignedChains, def mpnnWeightsNoise = null) {
    def resolved = resolveMpnnCheckpoint(params, mpnnWeightsNoise)
    def model_type = resolved.model_type
    def legacy_weights = resolved.legacy_weights
    def checkpoint_path = resolved.checkpoint_path

    def batch_size = resolveParam(params.pmpnn_seqs_per_struct, params.mpnn_batch_size, 1)
    def temperature = resolveParam(params.pmpnn_temperature, params.mpnn_temperature, 0.1)
    def structure_noise = resolveParam(params.pmpnn_augment_eps, params.mpnn_structure_noise, 0)

    // Resolve omit: new mpnn_omit can be either 1-letter (auto-converted) or already 3-letter JSON list
    def omit_raw = resolveParam(params.pmpnn_omit_aas, params.mpnn_omit, 'CX')

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
