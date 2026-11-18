#!/usr/bin/env nextflow

/*
Assertions for modules/local/common/refold_select_utils.nf.

A runnable script rather than an nf-test case so that it can be run with nothing but
Nextflow. From the repository root:

    nextflow run tests/lib/refold_selection.nf -lib lib

It logs one line and exits 0 when every assertion holds, and names the failing case
otherwise.
*/

nextflow.enable.dsl = 2

include {
    parseRefoldFilters
    refoldNumeric
    refoldPasses
    parseRefoldSortKey
    rankRefoldDesigns
    refoldMaxValue
    refoldNeedsCompleteSet
    refoldReportHeader
    refoldReportRow
    refoldSortKeySuggestions
    refoldStripModelSuffix
    refoldEngines
} from '../../modules/local/common/refold_select_utils'

def check(String name, Closure body) {
    try {
        body.call()
    }
    catch (Throwable t) {
        def cause = (t instanceof java.lang.reflect.InvocationTargetException) ? t.cause : t
        throw new RuntimeException("FAILED: ${name}\n  ${cause.class.simpleName}: ${cause.message}", cause)
    }
}

/*
Run something that must fail, and hand back what it threw.

Nextflow invokes an included function reflectively, so an exception raised inside one
reaches a caller that catches it wrapped in InvocationTargetException with a null
message. Nextflow unwraps it again when nothing catches it - a bad --refold_filter
does reach the user as its own message - but a test that inspects the exception has to
unwrap it itself.
*/
def expectFailure(String name, Closure body) {
    try {
        body.call()
    }
    catch (Throwable t) {
        return (t instanceof java.lang.reflect.InvocationTargetException) ? t.cause : t
    }
    throw new RuntimeException("FAILED: ${name}\n  expected a failure, got none")
}

workflow {

    // --- parsing filter expressions ---

    check('an empty spec is no filters') {
        assert parseRefoldFilters(null) == []
        assert parseRefoldFilters('') == []
        assert parseRefoldFilters(false) == []
        assert parseRefoldFilters('  ;  ') == []
    }

    check('filters split on semicolons and tolerate whitespace') {
        def parsed = parseRefoldFilters(' ipsae >= 0.5 ; pae_interaction<10 ')
        assert parsed.size() == 2
        assert parsed[0] == [metric: 'ipsae', op: '>=', value: '0.5']
        assert parsed[1] == [metric: 'pae_interaction', op: '<', value: '10']
    }

    check('a two-character operator wins over its first character') {
        assert parseRefoldFilters('x>=1')[0].op == '>='
        assert parseRefoldFilters('x<=1')[0].op == '<='
        assert parseRefoldFilters('x!=1')[0].op == '!='
        assert parseRefoldFilters('x>1')[0].op == '>'
    }

    check('an unparseable expression names itself') {
        def e = expectFailure('filter with no operator') { parseRefoldFilters('plddt') }
        assert e instanceof IllegalArgumentException
        assert e.message.contains("'plddt'")
    }

    // --- numeric coercion ---

    check('the ways a score table writes "missing" all read as null') {
        assert refoldNumeric(null) == null
        assert refoldNumeric('') == null
        assert refoldNumeric('  ') == null
        assert refoldNumeric('None') == null
        assert refoldNumeric('NA') == null
        assert refoldNumeric('nan') == null
        assert refoldNumeric(Double.NaN) == null
        assert refoldNumeric('not_a_number') == null
    }

    check('numbers read as numbers whatever their type') {
        assert refoldNumeric('0.5') == 0.5d
        assert refoldNumeric(' 12 ') == 12.0d
        assert refoldNumeric(3) == 3.0d
        assert refoldNumeric('-1e-3') == -0.001d
    }

    // --- evaluating filters ---

    def filters = parseRefoldFilters('ipsae>=0.5;pae_interaction<=10')

    check('a design passes only when every filter holds') {
        assert refoldPasses([ipsae: 0.7, pae_interaction: 8], filters)
        assert !refoldPasses([ipsae: 0.4, pae_interaction: 8], filters)
        assert !refoldPasses([ipsae: 0.7, pae_interaction: 12], filters)
    }

    check('boundaries are inclusive for >= and <=') {
        assert refoldPasses([ipsae: 0.5, pae_interaction: 10], filters)
    }

    check('no filters admits everything, including a design with no scores') {
        assert refoldPasses([:], [])
        assert refoldPasses([ipsae: 0.1], [])
    }

    check('a missing value fails a numeric comparison rather than passing it') {
        // '!=' included: a design with no ipsae is not a design whose ipsae differs.
        assert !refoldPasses([ipsae: 'None', pae_interaction: 5], filters)
        assert !refoldPasses([ipsae: ''], parseRefoldFilters('ipsae!=0.5'))
    }

    check('a filter on a metric the table lacks stops the run and lists what it has') {
        def e = expectFailure('unknown metric') {
            refoldPasses([plddt: 90, ptm: 0.8], parseRefoldFilters('iptm>0.5'), 'design_7')
        }
        assert e instanceof IllegalArgumentException
        assert e.message.contains("'iptm'")
        assert e.message.contains('design_7')
        assert e.message.contains('plddt')
        assert e.message.contains('ptm')
    }

    check('booleans and strings compare by equality') {
        assert refoldPasses([ok: 'True'], parseRefoldFilters('ok==true'))
        assert refoldPasses([ok: false], parseRefoldFilters('ok==false'))
        assert refoldPasses([ok: true], parseRefoldFilters('ok!=false'))
        assert refoldPasses([tag: 'keep'], parseRefoldFilters('tag==keep'))
        assert !refoldPasses([tag: 'drop'], parseRefoldFilters('tag==keep'))
    }

    // --- sort keys ---

    check('a leading dash means descending') {
        assert parseRefoldSortKey('ipsae') == [key: 'ipsae', descending: false]
        assert parseRefoldSortKey('-ipsae') == [key: 'ipsae', descending: true]
        assert parseRefoldSortKey(' -ipsae ') == [key: 'ipsae', descending: true]
    }

    check('no sort key is null, not a sort on nothing') {
        assert parseRefoldSortKey(null) == null
        assert parseRefoldSortKey('') == null
        assert parseRefoldSortKey(false) == null
    }

    check('a dash with no metric is rejected') {
        def e = expectFailure('bare dash') { parseRefoldSortKey('-') }
        assert e instanceof IllegalArgumentException
    }

    // --- ranking and the cap ---

    def entries = [
        ['a', 'a.cif', [ipsae: 0.3]],
        ['b', 'b.cif', [ipsae: 0.9]],
        ['c', 'c.cif', [ipsae: 0.6]],
        ['d', 'd.cif', [ipsae: 'None']],
    ]

    check('descending rank puts the best first') {
        assert rankRefoldDesigns(entries, parseRefoldSortKey('-ipsae'), null).collect { it[0] } == ['b', 'c', 'a', 'd']
    }

    check('ascending rank puts the lowest first') {
        assert rankRefoldDesigns(entries, parseRefoldSortKey('ipsae'), null).collect { it[0] } == ['a', 'c', 'b', 'd']
    }

    check('a design with no value for the metric ranks last either way') {
        assert rankRefoldDesigns(entries, parseRefoldSortKey('-ipsae'), null).last()[0] == 'd'
        assert rankRefoldDesigns(entries, parseRefoldSortKey('ipsae'), null).last()[0] == 'd'
    }

    check('the cap takes the best N once ranked') {
        assert rankRefoldDesigns(entries, parseRefoldSortKey('-ipsae'), 2).collect { it[0] } == ['b', 'c']
    }

    check('the cap without a sort key keeps arrival order') {
        assert rankRefoldDesigns(entries, null, 2).collect { it[0] } == ['a', 'b']
    }

    check('no cap keeps everything') {
        assert rankRefoldDesigns(entries, null, null).size() == 4
    }

    // --- the cap value itself ---

    check('--refold_max reads from the several ways it can be written') {
        assert refoldMaxValue(false) == null
        assert refoldMaxValue(null) == null
        assert refoldMaxValue(0) == null
        assert refoldMaxValue(-1) == null
        assert refoldMaxValue(5) == 5
        assert refoldMaxValue('5') == 5
        assert refoldMaxValue(' 5 ') == 5
    }

    check('the cap is what needs the whole set in hand, with or without a rank') {
        // A threshold gate admits designs one at a time, so refolding overlaps with
        // scoring. A cap cannot: "the best N" needs the last design scored and "the
        // first N" needs a count shared across the channel. A rank with no cap selects
        // nothing, so it does not force a wait either.
        assert refoldNeedsCompleteSet(null, 10)
        assert refoldNeedsCompleteSet('-ipsae', 10)
        assert !refoldNeedsCompleteSet('-ipsae', false)
        assert !refoldNeedsCompleteSet(null, false)
    }

    // --- the report ---

    check('the report row lines up with its header') {
        def header = refoldReportHeader().split('\t')
        def row = refoldReportRow('design_1', [ipsae: 0.72], parseRefoldSortKey('-ipsae'), true, 1).split('\t', -1)
        assert header.size() == row.size()
        assert row[0] == 'design_1'
        assert row[1] == 'selected'
        assert row[2] == '1'
        assert row[3] == 'ipsae'
        assert row[4] == '0.72'
    }

    check('a rejected design reports no rank') {
        def row = refoldReportRow('design_2', [ipsae: 0.1], parseRefoldSortKey('-ipsae'), false, null).split('\t', -1)
        assert row[1] == 'rejected'
        assert row[2] == ''
    }

    // --- suggesting sort keys when the named one is absent ---

    check('suggestions list only the numeric columns') {
        def suggest_entries = [
            [[id: 'design_1'], 'a.pdb', [Average_i_pTM: 0.81, Sequence: 'MKV', plddt_binder: 88.2]],
            [[id: 'design_2'], 'b.pdb', [Average_i_pTM: 0.44, Sequence: 'MKW', plddt_binder: 61.0]],
        ]
        def text = refoldSortKeySuggestions(suggest_entries)
        assert text.contains('Average_i_pTM')
        assert text.contains('plddt_binder')
        assert !text.contains('Sequence')
    }

    check('suggestions are capped and say how many were left out') {
        def scores = (1..20).collectEntries { i -> ["metric_${i}".toString(), i as Double] }
        def text = refoldSortKeySuggestions([[[id: 'design_1'], 'a.pdb', scores]], 5)
        assert text.contains('and 15 more')
    }

    check('suggestions cope with a set that has no numeric score at all') {
        def text = refoldSortKeySuggestions([[[id: 'design_1'], 'a.pdb', [Sequence: 'MKV']]])
        assert text.contains('nothing here can be ranked')
    }

    check('a sort key no design carries leaves every pair tied') {
        def tied_entries = [
            [[id: 'design_1'], 'a.pdb', [Average_i_pTM: 0.44]],
            [[id: 'design_2'], 'b.pdb', [Average_i_pTM: 0.81]],
        ]
        // The bare i_pTM that BindCraft never writes: ranking must not reorder, so the
        // cap falls back to input order rather than inventing one.
        def ranked = rankRefoldDesigns(tied_entries, parseRefoldSortKey('-i_pTM'), 1)
        assert ranked.size() == 1
        assert ranked[0][0].id == 'design_1'
    }

    // --- BindCraft's _model<N> suffix ---

    check('the model suffix is dropped from a BindCraft design name') {
        assert refoldStripModelSuffix('bindcraft_design_PDL1_0_l107_s306040_mpnn1_model2') ==
            'bindcraft_design_PDL1_0_l107_s306040_mpnn1'
    }

    check('a name without the suffix is left alone') {
        assert refoldStripModelSuffix('design_ppi_Tzn1xYg_1_dldesign_0_cycle3_mpnn0_af2pred') ==
            'design_ppi_Tzn1xYg_1_dldesign_0_cycle3_mpnn0_af2pred'
        // Only a trailing _model<digits> goes: _model2 in the middle is part of the name.
        assert refoldStripModelSuffix('a_model2_b') == 'a_model2_b'
        assert refoldStripModelSuffix('design_model') == 'design_model'
    }

    // --- engine list parsing ---

    // ''.split(',') is [''], and a one-element list is truthy: --refold_with unset used to
    // look like a request to refold with an engine named '' and failed validation.
    def quietLog = [warn: { String m -> }]
    [false, null, '', ' ', ',', ' , '].each { unset ->
        check("--refold_with=${unset == null ? 'null' : "'${unset}'"} means no engines") {
            assert refoldEngines([refold_with: unset], 'test', quietLog) == []
            assert !refoldEngines([refold_with: unset], 'test', quietLog)
        }
    }
    // The default in conf/refold.config is the empty string.
    check('the empty refold_with default gives no engines') {
        assert refoldEngines([refold_with: ''], 'test', quietLog) == []
        assert refoldEngines([:], 'test', quietLog) == []
    }

    check('engine lists are trimmed, lower-cased and tolerate stray commas') {
        assert refoldEngines([refold_with: 'Boltz, rf3,'], 'test', quietLog) == ['boltz', 'rf3']
    }

    check('deprecated --refold_oracles is merged without duplicates') {
        assert refoldEngines([refold_with: 'boltz', refold_oracles: 'boltz,rf3'], 'test', quietLog) == ['boltz', 'rf3']
        assert refoldEngines([refold_with: '', refold_oracles: 'rf3'], 'test', quietLog) == ['rf3']
    }

    log.info('refold_select_utils: all assertions passed')
}
