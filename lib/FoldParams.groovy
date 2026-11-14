/**
 * Parameter lookup for the fold stack's engine flags.
 *
 * The fold stack doubles as the refolding stack for binder-design workflows, and some
 * of those already declare parameters with the same names and different meanings:
 * rfd3's --rf3_batch_size is MPNN designs per RF3 task, not samples per fold job. A
 * workflow in that position declares the refold's copies with a full_refold_ prefix
 * instead (see workflows/rfd4_proteina.nf).
 *
 * Where the prefixed name is declared it is the only one read, so a design-time
 * parameter of the same name cannot leak into the refold - including when the
 * prefixed one is left at a falsy default.
 */
class FoldParams {

    /** The name in effect for an engine parameter: prefixed where that one is declared. */
    static String nameFor(Map params, String param) {
        def prefixed = "full_refold_${param}".toString()
        return params.containsKey(prefixed) ? prefixed : param
    }

    /** Its value, or null when neither name is declared (no undefined-parameter warning). */
    static Object engine(Map params, String param) {
        def key = nameFor(params, param)
        return params.containsKey(key) ? params[key] : null
    }
}
