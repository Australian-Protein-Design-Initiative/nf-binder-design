import java.io.File
import java.nio.file.Path

/*
Shared workflow.onComplete params.json helper, de-duplicated from the
identical paramsToMap() copies previously living in main.nf and af2.nf (see
plans/fold-nf-multi-method-folding.md). fold.nf and the refactored af2.nf both
use this.
*/
class ParamsHelper {
    static def paramsToMap(params) {
        def map = [:]
        params.each { key, value ->
            if (value instanceof Path || value instanceof File) {
                map[key] = value.toString()
            }
            else if (!(value instanceof Closure) && !(key in [
                'class',
                'launchDir',
                'projectDir',
                'workDir',
            ])) {
                map[key] = value
            }
        }
        return map
    }
}
