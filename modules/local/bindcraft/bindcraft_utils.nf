// Utility functions for BindCraft workflows

// --input_pdb may be a single file, a directory of PDBs, or a glob
def resolveInputPdbs(input_pdb) {
    def pattern = input_pdb.toString()
    if (!(pattern =~ /[*?\[]/) && file(pattern).isDirectory()) {
        pattern = "${pattern.replaceAll(/\/+$/, '')}/*.pdb"
    }
    def pdbs = file(pattern)
    if (!(pdbs instanceof List)) {
        pdbs = pdbs.exists() ? [pdbs] : []
    }
    if (pdbs.isEmpty()) {
        error "No PDB files found for --input_pdb '${input_pdb}' (resolved pattern: '${pattern}')"
    }
    def dupes = pdbs.groupBy { it.simpleName }.findAll { _k, v -> v.size() > 1 }.keySet()
    if (dupes) {
        error "Duplicate input PDB basenames are not supported (batch/design names would collide): ${dupes.join(', ')}"
    }
    return pdbs.sort { it.name }
}

def makeBindcraftBatchId(pdb, index) {
    def key = pdb.simpleName.replaceAll(/[^A-Za-z0-9]+/, '_')
    return "${key}_${index}"
}
