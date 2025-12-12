process BOLTZGEN_CHECK {

    container 'ghcr.io/australian-protein-design-initiative/containers/boltzgen:0.2.0'

    input:
    path 'input/*'
    val design_name

    output:
    path "${design_name}.cif"

    script:
    """
    set -euo pipefail

    boltzgen check ${design_name}.yaml
    """
}
