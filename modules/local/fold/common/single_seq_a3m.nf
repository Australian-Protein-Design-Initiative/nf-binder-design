// Write a query-only a3m (header + sequence) from a single-record FASTA.
// Used by fold_pulldown when --create_target_msa / --create_binder_msa is off,
// so every downstream engine still receives a valid a3m / rendered MSA.
process SINGLE_SEQ_A3M {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("${meta.id}.a3m"), emit: a3m

    script:
    """
    set -euo pipefail
    python3 - <<'PY'
from pathlib import Path

fasta = Path("${fasta}")
header = None
seq_parts = []
for line in fasta.read_text().splitlines():
    line = line.strip()
    if not line:
        continue
    if line.startswith(">"):
        if header is not None:
            break
        header = line[1:].strip() or "${meta.id}"
    else:
        seq_parts.append(line)

if header is None or not seq_parts:
    raise SystemExit(f"no FASTA record in {fasta}")

seq = "".join(seq_parts)
Path("${meta.id}.a3m").write_text(f">{header}\\n{seq}\\n")
PY
    """
}
