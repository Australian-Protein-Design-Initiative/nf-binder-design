#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Generic RF3 (RosettaFold3) input JSON for the standalone fold.nf
ROSETTAFOLD3_FOLD subworkflow.

NOT to be confused with bin/rfd3/make_rf3_input_spec.py, which is hard-wired
to a 2-body target+binder complex for the binder-design pipeline (rfd3.nf) and
always emits target+binder roles plus optional template handling.

Builds one {seq, chain_id, msa_path?} component per input FASTA record (chain
IDs A, B, C, ... in file order), with no target/binder roles, no template, and
no binder postprocessing - this is exactly what rf3's generic
InferenceInput.from_json_dict() / components_to_atom_array() accepts (ground
truth confirmed by inspecting rc-foundry:0.2.0-weights' rf3/utils/inference.py
on 2026-07-17: it builds an atom array generically from `components`, with no
target/binder-specific logic at all).

Phase 1 scope: fold.nf rejects multi-record FASTA before this script ever
runs (see fold.nf's monomer-only guard), so today this always emits a
single-chain (chain A) component list. The multi-chain per-record logic below
is deliberately generic so multimer support (Phase 2) only needs a caller-side
change, not a rewrite of this script.
"""

import argparse
import json
import logging
import string
import sys
from pathlib import Path
from typing import List, Optional

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def parse_fasta_records(fasta_path: Path) -> List[str]:
    """Return one sequence string per FASTA record, in file order (no external deps)."""
    records: List[str] = []
    current: List[str] = []
    for line in fasta_path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                records.append("".join(current))
                current = []
        else:
            current.append(line)
    if current:
        records.append("".join(current))
    return records


def make_rf3_fold_spec(fasta_path: Path, name: str, a3m_path: Optional[Path] = None) -> dict:
    sequences = parse_fasta_records(fasta_path)
    if not sequences:
        raise ValueError(f"No FASTA records found in {fasta_path}")

    chain_ids = list(string.ascii_uppercase)
    if len(sequences) > len(chain_ids):
        raise ValueError(f"{fasta_path} has more than {len(chain_ids)} chains; not supported")

    components = []
    for chain_id, seq in zip(chain_ids, sequences):
        comp = {"seq": seq, "chain_id": chain_id}
        # Phase 1 is monomer-only, so a single provided a3m always belongs to
        # the sole chain. Multimer per-chain MSA routing is Phase 2 work (see
        # module docstring).
        if a3m_path is not None and len(sequences) == 1:
            comp["msa_path"] = str(a3m_path.resolve())
        components.append(comp)

    return {"name": name, "components": components}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="FASTA file (one record per chain)")
    parser.add_argument("--name", required=True, help="RF3 'name' field")
    parser.add_argument("--a3m", default=None, help="Optional a3m MSA (monomer only in Phase 1)")
    parser.add_argument("-o", "--output", required=True, help="Output JSON path")
    args = parser.parse_args()

    spec = make_rf3_fold_spec(
        fasta_path=Path(args.fasta),
        name=args.name,
        a3m_path=Path(args.a3m) if args.a3m else None,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # rf3 fold's `inputs=` argument is a JSON list of examples (see
    # bin/rfd3/make_rf3_input_spec.py's json-batch mode for the same
    # convention), even for a single example.
    with open(out_path, "w") as f:
        json.dump([spec], f, indent=2)
    log.info("Wrote RF3 fold input JSON (1 example, %d chain(s)) to %s", len(spec["components"]), out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
