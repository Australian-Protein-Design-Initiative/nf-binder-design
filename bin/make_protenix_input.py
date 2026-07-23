#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Generic Protenix (AF3-style) input JSON for the standalone fold.nf
PROTENIX_FOLD subworkflow.

Builds one {"proteinChain": {"sequence": ..., "count": 1}} entry per input
FASTA record, wrapped in the top-level {"name": ..., "sequences": [...]}
list-of-jobs format `protenix pred -i` expects (ground truth confirmed by
inspecting protenix:v2.0.0-weights' protenix/data/inference/json_parser.py
and runner/inference.py on 2026-07-17: `sequences` entries are keyed by
entity type - "proteinChain"/"dnaSequence"/"rnaSequence"/"ligand"/"ion" - and
the top-level "name" field becomes the output sample_name/directory).

Precomputed MSA (monomer contract, confirmed against
runner/msa_search.py's need_msa_search() and
protenix/data/msa/msa_featurizer.py on 2026-07-17): a per-chain a3m is fed via
the proteinChain entry's "unpairedMsaPath" field (a plain a3m file path,
query sequence first) - this is the *current* field name; the older
`{"msa": {"precomputed_msa_dir": ...}}` form (with `pairing.a3m`/
`non_pairing.a3m` inside) still works but only logs a deprecation warning, so
we always emit the new field. If present and the path exists, Protenix skips
its own MSA search entirely (need_msa_search() returns False), matching how
boltz/rf3 consume our shared a3m via `msa:`/`msa_path`.

Multimer: pass one --a3m (unpairedMsaPath) per chain in record order, and one
--paired-a3m (pairedMsaPath) per chain. Protenix pairs chains internally by
species *mnemonic* (_HUMAN, _9BETA via MSAPairingEngine.get_species_ids), NOT
numeric TaxID=, so the paired a3m headers must carry the mnemonic
(bin/msa_taxonomy.py --tool protenix renders both files). A single --a3m with a
single-record FASTA is the monomer case (unpaired only; single chains don't
pair).
"""

import argparse
import json
import logging
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


def _match_per_chain(paths: Optional[List[Path]], n_seq: int, flag: str) -> Optional[List[Path]]:
    if not paths:
        return None
    if len(paths) != n_seq:
        raise ValueError(
            f"{flag} got {len(paths)} file(s) for {n_seq} chain(s); "
            f"pass exactly one per chain in record order"
        )
    return paths


def make_protenix_input(
    fasta_path: Path,
    name: str,
    unpaired_a3m_paths: Optional[List[Path]] = None,
    paired_a3m_paths: Optional[List[Path]] = None,
) -> dict:
    sequences = parse_fasta_records(fasta_path)
    if not sequences:
        raise ValueError(f"No FASTA records found in {fasta_path}")

    n = len(sequences)
    unpaired = _match_per_chain(unpaired_a3m_paths, n, "--a3m")
    paired = _match_per_chain(paired_a3m_paths, n, "--paired-a3m")

    entries = []
    for i, seq in enumerate(sequences):
        protein_chain = {"sequence": seq, "count": 1}
        # Basename so fold.nf can overwrite the staged a3m in-task (MSA
        # subsample) without rewriting this JSON. a3m paths match chains by
        # position (record order); a single file only for a monomer.
        if unpaired:
            protein_chain["unpairedMsaPath"] = unpaired[i].name
        if paired:
            protein_chain["pairedMsaPath"] = paired[i].name
        entries.append({"proteinChain": protein_chain})

    return {"name": name, "sequences": entries}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="FASTA file (one record per chain)")
    parser.add_argument("--name", required=True, help="Protenix 'name' field (also the output sample_name)")
    parser.add_argument(
        "--a3m",
        nargs="+",
        default=None,
        help="Optional unpairedMsaPath a3m(s): one per chain in record order (or single for a monomer)",
    )
    parser.add_argument(
        "--paired-a3m",
        dest="paired_a3m",
        nargs="+",
        default=None,
        help="Optional pairedMsaPath a3m(s) for multimer: one per chain in record order",
    )
    parser.add_argument("-o", "--output", required=True, help="Output JSON path")
    args = parser.parse_args()

    spec = make_protenix_input(
        fasta_path=Path(args.fasta),
        name=args.name,
        unpaired_a3m_paths=[Path(p) for p in args.a3m] if args.a3m else None,
        paired_a3m_paths=[Path(p) for p in args.paired_a3m] if args.paired_a3m else None,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # protenix pred's `-i` input is a JSON list of jobs, even for a single job
    # (see runner/inference.py's infer_predict: `if not isinstance(json_data, list)...`).
    with open(out_path, "w") as f:
        json.dump([spec], f, indent=2)
    log.info("Wrote Protenix input JSON (1 job, %d chain(s)) to %s", len(spec["sequences"]), out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
