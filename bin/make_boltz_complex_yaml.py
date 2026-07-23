#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "pyyaml",
# ]
# ///

"""
Generic N-chain Boltz-2 YAML for the standalone fold.nf BOLTZ_FOLD subworkflow
(multimer). One `protein:` entry per input FASTA record (chain IDs A, B, C, ...
in file order), each with its own `msa:` path.

Kept separate from bin/create_boltz_yaml.py (which is hard-wired to the
target/binder two-body layout of the binder-design pipeline) so the fold.nf
complex path stays a clean N-record loop.

MSA pairing (ground truth, plans/fold-nf-multimer-paired-msa.md §0b): Boltz
dispatches the `msa:` file by extension (boltz/main.py:615) - a `.csv` with
columns exactly `key,sequence` (key == taxonomy id) is the offline-pairable
form, so bin/msa_taxonomy.py --tool boltz renders one CSV per chain and we point
each chain's `msa:` at its CSV. With --use_msa_server, `msa:` is omitted
entirely and Boltz fetches + pairs its own MSA.

Homo-oligomers: pass the sequence as repeated FASTA records (one `protein:`
entry per copy, distinct chain IDs). A `count:`/id-list shorthand is deferred.
"""

import argparse
import glob
import os
import string
import sys
from pathlib import Path
from typing import List, Optional

import yaml  # type: ignore


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


def make_boltz_complex_yaml(
    fasta_path: Path,
    msa_paths: Optional[List[Path]] = None,
    use_msa_server: bool = False,
    templates_dir: Optional[str] = None,
) -> dict:
    sequences = parse_fasta_records(fasta_path)
    if not sequences:
        raise ValueError(f"No FASTA records found in {fasta_path}")

    chain_ids = list(string.ascii_uppercase)
    if len(sequences) > len(chain_ids):
        raise ValueError(f"{fasta_path} has more than {len(chain_ids)} chains; not supported")

    if msa_paths and len(msa_paths) != len(sequences):
        raise ValueError(
            f"--msa got {len(msa_paths)} file(s) for {len(sequences)} chain(s); "
            f"pass exactly one MSA per chain in record order"
        )

    entries = []
    for i, (chain_id, seq) in enumerate(zip(chain_ids, sequences)):
        protein = {"id": [chain_id], "sequence": seq}
        # Omit msa when using the Boltz MSA server (it fetches + pairs its own).
        if not use_msa_server and msa_paths:
            # Basename so fold.nf can overwrite the staged MSA in-task without
            # rewriting this YAML.
            protein["msa"] = os.path.basename(str(msa_paths[i]))
        entries.append({"protein": protein})

    data = {"version": 1, "sequences": entries}

    if templates_dir:
        cif_files = sorted(glob.glob(os.path.join(templates_dir, "*.cif")))
        if cif_files:
            data["templates"] = [{"cif": c} for c in cif_files]

    return data


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="FASTA file (one record per chain)")
    parser.add_argument(
        "--msa",
        nargs="+",
        default=None,
        help="Per-chain MSA file(s) (.csv/.a3m) in record order (or a single file for a monomer)",
    )
    parser.add_argument("--templates", default=None, help="Optional templates directory of .cif files")
    parser.add_argument("--use_msa_server", action="store_true", help="Omit msa: so Boltz fetches its own")
    parser.add_argument("--output_yaml", required=True, help="Output YAML path")
    args = parser.parse_args()

    data = make_boltz_complex_yaml(
        fasta_path=Path(args.fasta),
        msa_paths=[Path(p) for p in args.msa] if args.msa else None,
        use_msa_server=args.use_msa_server,
        templates_dir=args.templates,
    )

    out = Path(args.output_yaml)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        yaml.dump(data, f, sort_keys=False)
    print(f"Wrote Boltz complex YAML ({len(data['sequences'])} chain(s)) to {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
