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
form, so bin/fold/msa_taxonomy.py --tool boltz renders one CSV per chain and we point
each chain's `msa:` at its CSV. With --use_msa_server, `msa:` is omitted
entirely and Boltz fetches + pairs its own MSA.

Homo-oligomers: pass the sequence as repeated FASTA records (one `protein:`
entry per copy, distinct chain IDs). A `count:`/id-list shorthand is deferred.

Templates: --templates is the directory written by bin/fold/match_templates.py.
Each chain's templates are looked up in its index.json by sequence md5 and pinned
to that chain (`chain_id`), so Boltz never assigns a target template to a binder.
"""

import argparse
import csv
import hashlib
import json
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


def csv_query(csv_path: Path) -> Optional[str]:
    """Row 0's sequence column - Boltz always treats a chain MSA's first row as
    the query (render_boltz_csv in bin/fold/msa_taxonomy.py), so this is the
    sequence Boltz will fold against for that chain.
    """
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            return row.get("sequence")
    return None


def _check_query(path: Path, seq: str, chain_index: int) -> None:
    # A mis-ordered bundle otherwise folds each chain with another chain's MSA
    # without any error from Boltz (this is the bug class fixed in commit 1ae37a8).
    if path.suffix.lower() != ".csv":
        return
    query = csv_query(path)
    if query is not None and query.upper() != seq.upper():
        raise ValueError(
            f"chain {chain_index}: first sequence in {path.name} does not match the FASTA "
            f"record - MSA files are not in chain order"
        )


def make_boltz_complex_yaml(
    fasta_path: Path,
    msa_paths: Optional[List[Path]] = None,
    use_msa_server: bool = False,
    templates_dir: Optional[str] = None,
    query_only_chains: Optional[List[str]] = None,
    template_chains: Optional[List[str]] = None,
    template_force: bool = False,
    template_threshold: float = 1.0,
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

    query_only = set(query_only_chains or [])

    entries = []
    for i, (chain_id, seq) in enumerate(zip(chain_ids, sequences)):
        protein = {"id": [chain_id], "sequence": seq}
        if chain_id in query_only:
            # Keep this chain query-only even under --use_msa_server, e.g. a
            # binder chain the caller does not want the MSA server to fetch
            # for (--create_binder_msa false).
            protein["msa"] = "empty"
        elif not use_msa_server and msa_paths:
            _check_query(msa_paths[i], seq, i)
            # Basename so fold.nf can overwrite the staged MSA in-task without
            # rewriting this YAML.
            protein["msa"] = os.path.basename(str(msa_paths[i]))
        entries.append({"protein": protein})

    data = {"version": 1, "sequences": entries}

    if templates_dir:
        templates = chain_templates(
            templates_dir, chain_ids[: len(sequences)], sequences,
            template_chains, template_force, template_threshold,
        )
        if templates:
            data["templates"] = templates

    return data


def chain_templates(
    templates_dir: str,
    chain_ids: List[str],
    sequences: List[str],
    template_chains: Optional[List[str]],
    force: bool,
    threshold: float,
) -> List[dict]:
    index_path = Path(templates_dir) / "index.json"
    if not index_path.exists():
        return []
    index = json.loads(index_path.read_text())
    allowed = set(template_chains) if template_chains else set(chain_ids)
    entries = []
    for chain_id, seq in zip(chain_ids, sequences):
        if chain_id not in allowed:
            continue
        for hit in index.get(hashlib.md5(seq.upper().encode()).hexdigest(), []):
            entry = {
                "cif": os.path.join(templates_dir, hit["cif"]),
                "chain_id": chain_id,
                "template_id": hit["template_chain"],
            }
            if force:
                entry["force"] = True
                entry["threshold"] = threshold
            entries.append(entry)
    return entries


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="FASTA file (one record per chain)")
    parser.add_argument(
        "--msa",
        nargs="+",
        default=None,
        help="Per-chain MSA file(s) (.csv/.a3m) in record order (or a single file for a monomer)",
    )
    parser.add_argument("--templates", default=None, help="Optional matched-templates directory (bin/fold/match_templates.py)")
    parser.add_argument(
        "--template_chains",
        nargs="+",
        default=None,
        help="Chain ID(s) that may be templated (default: all)",
    )
    parser.add_argument("--template_force", action="store_true", help="Hold templated chains close to the template (Boltz force)")
    parser.add_argument("--template_threshold", type=float, default=1.0, help="Boltz force threshold in Angstrom (default: 1.0)")
    parser.add_argument("--use_msa_server", action="store_true", help="Omit msa: so Boltz fetches its own")
    parser.add_argument(
        "--query_only_chains",
        nargs="+",
        default=None,
        help="Chain ID(s) (e.g. B) to force msa: empty for, even under --use_msa_server",
    )
    parser.add_argument("--output_yaml", required=True, help="Output YAML path")
    args = parser.parse_args()

    data = make_boltz_complex_yaml(
        fasta_path=Path(args.fasta),
        msa_paths=[Path(p) for p in args.msa] if args.msa else None,
        use_msa_server=args.use_msa_server,
        templates_dir=args.templates,
        query_only_chains=args.query_only_chains,
        template_chains=args.template_chains,
        template_force=args.template_force,
        template_threshold=args.template_threshold,
    )

    out = Path(args.output_yaml)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        yaml.dump(data, f, sort_keys=False)
    print(f"Wrote Boltz complex YAML ({len(data['sequences'])} chain(s)) to {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
