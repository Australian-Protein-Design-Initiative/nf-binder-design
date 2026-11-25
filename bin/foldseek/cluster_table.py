#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Turn foldseek's cluster output into a per-design table.

`foldseek easy-cluster` writes a two-column TSV of representative<TAB>member, with one
row per design. This rewrites it as one row per design keyed on the design id, so it can
be joined onto combined_scores.tsv with bin/merge_scores.py --keys id.

The id is the structure's basename with its extension dropped, which is what foldseek
uses as the entry name, and what the design tables already call `id`.
"""

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

# Structure suffixes foldseek may leave on an entry name
SUFFIXES = (".pdb.gz", ".cif.gz", ".pdb", ".cif", ".bcif", ".ent")


def strip_suffix(name: str) -> str:
    for suffix in SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def read_clusters(path: Path) -> List[Tuple[str, str]]:
    """(representative, member) pairs, ids normalised."""
    pairs = []
    with open(path, newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 2 or not row[0].strip():
                continue
            pairs.append((strip_suffix(row[0].strip()), strip_suffix(row[1].strip())))
    return pairs


def build_rows(pairs: List[Tuple[str, str]]) -> List[Dict[str, str]]:
    sizes = Counter(rep for rep, _ in pairs)
    # Rank clusters largest first so cluster1 is the most populated, which is the
    # ordering a diversity summary is usually read in
    order = {
        rep: index + 1
        for index, (rep, _) in enumerate(sorted(sizes.items(), key=lambda kv: (-kv[1], kv[0])))
    }
    return [
        {
            "id": member,
            "foldseek_cluster_rep": rep,
            "foldseek_cluster_id": str(order[rep]),
            "foldseek_cluster_size": str(sizes[rep]),
            "foldseek_cluster_is_rep": "1" if member == rep else "0",
        }
        for rep, member in sorted(pairs, key=lambda p: (order[p[0]], p[1]))
    ]


FIELDS = [
    "id",
    "foldseek_cluster_rep",
    "foldseek_cluster_id",
    "foldseek_cluster_size",
    "foldseek_cluster_is_rep",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rewrite a foldseek cluster TSV as a per-design table", allow_abbrev=False
    )
    parser.add_argument("clusters", help="foldseek *_cluster.tsv (representative, member)")
    parser.add_argument("-o", "--output", default="-", help="Output TSV, or - for stdout")
    args = parser.parse_args()

    rows = build_rows(read_clusters(Path(args.clusters)))
    handle = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    try:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if handle is not sys.stdout:
            handle.close()

    clusters = len({r["foldseek_cluster_rep"] for r in rows})
    log.info("%d designs in %d structural clusters", len(rows), clusters)


if __name__ == "__main__":
    main()
