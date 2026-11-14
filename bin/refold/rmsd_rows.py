#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Turn the rmsd4all tables for one refolded design into long rows keyed by
(id, tool, batch, model).

rmsd4all writes one row per (fixed, mobile) structure pair, with its own column
names. A refold run compares one design against every model the engine predicted
for it, and the comparisons are run separately (superimposed on different chains),
so this joins them on the mobile structure and labels the columns by what was
compared rather than by which rmsd4all run produced them:

    id        tool   batch  model  target_aligned_binder_rmsd_pruned  complex_rmsd_all
    design_0  boltz  1      0      0.94                               1.21

The model index comes from the mobile filename, which every engine ends with a
sample or model number (see FoldNaming.isModelStructure). Where it cannot be read,
the rows are numbered in filename order so a design with several samples still
gives distinct rows.

The batch is the fold job the model came from, and is needed as well as the model
index because the refold runs --n_predictions as one job per seed: each of those
jobs writes its single sample under the same model index, so (id, tool, model)
alone would collapse every seed's sample onto one key.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from typing import Dict, List

# rmsd4all column -> suffix in the output. Anything else in the table is dropped:
# pair counts and residue counts describe the alignment, not the design.
VALUE_COLUMNS = {
    "rmsd_pruned": "rmsd_pruned",
    "rmsd_all": "rmsd_all",
    "tm_score_pruned": "tm_score_pruned",
    "tm_score_all": "tm_score_all",
}

MODEL_PATTERNS = (
    re.compile(r"_model_(\d+)\."),      # boltz
    re.compile(r"_sample[_-](\d+)[._]"),  # protenix, rf3, af3, openfold3, chai1, esmfold2
)


def model_index(filename: str, fallback: int) -> str:
    for pattern in MODEL_PATTERNS:
        match = pattern.search(filename)
        if match:
            return match.group(1)
    return str(fallback)


def read_table(path: str) -> Dict[str, dict]:
    with open(path, newline="") as handle:
        return {row["structure2"]: row for row in csv.DictReader(handle, delimiter="\t")}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--id", required=True, help="Design id")
    parser.add_argument("--tool", required=True, help="Folding engine that produced the models")
    parser.add_argument(
        "--batch", default="1",
        help="Fold job the models came from, written out as the batch column [default: 1]",
    )
    parser.add_argument(
        "--model",
        help="Sample label for every row, as the engine's own score table writes it "
             "(FoldNaming.sampleLabel). Given explicitly, the label is not guessed from "
             "the filename, which keeps the two tables joinable.",
    )
    parser.add_argument(
        "--table",
        action="append",
        default=[],
        metavar="PREFIX=PATH",
        help="An rmsd4all TSV and the prefix its columns take, eg "
             "target_aligned_binder=rmsd_target_aligned.tsv. Repeatable.",
    )
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    tables = []
    for spec in args.table:
        if "=" not in spec:
            print(f"ERROR: --table needs PREFIX=PATH, got '{spec}'", file=sys.stderr)
            return 1
        prefix, path = spec.split("=", 1)
        tables.append((prefix, read_table(path)))

    if not tables:
        print("ERROR: at least one --table is required", file=sys.stderr)
        return 1

    # Every table compares the same set of models, so the first one fixes the row
    # order; a model missing from a later table leaves its columns empty rather
    # than dropping the row.
    models = sorted(tables[0][1])

    columns = ["id", "tool", "batch", "model"]
    for prefix, _ in tables:
        columns += [f"{prefix}_{suffix}" for suffix in VALUE_COLUMNS.values()]

    rows: List[dict] = []
    for position, model in enumerate(models):
        row = {
            "id": args.id,
            "tool": args.tool,
            "batch": args.batch,
            "model": args.model if args.model else model_index(model, position),
        }
        for prefix, table in tables:
            source = table.get(model, {})
            for column, suffix in VALUE_COLUMNS.items():
                row[f"{prefix}_{suffix}"] = source.get(column, "")
        rows.append(row)

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
