#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Pivot fold_pulldown_summary.tsv into one row per design.

The summary is long: one row per (target, binder, tool). A design workflow's
combined_scores.tsv is wide, one row per design, so each tool's statistics become
their own columns, prefixed with the tool name:

    id            chai1_iptm_max  chai1_ipsae_max  boltz_iptm_max  consensus_z  n_tools
    rfd4_x_0      0.81            0.62             0.78            1.42         2

consensus_z, consensus_z_metric and n_tools describe the (target, binder) pair
rather than one tool, so they are written once per design rather than prefixed.

Only one target is expected, since every design in a run shares it. If several
appear, the target id is prefixed onto the column names so the rows stay distinct.
"""

from __future__ import annotations

import argparse
import csv
import sys
from typing import Dict, List

# Per-pair, not per-tool: written once per design rather than prefixed per tool.
PAIR_COLUMNS = ("consensus_z", "consensus_z_metric", "n_tools")

DEFAULT_METRICS = ("iptm_max", "iptm_mean", "ipsae_max", "ipsae_mean", "n")


def pivot(rows: List[dict], metrics: List[str]) -> List[dict]:
    targets = {r["target"] for r in rows}
    multi_target = len(targets) > 1

    by_design: Dict[str, dict] = {}
    columns: List[str] = []

    for row in rows:
        design = row["binder"]
        out = by_design.setdefault(design, {"id": design})

        tool = row["tool"]
        prefix = f"{row['target']}_{tool}" if multi_target else tool
        for metric in metrics:
            if metric not in row:
                continue
            column = f"{prefix}_{metric}"
            out[column] = row[metric]
            if column not in columns:
                columns.append(column)

        for column in PAIR_COLUMNS:
            if column in row and not out.get(column):
                out[column] = row[column]
                if column not in columns:
                    columns.append(column)

    ordered = ["id"] + sorted(c for c in columns if c not in PAIR_COLUMNS)
    ordered += [c for c in PAIR_COLUMNS if c in columns]
    return [{c: design.get(c, "") for c in ordered} for design in by_design.values()], ordered


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", help="fold_pulldown_summary.tsv")
    parser.add_argument("-o", "--output", default="refold_oracle_scores.tsv")
    parser.add_argument(
        "--metrics",
        default=",".join(DEFAULT_METRICS),
        help=f"Comma-separated summary columns to carry over [default: {','.join(DEFAULT_METRICS)}]",
    )
    args = parser.parse_args()

    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]

    with open(args.summary, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        # An empty summary still needs a readable header, so a downstream merge
        # that globs for this file does not fail on a zero-byte input.
        with open(args.output, "w", newline="") as handle:
            handle.write("id\n")
        return 0

    pivoted, columns = pivot(rows, metrics)

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(pivoted, key=lambda r: r["id"]):
            writer.writerow(row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
