#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Pivot refold scoring into one row per design.

Everything the refold produces is long - one row per (target, binder, tool) for the
confidence summary, one row per (design, tool, model) for the RMSD tables, one row
per refolded complex for the interface scores - while combined_scores.tsv is wide,
one row per design. Each engine's values become their own columns, prefixed with the
engine name:

    id        boltz_iptm_max  boltz_target_aligned_binder_rmsd_all  esmfold2_iptm_max
    rfd4_x_0  0.78            1.75                                  0.71

consensus_z, consensus_z_metric and n_tools describe the (target, binder) pair rather
than one engine, so they are written once per design rather than prefixed.

Several models per design (--n_predictions > 1) are reduced to one value per engine.
Confidence metrics arrive already aggregated in the summary, where iptm_max and
ipsae_max are the best of the samples. The structural metrics - RMSD against the
design, and the BindCraft interface scores - are not aggregated: they all describe one
set of coordinates, so mixing samples would produce a row that no predicted structure
matches. They are read from the single sample --samples selects, which is the one the
confidence ranking already picked out as the best. --structure-stat min / mean / max /
first override that with a statistic over the samples, applied the same way to every
structural column, which suits RMSD and does not suit the interface metrics.

Every sample's own values are kept, unreduced, in the per-sample table that
bin/refold/sample_scores.py writes.

Only one target is expected, since every design in a run shares it. If several
appear, the target id is prefixed onto the column names so the rows stay distinct.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from typing import Dict, List, Optional

# Per-pair, not per-tool: written once per design rather than prefixed per tool.
PAIR_COLUMNS = ("consensus_z", "consensus_z_metric", "n_tools")

DEFAULT_METRICS = (
    "iptm_max", "iptm_mean", "ipsae_max", "ipsae_mean",
    "plddt_max", "plddt_mean", "ptm_max", "pae_mean", "n",
)

# Columns in the long tables that identify the row rather than measure anything.
KEY_COLUMNS = ("id", "tool", "batch", "model", "filename")


def _f(value: str) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def reduce_values(values: List[str], stat: str) -> str:
    """Reduce one metric's values across a design's models to a single cell.

    Non-numeric values (an interface residue list, say) have no statistic, so the
    first non-empty one is kept.
    """
    numbers = [v for v in (_f(v) for v in values) if v is not None]
    if not numbers:
        for value in values:
            if value not in ("", None):
                return value
        return ""
    if stat == "first":
        return f"{numbers[0]:.6g}"
    if stat == "mean":
        return f"{statistics.mean(numbers):.6g}"
    if stat == "max":
        return f"{max(numbers):.6g}"
    return f"{min(numbers):.6g}"


def read_long(path: str) -> List[dict]:
    with open(path, newline="") as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t")]
    return rows


def _parse_interface_filename(filename: str):
    """(design, tool, batch, model) from the PDB name REFOLD_STRUCTURE_TO_PDB wrote.

    The current shape is <design>__<tool>__b<batch>__<model>; the two earlier ones,
    <design>__<tool>__b<batch> and <design>__<tool>, are still read so a score table
    from an older run can be merged.
    """
    name = (filename or "").rsplit(".", 1)[0]
    parts = name.split("__")
    if len(parts) == 4 and parts[2].startswith("b"):
        return parts[0], parts[1], parts[2][1:], parts[3]
    if len(parts) == 3 and parts[2].startswith("b"):
        return parts[0], parts[1], parts[2][1:], ""
    if len(parts) == 2:
        return parts[0], parts[1], "1", ""
    return None


def rows_from_bindcraft(path: str) -> List[dict]:
    """Interface score rows, attributed by the PDB filename they scored."""
    out = []
    for row in read_long(path):
        parsed = _parse_interface_filename(row.get("filename"))
        if parsed is None:
            print(
                f"ERROR: interface score row '{row.get('filename')}' is not named "
                "<design>__<tool>__b<batch>__<model>.pdb, so it cannot be attributed "
                "to a design.",
                file=sys.stderr,
            )
            raise SystemExit(1)
        design, tool, batch, model = parsed
        new = {k: v for k, v in row.items() if k != "filename"}
        new["id"] = design
        new["tool"] = tool
        new["batch"] = batch
        new["model"] = model
        out.append(new)
    return out


def _sample_order(row: dict) -> tuple:
    """Sort key putting a design's samples in (batch, model) order, so `first` means
    the first seed's sample rather than whichever row the collectFile happened to
    write first."""
    def as_int(value: str) -> int:
        try:
            return int(value)
        except (TypeError, ValueError):
            return 0

    return (as_int(row.get("batch", "")), as_int(row.get("model", "")))


def selected_samples(path: str, rank_metric: str, fallback_metric: str = "iptm") -> Dict[tuple, tuple]:
    """The (batch, model) of the best-ranked sample of each (design, engine), read
    from the per-sample table.

    Both parts are needed. An engine run once per seed writes one sample per batch,
    so the batch alone identifies it; an engine handed the whole seed list returns
    every sample from batch 1, and only the model label tells them apart.

    Best means the highest rank metric, the same quantity the confidence summary
    reports as <metric>_max and ranks designs by. Picking one sample, rather than a
    statistic per metric, is what keeps the structural columns describing a single
    set of coordinates. Ties and a missing metric both fall back to the first sample
    in (batch, model) order, so the choice is always defined and always a real
    sample.
    """
    grouped: Dict[tuple, List[dict]] = {}
    for row in read_long(path):
        grouped.setdefault((row.get("id", ""), row.get("tool", "")), []).append(row)

    chosen: Dict[tuple, tuple] = {}
    for key, group in grouped.items():
        group.sort(key=_sample_order)
        metric = rank_metric if any(_f(r.get(rank_metric, "")) is not None for r in group) \
            else fallback_metric
        best = group[0]
        best_value = _f(best.get(metric, ""))
        for row in group[1:]:
            value = _f(row.get(metric, ""))
            if value is not None and (best_value is None or value > best_value):
                best, best_value = row, value
        chosen[key] = (best.get("batch", "1"), best.get("model", ""))
    return chosen


def add_long_rows(
    by_design: Dict[str, dict],
    columns: List[str],
    rows: List[dict],
    stat: str,
    selected: Optional[Dict[tuple, str]] = None,
) -> None:
    grouped: Dict[tuple, List[dict]] = {}
    for row in rows:
        grouped.setdefault((row["id"], row["tool"]), []).append(row)

    for (design, tool), group in grouped.items():
        group.sort(key=_sample_order)
        row_stat = stat
        if stat == "selected":
            # Every structural column then comes from one sample's coordinates. Where
            # no selection was given for this (design, engine) - no per-sample table,
            # or no row for it in one - the first sample stands in, which is still a
            # single real structure rather than a mixture.
            choice = (selected or {}).get((design, tool))
            if choice is not None:
                batch, model = choice
                in_batch = [r for r in group if r.get("batch", "1") == batch]
                # On the model too where the table carries one. A table written before
                # the sample label was recorded - an older interface score table, say -
                # has none, and then the batch is as far as the match can go.
                picked = [r for r in in_batch if r.get("model", "") == model] if model else []
                group = picked or in_batch or group
            row_stat = "first"
        out = by_design.setdefault(design, {"id": design})
        metrics = [c for c in group[0] if c not in KEY_COLUMNS]
        for metric in metrics:
            column = f"{tool}_{metric}"
            out[column] = reduce_values([r.get(metric, "") for r in group], row_stat)
            if column not in columns:
                columns.append(column)


def pivot(summary_rows: List[dict], metrics: List[str]):
    targets = {r["target"] for r in summary_rows}
    multi_target = len(targets) > 1

    by_design: Dict[str, dict] = {}
    columns: List[str] = []

    for row in summary_rows:
        design = row["binder"]
        out = by_design.setdefault(design, {"id": design})

        prefix = f"{row['target']}_{row['tool']}" if multi_target else row["tool"]
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

    return by_design, columns


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("summary", help="fold_pulldown_summary.tsv from the refold")
    parser.add_argument("-o", "--output", default="refold_scores.tsv")
    parser.add_argument(
        "--metrics",
        default=",".join(DEFAULT_METRICS),
        help=f"Comma-separated summary columns to carry over [default: {','.join(DEFAULT_METRICS)}]",
    )
    parser.add_argument(
        "--long",
        action="append",
        default=[],
        metavar="PATH",
        help="A long table with id and tool columns (the RMSD tables). Repeatable; "
             "empty or header-only files are skipped.",
    )
    parser.add_argument(
        "--bindcraft",
        action="append",
        default=[],
        metavar="PATH",
        help="Interface scores whose filename column is <design>__<tool>.pdb. Repeatable.",
    )
    parser.add_argument(
        "--samples",
        metavar="PATH",
        help="refold_sample_scores.tsv, read to pick the sample the structural columns "
             "come from under --structure-stat selected.",
    )
    parser.add_argument(
        "--rank-metric",
        default="ipsae",
        help="Per-sample column the selection ranks on, falling back to iptm where a "
             "(design, engine) has none of it [default: ipsae]",
    )
    parser.add_argument(
        "--structure-stat",
        choices=("selected", "min", "first", "mean", "max"),
        default="selected",
        help="How several samples per design reduce to one value in the RMSD and "
             "interface columns [default: selected, every column read from the "
             "best-ranked sample]. The structural metrics all describe one set of "
             "coordinates, and they have no shared direction - interface_dG is better "
             "lower, interface_dSASA and interface_sc better higher - so a statistic "
             "per metric would describe a structure that was never predicted. The "
             "other choices apply one statistic to every structural column, which "
             "suits RMSD and does not suit the interface metrics.",
    )
    args = parser.parse_args()

    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    summary_rows = read_long(args.summary)

    if not summary_rows:
        # An empty summary still needs a readable header, so a downstream merge that
        # globs for this file does not fail on a zero-byte input.
        with open(args.output, "w", newline="") as handle:
            handle.write("id\n")
        return 0

    by_design, columns = pivot(summary_rows, metrics)

    selected = selected_samples(args.samples, args.rank_metric) if args.samples else {}

    for path in args.long:
        rows = read_long(path)
        if rows:
            add_long_rows(by_design, columns, rows, args.structure_stat, selected)
    for path in args.bindcraft:
        rows = rows_from_bindcraft(path)
        if rows:
            add_long_rows(by_design, columns, rows, args.structure_stat, selected)

    # Which sample the structural columns describe, so a row can be traced back to a
    # structure on disk and to its own line in the per-sample table.
    if args.structure_stat == "selected":
        for (design, tool), (batch, model) in selected.items():
            if design not in by_design:
                continue
            for suffix, value in (("selected_batch", batch), ("selected_model", model)):
                column = f"{tool}_{suffix}"
                by_design[design][column] = value
                if column not in columns:
                    columns.append(column)

    ordered = ["id"] + sorted(c for c in columns if c not in PAIR_COLUMNS)
    ordered += [c for c in PAIR_COLUMNS if c in columns]

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for design in sorted(by_design.values(), key=lambda r: r["id"]):
            writer.writerow({c: design.get(c, "") for c in ordered})
    return 0


if __name__ == "__main__":
    sys.exit(main())
