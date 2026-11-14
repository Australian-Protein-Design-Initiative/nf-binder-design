#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
One row per refolded sample, with every score that sample has.

combined_scores.tsv is one row per design, which means every number in it has been
reduced across the samples the engines produced: confidence metrics as aggregates over
the samples, and the structural metrics - RMSD and interface scores - read from the
single best-ranked sample, so that they all describe one set of coordinates. That
reduction is what makes the table rankable, and it is also what throws away the
spread. This table keeps the unreduced values, so any summary statistic - min, max,
mean, median, standard deviation - can be computed afterwards from the samples
themselves. It is also what the wide table reads to make that choice.

    id        tool   batch  model  iptm    ipsae   complex_rmsd_all  interface_dG
    design_0  boltz  1      0      0.9590  0.8965  1.27              -43.55
    design_0  boltz  2      0      0.9123  0.8411  1.44              -41.02

The key is (id, tool, batch, model): the refold runs --n_predictions as one fold job
per seed, and every job writes its sample under the same model index, so the batch is
what tells one seed's sample from another's.

Inputs are the long tables the refold already writes:

  fold_pulldown_scores.tsv  confidence, one row per predicted structure
  refold_rmsd.tsv           pose and complex RMSD, one row per sample  (--long)
  refold_monomer_rmsd.tsv   binder-alone RMSD, one row per sample      (--long)
  refold_interface_scores.tsv  BindCraft metrics, keyed by filename    (--bindcraft)

A column already present on a row is not overwritten by a later table; the clash is
reported on stderr rather than silently resolved.
"""

from __future__ import annotations

import argparse
import csv
import sys
from typing import Dict, List, Optional, Tuple

# Confidence columns that identify the row rather than measure anything. `id` and
# `binder` are handled separately: the score tables call the complex `id` and the
# design `binder`, while every other refold table calls the design `id`.
SCORE_DROP = ("original_file",)


def read_tsv(path: str) -> List[dict]:
    with open(path, newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t")]


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


def interface_rows(path: str) -> List[dict]:
    """Interface scores, attributed by the PDB filename they scored."""
    out = []
    for row in read_tsv(path):
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


class Index:
    """Long-table rows looked up by (id, tool, batch, model), falling back to
    (id, tool, batch) where that is unambiguous.

    The refold labels its rows with the engine's own sample label, the one the
    confidence table uses, so the four-part key is the normal path. The fallback
    covers a table written before that - or by an engine whose label cannot be
    recovered - and is only taken when the batch holds a single row, so it can
    never pair a sample with another sample's scores."""

    def __init__(self, rows: List[dict], keyed_on_model: bool):
        self.by_model: Dict[Tuple[str, str, str, str], dict] = {}
        self.by_batch: Dict[Tuple[str, str, str], List[dict]] = {}
        for row in rows:
            design = row.get("id", "")
            tool = row.get("tool", "")
            batch = row.get("batch", "1") or "1"
            if keyed_on_model:
                self.by_model[(design, tool, batch, row.get("model", ""))] = row
            self.by_batch.setdefault((design, tool, batch), []).append(row)

    def get(self, design: str, tool: str, batch: str, model: str) -> Optional[dict]:
        row = self.by_model.get((design, tool, batch, model))
        if row is not None:
            return row
        candidates = self.by_batch.get((design, tool, batch), [])
        return candidates[0] if len(candidates) == 1 else None


def merge(target: dict, columns: List[str], source: Optional[dict], where: str,
          clashes: set) -> None:
    if not source:
        return
    for key, value in source.items():
        if key in ("id", "tool", "batch", "model"):
            continue
        if key in target:
            if (where, key) not in clashes:
                clashes.add((where, key))
                print(
                    f"Warning: column '{key}' from {where} is already set on the row; "
                    "keeping the first value.",
                    file=sys.stderr,
                )
            continue
        target[key] = value
        if key not in columns:
            columns.append(key)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("scores", help="fold_pulldown_scores.tsv from the refold")
    parser.add_argument("-o", "--output", default="refold_sample_scores.tsv")
    parser.add_argument(
        "--long", action="append", default=[], metavar="PATH",
        help="A long table keyed by id, tool, batch and model (the RMSD tables). "
             "Repeatable; empty or header-only files are skipped.",
    )
    parser.add_argument(
        "--bindcraft", action="append", default=[], metavar="PATH",
        help="Interface scores whose filename column is <design>__<tool>__b<batch>.pdb. "
             "Repeatable.",
    )
    args = parser.parse_args()

    score_rows = read_tsv(args.scores)
    if not score_rows:
        # An empty input still needs a readable header, so a downstream read of this
        # file does not fail on a zero-byte table.
        with open(args.output, "w", newline="") as handle:
            handle.write("id\ttool\tbatch\tmodel\n")
        return 0

    long_indexes = [(path, Index(read_tsv(path), keyed_on_model=True)) for path in args.long]
    interface_indexes = [
        (path, Index(interface_rows(path), keyed_on_model=True)) for path in args.bindcraft
    ]

    columns: List[str] = ["id", "tool", "batch", "model", "pair_id", "target"]
    out_rows: List[dict] = []
    clashes: set = set()

    for row in score_rows:
        design = row.get("binder", "")
        tool = row.get("tool", "")
        batch = row.get("batch", "1") or "1"
        model = row.get("model", "")

        out = {
            "id": design,
            "tool": tool,
            "batch": batch,
            "model": model,
            "pair_id": row.get("id", ""),
            "target": row.get("target", ""),
        }
        for key, value in row.items():
            if key in ("id", "tool", "batch", "model", "binder", "target") or key in SCORE_DROP:
                continue
            out[key] = value
            if key not in columns:
                columns.append(key)

        for path, index in long_indexes:
            merge(out, columns, index.get(design, tool, batch, model), path, clashes)
        for path, index in interface_indexes:
            merge(out, columns, index.get(design, tool, batch, model), path, clashes)

        out_rows.append(out)

    out_rows.sort(key=lambda r: (r["id"], r["tool"], r["batch"], r["model"]))

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in out_rows:
            writer.writerow({c: row.get(c, "") for c in columns})
    return 0


if __name__ == "__main__":
    sys.exit(main())
