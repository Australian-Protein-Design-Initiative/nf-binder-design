#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///
"""Add a Target column to BindCraft per-batch stats CSVs and write target.txt."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
logger = logging.getLogger(__name__)

DEFAULT_CSVS = (
    "final_design_stats.csv",
    "trajectory_stats.csv",
    "mpnn_design_stats.csv",
)


def add_target_column(csv_path: Path, target: str) -> bool:
    if not csv_path.exists():
        logger.info("Skipping missing CSV: %s", csv_path)
        return False

    with csv_path.open(newline="") as handle:
        rows = list(csv.reader(handle))

    if not rows:
        logger.info("Skipping empty CSV: %s", csv_path)
        return False

    header = rows[0]
    if "Target" in header:
        logger.info("Target column already present in %s", csv_path)
        return False

    new_header = ["Target"] + header
    new_rows = [new_header]
    for row in rows[1:]:
        new_rows.append([target] + row)

    with csv_path.open("w", newline="") as handle:
        csv.writer(handle).writerows(new_rows)

    logger.info("Added Target=%s to %s (%d data rows)", target, csv_path, len(rows) - 1)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target",
        required=True,
        help="Target structure filename including extension, e.g. PDL1.pdb or PDL1-2.cif",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results"),
        help="BindCraft results directory containing stats CSVs (default: results)",
    )
    parser.add_argument(
        "--csv",
        action="append",
        dest="csvs",
        help="CSV filename under results-dir (repeatable). Default: final/trajectory/mpnn stats",
    )
    args = parser.parse_args()

    target = args.target.strip()
    if not target:
        logger.error("--target must be non-empty")
        return 1

    results_dir = args.results_dir
    results_dir.mkdir(parents=True, exist_ok=True)

    csv_names = args.csvs if args.csvs else list(DEFAULT_CSVS)
    for name in csv_names:
        add_target_column(results_dir / name, target)

    target_file = results_dir / "target.txt"
    target_file.write_text(f"{target}\n", encoding="utf-8")
    logger.info("Wrote %s", target_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
