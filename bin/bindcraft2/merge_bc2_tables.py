#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Concatenate BindCraft2 per-stage CSV tables across independently run campaigns.

`bindcraft campaign_output` rewrites each campaign's ranking and writes a
combined summary.csv, but it does not concatenate the stage tables, so a
parallel run has no single view of its trajectories, candidates or accepted
designs. This fills that gap.

A `campaign` column is prepended to every row so rows stay attributable to the
batch that produced them. Column sets are unioned rather than assumed identical,
because a campaign that never reached a stage, or ran against a different number
of targets, legitimately carries different columns.
"""

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

# BC2 stage tables, as written under each campaign's project_folder.
STAGE_TABLES = {
    "trajectories": Path("1_Trajectories") / "!_Trajectories.csv",
    "refolded": Path("2_Refolded") / "!_Refolded.csv",
    "ranked": Path("3_Ranked") / "!_Ranked.csv",
}

CAMPAIGN_COLUMN = "campaign"


def campaign_dirs(parent: Path) -> List[Path]:
    return sorted(p for p in parent.iterdir() if p.is_dir())


def read_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        return fieldnames, [row for row in reader]


def merge_stage(parent: Path, relative: Path, output: Path) -> int:
    columns: List[str] = [CAMPAIGN_COLUMN]
    merged: List[Dict[str, str]] = []

    for campaign in campaign_dirs(parent):
        table = campaign / relative
        if not table.is_file():
            continue
        try:
            fieldnames, rows = read_rows(table)
        except OSError as e:
            # One unreadable campaign must not cost the other several hundred.
            log.warning("Skipping unreadable table %s: %s", table, e)
            continue
        for name in fieldnames:
            if name and name not in columns:
                columns.append(name)
        for row in rows:
            row.pop(None, None)
            merged.append({CAMPAIGN_COLUMN: campaign.name, **row})

    if not merged:
        log.info("No %s rows found under %s, not writing %s", relative, parent, output)
        return 0

    with open(output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, restval="", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged)

    log.info("Wrote %d rows from %s to %s", len(merged), relative, output)
    return len(merged)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Concatenate BindCraft2 stage tables across campaign folders"
    )
    parser.add_argument("parent", help="Directory whose subdirectories are campaign folders")
    parser.add_argument("--trajectories", default=None, help="Output CSV for 1_Trajectories")
    parser.add_argument("--refolded", default=None, help="Output CSV for 2_Refolded")
    parser.add_argument("--ranked", default=None, help="Output CSV for 3_Ranked")
    args = parser.parse_args()

    parent = Path(args.parent)
    if not parent.is_dir():
        log.error("Not a directory: %s", parent)
        return 1

    outputs: Dict[str, Optional[str]] = {
        "trajectories": args.trajectories,
        "refolded": args.refolded,
        "ranked": args.ranked,
    }
    for stage, output in outputs.items():
        if output:
            merge_stage(parent, STAGE_TABLES[stage], Path(output))

    return 0


if __name__ == "__main__":
    sys.exit(main())
