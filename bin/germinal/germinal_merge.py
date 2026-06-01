#!/usr/bin/env python
"""Merge Germinal batch output directories into a single results tree."""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable


logger = logging.getLogger(__name__)

DESIGN_CSVS = {
    "all_trajectories.csv": "all_trajectories.csv",
    "accepted/designs.csv": "accepted_designs.csv",
    "trajectories/designs.csv": "trajectories/designs.csv",
    "redesign_candidates/designs.csv": "redesign_candidates/designs.csv",
}

STRUCTURE_DIRS = [
    "accepted/structures",
    "accepted/plots",
    "trajectories/structures",
    "trajectories/plots",
    "redesign_candidates/structures",
    "redesign_candidates/plots",
]


def parse_args() -> argparse.Namespace:
    if len(sys.argv) > 1 and sys.argv[1] in {
        "merge-batches",
        "merge-csvs",
        "merge-failure-counts",
    }:
        return _parse_subcommand_args()

    legacy_parser = argparse.ArgumentParser(
        description="Merge Germinal per-batch output directories."
    )
    legacy_parser.add_argument(
        "batch_dirs",
        nargs="+",
        type=Path,
        help="Batch result directories (each contains Germinal experiment output)",
    )
    legacy_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("merged"),
        help="Output directory for merged results (default: merged)",
    )
    args = legacy_parser.parse_args()
    args.command = "merge-batches"
    return args


def _parse_subcommand_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge Germinal per-batch output directories or CSV files."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    batch_parser = subparsers.add_parser(
        "merge-batches",
        help="Merge Germinal batch output directories",
    )
    batch_parser.add_argument(
        "batch_dirs",
        nargs="+",
        type=Path,
        help="Batch result directories (each contains Germinal experiment output)",
    )
    batch_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("merged"),
        help="Output directory for merged results (default: merged)",
    )

    csv_parser = subparsers.add_parser(
        "merge-csvs",
        help="Merge Germinal design CSV files with aligned columns",
    )
    csv_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output CSV path",
    )
    csv_parser.add_argument(
        "csv_files",
        nargs="+",
        type=Path,
        help="Input CSV files to merge",
    )
    csv_parser.add_argument(
        "--dedupe-column",
        default="design_name",
        help="Column used to drop duplicate rows (default: design_name)",
    )

    failure_parser = subparsers.add_parser(
        "merge-failure-counts",
        help="Sum Germinal failure_counts.csv files",
    )
    failure_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output CSV path",
    )
    failure_parser.add_argument(
        "csv_files",
        nargs="+",
        type=Path,
        help="Input failure_counts.csv files to merge",
    )

    return parser.parse_args()


def find_csv(batch_dirs: Iterable[Path], relative_path: str) -> list[Path]:
    files: list[Path] = []
    for batch_dir in batch_dirs:
        matches = sorted(batch_dir.rglob(relative_path))
        if not matches:
            logger.warning("Missing %s in %s", relative_path, batch_dir)
            continue
        if len(matches) > 1:
            logger.warning(
                "Multiple matches for %s in %s; using %s",
                relative_path,
                batch_dir,
                matches[0],
            )
        files.append(matches[0])
    return files


def merge_design_csvs(
    input_files: list[Path],
    output_path: Path,
    dedupe_column: str = "design_name",
) -> None:
    if not input_files:
        logger.warning("No input files for %s; skipping", output_path)
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    columns: list[str] = []
    rows: list[dict[str, str]] = []

    for input_file in input_files:
        if not input_file.exists():
            logger.warning("Missing file %s", input_file)
            continue

        if input_file.stat().st_size == 0:
            logger.warning("Skipping empty file %s", input_file)
            continue

        with input_file.open(newline="") as in_fh:
            reader = csv.DictReader(in_fh)
            if not reader.fieldnames:
                logger.warning("Skipping file with no header: %s", input_file)
                continue

            file_columns = list(reader.fieldnames)
            if columns and file_columns != columns:
                added = [col for col in file_columns if col not in columns]
                if added:
                    logger.warning(
                        "Header mismatch in %s; adding columns: %s",
                        input_file,
                        ", ".join(added),
                    )
            for col in file_columns:
                if col not in columns:
                    columns.append(col)

            for row in reader:
                if not row or not any((value or "").strip() for value in row.values()):
                    continue
                if dedupe_column in row and row[dedupe_column]:
                    key = row[dedupe_column]
                    if key in seen:
                        continue
                    seen.add(key)
                rows.append(row)

    if not columns:
        logger.warning(
            "No CSV headers found for %s (%d input files, all empty); skipping",
            output_path,
            len(input_files),
        )
        return

    with output_path.open("w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})

    logger.info("Wrote %d rows to %s", len(rows), output_path)


def merge_failure_counts(input_files: list[Path], output_path: Path) -> None:
    if not input_files:
        logger.warning("No failure_counts.csv files found; skipping")
        return

    sums: dict[str, int | float] = defaultdict(int)
    columns: list[str] = []

    for input_file in input_files:
        with input_file.open(newline="") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None:
                continue
            for col in reader.fieldnames:
                if col not in columns:
                    columns.append(col)
            for row in reader:
                for col in reader.fieldnames:
                    value = row.get(col, "")
                    if value == "":
                        continue
                    try:
                        numeric = int(value)
                    except ValueError:
                        try:
                            numeric = float(value)
                        except ValueError:
                            logger.warning(
                                "Non-numeric failure count %s=%s in %s",
                                col,
                                value,
                                input_file,
                            )
                            continue
                    sums[col] += numeric

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        writer.writerow([sums.get(col, 0) for col in columns])

    logger.info("Merged failure counts to %s", output_path)


def copy_tree_union(batch_dirs: Iterable[Path], relative_dir: str, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    copied = 0
    for batch_dir in batch_dirs:
        for src_dir in batch_dir.rglob(relative_dir):
            if not src_dir.is_dir():
                continue
            for src_file in src_dir.iterdir():
                if not src_file.is_file():
                    continue
                dest_file = output_dir / src_file.name
                if dest_file.exists():
                    continue
                shutil.copy2(src_file, dest_file)
                copied += 1
    logger.info("Copied %d files into %s", copied, output_dir)


def copy_resolved_config(batch_dirs: Iterable[Path], output_dir: Path) -> None:
    """Copy one resolved final_config.yaml into config/."""
    configs_dir = output_dir / "config"
    dest = configs_dir / "final_config.yaml"
    if dest.exists():
        logger.info("Keeping existing %s", dest)
        return
    for batch_dir in batch_dirs:
        matches = sorted(batch_dir.rglob("final_config.yaml"))
        if not matches:
            continue
        configs_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(matches[0], dest)
        logger.info("Copied final_config.yaml from %s to %s", matches[0], dest)
        return
    logger.warning("No final_config.yaml found in batch directories")


def merge_batch_directories(batch_dirs: list[Path], output_dir: Path) -> None:
    for relative_path, output_name in DESIGN_CSVS.items():
        merge_design_csvs(
            find_csv(batch_dirs, relative_path),
            output_dir / output_name,
        )

    merge_failure_counts(
        find_csv(batch_dirs, "failure_counts.csv"),
        output_dir / "failure_counts.csv",
    )

    for relative_dir in STRUCTURE_DIRS:
        copy_tree_union(batch_dirs, relative_dir, output_dir / relative_dir)

    copy_resolved_config(batch_dirs, output_dir)


def main() -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        stream=sys.stderr,
    )
    args = parse_args()

    if args.command == "merge-csvs":
        merge_design_csvs(
            [path.resolve() for path in args.csv_files],
            args.output.resolve(),
            dedupe_column=args.dedupe_column,
        )
        return 0

    if args.command == "merge-failure-counts":
        merge_failure_counts(
            [path.resolve() for path in args.csv_files],
            args.output.resolve(),
        )
        return 0

    batch_dirs = [path.resolve() for path in args.batch_dirs if path.exists()]
    if not batch_dirs:
        logger.error("No valid batch directories provided")
        return 1

    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    merge_batch_directories(batch_dirs, output_dir)

    logger.info("Merged %d batch directories into %s", len(batch_dirs), output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
