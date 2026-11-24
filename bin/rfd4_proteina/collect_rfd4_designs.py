#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Flatten an rfd4-generate run directory into one directory of uniquely named designs.

rfd4-generate writes `<out_dir>/inference_samples/<example_id>/<example_id>_generated.cif`,
with one subdirectory per input example and per replicate. Each Nextflow task runs its own
chunk of replicates, so the example ids repeat across tasks; this prefixes them with the
batch index to keep them unique once they are collected, and records the mapping in a TSV.

Compressed outputs (.cif.zst / .cif.gz) are carried through with their suffix intact.

The run's inline metrics (`<out_dir>/inference_metrics/per_example_metrics.csv`, written when
`metrics_tags` is set) are merged into the same TSV, so the oracle refold scores land in
`combined_scores.tsv` keyed to the design ids. The join is on example_id within this one run,
which is why it happens here rather than downstream: example ids repeat across batches.
"""

import argparse
import csv
import logging
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

GENERATED_SUFFIXES = ("_generated.cif", "_generated.cif.gz", "_generated.cif.zst")

# Inline metric columns are prefixed so they cannot collide with refold columns, which are
# prefixed with their engine name.
METRIC_PREFIX = "rfd4_"
# Hydra config and run identity, repeated on every row and already recorded elsewhere.
DROP_COLUMN_PREFIXES = ("generation.", "ckpt_")
DROP_COLUMNS = {
    "work_hash",
    "example_id",
    "designed_structure_path",
    "task",
    "dataset_name",
    "chain_of_interest",
    "seed",
    "input_json_path",
    "output_base_dir",
    "metrics_tags",
}


def metric_columns(fieldnames: List[str]) -> List[str]:
    return [
        name
        for name in fieldnames
        if name not in DROP_COLUMNS and not name.startswith(DROP_COLUMN_PREFIXES)
    ]


def rewrite_run_path(value: str, batch: str) -> str:
    """Point a path inside the run directory at where publishDir puts that batch."""
    marker = "rundir/"
    index = value.find(marker)
    return f"batches/batch{batch}/{value[index + len(marker):]}" if index >= 0 else value


def load_metrics(path: Path, batch: str) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """Inline metrics keyed by example id, with the column order they should be written in.

    Every row in a batch gets the same columns, because the per-batch tables are concatenated
    downstream with `collectFile(keepHeader: true)`, which assumes one shared header.
    """
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "example_id" not in reader.fieldnames:
            raise SystemExit(f"{path} has no example_id column")
        keep = metric_columns(reader.fieldnames)
        metrics = {}
        for row in reader:
            metrics[row["example_id"]] = {
                # Dots are awkward in a column name for pandas/R consumers downstream
                f"{METRIC_PREFIX}{name.replace('.', '_')}": (
                    rewrite_run_path(row[name], batch) if name.endswith("_path") else row[name]
                )
                for name in keep
            }
    log.info("Read inline metrics for %d examples from %s", len(metrics), path)
    return metrics, [f"{METRIC_PREFIX}{name.replace('.', '_')}" for name in keep]


def split_generated_name(name: str) -> Optional[Tuple[str, str]]:
    """Split '<example>_generated.cif[.gz|.zst]' into (example_id, extension)."""
    for suffix in GENERATED_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)], suffix[len("_generated") :]
    return None


def find_designs(samples_dir: Path) -> List[Path]:
    return sorted(p for p in samples_dir.rglob("*") if p.is_file() and split_generated_name(p.name))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collect rfd4-generate designs into a flat, uniquely named directory",
        allow_abbrev=False,
    )
    parser.add_argument("run_dir", help="rfd4-generate out_dir (holding inference_samples/)")
    parser.add_argument("--design-name", required=True)
    parser.add_argument("--batch", required=True, help="Batch index, used to keep ids unique")
    parser.add_argument("--seed", default="", help="Seed used for this batch, recorded in the TSV")
    parser.add_argument("-o", "--output-dir", default="output")
    parser.add_argument("--tsv", default=None, help="Where to write the design table")
    parser.add_argument(
        "--metrics",
        default=None,
        help="Inline metrics CSV [default: <run_dir>/inference_metrics/per_example_metrics.csv]",
    )
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Exit 0 when the run produced no designs, instead of failing",
    )
    args = parser.parse_args()

    samples_dir = Path(args.run_dir) / "inference_samples"
    if not samples_dir.is_dir():
        raise SystemExit(f"No inference_samples/ under {args.run_dir}; the run produced nothing")

    designs = find_designs(samples_dir)
    if not designs and not args.allow_empty:
        raise SystemExit(f"No *_generated.cif found under {samples_dir}")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics_path = (
        Path(args.metrics)
        if args.metrics
        else Path(args.run_dir) / "inference_metrics" / "per_example_metrics.csv"
    )
    # Absent whenever metrics_tags is false, which is a supported way to run
    metrics: Dict[str, Dict[str, str]] = {}
    metric_fields: List[str] = []
    if metrics_path.is_file():
        metrics, metric_fields = load_metrics(metrics_path, args.batch)
    else:
        # Absent whenever metrics_tags is false, which is a supported way to run
        log.info("No inline metrics at %s; writing the design table without them", metrics_path)

    rows = []
    for path in designs:
        example_id, extension = split_generated_name(path.name)
        design_id = f"{args.design_name}_b{args.batch}_{example_id}"
        destination = out_dir / f"{design_id}{extension}"
        shutil.copy2(path, destination)
        row = {
            "id": design_id,
            "filename": destination.name,
            "example_id": example_id,
            "batch": args.batch,
            "seed": args.seed,
        }
        design_metrics = metrics.get(example_id)
        if design_metrics:
            row.update(design_metrics)
        elif metrics:
            log.warning("No inline metrics row for example %s", example_id)
        rows.append(row)
        log.info("Collected %s -> %s", path, destination)

    tsv_path = Path(args.tsv) if args.tsv else out_dir.parent / "rfd4_designs.tsv"
    with open(tsv_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["id", "filename", "example_id", "batch", "seed"] + metric_fields,
            delimiter="\t",
            restval="",
        )
        writer.writeheader()
        writer.writerows(rows)
    log.info("Wrote %d design rows to %s", len(rows), tsv_path)


if __name__ == "__main__":
    main()
