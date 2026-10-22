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
"""

import argparse
import csv
import logging
import shutil
import sys
from pathlib import Path
from typing import List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

GENERATED_SUFFIXES = ("_generated.cif", "_generated.cif.gz", "_generated.cif.zst")


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

    rows = []
    for path in designs:
        example_id, extension = split_generated_name(path.name)
        design_id = f"{args.design_name}_b{args.batch}_{example_id}"
        destination = out_dir / f"{design_id}{extension}"
        shutil.copy2(path, destination)
        rows.append(
            {
                "id": design_id,
                "filename": destination.name,
                "example_id": example_id,
                "batch": args.batch,
                "seed": args.seed,
            }
        )
        log.info("Collected %s -> %s", path, destination)

    tsv_path = Path(args.tsv) if args.tsv else out_dir.parent / "rfd4_designs.tsv"
    with open(tsv_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["id", "filename", "example_id", "batch", "seed"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    log.info("Wrote %d design rows to %s", len(rows), tsv_path)


if __name__ == "__main__":
    main()
