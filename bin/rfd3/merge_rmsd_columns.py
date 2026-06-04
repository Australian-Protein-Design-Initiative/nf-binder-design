#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "pandas",
# ]
# ///

"""
Merge one or more RMSD TSVs (each with a structure path column + rmsd value
columns) onto a base score table in a single pass, applying a per-file column
prefix. Replaces the chain of `csvtk cut | merge_scores.py --column-prefix
--drop-columns` blocks in combine_rfd3_scores.nf.

Each RMSD spec is `PATH=PREFIX`. Only the columns named in --value-columns are
kept from each RMSD TSV (prefixed); the structure key column is dropped. The
merge key is built with the same basename normalisation as merge_scores.py so
behaviour is identical to the previous sequential merges.
"""

import argparse
import logging
import os
import sys
from typing import List, Optional, Tuple

import pandas as pd

# Reuse the exact merge-key logic from the sibling bin/merge_scores.py so the
# join behaviour (path-column detection, basename substitution, suffix stripping)
# stays byte-for-byte identical to the previous sequential merge_scores calls.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from merge_scores import (  # noqa: E402
    _compile_basename_subs,
    _is_file_empty_or_whitespace,
    _prepare_df_for_merge,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def _parse_spec(spec: str) -> Tuple[str, str]:
    if "=" not in spec:
        raise argparse.ArgumentTypeError(
            f"RMSD spec must be PATH=PREFIX, got: {spec!r}"
        )
    path, prefix = spec.split("=", 1)
    if not path or not prefix:
        raise argparse.ArgumentTypeError(
            f"RMSD spec must be PATH=PREFIX with non-empty parts, got: {spec!r}"
        )
    return path, prefix


def merge_rmsd_columns(
    base_file: str,
    specs: List[Tuple[str, str]],
    keys: str,
    strip_suffix: str,
    value_columns: List[str],
    basename_subs: Optional[List[Tuple[str, str]]],
    first_column: str,
    output_file: str,
) -> pd.DataFrame:
    base_df = pd.read_csv(base_file, sep="\t")
    potential_keys = [k.strip() for k in keys.split(",") if k.strip()]
    subs = _compile_basename_subs(basename_subs)

    base_df, _ = _prepare_df_for_merge(base_df, potential_keys, strip_suffix, subs)

    for path, prefix in specs:
        if _is_file_empty_or_whitespace(path):
            log.warning("RMSD file %s is empty or whitespace-only. Skipping.", path)
            continue
        right = pd.read_csv(path, sep="\t")
        right, _ = _prepare_df_for_merge(right, potential_keys, strip_suffix, subs)

        keep = [c for c in value_columns if c in right.columns]
        missing = [c for c in value_columns if c not in right.columns]
        if missing:
            log.warning(
                "Value column(s) %s not found in %s. Available: %s",
                missing,
                path,
                [c for c in right.columns if c != "_merge_key"],
            )
        right = right[["_merge_key"] + keep].rename(
            columns={c: f"{prefix}{c}" for c in keep}
        )
        base_df = pd.merge(base_df, right, on="_merge_key", how="left")

    base_df = base_df.drop(columns=["_merge_key"])

    first_cols = [c.strip() for c in first_column.split(",") if c.strip()]
    existing_first = [c for c in first_cols if c in base_df.columns]
    if existing_first:
        remaining = [c for c in base_df.columns if c not in existing_first]
        base_df = base_df[existing_first + remaining]

    if output_file == "-":
        base_df.to_csv(sys.stdout, sep="\t", index=False, lineterminator="\n")
    else:
        base_df.to_csv(output_file, sep="\t", index=False, lineterminator="\n")
        log.info("Wrote %s", output_file)

    return base_df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge RMSD value columns onto a base score table with per-file prefixes."
    )
    parser.add_argument("base", help="Base score TSV to merge RMSD columns onto.")
    parser.add_argument(
        "rmsd_specs",
        nargs="+",
        type=_parse_spec,
        metavar="PATH=PREFIX",
        help="One or more RMSD TSVs as PATH=PREFIX (prefix applied to kept value columns).",
    )
    parser.add_argument(
        "--keys",
        default="filename,structure1",
        help="Comma-separated key columns to try, in order, on every table "
        '(default: "filename,structure1").',
    )
    parser.add_argument(
        "--strip-suffix",
        default="(\\.pdb|\\.cif)$",
        help='Regex stripped from path-like keys when building the merge key (default: "(\\.pdb|\\.cif)$").',
    )
    parser.add_argument(
        "--value-columns",
        default="rmsd_all",
        help='Comma-separated columns to keep from each RMSD TSV (default: "rmsd_all").',
    )
    parser.add_argument(
        "--merge-key-replace-basename",
        nargs=2,
        metavar=("REGEX", "REPL"),
        action="append",
        default=None,
        help="Apply re.sub(REGEX, REPL) to path-like key basenames before --strip-suffix. May be repeated.",
    )
    parser.add_argument(
        "--first-column",
        default="filename",
        help='Column(s) to move to the front of the output, comma-separated (default: "filename").',
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help='Output TSV path (default: "-" for stdout).',
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    value_columns = [c.strip() for c in args.value_columns.split(",") if c.strip()]
    merge_rmsd_columns(
        args.base,
        args.rmsd_specs,
        args.keys,
        args.strip_suffix,
        value_columns,
        args.merge_key_replace_basename,
        args.first_column,
        args.output,
    )
