#!/usr/bin/env python
# /// script
# requires-python = ">=3.7"
# dependencies = [
#     "pandas",
# ]
# ///

import argparse
import os
import sys
import pandas as pd
import logging
import re
from typing import List, Optional, TextIO, Union
from functools import reduce

# Set up logging to stderr
logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stderr,
)
logger = logging.getLogger(__name__)


def _prepare_df_for_merge(
    df: pd.DataFrame, potential_keys: List[str], strip_suffix: str
) -> pd.DataFrame:
    """Finds a key in the dataframe and creates a `_merge_key` column for merging."""
    for key in potential_keys:
        if key in df.columns:
            # Check if the column seems to contain PDB file paths by looking for '.pdb' suffix
            col_series = df[key].dropna()
            is_pdb_col = False
            if (
                pd.api.types.is_string_dtype(col_series)
                and col_series.str.endswith(".pdb").any()
            ):
                is_pdb_col = True

            if is_pdb_col:
                logger.info(
                    f"Found PDB column '{key}'. Creating merge key from basenames."
                )
                # Handle non-string values gracefully (e.g., NaN)
                df["_merge_key"] = df[key].apply(
                    lambda x: (
                        re.sub(strip_suffix, "", os.path.basename(x))
                        if isinstance(x, str)
                        else x
                    )
                )
            else:
                logger.info(f"Found key column '{key}'. Using it as the merge key.")
                df["_merge_key"] = df[key]

            return df

    raise ValueError(
        f"No potential merge key ({potential_keys}) found in one of the dataframes. "
        f"Columns available: {df.columns.tolist()}"
    )


def merge_scores(
    tsv_files: List[str],
    keys_str: str,
    strip_suffix: str,
    output_file: Union[str, TextIO, None] = None,
    sort_by: str = "pae_interaction",
    first_column: str = "filename",
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Merges multiple TSV score files into a single DataFrame.

    It uses a list of potential key columns to find a common identifier for merging.
    If a key column contains PDB file paths, it intelligently extracts the basename
    to match against other keys.
    """
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    potential_keys = [k.strip() for k in keys_str.split(",")]

    if not tsv_files:
        logger.warning("No TSV files provided to merge.")
        return pd.DataFrame()

    logger.info(f"Reading {len(tsv_files)} TSV files.")
    dataframes = [pd.read_csv(f, sep="\t") for f in tsv_files]

    logger.info(f"Preparing dataframes for merge using keys: {potential_keys}")
    prepared_dfs = [
        _prepare_df_for_merge(df.copy(), potential_keys, strip_suffix)
        for df in dataframes
    ]

    # Merge the dataframes sequentially
    logger.info("Merging dataframes...")
    merged_df = reduce(
        lambda left, right: pd.merge(left, right, on="_merge_key", how="left"),
        prepared_dfs,
    )

    # Drop the temporary matching column
    merged_df = merged_df.drop(columns=["_merge_key"])

    # --- Deduplicate columns with the same base name if their values are identical ---
    def deduplicate_columns(df):
        from collections import defaultdict

        # Find columns with suffixes (_x, _y, _z, ...)
        col_map = defaultdict(list)
        for col in df.columns:
            if col.endswith(
                tuple([f"_{chr(i)}" for i in range(120, 123)])
            ):  # _x, _y, _z
                base = col[:-2]
                col_map[base].append(col)
        # Also check for more than 3 (e.g., _w, _v, ...)
        for col in df.columns:
            if len(col) > 2 and col[-2] == "_" and col[-1].isalpha():
                base = col[:-2]
                col_map[base].append(col)
        # Remove duplicates in col_map
        for base in list(col_map.keys()):
            col_map[base] = list(set(col_map[base]))
            if len(col_map[base]) < 2:
                del col_map[base]
        # For each set of duplicate columns
        for base, cols in col_map.items():
            # Add the base column if it exists (no suffix)
            if base in df.columns:
                cols = [base] + cols
            # Compare all columns
            arrays = [df[c] for c in cols]
            all_equal = True
            for i in range(1, len(arrays)):
                # Use pandas equals for robust comparison (handles NaN)
                if not arrays[0].equals(arrays[i]):
                    all_equal = False
                    break
            if all_equal:
                # Keep only one column (the base if present, else the first suffixed)
                keep_col = base if base in df.columns else cols[0]
                df[base] = df[keep_col]
                for c in cols:
                    if c != base:
                        df.drop(columns=[c], inplace=True, errors="ignore")
            else:
                logger.warning(
                    f"Duplicate column name '{base}' found with differing values. Keeping all with suffixes."
                )
        return df

    merged_df = deduplicate_columns(merged_df)
    # --- End deduplication logic ---

    # Sort the merged dataframe
    if sort_by in merged_df.columns:
        logger.info(f"Sorting merged dataframe by {sort_by}")
        merged_df = merged_df.sort_values(by=sort_by, ascending=True)
        # Move the sort_by column to the first position
        cols = [sort_by] + [col for col in merged_df.columns if col != sort_by]
        merged_df = merged_df[cols]
    else:
        logger.warning(
            f"Sort key '{sort_by}' not found in merged dataframe. Skipping sorting."
        )

    # Ensure 'filename' is the first column if present
    if first_column in merged_df.columns:
        cols = [first_column] + [
            col for col in merged_df.columns if col != first_column
        ]
        merged_df = merged_df[cols]
    else:
        logger.warning(
            f"First column '{first_column}' not found in merged dataframe. Skipping."
        )

    # Write to output file if specified
    if output_file:
        if output_file == "-":
            logger.info("Writing combined scores to stdout")
            merged_df.to_csv(sys.stdout, sep="\t", index=False)
        else:
            logger.info(f"Writing combined scores to {output_file}")
            merged_df.to_csv(output_file, sep="\t", index=False)

    return pd.DataFrame(merged_df)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge multiple score TSV files based on a shared key."
    )
    parser.add_argument(
        "tsv_files", nargs="+", help="Paths to one or more TSV score files to merge."
    )
    parser.add_argument(
        "--keys",
        default="description,filename",
        help='Comma-separated list of column names to look for, in order of preference, to use as the merge key (default: "description,filename").',
    )
    parser.add_argument(
        "--strip-suffix",
        default="(\\.pdb)$",
        help='Regex to strip from filenames when creating the merge key (default: "(_af2pred\\.pdb|\\.pdb)$").',
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help='Output filename for combined scores (default: "-" for stdout)',
    )
    parser.add_argument(
        "--sort-by",
        default="pae_interaction",
        help='Column name to sort the merged DataFrame by (default: "pae_interaction")',
    )
    parser.add_argument(
        "--first-column",
        default="filename",
        help='Column name to move to the first position in the output (default: "filename")',
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    merge_scores(
        args.tsv_files,
        args.keys,
        args.strip_suffix,
        args.output,
        args.sort_by,
        args.first_column,
        args.verbose,
    )
    logger.info("Completed successfully")
