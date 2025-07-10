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
        args.verbose,
    )
    logger.info("Completed successfully")
