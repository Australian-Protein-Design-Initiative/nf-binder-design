#!/usr/bin/env python
# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "pandas",
# ]
# ///

"""
Parse Boltz confidence JSON files and output flattened TSV.

This script takes a Boltz confidence JSON file and flattens it into a TSV format,
adding metadata columns for design ID, target, and binder names.
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Parse Boltz confidence JSON to TSV")
    parser.add_argument("--json", required=True, help="Path to confidence JSON file")
    parser.add_argument("--id", required=True, help="Design ID")
    parser.add_argument("--target", help="Target name (optional)")
    parser.add_argument("--binder", help="Binder name (optional)")

    args = parser.parse_args()

    # Read the JSON file
    with open(args.json, "r") as f:
        data = json.load(f)

    # Flatten the nested dictionaries from the JSON data
    df_flat = pd.json_normalize(data, sep="_")

    # Add the metadata columns to the beginning of the DataFrame
    df_flat.insert(0, "id", args.id)
    if args.target:
        df_flat.insert(1, "target", args.target)
    if args.binder:
        df_flat.insert(2, "binder", args.binder)

    # Write the flattened data to stdout as a TSV
    df_flat.to_csv(sys.stdout, sep="\t", index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
