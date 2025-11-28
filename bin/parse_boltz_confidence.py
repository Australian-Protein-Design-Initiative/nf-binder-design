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
    parser.add_argument("--merge-ipsae", help="Path to *_ipsae.tsv file to merge (optional)")

    args = parser.parse_args()

    # Read the JSON file
    with open(args.json, "r") as f:
        data = json.load(f)

    # Flatten the nested dictionaries from the JSON data
    df_flat = pd.json_normalize(data, sep="_")

    # Merge IPSAE data if provided
    if args.merge_ipsae:
        try:
            # Read IPSAE TSV
            df_ipsae = pd.read_csv(args.merge_ipsae, sep="\t")
            
            # Filter for Type == 'min'
            if "Type" in df_ipsae.columns:
                row_min = df_ipsae[df_ipsae["Type"] == "min"]
                
                if not row_min.empty:
                    # Take the first match (should be unique for Type=min)
                    row_to_merge = row_min.iloc[[0]].copy()
                    
                    # Identify columns to exclude
                    cols_to_exclude = {"Chn1", "Chn2", "PAE", "Dist", "nres1", "nres2", 
                                       "dist1", "dist2", "Type", "Model"}
                    cols_to_drop = [c for c in row_to_merge.columns 
                                    if c in cols_to_exclude or c.startswith("d0") or c.startswith("n0")]
                    
                    row_to_merge.drop(columns=cols_to_drop, inplace=True, errors="ignore")
                    
                    # Rename ipSAE to ipSAE_min
                    if "ipSAE" in row_to_merge.columns:
                        row_to_merge.rename(columns={"ipSAE": "ipSAE_min"}, inplace=True)
                    
                    # Merge with df_flat
                    df_flat = pd.concat([df_flat.reset_index(drop=True), 
                                         row_to_merge.reset_index(drop=True)], axis=1)
                else:
                    print(f"Warning: No row with Type='min' found in {args.merge_ipsae}", file=sys.stderr)
            else:
                print(f"Warning: 'Type' column not found in {args.merge_ipsae}", file=sys.stderr)
                
        except Exception as e:
            print(f"Error merging IPSAE file: {e}", file=sys.stderr)
            sys.exit(1)

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
