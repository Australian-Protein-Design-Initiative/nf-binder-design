#!/usr/bin/env python3

import sys
from pathlib import Path
import argparse
import pandas as pd
import io


def parse_score_file(file_path):
    """Parse a single .cs file and return a pandas DataFrame."""
    try:
        with open(file_path, "r") as f:
            lines = f.readlines()

        # Find the header line (starts with SCORE:)
        header_line = None
        data_line = None
        for i, line in enumerate(lines):
            if line.startswith("SCORE:"):
                if header_line is None:
                    header_line = line.strip()
                else:
                    data_line = line.strip()
                    break

        if header_line is None or data_line is None:
            print(
                f"Warning: Could not find SCORE: lines in {file_path}", file=sys.stderr
            )
            return None

        # Extract column names and data
        columns = [x.strip() for x in header_line.split()[1:]]
        data = [x.strip() for x in data_line.split()[1:]]

        # Create DataFrame with a single row
        df = pd.DataFrame([data], columns=columns)

        # Add filename column
        df["source_file"] = file_path.name

        return df

    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}", file=sys.stderr)
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Combine multiple .cs score files into a single table"
    )
    parser.add_argument("path", help="Path to search for .cs files")
    parser.add_argument(
        "--output",
        "-o",
        help="Output file path (default: stdout, use - for stdout)",
        default="-",
    )
    parser.add_argument(
        "--sort-by",
        help="Column to sort results by (default: pae_interaction, ascending)",
        default="pae_interaction",
    )
    parser.add_argument(
        "--descending",
        action="store_true",
        help="Sort in descending order (default: ascending)",
    )

    args = parser.parse_args()

    # Convert path to Path object
    search_path = Path(args.path)

    if not search_path.exists():
        print(f"Error: Path {args.path} does not exist", file=sys.stderr)
        sys.exit(1)

    # Find all .cs files
    cs_files = list(search_path.glob("**/*.cs"))

    if not cs_files:
        print(f"No .cs files found in {args.path}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(cs_files)} .cs files", file=sys.stderr)

    # Process each file and combine into a single DataFrame
    dfs = []
    for file_path in cs_files:
        df = parse_score_file(file_path)
        if df is not None:
            dfs.append(df)

    if not dfs:
        print("No valid data found in any of the .cs files", file=sys.stderr)
        sys.exit(1)

    # Combine all DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)

    # Convert numeric columns to float
    numeric_columns = [
        "binder_aligned_rmsd",
        "pae_binder",
        "pae_interaction",
        "pae_target",
        "plddt_binder",
        "plddt_target",
        "plddt_total",
        "target_aligned_rmsd",
        "time",
    ]

    for col in numeric_columns:
        if col in combined_df.columns:
            combined_df[col] = pd.to_numeric(combined_df[col], errors="coerce")

    # Sort the dataframe
    if args.sort_by not in combined_df.columns:
        print(
            f"Warning: Sort column '{args.sort_by}' not found in data", file=sys.stderr
        )
    else:
        combined_df = combined_df.sort_values(
            by=args.sort_by, ascending=not args.descending
        )

    # Save to TSV
    if args.output == "-":
        combined_df.to_csv(sys.stdout, index=False, sep="\t")
    else:
        combined_df.to_csv(args.output, index=False, sep="\t")
        print(f"Combined scores saved to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
