#!/usr/bin/env python
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "mdanalysis",
#     "pandas",
# ]
# ///

import argparse
import sys
import pandas as pd
import MDAnalysis as mda
from pathlib import Path
import logging


def register_metrics() -> list[str]:
    """Returns a list of metrics that this plugin can calculate."""
    return ["rg"]


def calculate_rg(pdb_file: str, chain: str = "A") -> float:
    """Calculates the radius of gyration for a given chain in a PDB file."""
    u = mda.Universe(pdb_file)
    selection = u.select_atoms(f"chainID {chain}")
    if not selection:
        raise ValueError(f"Chain {chain} not found in {pdb_file}")
    return selection.radius_of_gyration()


def calculate_metrics(
    pdb_files: list[str], binder_chains: list[str] = ["A"]
) -> pd.DataFrame:
    """
    Calculates the radius of gyration for a list of PDB files.
    """
    # For now, we only handle a single binder chain for simplicity
    if len(binder_chains) > 1:
        logging.warning(
            "This plugin currently only calculates rg for the first specified binder chain."
        )

    chain = binder_chains[0]

    results = []
    for pdb_path in pdb_files:
        try:
            rg = calculate_rg(str(pdb_path), chain)
            results.append({"design_id": Path(pdb_path).stem, "rg": rg})
        except Exception as e:
            logging.error(f"Error processing {pdb_path}: {e}")
            results.append({"design_id": Path(pdb_path).stem, "rg": None})

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.set_index("design_id")
    return df


def main():
    """For standalone execution."""
    parser = argparse.ArgumentParser(
        description="Calculate Radius of Gyration (Rg) for a PDB file."
    )
    parser.add_argument("pdb_file", help="Path to the PDB file.")
    parser.add_argument(
        "--chain", default="A", help="Chain to calculate Rg for. Defaults to 'A'."
    )
    parser.add_argument(
        "--format", default="tsv", choices=["tsv", "json"], help="Output format."
    )

    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, stream=sys.stderr, format="%(levelname)s: %(message)s"
    )

    try:
        rg = calculate_rg(args.pdb_file, args.chain)
        design_id = Path(args.pdb_file).stem

        if args.format == "json":
            import json

            print(json.dumps({design_id: {"rg": rg}}))
        else:  # tsv
            print(f"design_id\\trg")
            print(f"{design_id}\\t{rg}")

    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
