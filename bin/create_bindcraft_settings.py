#!/usr/bin/env python

import argparse
import json
import logging
import math
import os
import random
import sys


def main():
    
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description="Create a BindCraft settings JSON file."
    )
    parser.add_argument(
        "--input_pdb", required=True, help="Input PDB file for the target."
    )
    parser.add_argument(
        "--hotspot_res",
        default="",
        required=False,
        help="Hotspot residues, e.g., 'A29,A32,A34'.",
    )
    parser.add_argument(
        "--target_chains",
        required=True,
        help="The target chain(s), e.g., 'A' or 'A,B'.",
    )
    parser.add_argument(
        "--binder_length_range",
        required=True,
        help="Dash-separated min and max binder lengths, e.g., '65-170'.",
    )
    parser.add_argument("--design_name", required=True, help="Name of the design.")
    parser.add_argument(
        "--n_designs",
        required=True,
        type=int,
        help="Number of final designs to generate.",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory for BindCraft to write results to.",
    )
    parser.add_argument(
        "--output_json",
        default="settings.json",
        help="Path for the output settings JSON file.",
    )
    parser.add_argument(
        "--hotspot_subsample",
        type=float,
        default=1.0,
        help="Fraction of hotspot residues to subsample (0.0-1.0).",
    )

    args = parser.parse_args()

    try:
        min_len, max_len = map(int, args.binder_length_range.split("-"))
    except ValueError:
        sys.exit("Error: --binder_length_range must be two dash-separated integers.")

    
    if not (0.0 <= args.hotspot_subsample <= 1.0):
        sys.exit("Error: --hotspot_subsample must be between 0.0 and 1.0.")

    
    # UPDATED SECTION START
    
    if args.hotspot_res and args.hotspot_res != "false":
        hotspot_list = [
            res.strip() for res in args.hotspot_res.split(",") if res.strip()
        ]

        
        if not hotspot_list:
            sys.exit(
                "Error: Hotspot residues are empty. Provide values like A123,B456."
            )

       
        for res in hotspot_list:
            if len(res) < 2 or not res[0].isalpha() or not res[1:].isdigit():
                sys.exit(
                    f"Error: Invalid hotspot residue '{res}'. "
                    "Expected format like A123 (chain + residue number)."
                )

        hotspot_residues = args.hotspot_res

        # Subsample logic
        if args.hotspot_subsample < 1.0:
            n_total = len(hotspot_list)
            n_keep = max(1, math.ceil(n_total * args.hotspot_subsample))

            random.seed()
            subsampled_hotspots = random.sample(hotspot_list, n_keep)
            hotspot_residues = ",".join(subsampled_hotspots)

            logging.info(
                f"Subsampled {n_keep}/{n_total} hotspot residues: {hotspot_residues}"
            )
    else:
        hotspot_residues = ""
    
    # UPDATED SECTION END
    

    settings = {
        "design_path": args.output_dir,
        "binder_name": args.design_name,
        "starting_pdb": args.input_pdb,
        "chains": args.target_chains,
        "target_hotspot_residues": hotspot_residues,
        "lengths": [min_len, max_len],
        "number_of_final_designs": args.n_designs,
    }

    with open(args.output_json, "w") as f:
        json.dump(settings, f, indent=4)


if __name__ == "__main__":
    main()