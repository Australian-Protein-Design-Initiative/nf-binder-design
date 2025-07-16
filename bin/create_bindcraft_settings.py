#!/usr/bin/env python

import argparse
import json
import os
import sys

# TODO: If we have a --contigs arg in the same format as used for RFDiffusion (main.nf), 
#       extract the target chain and binder lengths from that and trim the PDB to tne contigs

def main():
    parser = argparse.ArgumentParser(description="Create a BindCraft settings JSON file.")
    parser.add_argument("--input_pdb", required=True, help="Input PDB file for the target.")
    parser.add_argument("--hotspot_res", required=True, help="Hotspot residues, e.g., 'A29,A32,A34'.")
    parser.add_argument("--target_chains", required=True, help="The target chain(s), e.g., 'A' or 'A,B'.")
    parser.add_argument("--binder_length_range", required=True, help="Dash-separated min and max binder lengths, e.g., '65-170'.")
    parser.add_argument("--design_name", required=True, help="Name of the design.")
    parser.add_argument("--n_designs", required=True, type=int, help="Number of final designs to generate.")
    parser.add_argument("--output_dir", required=True, help="Directory for BindCraft to write results to.")
    parser.add_argument("--output_json", default="settings.json", help="Path for the output settings JSON file.")

    args = parser.parse_args()

    try:
        min_len, max_len = map(int, args.binder_length_range.split('-'))
    except ValueError:
        sys.exit("Error: --binder_length_range must be two dash-separated integers.")

    # Convert relative input_pdb to absolute path, as BindCraft may be run in a different working directory
    # within a container.
    settings = {
        "design_path": args.output_dir,
        "binder_name": args.design_name,
        "starting_pdb": args.input_pdb,
        "chains": args.target_chains,   # comma-separated list of chains
        "target_hotspot_residues": args.hotspot_res,
        "lengths": [min_len, max_len],
        "number_of_final_designs": args.n_designs,
    }

    with open(args.output_json, 'w') as f:
        json.dump(settings, f, indent=4)

if __name__ == "__main__":
    main()
