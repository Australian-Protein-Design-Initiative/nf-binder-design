#!/usr/bin/env python

import argparse
import yaml
import sys
import re
import os
import glob


def sanitize_for_filename(text: str) -> str:
    """Removes characters that are not filename-friendly."""
    return re.sub(r"[^a-zA-Z0-9_\-.]", "_", text)


def main():
    parser = argparse.ArgumentParser(description="Create a YAML file for BOLTZ input.")
    parser.add_argument("--target_id", required=False, help="ID of the target protein.")
    parser.add_argument(
        "--target_seq", required=False, help="Sequence of the target protein."
    )
    parser.add_argument("--binder_id", required=False, help="ID of the binder protein.")
    parser.add_argument(
        "--binder_seq", required=False, help="Sequence of the binder protein."
    )
    parser.add_argument(
        "--target_msa", required=False, help="Path to the target MSA file."
    )
    parser.add_argument(
        "--binder_msa", required=False, help="Path to the binder MSA file."
    )
    parser.add_argument(
        "--templates", required=False, help="Path to the templates directory."
    )
    parser.add_argument(
        "--output_yaml", required=True, help="Path to the output YAML file."
    )
    parser.add_argument(
        "--target_chain",
        required=False,
        help="Optional chain ID to use for target (default A).",
    )
    parser.add_argument(
        "--binder_chain",
        required=False,
        help="Optional chain ID to use for binder (default B).",
    )
    parser.add_argument(
        "--use_msa_server", action="store_true", help="Omit 'msa: empty' if set."
    )

    args = parser.parse_args()

    # Validate at least one of target or binder provided
    if not (
        (args.target_id and args.target_seq) or (args.binder_id and args.binder_seq)
    ):
        print(
            "Error: Must provide at least one of target_* or binder_* (id and seq)",
            file=sys.stderr,
        )
        sys.exit(2)

    sequences_list = []

    if args.target_id and args.target_seq:
        target_chain = args.target_chain if args.target_chain else "A"
        target_protein = {"id": [target_chain], "sequence": args.target_seq}
        if args.target_msa:
            target_protein["msa"] = args.target_msa
        sequences_list.append({"protein": target_protein})

    if args.binder_id and args.binder_seq:
        binder_chain = args.binder_chain if args.binder_chain else "B"
        binder_protein = {"id": [binder_chain], "sequence": args.binder_seq}
        if args.binder_msa:
            binder_protein["msa"] = args.binder_msa
        sequences_list.append({"protein": binder_protein})

    data = {
        "version": 1,
        "sequences": sequences_list,
    }

    if args.templates:
        # Find all .cif and .pdb files in the templates directory
        cif_files = glob.glob(os.path.join(args.templates, "*.cif"))
        pdb_files = glob.glob(os.path.join(args.templates, "*.pdb"))

        templates_list = []

        for cif_file in sorted(cif_files):
            templates_list.append({"cif": cif_file})

        for pdb_file in sorted(pdb_files):
            templates_list.append({"pdb": pdb_file})

        if templates_list:
            data["templates"] = templates_list

    with open(args.output_yaml, "w") as f:
        yaml.dump(data, f, sort_keys=False)


if __name__ == "__main__":
    main()
