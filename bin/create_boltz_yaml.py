#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "biopython",
#     "pyyaml",
# ]
# ///

import argparse
import yaml
import sys
import re
import os
import glob
from typing import Optional
from Bio import PDB  # type: ignore
from Bio.PDB.Polypeptide import protein_letters_3to1, is_aa  # type: ignore
from Bio import SeqIO  # type: ignore


def sanitize_for_filename(text: str) -> str:
    """Removes characters that are not filename-friendly."""
    return re.sub(r"[^a-zA-Z0-9_\-.]", "_", text)


def get_chain_sequences(pdb_file: str, chain_ids: Optional[list] = None) -> dict:
    """Extract one-letter amino acid sequences for specified chains in a PDB file.

    Args:
        pdb_file: Path to the PDB file
        chain_ids: List of chain IDs to extract sequences from. If None, extract all chains.

    Returns:
        Dictionary mapping chain IDs to their amino acid sequences
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    sequences = {}

    # Process each chain in the structure
    for chain in structure.get_chains():
        chain_id = chain.id.strip()
        if not chain_id:
            continue

        # Skip if specific chains requested and this one isn't in the list
        if chain_ids and chain_id not in chain_ids:
            continue

        # Extract amino acid sequence
        sequence = ""
        for residue in chain:
            if is_aa(residue):
                try:
                    # Use the protein_letters_3to1 dictionary
                    res_name = residue.get_resname().upper()
                    one_letter = protein_letters_3to1.get(res_name, "X")
                    sequence += one_letter
                except (KeyError, AttributeError):
                    # Skip non-standard amino acids
                    sequence += "X"

        if sequence:
            sequences[chain_id] = sequence

    return sequences


def get_sequence_from_fasta(fasta_file: str, sequence_id: Optional[str] = None) -> str:
    """Extract a sequence from a FASTA file by ID.

    Args:
        fasta_file: Path to the FASTA file
        sequence_id: ID of the sequence to extract. If None, returns the first sequence.

    Returns:
        The sequence string
    """
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file {fasta_file} not found")

    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {fasta_file}")

    if sequence_id is None:
        print(
            f"No {fasta_file.split('/')[-1].replace('.fasta', '')}_fasta_id specified, taking first sequence",
            file=sys.stderr,
        )
        return str(records[0].seq)

    for record in records:
        if record.id == sequence_id or record.description == sequence_id:
            return str(record.seq)

    raise ValueError(f"Sequence with ID '{sequence_id}' not found in {fasta_file}")


def main():
    parser = argparse.ArgumentParser(description="Create a YAML file for BOLTZ input.")
    parser.add_argument("--target_id", required=False, help="ID of the target protein.")
    parser.add_argument(
        "--target_seq", required=False, help="Sequence of the target protein."
    )
    parser.add_argument("--binder_id", required=True, help="ID of the binder protein.")
    parser.add_argument(
        "--binder_seq", required=False, help="Sequence of the binder protein."
    )
    parser.add_argument(
        "--target_from_pdb",
        required=False,
        help="Path to PDB file to extract target sequence from.",
    )
    parser.add_argument(
        "--binder_from_pdb",
        required=False,
        help="Path to PDB file to extract binder sequence from.",
    )
    parser.add_argument(
        "--target_from_fasta",
        required=False,
        help="Path to FASTA file to extract target sequence from.",
    )
    parser.add_argument(
        "--target_fasta_id",
        required=False,
        help="Sequence ID to extract from FASTA file.",
    )
    parser.add_argument(
        "--binder_from_fasta",
        required=False,
        help="Path to FASTA file to extract binder sequence from.",
    )
    parser.add_argument(
        "--binder_fasta_id",
        required=False,
        help="Sequence ID to extract from FASTA file.",
    )
    parser.add_argument(
        "--target_chains",
        required=False,
        help="Comma-separated list of chain IDs for target (default: all non-A chains).",
    )
    parser.add_argument(
        "--binder_chains",
        required=False,
        help="Comma-separated list of chain IDs for binder (default: A).",
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
        "--use_msa_server", action="store_true", help="Omit 'msa: empty' if set."
    )

    args = parser.parse_args()

    # Determine if target is provided
    target_sources = [args.target_seq, args.target_from_pdb, args.target_from_fasta]
    target_sources = [s for s in target_sources if s]  # Remove None values
    has_target = len(target_sources) > 0

    # Validate that only one target source is specified (if any)
    if len(target_sources) > 1:
        print(
            "Error: Only one of --target_seq, --target_from_pdb, or --target_from_fasta can be specified",
            file=sys.stderr,
        )
        sys.exit(1)

    # Validate that only one binder source is specified
    binder_sources = [args.binder_seq, args.binder_from_pdb, args.binder_from_fasta]
    binder_sources = [s for s in binder_sources if s]  # Remove None values
    if len(binder_sources) > 1:
        print(
            "Error: Only one of --binder_seq, --binder_from_pdb, or --binder_from_fasta can be specified",
            file=sys.stderr,
        )
        sys.exit(1)

    # Extract sequences from PDB files if provided
    target_seq = args.target_seq
    binder_seq = args.binder_seq

    # Only extract target sequence if target is provided
    if has_target:
        if args.target_from_fasta:
            try:
                target_seq = get_sequence_from_fasta(
                    args.target_from_fasta, args.target_fasta_id
                )
            except (FileNotFoundError, ValueError) as e:
                print(f"Error: {e}", file=sys.stderr)
                sys.exit(1)

        elif args.target_from_pdb:
            if not os.path.exists(args.target_from_pdb):
                print(
                    f"Error: Target PDB file {args.target_from_pdb} not found",
                    file=sys.stderr,
                )
                sys.exit(1)

            # Determine target chains
            if args.target_chains:
                target_chain_ids = [c.strip() for c in args.target_chains.split(",")]
            else:
                # Default: all chains except A
                all_sequences = get_chain_sequences(args.target_from_pdb)
                target_chain_ids = [c for c in all_sequences.keys() if c != "A"]

            target_sequences = get_chain_sequences(
                args.target_from_pdb, target_chain_ids
            )
            if not target_sequences:
                print(
                    f"Error: No target sequences found in {args.target_from_pdb} for chains {target_chain_ids}",
                    file=sys.stderr,
                )
                sys.exit(1)

            # Combine multiple chains with * separator
            target_seq = "*".join(
                [target_sequences[c] for c in sorted(target_sequences.keys())]
            )

    if args.binder_from_pdb:
        if not os.path.exists(args.binder_from_pdb):
            print(
                f"Error: Binder PDB file {args.binder_from_pdb} not found",
                file=sys.stderr,
            )
            sys.exit(1)

        # Determine binder chains
        if args.binder_chains:
            binder_chain_ids = [c.strip() for c in args.binder_chains.split(",")]
        else:
            binder_chain_ids = ["A"]  # Default: chain A

        binder_sequences = get_chain_sequences(args.binder_from_pdb, binder_chain_ids)
        if not binder_sequences:
            print(
                f"Error: No binder sequences found in {args.binder_from_pdb} for chains {binder_chain_ids}",
                file=sys.stderr,
            )
            sys.exit(1)

        # Combine multiple chains with * separator
        binder_seq = "*".join(
            [binder_sequences[c] for c in sorted(binder_sequences.keys())]
        )

    if args.binder_from_fasta:
        try:
            binder_seq = get_sequence_from_fasta(
                args.binder_from_fasta, args.binder_fasta_id
            )
        except (FileNotFoundError, ValueError) as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    # Validate that we have sequences
    if has_target:
        if not target_seq:
            print(
                "Error: Failed to extract target sequence from PDB or FASTA",
                file=sys.stderr,
            )
            sys.exit(1)

    if not binder_seq:
        if not args.binder_from_pdb and not args.binder_from_fasta:
            print(
                "Error: Binder sequence is required. Provide --binder_seq, --binder_from_pdb, or --binder_from_fasta",
                file=sys.stderr,
            )
            sys.exit(1)
        else:
            print(
                "Error: Failed to extract binder sequence from PDB or FASTA",
                file=sys.stderr,
            )
            sys.exit(1)

    # Create YAML structure based on mode
    if not has_target:
        # Binder-only mode: single protein entry
        binder_protein = {"id": ["A"], "sequence": binder_seq}

        if args.binder_msa:
            binder_protein["msa"] = args.binder_msa

        data = {"version": 1, "sequences": [{"protein": binder_protein}]}
    else:
        # Complex mode: target and binder
        target_protein = {"id": ["A"], "sequence": target_seq}

        binder_protein = {"id": ["B"], "sequence": binder_seq}

        if args.target_msa:
            target_protein["msa"] = args.target_msa
        if args.binder_msa:
            binder_protein["msa"] = args.binder_msa

        data = {
            "version": 1,
            "sequences": [{"protein": target_protein}, {"protein": binder_protein}],
        }

    if args.templates:
        # Find all .cif files in the templates directory
        cif_files = glob.glob(os.path.join(args.templates, "*.cif"))
        if cif_files:
            data["templates"] = []
            for cif_file in sorted(cif_files):
                data["templates"].append({"cif": cif_file})

    with open(args.output_yaml, "w") as f:
        yaml.dump(data, f, sort_keys=False)


if __name__ == "__main__":
    main()
