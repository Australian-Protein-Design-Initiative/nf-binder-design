#!/usr/bin/env python
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "biopython>=1.75",
#     "numpy",
# ]
# ///

import sys
import logging
import argparse
import warnings
from pathlib import Path
from typing import List, Optional, TextIO
import numpy as np

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Structure, Residue
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.Polypeptide import is_aa, PPBuilder

logger = logging.getLogger(__name__)

# Suppress PDBConstructionWarning, typically about discontinuous chains or missing atoms.
warnings.filterwarnings("ignore", category=PDBConstructionWarning)
# Also suppress UserWarning from ShrakeRupley for atoms with unknown radii (e.g. H, or other non C,N,O,S atoms)
warnings.filterwarnings("ignore", message="WARNING: Unrecognized atom type")
warnings.filterwarnings(
    "ignore", message="WARNING: Negative sasa result!"
)  # Can happen with odd geometries

# Maximum residue accessible surface area (theoretical) from Tien et al., 2013
# https://doi.org/10.1371/journal.pone.0080635
TIEN_2023_THEORETICAL: dict[str, float] = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLU": 223.0,
    "GLN": 225.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0,
}


def calculate_per_residue_sasa(
    structure: Structure.Structure,
    probe_radius: float = 1.4,
    n_points: int = 100,
) -> List[dict]:
    """
    Calculate SASA for each residue in all chains of the structure.

    Args:
        structure: The PDB structure
        probe_radius: Radius of the probe for SASA calculation (in Angstroms)
        n_points: Number of points for SASA sphere resolution

    Returns:
        List of dictionaries with residue SASA data
    """
    pdb_id = structure.id
    sasa_data = []

    try:
        model = structure[0]

        # Process all chains in the structure
        for chain in model:
            chain_id = chain.get_id()

            # Get all standard amino acid residues
            aa_residues: List[Residue.Residue] = [
                res for res in chain.get_residues() if is_aa(res, standard=True)
            ]

            if not aa_residues:
                logger.debug(
                    f"No standard amino acid residues found in chain '{chain_id}' of PDB ID {pdb_id}."
                )
                continue

            for residue in aa_residues:
                res_id = residue.get_id()
                resnum = res_id[1]  # Residue number
                aa_type = residue.get_resname().upper()

                # Get SASA value
                sasa_angstrom = None
                sasa_percent = None

                if hasattr(residue, "sasa"):
                    sasa_angstrom = round(float(getattr(residue, "sasa")), 2)

                    # Calculate percent SASA
                    standard_sasa = TIEN_2023_THEORETICAL.get(aa_type)
                    if standard_sasa and standard_sasa > 0:
                        sasa_percent = round((sasa_angstrom / standard_sasa) * 100, 2)

                sasa_data.append(
                    {
                        "resnum": resnum,
                        "chain": chain_id,
                        "aa": aa_type,
                        "sasa_angstrom": sasa_angstrom,
                        "sasa_percent": sasa_percent,
                        "probe_radius": probe_radius,
                    }
                )

    except Exception as e:
        logger.error(
            f"Error processing structure {pdb_id}: {e}",
            exc_info=True,
        )

    return sasa_data


def main():
    parser = argparse.ArgumentParser(
        description="Calculate Solvent Accessible Surface Area (SASA) per residue for PDB files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "pdb_files",
        nargs="*",
        type=str,
        help="PDB files to process. Can use shell globs like *.pdb",
    )
    parser.add_argument(
        "--pdbdir",
        type=str,
        help="Directory containing PDB files (scanned recursively). If provided, positional PDB files are ignored.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="-",
        help="Path to the output TSV file. Use '-' for stdout.",
    )
    parser.add_argument(
        "--probe-radius",
        type=float,
        default=1.4,
        help="Radius of the probe for SASA calculation (in Angstroms).",
    )
    parser.add_argument(
        "--n-points",
        type=int,
        default=100,
        help="Number of points for SASA sphere resolution (higher is more precise but slower).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging output.",
    )

    args = parser.parse_args()

    # Set up logging based on verbose flag
    if args.verbose:
        logging.basicConfig(
            stream=sys.stderr,
            level=logging.INFO,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        )
    else:
        logging.basicConfig(
            stream=sys.stderr,
            level=logging.WARNING,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        )

    pdb_files = []

    if args.pdbdir:
        # Process directory if --pdbdir is specified
        input_dir = Path(args.pdbdir)
        if not input_dir.is_dir():
            logger.error(
                f"Input directory '{args.pdbdir}' does not exist or is not a directory."
            )
            sys.exit(1)

        # Find all PDB files recursively
        pdb_files = list(input_dir.rglob("*.pdb"))
        if not pdb_files:
            logger.warning(f"No PDB files found in directory '{args.pdbdir}'.")
            sys.exit(0)
    elif args.pdb_files:
        # Process individual PDB files
        for pdb_path_str in args.pdb_files:
            pdb_path = Path(pdb_path_str)
            if pdb_path.is_file() and pdb_path.suffix.lower() == ".pdb":
                pdb_files.append(pdb_path.resolve())
            elif pdb_path.is_file():
                logger.warning(
                    f"File '{pdb_path_str}' is not a PDB file (expected .pdb extension). Skipping."
                )
            else:
                logger.warning(
                    f"File '{pdb_path_str}' does not exist or is not a file. Skipping."
                )

        if not pdb_files:
            logger.warning("No valid PDB files found in the provided arguments.")
            sys.exit(0)
    else:
        logger.error(
            "Either provide PDB files as positional arguments or use --pdbdir option."
        )
        parser.print_help()
        sys.exit(1)

    logger.info(f"Found {len(pdb_files)} PDB file(s) to process.")

    pdb_parser = PDBParser(QUIET=True)
    sasa_calculator = ShrakeRupley(
        probe_radius=args.probe_radius, n_points=args.n_points
    )

    # Set up output
    output_file: Optional[TextIO] = None
    if args.output != "-":
        output_file_path = Path(args.output)
        try:
            output_file_path.parent.mkdir(parents=True, exist_ok=True)
            output_file = open(output_file_path, "w")
        except IOError as e:
            logger.error(f"Could not open output file {output_file_path}: {e}")
            sys.exit(1)

    # Write header
    header = "pdb_file\tresnum\tchain\taa\tsasa_angstrom\tsasa_percent\tprobe_radius"
    if args.output == "-":
        print(header)
    else:
        if output_file is not None:
            output_file.write(header + "\n")

    processed_count = 0
    for pdb_file in pdb_files:
        pdb_name = pdb_file.name
        logger.info(f"Processing {pdb_name}...")

        try:
            structure = pdb_parser.get_structure(pdb_name, str(pdb_file))
        except Exception as e:
            logger.error(f"Could not parse PDB file {pdb_file}: {e}")
            continue

        if structure is None:
            logger.error(f"Failed to parse structure from {pdb_file}")
            continue

        try:
            sasa_calculator.compute(structure, level="R")
        except Exception as e:
            logger.error(f"SASA computation failed for {pdb_file}: {e}. Skipping.")
            continue

        sasa_data = calculate_per_residue_sasa(
            structure, args.probe_radius, args.n_points
        )

        # Output rows immediately as they are calculated
        for data in sasa_data:
            sasa_angstrom_str = (
                f"{data['sasa_angstrom']:.2f}"
                if data["sasa_angstrom"] is not None
                else ""
            )
            sasa_percent_str = (
                f"{data['sasa_percent']:.2f}"
                if data["sasa_percent"] is not None
                else ""
            )

            output_line = (
                f"{pdb_name}\t{data['resnum']}\t{data['chain']}\t{data['aa']}\t"
                f"{sasa_angstrom_str}\t{sasa_percent_str}\t{data['probe_radius']}"
            )

            if args.output == "-":
                print(output_line)
            else:
                if output_file is not None:
                    output_file.write(output_line + "\n")

        processed_count += 1

        # Clear the structure from memory after processing
        del structure
        del sasa_data

    logger.info(
        f"Processed {processed_count} PDB files successfully out of {len(pdb_files)} found."
    )

    # Close output file if writing to file
    if output_file is not None:
        output_file.close()
        logger.info(f"Output written to {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()
