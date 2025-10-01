#!/usr/bin/env python
# /// script
# dependencies = [
#     "biotite",
#     "numpy",
# ]
# ///
#
# Usage: uv run rmsd4all.py <directory> [--all-atom] [--method METHOD] [--tm-score] [--processes N]
# Example: uv run rmsd4all.py ./pdbs --all-atom --method 3di --tm-score --processes 4

import argparse
import gzip
import logging
import sys
from itertools import combinations
from pathlib import Path
from typing import List, Tuple, Union
from multiprocessing import Pool

import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    stream=sys.stderr,
)


def read_structure(file_path: Union[str, Path]):
    """Read structure from PDB or CIF file, handling compression."""
    file_path = Path(file_path)

    # Handle compressed files
    if file_path.suffix == ".gz":
        opener = gzip.open
        actual_suffix = file_path.suffixes[-2] if len(file_path.suffixes) >= 2 else ""
    else:
        opener = open
        actual_suffix = file_path.suffix

    # Read based on file type
    if actual_suffix in [".pdb", ".pdb1"]:
        with opener(file_path, "rt") as f:
            pdb_file = pdb.PDBFile.read(f)
            structure = pdb.get_structure(pdb_file)
    elif actual_suffix in [".cif", ".mmcif"]:
        with opener(file_path, "rt") as f:
            cif_file = pdbx.CIFFile.read(f)
            structure = pdbx.get_structure(cif_file)
    else:
        raise ValueError(f"Unsupported file format: {file_path}")

    # Return first model if multiple models exist
    if isinstance(structure, struc.AtomArrayStack):
        structure = structure[0]

    # Ensure we have an AtomArray
    if not isinstance(structure, struc.AtomArray):
        raise ValueError(f"Expected AtomArray, got {type(structure)}")

    return structure


def rmsd_pair(
    pdb_path1: Union[str, Path],
    pdb_path2: Union[str, Path],
    all_atom: bool = False,
    method: str = "blosum62",
    tm_score: bool = False,
) -> Tuple[str, str, float]:
    """
    Calculate RMSD or TM-score between two protein structures.

    Args:
        pdb_path1: Path to first PDB file
        pdb_path2: Path to second PDB file
        all_atom: If True, use all atoms; if False, use only CA atoms
        method: Alignment method - 'blosum62' for sequence-based alignment,
                '3di' or 'pb' for structural alignment with 3Di or Protein Blocks
        tm_score: If True, calculate TM-score instead of RMSD

    Returns:
        Tuple of (pdb1_basename, pdb2_basename, score)
    """
    try:
        # Read structures
        fixed = read_structure(pdb_path1)
        mobile = read_structure(pdb_path2)

        # Filter to CA atoms if requested
        if not all_atom:
            fixed_ca = fixed[fixed.atom_name == "CA"]
            mobile_ca = mobile[mobile.atom_name == "CA"]
        else:
            fixed_ca = fixed
            mobile_ca = mobile

        # Superimpose structures based on method
        if method in ["3di", "pb"]:
            # Use structural alignment with specified structural alphabet
            if tm_score:
                logging.info(
                    f"Using superimpose_structural_homologs with {method} for TM-score calculation"
                )
            fitted, transform, fixed_indices, mobile_indices = (
                struc.superimpose_structural_homologs(
                    fixed_ca, mobile_ca, structural_alphabet=method
                )
            )
        elif method == "blosum62":
            # Use sequence-based alignment with BLOSUM62
            fitted, transform, fixed_indices, mobile_indices = (
                struc.superimpose_homologs(fixed_ca, mobile_ca)
            )
        else:
            raise ValueError(
                f"Unknown method: {method}. Use 'blosum62', '3di', or 'pb'"
            )

        # Calculate RMSD or TM-score only on the aligned atoms
        if len(fixed_indices) > 0 and len(mobile_indices) > 0:
            # Ensure fitted is an AtomArray (not AtomArrayStack)
            if isinstance(fitted, struc.AtomArrayStack):
                fitted_array = fitted[0]
            else:
                fitted_array = fitted

            # Ensure we have AtomArray objects
            if isinstance(fixed_ca, struc.AtomArray) and isinstance(
                fitted_array, struc.AtomArray
            ):
                try:
                    if tm_score:
                        # Calculate TM-score using the aligned atoms
                        score_value = struc.tm_score(
                            fixed_ca, fitted_array, fixed_indices, mobile_indices
                        )
                    else:
                        # Calculate RMSD using the aligned atoms
                        score_value = struc.rmsd(
                            fixed_ca[fixed_indices], fitted_array[mobile_indices]
                        )
                except Exception as score_error:
                    score_type = "TM-score" if tm_score else "RMSD"
                    logging.warning(
                        f"{score_type} calculation failed for {Path(pdb_path1).name} vs {Path(pdb_path2).name}: {score_error}"
                    )
                    score_value = float("nan")
            else:
                logging.warning(
                    f"Type mismatch for {Path(pdb_path1).name} vs {Path(pdb_path2).name}: fixed_ca={type(fixed_ca)}, fitted_array={type(fitted_array)}"
                )
                score_value = float("nan")
        else:
            # No aligned atoms found
            logging.warning(
                f"No aligned atoms found for {Path(pdb_path1).name} vs {Path(pdb_path2).name}"
            )
            score_value = float("nan")

        return (Path(pdb_path1).name, Path(pdb_path2).name, score_value)

    except Exception as e:
        logging.error(f"Error processing {pdb_path1} vs {pdb_path2}: {e}")
        return (Path(pdb_path1).name, Path(pdb_path2).name, float("nan"))


def rmsd_all_pairs(
    pdb_paths: List[Union[str, Path]],
    all_atom: bool = False,
    n_processes: int = 1,
    method: str = "blosum62",
    tm_score: bool = False,
) -> List[Tuple[str, str, float]]:
    """
    Calculate RMSD or TM-score for all pairs of structures.

    Args:
        pdb_paths: List of paths to PDB files
        all_atom: If True, use all atoms; if False, use only CA atoms
        n_processes: Number of processes for multiprocessing
        method: Alignment method - 'blosum62' for sequence-based alignment,
                '3di' or 'pb' for structural alignment with 3Di or Protein Blocks
        tm_score: If True, calculate TM-score instead of RMSD

    Returns:
        List of tuples (fixed_pdb_path, mobile_pdb_path, score)
    """
    pairs = list(combinations(pdb_paths, 2))

    if n_processes > 1:
        with Pool(n_processes) as pool:
            args = [(p1, p2, all_atom, method, tm_score) for p1, p2 in pairs]
            results = pool.starmap(rmsd_pair, args)
    else:
        results = []
        for p1, p2 in pairs:
            logging.info(f"Processing pair: {Path(p1).name} vs {Path(p2).name}")
            result = rmsd_pair(p1, p2, all_atom, method, tm_score)
            results.append(result)

    return results


def create_rmsd_matrix(
    results: List[Tuple[str, str, float]], all_files: List[str]
) -> np.ndarray:
    """Create symmetric RMSD matrix from pairwise results."""
    n_files = len(all_files)
    file_to_idx = {filename: i for i, filename in enumerate(all_files)}

    # Initialize matrix with zeros on diagonal, NaN elsewhere
    matrix = np.full((n_files, n_files), np.nan)
    np.fill_diagonal(matrix, 0.0)

    # Fill in pairwise distances
    for file1, file2, rmsd_val in results:
        if file1 in file_to_idx and file2 in file_to_idx:
            i, j = file_to_idx[file1], file_to_idx[file2]
            matrix[i, j] = rmsd_val
            matrix[j, i] = rmsd_val  # Symmetric

    return matrix


def main():
    """Main CLI function."""
    import os

    cpu_count = max(1, (os.cpu_count() or 1) - 2)

    parser = argparse.ArgumentParser(
        description="Calculate pairwise RMSD between protein structures"
    )
    parser.add_argument(
        "directory", type=Path, help="Directory containing PDB/CIF files"
    )
    parser.add_argument(
        "--all-atom", action="store_true", help="Use all atoms instead of just CA atoms"
    )
    parser.add_argument(
        "--method",
        choices=["blosum62", "3di", "pb"],
        default="blosum62",
        help="Alignment method: 'blosum62' for sequence-based alignment (superimpose_homologs) with BLOSUM62, "
        "'3di' for TM-align inspired structural alignment (superimpose_structural_homologs) with 3Di alphabet, "
        "'pb' for TM-align inspired structural alignment (superimpose_structural_homologs) with Protein Blocks (default: blosum62)",
    )
    parser.add_argument(
        "--tm-score",
        action="store_true",
        help="Calculate TM-score instead of RMSD (defaults to --method 3di if not specified)",
    )

    parser.add_argument(
        "--processes",
        "-p",
        type=int,
        default=cpu_count,
        help=f"Number of parallel processes (default: number of CPUs: {cpu_count})",
    )

    args = parser.parse_args()

    # Set default method to 3di when tm-score is specified
    if args.tm_score and args.method == "blosum62":
        args.method = "3di"
        logging.info("TM-score specified: defaulting to --method 3di")

    # Find all supported structure files
    supported_extensions = [".pdb", ".pdb.gz", ".cif", ".cif.gz", ".mmcif", ".mmcif.gz"]
    pdb_files = []

    for ext in supported_extensions:
        pdb_files.extend(args.directory.glob(f"*{ext}"))

    if not pdb_files:
        logging.error(f"No structure files found in {args.directory}")
        sys.exit(1)

    logging.info(f"Found {len(pdb_files)} structure files")

    # Calculate pairwise RMSDs or TM-scores
    results = rmsd_all_pairs(
        pdb_files,
        all_atom=args.all_atom,
        n_processes=args.processes,
        method=args.method,
        tm_score=args.tm_score,
    )

    # Create matrix and output TSV
    filenames = [f.name for f in pdb_files]
    matrix = create_rmsd_matrix(results, filenames)

    # Output TSV to stdout
    print("\t" + "\t".join(filenames))
    for i, filename in enumerate(filenames):
        row_values = [f"{val:.2f}" if not np.isnan(val) else "NaN" for val in matrix[i]]
        print(f"{filename}\t" + "\t".join(row_values))


if __name__ == "__main__":
    main()
