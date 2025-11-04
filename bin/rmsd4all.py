#!/usr/bin/env python
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "biotite<=1.4",
#     "numpy",
# ]
# ///
#
# Usage: uv run rmsd4all.py <directory> [--all-atom] [--method-rmsd METHOD] [--method-tm METHOD] [--tm-score] [--threads N] [--matrix SCORE] [--output FILE]
# Example: uv run rmsd4all.py ./pdbs --all-atom --method-rmsd blosum62 --method-tm 3di --tm-score --threads 4
# Example: uv run rmsd4all.py ./pdbs --matrix rmsd_pruned

import argparse
import gzip
import logging
import sys
from itertools import combinations, product
from pathlib import Path
from typing import List, Tuple, Union, Optional, Iterable, cast, Dict, Any
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


def _parse_chains_arg(chains_arg: Optional[str]) -> Optional[List[str]]:
    """Parse a comma/space-separated chains argument into a list of chain IDs."""
    if chains_arg is None:
        return None
    if isinstance(chains_arg, str):
        parts = [
            c.strip() for c in chains_arg.replace(" ", ",").split(",") if c.strip()
        ]
        return parts if parts else None
    return None


def _filter_structure_by_chains(
    structure: struc.AtomArray, chains: Optional[Iterable[str]]
) -> struc.AtomArray:
    """Return structure filtered to specified chain IDs. If chains is None, return unchanged."""
    if chains is None:
        return structure
    chain_ids = getattr(structure, "chain_id", None)
    if chain_ids is None:
        return structure
    chain_ids_str = chain_ids.astype(str)
    mask = np.isin(chain_ids_str, np.array(list(chains), dtype=str))
    filtered = structure[mask]
    return cast(struc.AtomArray, filtered)


def _is_structure_file(path: Union[str, Path]) -> bool:
    """Check if a path is a file with a supported structure file extension."""
    path = Path(path)
    if not path.is_file():
        return False
    supported_extensions = [".pdb", ".pdb.gz", ".cif", ".cif.gz", ".mmcif", ".mmcif.gz"]
    return any(str(path).endswith(ext) for ext in supported_extensions)


def _get_aligned_filename(mobile_path: Union[str, Path]) -> str:
    """Generate output filename for aligned structure by removing extensions and adding '_aligned.pdb'."""
    mobile_path = Path(mobile_path)
    # Remove all known structure file extensions (.pdb, .pdb1, .cif, .mmcif, .gz)
    # Iteratively remove suffixes until we get to the base name
    base_name = mobile_path.name
    known_extensions = (".pdb", ".pdb1", ".cif", ".mmcif", ".gz")
    previous_name = None
    while base_name != previous_name:
        previous_name = base_name
        if base_name.endswith(known_extensions):
            base_name = Path(base_name).stem
        else:
            break
    return f"{base_name}_aligned.pdb"


def superimpose_structures(
    structure_path1: Union[str, Path],
    structure_path2: Union[str, Path],
    method: str = "blosum62",
    chains1: Optional[Iterable[str]] = None,
    chains2: Optional[Iterable[str]] = None,
) -> Tuple[struc.AtomArray, struc.AtomArray, Any, np.ndarray, np.ndarray]:
    """
    Superimpose two protein structures using the specified method.

    Args:
        structure_path1: Path to first PDB/CIF file
        structure_path2: Path to second PDB/CIF file
        method: Alignment method - 'blosum62' for sequence-based alignment,
                '3di' or 'pb' for structural alignment with 3Di or Protein Blocks
        chains1: Optional list of chain IDs to use from first structure
        chains2: Optional list of chain IDs to use from second structure

    Returns:
        Tuple of (fixed_ca, fitted_array, transform, fixed_indices, mobile_indices)
    """
    # Read structures
    fixed: struc.AtomArray = read_structure(structure_path1)
    mobile: struc.AtomArray = read_structure(structure_path2)

    # Select chains if specified
    if chains1:
        fixed = _filter_structure_by_chains(fixed, chains1)
    if chains2:
        mobile = _filter_structure_by_chains(mobile, chains2)

    # Select CA atoms
    fixed_ca = fixed[fixed.atom_name == "CA"]
    mobile_ca = mobile[mobile.atom_name == "CA"]

    if fixed_ca.array_length == 0 or mobile_ca.array_length == 0:
        raise ValueError(
            f"No CA atoms found for {Path(structure_path1).name} or {Path(structure_path2).name}"
        )

    # Superimpose structures based on method
    if method in ["3di", "pb"]:
        # Use structural alignment with specified structural alphabet
        fitted, transform, fixed_indices, mobile_indices = (
            struc.superimpose_structural_homologs(
                fixed_ca, mobile_ca, structural_alphabet=method
            )
        )
    elif method == "blosum62":
        # Use sequence-based alignment with BLOSUM62
        fitted, transform, fixed_indices, mobile_indices = struc.superimpose_homologs(
            fixed_ca, mobile_ca
        )
    else:
        raise ValueError(f"Unknown method: {method}. Use 'blosum62', '3di', or 'pb'")

    # Ensure fitted is an AtomArray (not AtomArrayStack)
    if isinstance(fitted, struc.AtomArrayStack):
        fitted_array = cast(struc.AtomArray, fitted[0])
    else:
        fitted_array = cast(struc.AtomArray, fitted)

    return (
        cast(struc.AtomArray, fixed_ca),
        fitted_array,
        transform,
        fixed_indices,
        mobile_indices,
    )


def calculate_scores(
    fixed_ca: struc.AtomArray,
    fitted_array: struc.AtomArray,
    fixed_indices: np.ndarray,
    mobile_indices: np.ndarray,
    structure_path1: Union[str, Path],
    structure_path2: Union[str, Path],
    tm_score: bool = False,
) -> Dict[str, Any]:
    """
    Calculate RMSD and optionally TM-score from superimposed structures.

    Args:
        fixed_ca: Fixed structure (CA atoms)
        fitted_array: Fitted mobile structure
        fixed_indices: Indices of aligned atoms in fixed structure
        mobile_indices: Indices of aligned atoms in mobile structure
        structure_path1: Path to first structure (for error messages)
        structure_path2: Path to second structure (for error messages)
        tm_score: If True, also calculate TM-score in addition to RMSD

    Returns:
        Dictionary with score results
    """
    # Initialize all values
    rmsd_pruned = float("nan")
    n_pairs_rmsd_pruned = 0
    rmsd_all = float("nan")
    n_pairs_rmsd_all = 0
    tm_score_pruned = float("nan")
    n_pairs_tm_score_pruned = 0
    tm_score_all = float("nan")
    n_pairs_tm_score_all = 0

    # Ensure we have AtomArray objects
    if isinstance(fixed_ca, struc.AtomArray) and isinstance(
        fitted_array, struc.AtomArray
    ):
        # Calculate RMSD using the aligned atoms
        try:
            # RMSD over the pruned pairs (outliers excluded)
            fixed_sel = cast(struc.AtomArray, fixed_ca[fixed_indices])
            mobile_sel = cast(struc.AtomArray, fitted_array[mobile_indices])
            rmsd_pruned = struc.rmsd(fixed_sel, mobile_sel)
            n_pairs_rmsd_pruned = len(fixed_indices)
        except Exception as rmsd_pruned_error:
            logging.warning(
                f"rmsd_pruned calculation failed for {Path(structure_path1).name} vs {Path(structure_path2).name}: {rmsd_pruned_error}"
            )

        # Check if structures have the same length for RMSD all calculation
        if len(fixed_ca) != len(fitted_array):
            logging.warning(
                f"Skipping rmsd_all for {Path(structure_path1).name} and {Path(structure_path2).name}, chains are not the same length"
            )
        else:
            try:
                # RMSD over all pairs (including outlier loops)
                rmsd_all = struc.rmsd(cast(struc.AtomArray, fixed_ca), fitted_array)
                n_pairs_rmsd_all = len(fixed_ca)
            except Exception as rmsd_all_error:
                logging.warning(
                    f"rmsd_all calculation failed for {Path(structure_path1).name} vs {Path(structure_path2).name}: {rmsd_all_error}"
                )

        # Optionally calculate TM-score using the aligned atoms
        if tm_score:
            try:
                tm_score_pruned = struc.tm_score(
                    fixed_ca, fitted_array, fixed_indices, mobile_indices
                )
                n_pairs_tm_score_pruned = len(fixed_indices)
            except Exception as tm_score_pruned_error:
                logging.warning(
                    f"tm_score_pruned calculation failed for {Path(structure_path1).name} vs {Path(structure_path2).name}: {tm_score_pruned_error}"
                )

            # Check if structures have the same length for TM-score all calculation
            if len(fixed_ca) != len(fitted_array):
                logging.warning(
                    f"Skipping tm_score_all for {Path(structure_path1).name} and {Path(structure_path2).name}, chains are not the same length"
                )
            else:
                try:
                    tm_score_all = struc.tm_score(
                        fixed_ca,
                        fitted_array,
                        np.arange(len(fixed_ca)),
                        np.arange(len(fitted_array)),
                    )
                    n_pairs_tm_score_all = len(fixed_ca)
                except Exception as tm_score_all_error:
                    logging.warning(
                        f"tm_score_all calculation failed for {Path(structure_path1).name} vs {Path(structure_path2).name}: {tm_score_all_error}"
                    )
    else:
        logging.warning(
            f"Type mismatch for {Path(structure_path1).name} vs {Path(structure_path2).name}: fixed_ca={type(fixed_ca)}, fitted_array={type(fitted_array)}"
        )

    return {
        "rmsd_pruned": rmsd_pruned,
        "n_pairs_rmsd_pruned": n_pairs_rmsd_pruned,
        "rmsd_all": rmsd_all,
        "n_pairs_rmsd_all": n_pairs_rmsd_all,
        "tm_score_pruned": tm_score_pruned,
        "n_pairs_tm_score_pruned": n_pairs_tm_score_pruned,
        "tm_score_all": tm_score_all,
        "n_pairs_tm_score_all": n_pairs_tm_score_all,
    }


def rmsd_pair(
    structure_path1: Union[str, Path],
    structure_path2: Union[str, Path],
    method_rmsd: str = "blosum62",
    method_tm: str = "3di",
    tm_score: bool = False,
    superimpose_chains1: Optional[Iterable[str]] = None,
    superimpose_chains2: Optional[Iterable[str]] = None,
    score_chains1: Optional[Iterable[str]] = None,
    score_chains2: Optional[Iterable[str]] = None,
    output_transformed_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Calculate RMSD and optionally TM-score between two protein structures.

    Args:
        structure_path1: Path to first PDB/CIF file
        structure_path2: Path to second PDB/CIF file
        method_rmsd: Alignment method for RMSD calculation - 'blosum62' for sequence-based alignment,
                    '3di' or 'pb' for structural alignment with 3Di or Protein Blocks
        method_tm: Alignment method for TM-score calculation - 'blosum62' for sequence-based alignment,
                  '3di' or 'pb' for structural alignment with 3Di or Protein Blocks
        tm_score: If True, also calculate TM-score in addition to RMSD
        superimpose_chains1: Optional list of chain IDs to use for superimposition from first structure
        superimpose_chains2: Optional list of chain IDs to use for superimposition from second structure
        score_chains1: Optional list of chain IDs to use for scoring from first structure
        score_chains2: Optional list of chain IDs to use for scoring from second structure
        output_transformed_dir: Optional directory or .pdb file to save transformed mobile structures.
            If a .pdb file, writes directly to that file. If a directory, files are saved with '_aligned.pdb' suffix.

    Returns:
        Dictionary with keys: 'structure1', 'structure2', 'rmsd_pruned', 'n_pairs_rmsd_pruned',
        'rmsd_all', 'n_pairs_rmsd_all', 'tm_score_pruned', 'n_pairs_tm_score_pruned',
        'tm_score_all', 'n_pairs_tm_score_all', 'n_residues_structure1', 'n_residues_structure2'
        (TM-score values are NaN if tm_score=False)
    """
    try:
        # Initialize result dictionary
        result = {
            "structure1": Path(structure_path1).name,
            "structure2": Path(structure_path2).name,
            "superimposed_chains_fixed": (
                ",".join(superimpose_chains1) if superimpose_chains1 else "all"
            ),
            "superimposed_chains_mobile": (
                ",".join(superimpose_chains2) if superimpose_chains2 else "all"
            ),
            "scored_chains_fixed": ",".join(score_chains1) if score_chains1 else "all",
            "scored_chains_mobile": ",".join(score_chains2) if score_chains2 else "all",
            "rmsd_pruned": float("nan"),
            "n_pairs_rmsd_pruned": 0,
            "rmsd_all": float("nan"),
            "n_pairs_rmsd_all": 0,
            "tm_score_pruned": float("nan"),
            "n_pairs_tm_score_pruned": 0,
            "tm_score_all": float("nan"),
            "n_pairs_tm_score_all": 0,
            "n_residues_structure1": 0,
            "n_residues_structure2": 0,
        }

        # Read full structures
        fixed_full: struc.AtomArray = read_structure(structure_path1)
        mobile_full: struc.AtomArray = read_structure(structure_path2)

        # Step 1: Superimpose using superimpose-chains
        try:
            (
                fixed_superimpose_ca,
                mobile_superimpose_ca,
                transform,
                _,
                _,
            ) = superimpose_structures(
                structure_path1,
                structure_path2,
                method_rmsd,
                superimpose_chains1,
                superimpose_chains2,
            )

            # Apply transformation to full mobile structure
            mobile_transformed = transform.apply(mobile_full)

            # Save transformed structure if output directory/file is specified
            if output_transformed_dir is not None:
                try:
                    output_path_obj = Path(output_transformed_dir)
                    # If output path ends with .pdb, write directly to that file
                    if str(output_path_obj).endswith(".pdb"):
                        output_path = output_path_obj
                        # Ensure parent directory exists
                        output_path.parent.mkdir(parents=True, exist_ok=True)
                    else:
                        # Otherwise, treat as directory and create filename
                        output_dir = output_path_obj
                        output_dir.mkdir(parents=True, exist_ok=True)
                        output_filename = _get_aligned_filename(structure_path2)
                        output_path = output_dir / output_filename

                    # Write as PDB file
                    pdb_file = pdb.PDBFile()
                    pdb_file.set_structure(mobile_transformed)
                    with open(output_path, "w") as f:
                        pdb_file.write(f)
                    logging.debug(f"Saved transformed structure to {output_path}")
                except Exception as save_error:
                    logging.warning(
                        f"Failed to save transformed structure for {Path(structure_path2).name}: {save_error}"
                    )

        except Exception as superimpose_error:
            logging.warning(
                f"Superimposition failed for {Path(structure_path1).name} vs {Path(structure_path2).name}: {superimpose_error}"
            )
            return result

        # Step 2: Filter to score-chains and calculate scores
        try:
            # Filter structures to score-chains
            fixed_score = _filter_structure_by_chains(fixed_full, score_chains1)
            mobile_score = _filter_structure_by_chains(
                mobile_transformed, score_chains2
            )

            # Get CA atoms for scoring
            fixed_score_ca = cast(
                struc.AtomArray, fixed_score[fixed_score.atom_name == "CA"]
            )
            mobile_score_ca = cast(
                struc.AtomArray, mobile_score[mobile_score.atom_name == "CA"]
            )

            if fixed_score_ca.array_length == 0 or mobile_score_ca.array_length == 0:
                raise ValueError(
                    f"No CA atoms found in score-chains for {Path(structure_path1).name} or {Path(structure_path2).name}"
                )

            # Calculate residue counts for the structures being compared
            result["n_residues_structure1"] = len(fixed_score_ca)
            result["n_residues_structure2"] = len(mobile_score_ca)

            # Perform alignment on score-chains to get pruned pairs for RMSD
            # We only need the indices, not the transformation (we use the superimpose-chain transformation)
            if method_rmsd in ["3di", "pb"]:
                _, _, fixed_indices_score, mobile_indices_score = (
                    struc.superimpose_structural_homologs(
                        fixed_score_ca, mobile_score_ca, structural_alphabet=method_rmsd
                    )
                )
            elif method_rmsd == "blosum62":
                _, _, fixed_indices_score, mobile_indices_score = (
                    struc.superimpose_homologs(fixed_score_ca, mobile_score_ca)
                )
            else:
                raise ValueError(f"Unknown method: {method_rmsd}")

            # Calculate RMSD scores using mobile_score_ca (already transformed via superimpose-chains)
            # NOT using the fitted_score from the alignment above
            rmsd_scores = calculate_scores(
                fixed_score_ca,
                mobile_score_ca,
                fixed_indices_score,
                mobile_indices_score,
                structure_path1,
                structure_path2,
                tm_score=False,
            )
            result.update(rmsd_scores)

        except Exception as rmsd_error:
            logging.warning(
                f"RMSD calculation failed for {Path(structure_path1).name} vs {Path(structure_path2).name}: {rmsd_error}"
            )

        # Step 3: Calculate TM-score if requested
        if tm_score:
            try:
                # Use same alignment if methods are the same
                if method_rmsd == method_tm:
                    tm_scores = calculate_scores(
                        fixed_score_ca,
                        mobile_score_ca,
                        fixed_indices_score,
                        mobile_indices_score,
                        structure_path1,
                        structure_path2,
                        tm_score=True,
                    )
                else:
                    # Different method - need separate alignment for TM-score to get different pairs
                    # We only need the indices, not the transformation (we use the superimpose-chain transformation)
                    if method_tm in ["3di", "pb"]:
                        _, _, fixed_indices_tm, mobile_indices_tm = (
                            struc.superimpose_structural_homologs(
                                fixed_score_ca,
                                mobile_score_ca,
                                structural_alphabet=method_tm,
                            )
                        )
                    elif method_tm == "blosum62":
                        _, _, fixed_indices_tm, mobile_indices_tm = (
                            struc.superimpose_homologs(fixed_score_ca, mobile_score_ca)
                        )
                    else:
                        raise ValueError(f"Unknown method: {method_tm}")

                    tm_scores = calculate_scores(
                        fixed_score_ca,
                        mobile_score_ca,
                        fixed_indices_tm,
                        mobile_indices_tm,
                        structure_path1,
                        structure_path2,
                        tm_score=True,
                    )

                # Update only TM-score values
                result["tm_score_pruned"] = tm_scores["tm_score_pruned"]
                result["n_pairs_tm_score_pruned"] = tm_scores["n_pairs_tm_score_pruned"]
                result["tm_score_all"] = tm_scores["tm_score_all"]
                result["n_pairs_tm_score_all"] = tm_scores["n_pairs_tm_score_all"]

            except Exception as tm_error:
                logging.warning(
                    f"TM-score calculation failed for {Path(structure_path1).name} vs {Path(structure_path2).name}: {tm_error}"
                )

    except Exception as e:
        logging.warning(f"Error processing {structure_path1} vs {structure_path2}: {e}")

    return result


def rmsd_all_pairs(
    structure_paths_a: List[Union[str, Path]],
    n_processes: int = 1,
    method_rmsd: str = "blosum62",
    method_tm: str = "3di",
    tm_score: bool = False,
    structure_paths_b: Union[None, List[Union[str, Path]]] = None,
    superimpose_chains_a: Optional[Iterable[str]] = None,
    superimpose_chains_b: Optional[Iterable[str]] = None,
    score_chains_a: Optional[Iterable[str]] = None,
    score_chains_b: Optional[Iterable[str]] = None,
    output_transformed_dir: Optional[Union[str, Path]] = None,
) -> List[Dict[str, Any]]:
    """
    Calculate RMSD and optionally TM-score for all pairs of structures.

    Args:
        structure_paths_a: List of paths to PDB/CIF files for rows
        n_processes: Number of processes for multiprocessing
        method_rmsd: Alignment method for RMSD calculation - 'blosum62' for sequence-based alignment,
                    '3di' or 'pb' for structural alignment with 3Di or Protein Blocks
        method_tm: Alignment method for TM-score calculation - 'blosum62' for sequence-based alignment,
                  '3di' or 'pb' for structural alignment with 3Di or Protein Blocks
        tm_score: If True, also calculate TM-score in addition to RMSD
        structure_paths_b: Optional list of paths to PDB/CIF files for columns. If provided,
            perform all-vs-all between A and B (Cartesian product). If None, do all-vs-all within A (combinations).
        superimpose_chains_a: Optional list of chain IDs to use for superimposition from structures in A
        superimpose_chains_b: Optional list of chain IDs to use for superimposition from structures in B
        score_chains_a: Optional list of chain IDs to use for scoring from structures in A
        score_chains_b: Optional list of chain IDs to use for scoring from structures in B
        output_transformed_dir: Optional directory or .pdb file to save transformed mobile structures.
            If a .pdb file, writes directly to that file. If a directory, files are saved with '_aligned.pdb' suffix.

    Returns:
        List of dictionaries with RMSD and/or TM-score results
    """
    if structure_paths_b is None:
        pairs = list(combinations(structure_paths_a, 2))
    else:
        pairs = list(product(structure_paths_a, structure_paths_b))

    if n_processes > 1:
        with Pool(n_processes) as pool:
            args = [
                (
                    p1,
                    p2,
                    method_rmsd,
                    method_tm,
                    tm_score,
                    superimpose_chains_a,
                    (
                        superimpose_chains_b
                        if structure_paths_b is not None
                        else superimpose_chains_a
                    ),
                    score_chains_a,
                    (
                        score_chains_b
                        if structure_paths_b is not None
                        else score_chains_a
                    ),
                    output_transformed_dir,
                )
                for p1, p2 in pairs
            ]
            results = pool.starmap(rmsd_pair, args)
    else:
        results = []
        for p1, p2 in pairs:
            logging.info(f"Processing pair: {Path(p1).name} vs {Path(p2).name}")
            result = rmsd_pair(
                p1,
                p2,
                method_rmsd,
                method_tm,
                tm_score,
                superimpose_chains_a,
                (
                    superimpose_chains_b
                    if structure_paths_b is not None
                    else superimpose_chains_a
                ),
                score_chains_a,
                (score_chains_b if structure_paths_b is not None else score_chains_a),
                output_transformed_dir,
            )
            results.append(result)

    return results


def create_rmsd_matrix(
    results: List[Dict[str, Any]],
    row_files: List[str],
    col_files: Union[None, List[str]] = None,
    score_key: str = "rmsd_pruned",
) -> np.ndarray:
    """Create RMSD/TM-score matrix from pairwise results.

    If col_files is None, create a symmetric square matrix for all-vs-all within row_files.
    If col_files is provided, create a rectangular matrix with rows=row_files and cols=col_files.

    Args:
        results: List of result dictionaries from rmsd_pair
        row_files: List of filenames for rows
        col_files: List of filenames for columns (None for symmetric matrix)
        score_key: Key in result dict to use for matrix values
    """
    if col_files is None:
        n_rows = len(row_files)
        row_to_idx = {filename: i for i, filename in enumerate(row_files)}
        matrix = np.full((n_rows, n_rows), np.nan)
        np.fill_diagonal(matrix, 0.0)
        for result in results:
            file1, file2 = result["structure1"], result["structure2"]
            score_val = result[score_key]
            if file1 in row_to_idx and file2 in row_to_idx:
                i, j = row_to_idx[file1], row_to_idx[file2]
                matrix[i, j] = score_val
                matrix[j, i] = score_val
        return matrix
    else:
        n_rows = len(row_files)
        n_cols = len(col_files)
        row_to_idx = {filename: i for i, filename in enumerate(row_files)}
        col_to_idx = {filename: j for j, filename in enumerate(col_files)}
        matrix = np.full((n_rows, n_cols), np.nan)
        for result in results:
            file1, file2 = result["structure1"], result["structure2"]
            score_val = result[score_key]
            if file1 in row_to_idx and file2 in col_to_idx:
                i = row_to_idx[file1]
                j = col_to_idx[file2]
                matrix[i, j] = score_val
        return matrix


def main():
    """Main CLI function."""
    import os

    cpu_count = max(1, (os.cpu_count() or 1) - 2)

    parser = argparse.ArgumentParser(
        description="Calculate pairwise RMSD between protein structures"
    )
    parser.add_argument(
        "fixed_directory",
        type=Path,
        help="Directory containing fixed/reference PDB/CIF files, or a single PDB/CIF file",
    )
    parser.add_argument(
        "mobile_directory",
        type=Path,
        nargs="?",
        default=None,
        help=(
            "Optional directory containing mobile PDB/CIF files to compare against, or a single PDB/CIF file. "
            "If provided, compute all files in {fixed_directory} vs all files in {mobile_directory}. "
            "If omitted, perform all-against-all comparisons within {fixed_directory}."
        ),
    )
    parser.add_argument(
        "-a",
        "--superimpose-chains",
        type=str,
        default=None,
        help=(
            "Chain IDs to use for superimposition from fixed structures. "
            "Accepts comma or space-separated list, e.g. 'A' or 'A,B'."
        ),
    )
    parser.add_argument(
        "-A",
        "--mobile-superimpose-chains",
        type=str,
        default=None,
        help=(
            "Chain IDs to use for superimposition from mobile structures. "
            "If not provided but mobile_directory is set, defaults to --superimpose-chains."
        ),
    )
    parser.add_argument(
        "-s",
        "--score-chains",
        type=str,
        default=None,
        help=(
            "Chain IDs to use for scoring from fixed structures. "
            "Accepts comma or space-separated list, e.g. 'A' or 'A,B'."
        ),
    )
    parser.add_argument(
        "-S",
        "--mobile-score-chains",
        type=str,
        default=None,
        help=(
            "Chain IDs to use for scoring from mobile structures. "
            "If not provided but mobile_directory is set, defaults to --score-chains."
        ),
    )
    parser.add_argument(
        "--method-rmsd",
        choices=["blosum62", "3di", "pb"],
        default="blosum62",
        help="Alignment method for RMSD calculation: 'blosum62' for sequence-based alignment (superimpose_homologs) with BLOSUM62, "
        "'3di' for TM-align inspired structural alignment (superimpose_structural_homologs) with 3Di alphabet, "
        "'pb' for TM-align inspired structural alignment (superimpose_structural_homologs) with Protein Blocks (default: blosum62)",
    )
    parser.add_argument(
        "--method-tm",
        choices=["blosum62", "3di", "pb"],
        default="3di",
        help="Alignment method for TM-score calculation: 'blosum62' for sequence-based alignment (superimpose_homologs) with BLOSUM62, "
        "'3di' for TM-align inspired structural alignment (superimpose_structural_homologs) with 3Di alphabet, "
        "'pb' for TM-align inspired structural alignment (superimpose_structural_homologs) with Protein Blocks (default: 3di)",
    )
    parser.add_argument(
        "--tm-score",
        action="store_true",
        help="Also calculate TM-score in addition to RMSD using --method-tm alignment",
    )

    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=cpu_count,
        help=f"Number of worker processes (default: {cpu_count})",
    )
    parser.add_argument(
        "--matrix",
        type=str,
        metavar="SCORE",
        help="Output results as a symmetric matrix format using specified score (e.g., rmsd_pruned, rmsd_all, tm_score_pruned, tm_score_all). Default: detailed pairwise TSV",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="-",
        help="Output file (default: stdout)",
    )
    parser.add_argument(
        "--output-transformed",
        type=Path,
        metavar="PATH",
        default=None,
        help=(
            "Directory or .pdb file to save transformed mobile structures after superimposition. "
            "If a .pdb file, writes directly to that file. If a directory, files are saved "
            "with '_aligned.pdb' suffix (e.g., 'structure_aligned.pdb')."
        ),
    )

    args = parser.parse_args()

    # Warn if using blosum62 for TM-score but allow it
    if args.tm_score and args.method_tm == "blosum62":
        logging.warning(
            "Using --method-tm blosum62 and --tm-score together is not conventional, "
            "TM-scores will not reflect common practice due to sequence-based superposition."
        )

    # Log the methods being used
    logging.info(f"Using --method-rmsd {args.method_rmsd}")
    if args.tm_score:
        logging.info(f"Using --method-tm {args.method_tm}")

    # Validate matrix score if specified
    if args.matrix is not None:
        valid_scores = ["rmsd_pruned", "rmsd_all", "tm_score_pruned", "tm_score_all"]
        if args.matrix not in valid_scores:
            logging.error(
                f"Invalid matrix score '{args.matrix}'. Valid options: {', '.join(valid_scores)}"
            )
            sys.exit(1)

        # Check if TM-score is requested but not calculated
        if args.matrix.startswith("tm_score") and not args.tm_score:
            logging.error(f"Matrix score '{args.matrix}' requires --tm-score flag")
            sys.exit(1)

    # If doing TM-score with structural alignment, emit a single informative log once
    if args.tm_score and args.method_tm in ("3di", "pb"):
        logging.info(
            f"Using superimpose_structural_homologs with {args.method_tm} for TM-score calculation"
        )

    # Log thread/process worker count
    logging.info(f"Using {args.threads} worker processes")

    # Parse chain arguments
    superimpose_chains_fixed = _parse_chains_arg(args.superimpose_chains)
    superimpose_chains_mobile = _parse_chains_arg(args.mobile_superimpose_chains)
    score_chains_fixed = _parse_chains_arg(args.score_chains)
    score_chains_mobile = _parse_chains_arg(args.mobile_score_chains)

    # Set defaults for mobile chains if not provided
    if args.mobile_directory is not None:
        if superimpose_chains_mobile is None:
            superimpose_chains_mobile = superimpose_chains_fixed
        if score_chains_mobile is None:
            score_chains_mobile = score_chains_fixed

    # Log chain usage
    if args.mobile_directory is not None:
        logging.info(
            f"Superimposing chains {superimpose_chains_fixed} from {args.fixed_directory} onto chains {superimpose_chains_mobile} from {args.mobile_directory}"
        )
        logging.info(
            f"Scoring using chains {score_chains_fixed} from {args.fixed_directory} against chains {score_chains_mobile} from {args.mobile_directory}"
        )
    else:
        logging.info(
            f"Superimposing using chains {superimpose_chains_fixed} from {args.fixed_directory}"
        )
        logging.info(
            f"Scoring using chains {score_chains_fixed} from {args.fixed_directory}"
        )

    # Find all supported structure files
    supported_extensions = [".pdb", ".pdb.gz", ".cif", ".cif.gz", ".mmcif", ".mmcif.gz"]

    # Check if fixed_directory is a file or directory
    if _is_structure_file(args.fixed_directory):
        pdb_files_fixed = [args.fixed_directory]
        logging.info(f"Using single file: {args.fixed_directory}")
    elif args.fixed_directory.is_dir():
        pdb_files_fixed = []
        for ext in supported_extensions:
            pdb_files_fixed.extend(args.fixed_directory.glob(f"*{ext}"))
        if not pdb_files_fixed:
            logging.error(f"No structure files found in {args.fixed_directory}")
            sys.exit(1)
        logging.info(
            f"Found {len(pdb_files_fixed)} structure files in {args.fixed_directory}"
        )
    else:
        logging.error(
            f"Path {args.fixed_directory} is neither a valid structure file nor a directory"
        )
        sys.exit(1)

    if args.mobile_directory is not None:
        # Check if mobile_directory is a file or directory
        if _is_structure_file(args.mobile_directory):
            pdb_files_mobile = [args.mobile_directory]
            logging.info(f"Using single file: {args.mobile_directory}")
        elif args.mobile_directory.is_dir():
            pdb_files_mobile = []
            for ext in supported_extensions:
                pdb_files_mobile.extend(args.mobile_directory.glob(f"*{ext}"))
            if not pdb_files_mobile:
                logging.error(f"No structure files found in {args.mobile_directory}")
                sys.exit(1)
            logging.info(
                f"Found {len(pdb_files_mobile)} structure files in {args.mobile_directory}"
            )
        else:
            logging.error(
                f"Path {args.mobile_directory} is neither a valid structure file nor a directory"
            )
            sys.exit(1)
        logging.info(
            f"Comparing {len(pdb_files_fixed)} fixed structure(s) against {len(pdb_files_mobile)} mobile structure(s)"
        )
    else:
        pdb_files_mobile = None

    # Calculate pairwise RMSDs or TM-scores
    results = rmsd_all_pairs(
        pdb_files_fixed,
        n_processes=args.threads,
        method_rmsd=args.method_rmsd,
        method_tm=args.method_tm,
        tm_score=args.tm_score,
        structure_paths_b=pdb_files_mobile,
        superimpose_chains_a=superimpose_chains_fixed,
        superimpose_chains_b=superimpose_chains_mobile,
        score_chains_a=score_chains_fixed,
        score_chains_b=score_chains_mobile,
        output_transformed_dir=args.output_transformed,
    )

    # Handle output
    if args.output == "-":
        output_file = sys.stdout
    else:
        output_file = open(args.output, "w")

    try:
        if args.matrix is None:
            # Output detailed results as TSV
            import csv

            writer = csv.writer(output_file, delimiter="\t")

            # Create header - always include RMSD, optionally include TM-score
            header = [
                "structure1",
                "structure2",
                "rmsd_pruned",
                "n_pairs_rmsd_pruned",
                "rmsd_all",
                "n_pairs_rmsd_all",
            ]
            if args.tm_score:
                header.extend(
                    [
                        "tm_score_pruned",
                        "n_pairs_tm_score_pruned",
                        "tm_score_all",
                        "n_pairs_tm_score_all",
                    ]
                )

            # Add residue count columns
            header.extend(
                [
                    "n_residues_structure1",
                    "n_residues_structure2",
                ]
            )

            # Add chain information columns at the end
            header.extend(
                [
                    "superimposed_chains_fixed",
                    "superimposed_chains_mobile",
                    "scored_chains_fixed",
                    "scored_chains_mobile",
                ]
            )

            writer.writerow(header)

            for result in results:
                # Always include RMSD values
                row = [
                    result["structure1"],
                    result["structure2"],
                    (
                        f"{result['rmsd_pruned']:.2f}"
                        if not np.isnan(result["rmsd_pruned"])
                        else "NaN"
                    ),
                    result["n_pairs_rmsd_pruned"],
                    (
                        f"{result['rmsd_all']:.2f}"
                        if not np.isnan(result["rmsd_all"])
                        else "NaN"
                    ),
                    result["n_pairs_rmsd_all"],
                ]

                # Optionally include TM-score values
                if args.tm_score:
                    row.extend(
                        [
                            (
                                f"{result['tm_score_pruned']:.4f}"
                                if not np.isnan(result["tm_score_pruned"])
                                else "NaN"
                            ),
                            result["n_pairs_tm_score_pruned"],
                            (
                                f"{result['tm_score_all']:.4f}"
                                if not np.isnan(result["tm_score_all"])
                                else "NaN"
                            ),
                            result["n_pairs_tm_score_all"],
                        ]
                    )

                # Add residue count columns
                row.extend(
                    [
                        result["n_residues_structure1"],
                        result["n_residues_structure2"],
                    ]
                )

                # Add chain information columns at the end
                row.extend(
                    [
                        result["superimposed_chains_fixed"],
                        result["superimposed_chains_mobile"],
                        result["scored_chains_fixed"],
                        result["scored_chains_mobile"],
                    ]
                )

                writer.writerow(row)
        else:
            # Output matrix format
            row_filenames = sorted([f.name for f in pdb_files_fixed])
            col_filenames = (
                sorted([f.name for f in pdb_files_mobile])
                if pdb_files_mobile is not None
                else None
            )

            # Use the specified score for matrix
            score_key = args.matrix
            matrix = create_rmsd_matrix(
                results, row_filenames, col_filenames, score_key=score_key
            )

            # Output TSV to stdout
            header_cols = col_filenames if col_filenames is not None else row_filenames
            print("\t" + "\t".join(header_cols), file=output_file)
            for i, filename in enumerate(row_filenames):
                row_values = [
                    f"{val:.2f}" if not np.isnan(val) else "NaN" for val in matrix[i]
                ]
                print(f"{filename}\t" + "\t".join(row_values), file=output_file)
    finally:
        if args.output != "-":
            output_file.close()


if __name__ == "__main__":
    main()
