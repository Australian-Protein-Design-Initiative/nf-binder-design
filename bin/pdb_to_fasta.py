#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "biopython",
# ]
# ///

import sys
import os
import argparse
import logging
import csv
from typing import List, Optional, Dict, Any
from Bio import PDB
from Bio.PDB.Polypeptide import protein_letters_3to1, is_aa


def get_chain_sequences(pdb_file: str, chain_ids: Optional[List[str]] = None) -> Dict[str, str]:
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


def load_scores_table(scores_file: str) -> Dict[str, Dict[str, Any]]:
    """Load scores from a TSV file.
    
    Args:
        scores_file: Path to the TSV file containing scores
        
    Returns:
        Dictionary mapping filenames to dictionaries of scores
    """
    scores_data = {}
    
    try:
        with open(scores_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            # Check if 'filename' column exists
            if reader.fieldnames is None or 'filename' not in reader.fieldnames:
                logging.error(f"Error: 'filename' column not found in {scores_file}")
                return {}
                
            for row in reader:
                filename = row['filename']
                # Extract just the basename without extension to match PDB files
                basename = os.path.basename(filename)
                basename_no_ext = os.path.splitext(basename)[0]
                scores_data[basename_no_ext] = row
                
    except Exception as e:
        logging.error(f"Error reading scores file {scores_file}: {e}")
        return {}
        
    return scores_data


def format_score_string(scores_data: Dict[str, Any], score_columns: List[str]) -> str:
    """Format scores as a string for inclusion in FASTA headers.
    
    Args:
        scores_data: Dictionary of scores for a specific file
        score_columns: List of column names to include
        
    Returns:
        Formatted string of scores
    """
    score_parts = []
    
    for col in score_columns:
        if col in scores_data and scores_data[col]:
            score_parts.append(f"{col}={scores_data[col]}")
            
    if score_parts:
        return " " + "|".join(score_parts)
    else:
        return ""


def write_fasta(pdb_files: List[str], chain_ids: Optional[List[str]] = None, 
                output_file: Optional[str] = None, scores_file: Optional[str] = None,
                score_columns: Optional[List[str]] = None):
    """Extract sequences from PDB files and write them in FASTA format.
    
    Args:
        pdb_files: List of PDB files to process
        chain_ids: List of chain IDs to extract sequences from. If None, extract all chains.
        output_file: Path to output file. If None or "-", write to stdout.
        scores_file: Path to TSV file containing scores
        score_columns: List of column names to include in the FASTA headers
    """
    output_lines = []
    
    # Load scores if a scores file is provided
    scores_data = {}
    if scores_file:
        scores_data = load_scores_table(scores_file)
        logging.info(f"Loaded scores for {len(scores_data)} structures")
    
    for pdb_file in pdb_files:
        logging.info(f"Processing {pdb_file}")
        
        # Get the basename without extension for the sequence ID
        basename = os.path.basename(pdb_file)
        seq_id_base = os.path.splitext(basename)[0]
        
        # Get scores string if available
        score_string = ""
        if scores_file and seq_id_base in scores_data and score_columns:
            score_string = format_score_string(scores_data[seq_id_base], score_columns)
        
        # Get sequences for the chains
        sequences = get_chain_sequences(pdb_file, chain_ids)
        
        if not sequences:
            logging.warning(f"No sequences found in {pdb_file}")
            continue
            
        if chain_ids:
            # Output only requested chains in the order specified
            for chain_id in chain_ids:
                if chain_id in sequences:
                    seq_id = f"{seq_id_base}"
                    output_lines.append(f">{seq_id}{score_string}")
                    output_lines.append(f"{sequences[chain_id]}")
                else:
                    logging.warning(f"Chain {chain_id} not found in {pdb_file}")
        else:
            # Output all chains, either combined or separate
            chain_ids_sorted = sorted(sequences.keys())
            
            if len(chain_ids_sorted) == 1:
                # Single chain, output without chain ID suffix
                chain_id = chain_ids_sorted[0]
                output_lines.append(f">{seq_id_base}{score_string}|chain={chain_id}")
                output_lines.append(f"{sequences[chain_id]}")
            else:
                # Multiple chains, output combined with * separator
                combined_seq = "*".join([sequences[c] for c in chain_ids_sorted])
                output_lines.append(f">{seq_id_base}{score_string}")
                output_lines.append(f"{combined_seq}")
                
                # Also output individual chains
                for chain_id in chain_ids_sorted:
                    seq_id = f"{seq_id_base}_{chain_id}"
                    output_lines.append(f">{seq_id}{score_string}")
                    output_lines.append(f"{sequences[chain_id]}")
    
    # Write output
    output_text = "\n".join(output_lines) + "\n"
    
    if not output_file or output_file == "-":
        sys.stdout.write(output_text)
    else:
        with open(output_file, "w") as f:
            f.write(output_text)


def main():
    parser = argparse.ArgumentParser(
        description="Extract amino acid sequences from PDB files and output in FASTA format"
    )
    parser.add_argument("files", nargs="+", help="PDB files to process")
    parser.add_argument(
        "-o", "--output", help="Output file (default: stdout)", default="-"
    )
    parser.add_argument(
        "--chains", help="Chain ID(s) to extract (comma-separated for multiple chains)", default=None
    )
    parser.add_argument(
        "--scores-table", help="TSV file containing scores to include in FASTA headers", default=None
    )
    parser.add_argument(
        "--scores", help="Comma-separated list of score columns to include (default: pae_interaction,rg,length)",
        default="pae_interaction,rg,length"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", 
        help="Enable verbose output"
    )
    args = parser.parse_args()

    # Set logging level based on verbose flag
    log_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=log_level, format="%(message)s", stream=sys.stderr)

    # Parse chain IDs if provided
    chain_ids = None
    if args.chains:
        chain_ids = [c.strip() for c in args.chains.split(",")]

    # Parse score columns if scores table is provided
    score_columns = None
    if args.scores_table:
        score_columns = [c.strip() for c in args.scores.split(",")]

    write_fasta(args.files, chain_ids, args.output, args.scores_table, score_columns)


if __name__ == "__main__":
    main() 