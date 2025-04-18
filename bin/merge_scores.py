#!/usr/bin/env python
# /// script
# requires-python = ">=3.7"
# dependencies = [
#     "pandas",
# ]
# ///

import argparse
import os
import sys
import pandas as pd
import logging
from typing import List, Optional, TextIO, Union

# Set up logging to stderr
logging.basicConfig(
    level=logging.WARNING,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stderr
)
logger = logging.getLogger(__name__)

def merge_scores(
    af2_scores_file: str, 
    shape_scores_file: str, 
    output_file: Union[str, TextIO, None] = None,
    sort_by: str = 'pae_interaction',
    verbose: bool = False
) -> pd.DataFrame:
    """
    Merge AlphaFold2 scores with shape scores into a single DataFrame.
    
    Args:
        af2_scores_file: Path to the AlphaFold2 scores TSV file
        shape_scores_file: Path to the shape scores TSV file
        output_file: Path to write the combined TSV file or '-' for stdout
        sort_by: Column name to sort the merged DataFrame by (default: 'pae_interaction')
        verbose: Whether to enable verbose logging (default: False)
        
    Returns:
        The merged DataFrame
    """
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    
    # Read the input files
    logger.info(f"Reading AF2 scores from {af2_scores_file}")
    af2_scores = pd.read_csv(af2_scores_file, sep='\t')
    
    logger.info(f"Reading shape scores from {shape_scores_file}")
    shape_scores = pd.read_csv(shape_scores_file, sep='\t')
    
    # Remove path and .pdb extension from filenames
    logger.info("Extracting basename and removing .pdb extension from filenames for matching")
    shape_scores['match_key'] = shape_scores['filename'].apply(
        lambda x: os.path.splitext(os.path.basename(x))[0]
    )
    
    # Merge the dataframes
    logger.info("Merging dataframes on description and filename (without path or .pdb)")
    merged_df = pd.merge(
        af2_scores,
        shape_scores,
        left_on='description',
        right_on='match_key',
        how='left'
    )
    
    # Drop the temporary matching column
    merged_df = merged_df.drop(columns=['match_key'])
    
    # Sort the merged dataframe
    logger.info(f"Sorting merged dataframe by {sort_by}")
    merged_df = merged_df.sort_values(by=sort_by, ascending=True)
    
    # Move the sort_by column to the first position
    cols = [sort_by] + [col for col in merged_df.columns if col != sort_by]
    merged_df = merged_df[cols]
    
    # Write to output file if specified
    if output_file:
        if output_file == '-':
            logger.info("Writing combined scores to stdout")
            merged_df.to_csv(sys.stdout, sep='\t', index=False)
        else:
            logger.info(f"Writing combined scores to {output_file}")
            merged_df.to_csv(output_file, sep='\t', index=False)
    
    return merged_df

def parse_args():
    parser = argparse.ArgumentParser(description='Merge AF2 scores with shape scores')
    parser.add_argument('af2_scores', help='Path to AlphaFold2 scores TSV file')
    parser.add_argument('shape_scores', help='Path to shape scores TSV file')
    parser.add_argument('-o', '--output', default='-',
                        help='Output filename for combined scores (default: "-" for stdout)')
    parser.add_argument('--sort-by', default='pae_interaction',
                        help='Column name to sort the merged DataFrame by (default: "pae_interaction")')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose logging')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    merge_scores(args.af2_scores, args.shape_scores, args.output, args.sort_by, args.verbose)
    logger.info("Completed successfully") 