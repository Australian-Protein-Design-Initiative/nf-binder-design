#!/usr/bin/env python

import argparse
import yaml
import sys
import re
import os
import glob

def sanitize_for_filename(text: str) -> str:
    """Removes characters that are not filename-friendly."""
    return re.sub(r'[^a-zA-Z0-9_\-.]', '_', text)

def main():
    parser = argparse.ArgumentParser(description="Create a YAML file for BOLTZ input.")
    parser.add_argument("--target_id", required=True, help="ID of the target protein.")
    parser.add_argument("--target_seq", required=True, help="Sequence of the target protein.")
    parser.add_argument("--binder_id", required=True, help="ID of the binder protein.")
    parser.add_argument("--binder_seq", required=True, help="Sequence of the binder protein.")
    parser.add_argument("--target_msa", required=False, help="Path to the target MSA file.")
    parser.add_argument("--binder_msa", required=False, help="Path to the binder MSA file.")
    parser.add_argument("--templates", required=False, help="Path to the templates directory.")
    parser.add_argument("--output_yaml", required=True, help="Path to the output YAML file.")
    parser.add_argument("--use_msa_server", action="store_true", help="Omit 'msa: empty' if set.")

    args = parser.parse_args()

    target_protein = {
        'id': ['A'],
        'sequence': args.target_seq
    }

    binder_protein = {
        'id': ['B'],
        'sequence': args.binder_seq
    }
    
    # if not args.use_msa_server:
    #     target_protein['msa'] = 'empty'
    #     binder_protein['msa'] = 'empty'
    if args.target_msa:
        target_protein['msa'] = args.target_msa
    if args.binder_msa:
        binder_protein['msa'] = args.binder_msa

    data = {
        'version': 1,
        'sequences': [
            {'protein': target_protein},
            {'protein': binder_protein}
        ]
    }

    if args.templates:
        # Find all .cif files in the templates directory
        cif_files = glob.glob(os.path.join(args.templates, '*.cif'))
        if cif_files:
            data['templates'] = []
            for cif_file in sorted(cif_files):
                data['templates'].append({'cif': cif_file})

    with open(args.output_yaml, 'w') as f:
        yaml.dump(data, f, sort_keys=False)

if __name__ == "__main__":
    main() 