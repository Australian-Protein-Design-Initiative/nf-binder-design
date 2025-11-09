#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "pyyaml",
# ]
# ///

import argparse
import yaml
import os
import shutil
import sys
from pathlib import Path


def stage_input_files(config_path: str, input_files_dir: str, config_dir: str = None):
    """Stage input files from config.yaml at correct relative paths.
    
    Args:
        config_path: Path to the config YAML file
        input_files_dir: Directory containing staged input files (from Nextflow)
        config_dir: Directory containing the config file (for resolving relative paths)
                   If None, uses the directory of config_path
    """
    config_path_obj = Path(config_path)
    
    if config_dir is None:
        config_dir = str(config_path_obj.parent)
    else:
        config_dir = str(Path(config_dir).resolve())
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    if 'entities' not in config:
        return
    
    # Create a mapping of basename to full path in input_files_dir
    input_files_map = {}
    if os.path.exists(input_files_dir):
        for root, dirs, files in os.walk(input_files_dir):
            for file in files:
                basename = os.path.basename(file)
                if basename not in input_files_map:
                    input_files_map[basename] = os.path.join(root, file)
    
    for entity in config['entities']:
        if 'file' in entity and 'path' in entity['file']:
            file_path = entity['file']['path']
            
            # Resolve relative to config directory
            if not os.path.isabs(file_path):
                resolved_path = os.path.join(config_dir, file_path)
            else:
                resolved_path = file_path
            
            # Get relative path from config location
            rel_path = os.path.relpath(resolved_path, config_dir) if not os.path.isabs(file_path) else file_path
            
            # Create directory structure if needed
            rel_dir = os.path.dirname(rel_path)
            if rel_dir:
                os.makedirs(rel_dir, exist_ok=True)
            
            # Find the input file by basename
            file_basename = os.path.basename(file_path)
            if file_basename in input_files_map:
                # Copy to correct relative path
                shutil.copy2(input_files_map[file_basename], rel_path)
            elif os.path.exists(file_basename):
                # Fallback: file might already be in current directory
                shutil.copy2(file_basename, rel_path)
            else:
                # Check if file already exists at the target location
                if not os.path.exists(rel_path):
                    print(f"Warning: Input file not found: {file_path} (basename: {file_basename})", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Stage input files from BoltzGen config YAML at correct relative paths'
    )
    parser.add_argument(
        'config_yaml',
        type=str,
        help='Path to BoltzGen config YAML file'
    )
    parser.add_argument(
        'input_files_dir',
        type=str,
        help='Directory containing staged input files from Nextflow'
    )
    parser.add_argument(
        '--config-dir',
        type=str,
        default=None,
        help='Directory containing config file (for resolving relative paths). Default: directory of config_yaml'
    )
    
    args = parser.parse_args()
    
    try:
        stage_input_files(args.config_yaml, args.input_files_dir, args.config_dir)
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

