#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "pyyaml",
# ]
# ///

import argparse
import yaml
import sys
import os
from pathlib import Path
from typing import List


def extract_file_paths(config_path: str, config_dir: str = None) -> List[str]:
    """Extract all file paths from entities[].file.path in BoltzGen config YAML.
    
    Args:
        config_path: Path to the config YAML file
        config_dir: Directory containing the config file (for resolving relative paths)
                   If None, uses the directory of config_path
    
    Returns:
        List of file paths found in entities[].file.path
    """
    config_path_obj = Path(config_path)
    
    if config_dir is None:
        config_dir = str(config_path_obj.parent)
    else:
        config_dir = str(Path(config_dir).resolve())
    
    if not config_path_obj.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    file_paths = []
    
    if 'entities' not in config:
        return file_paths
    
    for entity in config['entities']:
        if 'file' in entity and 'path' in entity['file']:
            file_path = entity['file']['path']
            # Resolve relative to config directory
            if not os.path.isabs(file_path):
                resolved_path = os.path.join(config_dir, file_path)
            else:
                resolved_path = file_path
            file_paths.append(resolved_path)
    
    return file_paths


def main():
    parser = argparse.ArgumentParser(
        description='Extract file paths from BoltzGen config YAML entities[].file.path'
    )
    parser.add_argument(
        'config_yaml',
        type=str,
        help='Path to BoltzGen config YAML file'
    )
    parser.add_argument(
        '--config-dir',
        type=str,
        default=None,
        help='Directory containing config file (for resolving relative paths). Default: directory of config_yaml'
    )
    parser.add_argument(
        '--format',
        type=str,
        choices=['list', 'json', 'tsv'],
        default='list',
        help='Output format: list (one path per line), json, or tsv'
    )
    
    args = parser.parse_args()
    
    try:
        file_paths = extract_file_paths(args.config_yaml, args.config_dir)
        
        if args.format == 'list':
            for path in file_paths:
                print(path)
        elif args.format == 'json':
            import json
            print(json.dumps(file_paths, indent=2))
        elif args.format == 'tsv':
            print('\n'.join(file_paths))
        
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

