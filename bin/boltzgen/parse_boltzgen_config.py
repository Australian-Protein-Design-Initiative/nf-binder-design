#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "pyyaml",
# ]
# ///

import argparse
import os
import sys
from pathlib import Path
from typing import List, Optional, Tuple, Union

import yaml


def _path_to_list(path_val: Union[str, List[str]]) -> List[str]:
    """Normalise path to a list of strings."""
    if isinstance(path_val, str):
        return [path_val]
    return list(path_val)


def _entity_file_path(entity: dict) -> Optional[Union[str, List[str]]]:
    """Return the path value handling both entity['file']['path'] and misindented entity['path']."""
    file_val = entity.get('file')
    if isinstance(file_val, dict) and 'path' in file_val:
        return file_val['path']
    if 'path' in entity:
        return entity['path']
    return None


def _referenced_config_inner_path(referenced_config_path: str) -> Optional[str]:
    """Return the resolved absolute path of the file referenced by a referenced config YAML's top-level 'path' key."""
    try:
        with open(referenced_config_path, 'r') as f:
            content = yaml.safe_load(f)
    except Exception:
        return None
    if not isinstance(content, dict):
        return None
    inner = content.get('path')
    if not inner or not isinstance(inner, str):
        return None
    ref_config_dir = os.path.dirname(os.path.abspath(referenced_config_path))
    if os.path.isabs(inner):
        return inner
    return os.path.normpath(os.path.join(ref_config_dir, inner))


def extract_file_paths_from_config(config: dict, config_dir: str) -> List[str]:
    """Extract all file paths that need to be staged for a BoltzGen config.

    For each entity, collects:
    - The path(s) listed under entities[].file.path (str or list)
    - For any .yaml/.yml listed there: also the inner top-level 'path' from that referenced config
    """
    file_paths: List[str] = []
    if 'entities' not in config:
        return file_paths
    config_dir = str(Path(config_dir).resolve())
    for entity in config['entities']:
        path_val = _entity_file_path(entity)
        if path_val is None:
            continue
        for file_path in _path_to_list(path_val):
            if not isinstance(file_path, str):
                continue
            if not os.path.isabs(file_path):
                resolved = os.path.normpath(os.path.join(config_dir, file_path))
            else:
                resolved = file_path
            file_paths.append(resolved)
            if file_path.lower().endswith(('.yaml', '.yml')) and os.path.exists(resolved):
                inner = _referenced_config_inner_path(resolved)
                if inner:
                    file_paths.append(inner)
    return file_paths


def extract_file_paths(
    config_path: str,
    config_dir: Optional[str] = None,
) -> Tuple[dict, List[str]]:
    """Load config and extract all file paths that need to be staged.

    Returns:
        (config dict, list of resolved absolute file paths)
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
    if config is None:
        config = {}
    file_paths = extract_file_paths_from_config(config, config_dir)
    return config, file_paths


def main() -> int:
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
        _, file_paths = extract_file_paths(args.config_yaml, args.config_dir)
        unique_paths = list(dict.fromkeys(file_paths))
        if args.format == 'list':
            for path in unique_paths:
                print(path)
        elif args.format == 'json':
            import json
            print(json.dumps(unique_paths, indent=2))
        elif args.format == 'tsv':
            print('\n'.join(unique_paths))
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
