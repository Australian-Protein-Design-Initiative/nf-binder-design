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
import yaml
from pathlib import Path
from typing import Dict, List, Optional, Union


def _path_to_list(path_val: Union[str, List[str]]) -> List[str]:
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


def _make_symlink(target_path: Path, source_abs: Path) -> None:
    """Create a relative symlink at target_path pointing to source_abs."""
    target_path.parent.mkdir(parents=True, exist_ok=True)
    rel_source = os.path.relpath(source_abs, target_path.parent)
    if target_path.exists() or target_path.is_symlink():
        target_path.unlink()
    target_path.symlink_to(rel_source)


def _stage_one(
    file_path: str,
    config_dir: str,
    input_files_map: Dict[str, str],
) -> None:
    """Stage a single file at its relative path under config_dir."""
    if not os.path.isabs(file_path):
        resolved = os.path.join(config_dir, file_path)
    else:
        resolved = file_path
    rel_path = os.path.relpath(resolved, config_dir)
    file_basename = os.path.basename(file_path)
    target_path = Path(config_dir) / rel_path

    if file_basename in input_files_map:
        source_path = Path(input_files_map[file_basename]).resolve()
        _make_symlink(target_path, source_path)
    elif Path(file_basename).exists():
        _make_symlink(target_path, Path(file_basename).resolve())
    else:
        if not target_path.exists():
            raise FileNotFoundError(
                f"Input file not found: {file_path} (basename: {file_basename})"
            )


def _stage_referenced_config_inner(
    staged_yaml_path: Path,
    config_dir: str,
    input_files_map: Dict[str, str],
) -> None:
    """After staging a referenced config YAML, also stage the PDB/CIF it references internally.

    Referenced configs have a top-level 'path: foo.cif' that is resolved relative to
    the YAML's own location. In the work dir the CIF must sit next to the YAML.
    """
    try:
        with open(staged_yaml_path, 'r') as f:
            content = yaml.safe_load(f)
    except Exception:
        return
    if not isinstance(content, dict):
        return
    inner = content.get('path')
    if not inner or not isinstance(inner, str):
        return
    inner_basename = os.path.basename(inner)
    target_path = staged_yaml_path.parent / (
        inner if not os.path.isabs(inner) else inner_basename
    )
    if target_path.exists() or target_path.is_symlink():
        return
    if inner_basename in input_files_map:
        source_path = Path(input_files_map[inner_basename]).resolve()
        _make_symlink(target_path, source_path)


def stage_input_files(config_path: str, input_files_dir: str, config_dir: Optional[str] = None):
    """Stage input files from config.yaml at correct relative paths using symlinks.

    Creates relative symlinks preserving the directory structure expected by BoltzGen.
    For referenced config YAML files, also stages the PDB/CIF they reference internally.

    Args:
        config_path: Path to the config YAML file
        input_files_dir: Directory containing staged input files (from Nextflow)
        config_dir: Directory containing the config file (for resolving relative paths).
                    If None, uses the directory of config_path.
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

    input_files_map: Dict[str, str] = {}
    if os.path.exists(input_files_dir):
        for root, _dirs, files in os.walk(input_files_dir):
            for fname in files:
                basename = os.path.basename(fname)
                if basename not in input_files_map:
                    input_files_map[basename] = os.path.join(root, fname)

    for entity in config['entities']:
        path_val = _entity_file_path(entity)
        if path_val is None:
            continue
        for file_path in _path_to_list(path_val):
            if not isinstance(file_path, str):
                continue
            _stage_one(file_path, config_dir, input_files_map)
            if file_path.lower().endswith(('.yaml', '.yml')):
                if not os.path.isabs(file_path):
                    staged_yaml = Path(config_dir) / os.path.relpath(
                        os.path.join(config_dir, file_path), config_dir
                    )
                else:
                    staged_yaml = Path(config_dir) / os.path.basename(file_path)
                _stage_referenced_config_inner(staged_yaml, config_dir, input_files_map)


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
        help='Directory containing config file. Default: directory of config_yaml'
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
