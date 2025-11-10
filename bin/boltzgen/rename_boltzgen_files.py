#!/usr/bin/env python
"""
Rename boltzgen design output files to add start_index offset to indices.

Handles two cases:
- {design}.{ext} -> {design}_{start_index}.{ext} (no index)
- {design}_{i}.{ext} -> {design}_{i+start_index}.{ext} (has index)

Only processes .cif and .npz files in intermediate_designs directory.
"""
import argparse
import logging
import re
from pathlib import Path
from typing import Optional

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Allowed extensions for renaming
ALLOWED_EXTS = {'.cif', '.npz'}


def parse_name(
    filename: str,
    design_name: str,
) -> Optional[tuple[Optional[int], str]]:
    """
    Parse filename into (index, ext) where:
    - index: design index if present (None if bare name)
    - ext: file extension including leading dot
    Returns None if the pattern doesn't match.
    """
    # Pattern: design_name, optional _i, dot, ext
    m = re.match(
        rf'^{re.escape(design_name)}'
        rf'(?:_(?P<idx>\d+))?'
        rf'\.(?P<ext>[^.]+)$',
        filename,
    )
    if not m:
        return None
    ext = '.' + m.group('ext')
    if ext not in ALLOWED_EXTS:
        return None
    idx_str = m.group('idx')
    idx = int(idx_str) if idx_str is not None else None
    return idx, ext


def format_name(
    design_name: str,
    index: Optional[int],
    ext: str,
) -> str:
    """Construct filename from components."""
    if index is not None:
        return f'{design_name}_{index}{ext}'
    return f'{design_name}{ext}'


def rename_file_with_offset(
    file_path: Path,
    design_name: str,
    start_index: int,
) -> bool:
    """
    Rename a file by applying start_index offset to design index.
    Returns True if file was renamed, False otherwise.
    """
    parsed = parse_name(file_path.name, design_name)
    if parsed is None:
        return False

    idx, ext = parsed

    if idx is None:
        # Bare name: add start_index
        new_idx = start_index
    else:
        # Has index: add start_index to it
        new_idx = idx + start_index

    new_name = format_name(design_name, new_idx, ext)
    if new_name == file_path.name:
        return False

    new_path = file_path.parent / new_name
    file_path.rename(new_path)
    logger.info(f"Renamed: {file_path.name} -> {new_name}")
    return True


def rename_files_in_directory(
    directory: Path,
    design_name: str,
    start_index: int,
) -> int:
    """Rename all matching files in directory."""
    renamed_count = 0

    # Only process files directly in the directory (not recursive)
    for file_path in directory.iterdir():
        if not file_path.is_file():
            continue
        if rename_file_with_offset(
            file_path=file_path,
            design_name=design_name,
            start_index=start_index,
        ):
            renamed_count += 1

    return renamed_count


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Rename boltzgen design output files to add start_index offset'
        )
    )
    parser.add_argument(
        'directory',
        type=Path,
        help='Directory containing files to rename',
    )
    parser.add_argument('design_name', help='Design name (e.g., "config")')
    parser.add_argument(
        'start_index',
        type=int,
        help='Start index offset to add',
    )
    parser.add_argument(
        '--num_designs',
        type=int,
        default=None,
        help='Number of designs (for validation, optional)',
    )

    args = parser.parse_args()

    if not args.directory.exists():
        logger.error(f"Directory does not exist: {args.directory}")
        return 1

    if not args.directory.is_dir():
        logger.error(f"Path is not a directory: {args.directory}")
        return 1

    renamed_count = rename_files_in_directory(
        directory=args.directory,
        design_name=args.design_name,
        start_index=args.start_index,
    )

    logger.info(f"Renamed {renamed_count} file(s)")
    return 0


if __name__ == '__main__':
    exit(main())
