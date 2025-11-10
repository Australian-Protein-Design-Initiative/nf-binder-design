#!/usr/bin/env python
"""
Rename boltzgen output files to add start_index offset to indices in filenames.

Covers cases:
- {design}_{i}.{ext} -> {design}_{i+start}.{ext}
- {design}_{i}_{j}.{ext} -> {design}_{i+start}_{j}.{ext}
- When num_designs == 1:
  - {design}.{ext} -> {design}_{start}.{ext}
  - {design}_{j}.{ext} (in inverse-folded dirs) -> {design}_{start}_{j}.{ext}
- metrics/data variants:
  - data_{design}.npz -> data_{design}_{start}.npz
  - data_{design}_{i|j}.npz and data_{design}_{i}_{j}.npz handled analogously
  - metrics_{design}.npz / metrics_{design}_{...}.npz handled analogously
"""
import argparse
import logging
import re
from pathlib import Path
from typing import Optional, Tuple

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Allowed extensions for renaming
ALLOWED_EXTS = {'.cif', '.npz', '.pdb'}


def detect_batch_root(directory: Path) -> Path:
    """
    Given a directory inside a batch, return the batch root
    (the directory that contains 'intermediate_designs' or
    'intermediate_designs_inverse_folded').
    """
    current = directory.resolve()
    for _ in range(4):
        if (
            (current / 'intermediate_designs').exists()
            or (current / 'intermediate_designs_inverse_folded').exists()
        ):
            return current
        if current.parent == current:
            break
        current = current.parent
    return directory.resolve()


def detect_num_designs(batch_root: Path, design_name: str) -> Optional[int]:
    """
    Infer number of designs from the intermediate_designs directory.
    Returns None if cannot infer.
    """
    inter_dir = batch_root / 'intermediate_designs'
    if not inter_dir.exists():
        return None

    # If a bare file exists -> single design
    if (
        (inter_dir / f'{design_name}.cif').exists()
        or (inter_dir / f'{design_name}.npz').exists()
    ):
        return 1

    # Count unique indices among files like design_i.(cif|npz)
    idx_set = set()
    pat = re.compile(rf'^{re.escape(design_name)}_(\d+)\.(?:cif|npz)$')
    for p in inter_dir.iterdir():
        if not p.is_file():
            continue
        m = pat.match(p.name)
        if m:
            idx_set.add(int(m.group(1)))
    if idx_set:
        return len(idx_set)
    return None


def determine_context_is_inverse(directory: Path) -> bool:
    """
    Determine if the directory belongs to inverse-folded outputs.
    """
    parts = {part for part in directory.resolve().parts}
    return 'intermediate_designs_inverse_folded' in parts


def parse_name(
    filename: str,
    design_name: str,
) -> Optional[Tuple[str, Optional[int], Optional[int], str]]:
    """
    Parse filename into (prefix, i, j, ext) where:
    - prefix: '', 'data_', or 'metrics_'
    - i: design index (if present, else None)
    - j: inverse-folding index (if present, else None)
    - ext: including leading dot
    Returns None if the pattern doesn't correspond to files we should rename.
    """
    # Unified pattern: optional (data_|metrics_), design_name, optional _i,
    # optional _j, dot, ext
    m = re.match(
        rf'^(?P<prefix>(?:data_|metrics_)?)'
        rf'{re.escape(design_name)}'
        rf'(?:_(?P<n>\d+))?'
        rf'(?:_(?P<m>\d+))?'
        rf'\.(?P<ext>[^.]+)$',
        filename,
    )
    if not m:
        return None
    ext = '.' + m.group('ext')
    if ext not in ALLOWED_EXTS:
        return None
    prefix = m.group('prefix') or ''
    n_str = m.group('n')
    m_str = m.group('m')
    i = int(n_str) if n_str is not None else None
    j = int(m_str) if m_str is not None else None
    return prefix, i, j, ext


def format_name(
    prefix: str,
    design_name: str,
    i: Optional[int],
    j: Optional[int],
    ext: str,
) -> str:
    """
    Construct filename from components.
    """
    base = f'{prefix}{design_name}'
    if i is not None and j is not None:
        return f'{base}_{i}_{j}{ext}'
    if i is not None:
        return f'{base}_{i}{ext}'
    return f'{base}{ext}'


def rename_file_with_offset(
    file_path: Path,
    directory_root: Path,
    design_name: str,
    start_index: int,
    num_designs_hint: Optional[int],
    invfold_hint: Optional[int],
) -> bool:
    """
    Rename a file by applying start_index to design index position,
    preserving any inv-folding index.
    """
    parsed = parse_name(file_path.name, design_name)
    if parsed is None:
        return False

    prefix, i, j, ext = parsed

    # Determine context
    is_inverse_ctx = determine_context_is_inverse(directory_root)

    # Decide how to treat single-index names
    # If two indices present -> i is design, j is invfold
    # If one index present:
    #   - In design context -> it's design index
    #   - In inverse context:
    #       - If num_designs == 1 -> it's invfold index
    #       - Else -> it's design index
    # If none present -> add design index (start_index)
    if i is None and j is None:
        # Bare name; add design index only
        new_i = start_index
        new_j = None
    elif i is not None and j is None:
        if not is_inverse_ctx:
            # Design context: single index is design
            new_i = i + start_index
            new_j = None
        else:
            # Inverse context: disambiguate with hint/detection
            # If inverse_fold_num_sequences is set, single index is invfold index
            if invfold_hint is not None and invfold_hint > 1:
                # Inverse folding multiplicity > 1: single index is invfold index
                new_i = start_index
                new_j = i
            else:
                # Check if single design run
                batch_root = detect_batch_root(directory_root)
                inferred_num_designs = (
                    num_designs_hint
                    or detect_num_designs(batch_root, design_name)
                    or None
                )
                if inferred_num_designs == 1:
                    # Single-design run: this index is invfold index
                    new_i = start_index
                    new_j = i
                else:
                    # Multi-design run (or unknown): treat as design index
                    new_i = i + start_index
                    new_j = None
    else:
        # Both i and j present
        new_i = (i + start_index) if i is not None else start_index
        new_j = j

    new_name = format_name(prefix, design_name, new_i, new_j, ext)
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
    recursive: bool = True,
    num_designs: Optional[int] = None,
    inverse_fold_num_sequences: Optional[int] = None,
) -> int:
    """
    Rename all matching files in directory.
    """
    renamed_count = 0

    # Snapshot file list first to avoid iterator issues while renaming
    if recursive:
        files = [p for p in directory.rglob('*') if p.is_file()]
    else:
        files = [p for p in directory.iterdir() if p.is_file()]

    for file_path in files:
        if rename_file_with_offset(
            file_path=file_path,
            directory_root=directory,
            design_name=design_name,
            start_index=start_index,
            num_designs_hint=num_designs,
            invfold_hint=inverse_fold_num_sequences,
        ):
            renamed_count += 1

    return renamed_count


def main():
    parser = argparse.ArgumentParser(
        description='Rename boltzgen output files to add start_index offset'
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
        help='Hint: number of designs (optional)',
    )
    parser.add_argument(
        '--inverse_fold_num_sequences',
        type=int,
        default=None,
        help='Hint: inverse folding multiplicity (optional)',
    )
    parser.add_argument(
        '--no-recursive',
        action='store_true',
        help='Only rename files in directory, not subdirectories',
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
        recursive=not args.no_recursive,
        num_designs=args.num_designs,
        inverse_fold_num_sequences=args.inverse_fold_num_sequences,
    )

    logger.info(f"Renamed {renamed_count} file(s)")
    return 0


if __name__ == '__main__':
    exit(main())
