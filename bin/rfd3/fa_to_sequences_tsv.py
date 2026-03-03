#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Build a TSV (filename, sequence, length, chain) from MPNN output CIF files.
Extracts the amino acid sequence for the specified chain (default: B, the binder in RFD3).
The filename column matches the CIF filename, which is the 'filename' key in combined_scores.tsv.
"""

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

AA3TO1: dict[str, str] = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O',
}


def read_cif_chain_sequence(path: Path, chain_id: str) -> str | None:
    """Extract amino acid sequence for a chain from an MPNN CIF file via atom_site records."""
    try:
        content = path.read_text()
    except OSError as e:
        log.warning("Cannot read %s: %s", path, e)
        return None

    loop_start = content.find('loop_\n_atom_site')
    if loop_start == -1:
        log.warning("No _atom_site loop in %s", path)
        return None

    atom_section = content[loop_start:]
    col_names: list[str] = []
    for line in atom_section.split('\n'):
        stripped = line.strip()
        if stripped.startswith('_atom_site.'):
            col_names.append(stripped[len('_atom_site.'):])
        elif col_names and stripped.startswith('ATOM'):
            break

    try:
        idx_chain = col_names.index('label_asym_id')
        idx_res = col_names.index('label_comp_id')
        idx_seq = col_names.index('label_seq_id')
    except ValueError as e:
        log.warning("Missing column in %s: %s", path, e)
        return None

    seen: set[tuple[str, str]] = set()
    residues: list[str] = []
    data_start = atom_section.find('\nATOM')
    for line in atom_section[data_start:].split('\n'):
        if not line.startswith('ATOM'):
            continue
        parts = line.split()
        if len(parts) <= max(idx_chain, idx_res, idx_seq):
            continue
        if parts[idx_chain] != chain_id:
            continue
        key = (parts[idx_chain], parts[idx_seq])
        if key in seen:
            continue
        seen.add(key)
        residues.append(AA3TO1.get(parts[idx_res], 'X'))

    if not residues:
        log.warning("Chain %s not found in %s", chain_id, path)
        return None
    return ''.join(residues)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build sequences.tsv from MPNN output CIF files")
    parser.add_argument("inputs", type=Path, nargs='+',
                        help="CIF files or a directory containing CIF files")
    parser.add_argument("-o", "--output", default="-",
                        help="Output TSV path (default: stdout)")
    parser.add_argument("--chain", default="B",
                        help="Chain ID to extract sequence for (default: B for RFD3 binder)")
    args = parser.parse_args()

    paths: list[Path] = []
    for p in args.inputs:
        if p.is_dir():
            paths.extend(sorted(p.glob("*.cif")))
        else:
            paths.append(p)

    rows: list[tuple[str, str, int]] = []
    for cif in paths:
        if cif.name in ("empty.cif", "empty"):
            continue
        seq = read_cif_chain_sequence(cif, args.chain)
        if seq is None:
            log.warning("No sequence extracted from %s (chain %s)", cif.name, args.chain)
            continue
        rows.append((cif.name, seq, len(seq)))

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        out.write("filename\tsequence\tlength\tchain\n")
        for filename, seq, length in rows:
            out.write(f"{filename}\t{seq}\t{length}\t{args.chain}\n")
    finally:
        if out is not sys.stdout:
            out.close()

    if args.output != "-":
        log.info("Wrote %d rows to %s", len(rows), args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
