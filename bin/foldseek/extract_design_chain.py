#!/usr/bin/env python3
"""
Extract the design (binder) chain from protein structure files.

Supports PDB (.pdb), mmCIF (.cif), and gzipped variants (.pdb.gz, .cif.gz).
Uses only the Python standard library - no external dependencies.

Usage: extract_design_chain.py <output_dir> <input1> [input2 ...]
       extract_design_chain.py <output_dir> --last-chains 2 <input1> [...]

By default the design chain is taken to be the polymer chain with the fewest
atoms, which holds for a one-target/one-binder complex. --last-chains N instead
keeps the last N chains in file order, for formats that write the binder last.

Output is always PDB format for maximum compatibility.
"""

import argparse
import gzip
import os
import sys


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

def select_design_chains(chain_order, chain_atoms, last_chains=0):
    """Which chain IDs make up the binder.

    With last_chains == 0, the smallest chain by atom count, which is right for
    a complex of one target chain plus one binder chain. Otherwise the last N
    chains in file order, for writers that put every target chain before the
    binder.
    """
    if last_chains > 0:
        if last_chains > len(chain_order):
            return list(chain_order)
        return list(chain_order[-last_chains:])
    return [min(chain_atoms, key=chain_atoms.get)]


def extract_from_pdb(lines, source_name, last_chains=0):
    """Parse PDB ATOM/HETATM lines, return design chain ID(s) + filtered lines."""
    chain_atoms = {}  # chain_id -> count
    chain_order = []  # chain ids, in order of first appearance
    atom_lines = []   # (chain_id, line)

    for line in lines:
        rec = line[:6].strip()
        if rec not in ('ATOM', 'HETATM'):
            continue
        if len(line) < 22:
            continue
        chain_id = line[21]
        if chain_id not in chain_atoms:
            chain_order.append(chain_id)
        chain_atoms[chain_id] = chain_atoms.get(chain_id, 0) + 1
        atom_lines.append((chain_id, line))

    if not chain_atoms:
        return None, [], 0

    design_chains = select_design_chains(chain_order, chain_atoms, last_chains)
    keep = set(design_chains)
    filtered = [line for cid, line in atom_lines if cid in keep]
    return ','.join(design_chains), filtered, len(chain_atoms)


def write_pdb(atom_lines, out_path):
    """Write ATOM lines as a PDB file."""
    with open(out_path, 'w') as f:
        f.write(f"HEADER    extracted design chain\n")
        for line in atom_lines:
            f.write(line)
            if not line.endswith('\n'):
                f.write('\n')
        f.write("END\n")


# ---------------------------------------------------------------------------
# mmCIF parsing (minimal, stdlib only)
# ---------------------------------------------------------------------------

def find_atom_site_loop(lines):
    """Find _atom_site loop: returns (header_start, col_map, data_start, data_end)."""
    i = 0
    while i < len(lines):
        if lines[i].strip() != 'loop_':
            i += 1
            continue
        j = i + 1
        cols = {}
        while j < len(lines) and lines[j].strip().startswith('_atom_site.'):
            cols[lines[j].strip()] = len(cols)
            j += 1
        if '_atom_site.label_asym_id' in cols:
            data_end = j
            for k in range(j, len(lines)):
                s = lines[k].strip()
                if not s or s.startswith('_') or s.startswith('loop_') or s.startswith('data_') or s.startswith('#'):
                    data_end = k
                    break
            return i, cols, j, data_end
        i += 1
    return None


def extract_from_cif(lines, source_name, last_chains=0):
    """Parse mmCIF, return design chain ID(s) + PDB-formatted atom lines."""
    result = find_atom_site_loop(lines)
    if result is None:
        return None, [], 0

    _, cols, data_start, data_end = result
    asym_idx = cols['_atom_site.label_asym_id']
    group_idx = cols.get('_atom_site.group_PDB')
    name_idx = cols.get('_atom_site.label_atom_id')
    altloc_idx = cols.get('_atom_site.label_alt_id')
    resname_idx = cols.get('_atom_site.label_comp_id')
    seqid_idx = cols.get('_atom_site.label_seq_id')
    x_idx = cols.get('_atom_site.Cartn_x')
    y_idx = cols.get('_atom_site.Cartn_y')
    z_idx = cols.get('_atom_site.Cartn_z')
    occ_idx = cols.get('_atom_site.occupancy')
    bfactor_idx = cols.get('_atom_site.B_iso_or_equiv')
    elem_idx = cols.get('_atom_site.type_symbol')
    max_col = max(cols.values())

    chain_atoms = {}
    chain_order = []   # chain ids, in order of first appearance
    atom_records = []  # (chain_id, fields)

    for i in range(data_start, data_end):
        parts = lines[i].split()
        if len(parts) <= max_col:
            continue
        chain_id = parts[asym_idx]
        if chain_id not in chain_atoms:
            chain_order.append(chain_id)
        chain_atoms[chain_id] = chain_atoms.get(chain_id, 0) + 1
        atom_records.append((chain_id, parts))

    if not chain_atoms:
        return None, [], 0

    design_chains = select_design_chains(chain_order, chain_atoms, last_chains)
    keep = set(design_chains)
    filtered = [(cid, p) for cid, p in atom_records if cid in keep]

    # Convert to PDB ATOM lines
    pdb_lines = []
    for idx, (cid, p) in enumerate(filtered, 1):
        rec = p[group_idx] if group_idx is not None else 'ATOM'
        if rec not in ('ATOM', 'HETATM'):
            rec = 'ATOM'
        atom_name = p[name_idx] if name_idx is not None else 'CA'
        altloc = p[altloc_idx] if altloc_idx is not None else ''
        resname = p[resname_idx] if resname_idx is not None else 'UNK'
        seqid = p[seqid_idx] if seqid_idx is not None else '.'
        x = p[x_idx] if x_idx is not None else '0.0'
        y = p[y_idx] if y_idx is not None else '0.0'
        z = p[z_idx] if z_idx is not None else '0.0'
        occ = p[occ_idx] if occ_idx is not None else '1.0'
        bfac = p[bfactor_idx] if bfactor_idx is not None else '0.0'
        elem = p[elem_idx] if elem_idx is not None else ''

        # PDB ATOM format (fixed width)
        # Atom name alignment: 1-char names -> " X  ", 2+ char -> " XX "
        if len(atom_name) < 4:
            atom_name_formatted = f" {atom_name:<3s}"
        else:
            atom_name_formatted = atom_name[:4]

        try:
            seqid_int = int(seqid) if seqid != '.' else 0
        except ValueError:
            seqid_int = 0

        pdb_line = (
            f"{rec:<6s}{idx:>5d} {atom_name_formatted}{altloc:>1s}"
            f"{resname:>3s} {cid:>1s}{seqid_int:>4d}    "
            f"{float(x):>8.3f}{float(y):>8.3f}{float(z):>8.3f}"
            f"{float(occ):>6.2f}{float(bfac):>6.2f}          "
            f"{elem:>2s}\n"
        )
        pdb_lines.append(pdb_line)

    return ','.join(design_chains), pdb_lines, len(chain_atoms)


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------

def read_lines(path):
    """Read lines from a file, handling .gz transparently."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        return f.read().split('\n')


def detect_format(path):
    """Detect structure format from extension (ignoring .gz)."""
    name = path[:-3] if path.endswith('.gz') else path
    if name.endswith('.pdb'):
        return 'pdb'
    elif name.endswith('.cif'):
        return 'cif'
    return None


def process_file(input_path, output_dir, last_chains=0):
    """Process one structure file, extract the design chain(s), write PDB."""
    fmt = detect_format(input_path)
    if fmt is None:
        print(f"WARNING: Unknown format for {input_path}, skipping", file=sys.stderr)
        return None

    lines = read_lines(input_path)
    base = os.path.basename(input_path)
    # Strip all extensions to get clean name
    stem = base
    for ext in ('.gz', '.cif', '.pdb'):
        if stem.endswith(ext):
            stem = stem[:-len(ext)]

    out_path = os.path.join(output_dir, f"{stem}.pdb")

    if fmt == 'pdb':
        chain, atom_lines, n_chains = extract_from_pdb(lines, base, last_chains)
    else:
        chain, atom_lines, n_chains = extract_from_cif(lines, base, last_chains)

    if chain is None:
        print(f"WARNING: No atoms found in {base}, skipping", file=sys.stderr)
        return None

    write_pdb(atom_lines, out_path)
    print(f"  {base}: extracted chain(s) {chain} ({len(atom_lines)} atoms, {n_chains} chains total)")
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Extract the design (binder) chain from structure files (.pdb, .cif, .gz)"
    )
    parser.add_argument("output_dir", help="Directory to write extracted PDB files into")
    parser.add_argument("inputs", nargs="+", help="Input structure files")
    parser.add_argument(
        "--last-chains",
        type=int,
        default=0,
        metavar="N",
        help="Keep the last N chains in file order instead of the single smallest chain. "
        "Use for writers that place every target chain before the binder (default: 0, "
        "meaning use the smallest-chain heuristic)",
    )
    args = parser.parse_args()

    if args.last_chains < 0:
        parser.error("--last-chains must not be negative")

    os.makedirs(args.output_dir, exist_ok=True)

    extracted = 0
    for path in sorted(args.inputs):
        if process_file(path, args.output_dir, args.last_chains):
            extracted += 1

    print(f"\nExtracted {extracted}/{len(args.inputs)} design chains -> {args.output_dir}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
