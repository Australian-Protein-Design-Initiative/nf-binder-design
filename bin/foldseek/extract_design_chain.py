#!/usr/bin/env python3
"""
Extract the shorter (design) chain from protein structure files.

Supports PDB (.pdb), mmCIF (.cif), and gzipped variants (.pdb.gz, .cif.gz).
Uses only the Python standard library — no external dependencies.

Usage: extract_design_chain.py <output_dir> <input1> [input2 ...]

For each input structure, identifies the polymer chain with fewer atoms
(the design/binder chain) and writes a new PDB file containing only that
chain's atoms. Output is always PDB format for maximum compatibility.
"""

import gzip
import os
import sys


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

def extract_from_pdb(lines, source_name):
    """Parse PDB ATOM/HETATM lines, return design chain ID + filtered lines."""
    chain_atoms = {}  # chain_id -> count
    atom_lines = []   # (chain_id, line)

    for line in lines:
        rec = line[:6].strip()
        if rec not in ('ATOM', 'HETATM'):
            continue
        if len(line) < 22:
            continue
        chain_id = line[21]
        chain_atoms[chain_id] = chain_atoms.get(chain_id, 0) + 1
        atom_lines.append((chain_id, line))

    if not chain_atoms:
        return None, [], 0

    design_chain = min(chain_atoms, key=chain_atoms.get)
    filtered = [line for cid, line in atom_lines if cid == design_chain]
    return design_chain, filtered, len(chain_atoms)


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


def extract_from_cif(lines, source_name):
    """Parse mmCIF, return design chain ID + PDB-formatted atom lines."""
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
    atom_records = []  # (chain_id, fields)

    for i in range(data_start, data_end):
        parts = lines[i].split()
        if len(parts) <= max_col:
            continue
        chain_id = parts[asym_idx]
        chain_atoms[chain_id] = chain_atoms.get(chain_id, 0) + 1
        atom_records.append((chain_id, parts))

    if not chain_atoms:
        return None, [], 0

    design_chain = min(chain_atoms, key=chain_atoms.get)
    filtered = [(cid, p) for cid, p in atom_records if cid == design_chain]

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

    return design_chain, pdb_lines, len(chain_atoms)


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


def process_file(input_path, output_dir):
    """Process one structure file, extract design chain, write PDB."""
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
        chain, atom_lines, n_chains = extract_from_pdb(lines, base)
    else:
        chain, atom_lines, n_chains = extract_from_cif(lines, base)

    if chain is None:
        print(f"WARNING: No atoms found in {base}, skipping", file=sys.stderr)
        return None

    write_pdb(atom_lines, out_path)
    print(f"  {base}: extracted chain {chain} ({len(atom_lines)} atoms, {n_chains} chains total)")
    return out_path


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <output_dir> <input> [input2 ...]", file=sys.stderr)
        print(f"  Supports: .pdb, .cif, .pdb.gz, .cif.gz", file=sys.stderr)
        sys.exit(1)

    output_dir = sys.argv[1]
    os.makedirs(output_dir, exist_ok=True)

    inputs = sys.argv[2:]
    extracted = 0

    for path in sorted(inputs):
        result = process_file(path, output_dir)
        if result:
            extracted += 1

    print(f"\nExtracted {extracted}/{len(inputs)} design chains -> {output_dir}")


if __name__ == '__main__':
    main()
