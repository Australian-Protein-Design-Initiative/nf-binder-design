#!/usr/bin/env python3
"""
Annotate FoldSeek results with CATH hierarchy descriptions.

Parses the FoldSeek target names (e.g. af_A0A096MJB5_3_128_1.25.40.10) to
extract CATH codes, then looks up descriptions at each level (Class,
Architecture, Topology, Homologous Superfamily) from the CATH names file.

Usage:
    annotate_cath.py --input foldseek_results.tsv --cath-names cath-names.txt
    annotate_cath.py -i foldseek_results.tsv -c cath-names.txt --gzip -o annotated.tsv.gz
"""

import argparse
import gzip
import re
import sys
from collections import OrderedDict

# Regex to extract CATH code from FoldSeek target names.
# Handles: af_{ACCESSION}_{START}_{END}_{CATHCODE} and PDB-like targets
CATH_CODE_RE = re.compile(r'(\d+\.\d+\.\d+\.\d+)$')


def parse_cath_names(path):
    """
    Parse CATH names file (cath-names.txt).

    Returns dict mapping CATH codes to their description strings.
    E.g. {'1': 'Mainly Alpha', '1.10': 'Orthogonal Bundle',
          '1.10.8': 'Helicase, Ruva Protein; domain 3',
          '1.10.8.10': 'DNA helicase RuvA subunit, C-terminal domain'}
    """
    names = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            # Format: NODE_NUM    REP_DOMAIN    :DESCRIPTION
            # Columns are whitespace-separated, description starts after ':'
            colon_idx = line.find(':')
            if colon_idx < 0:
                continue
            desc = line[colon_idx + 1:].strip()
            # The node number is the first whitespace-delimited field
            parts = line[:colon_idx].split()
            if not parts:
                continue
            cath_code = parts[0]
            names[cath_code] = desc
    return names


def extract_cath_code(target):
    """
    Extract CATH code from a FoldSeek target name.

    Examples:
        af_A0A096MJB5_3_128_1.25.40.10  -> 1.25.40.10
        af_A0A023GPK8_319_407_2.60.40.10 -> 2.60.40.10
        some_target_1.10.8.10             -> 1.10.8.10
        no_cath_code_here                 -> None
    """
    m = CATH_CODE_RE.search(target)
    return m.group(1) if m else None


def cath_code_to_levels(code):
    """
    Decompose CATH code into its four hierarchy levels.

    E.g. '1.25.40.10' -> {
        'C': '1',           Class
        'A': '1.25',        Architecture
        'T': '1.25.40',     Topology
        'H': '1.25.40.10',  Homologous Superfamily
    }
    """
    parts = code.split('.')
    if len(parts) != 4:
        return None
    return {
        'C': parts[0],
        'A': f'{parts[0]}.{parts[1]}',
        'T': f'{parts[0]}.{parts[1]}.{parts[2]}',
        'H': code,
    }


def annotate_row(fields, header, cath_names, out_columns):
    """
    Annotate a single TSV row with CATH descriptions.

    Returns list of annotation values (one per out_column).
    """
    target_idx = header.index('target')
    target = fields[target_idx]

    cath_code = extract_cath_code(target)
    if cath_code is None:
        return [''] * len(out_columns)

    levels = cath_code_to_levels(cath_code)
    if levels is None:
        return [''] * len(out_columns)

    result = []
    for col in out_columns:
        # col is like 'cath_class', 'cath_architecture', etc.
        level_key = col.split('_', 1)[1][0].upper()  # 'class' -> 'C'
        cath_level = levels.get(level_key, '')
        result.append(cath_names.get(cath_level, ''))

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Annotate FoldSeek results with CATH hierarchy descriptions.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output columns added:
  cath_class                    Class level (e.g. "Mainly Alpha")
  cath_architecture             Architecture level (e.g. "Orthogonal Bundle")
  cath_topology                 Topology level (e.g. "Rossmann fold")
  cath_homologous_superfamily   Homologous superfamily level (e.g. "UBA domain")
  cath_code                     Full CATH code (e.g. "1.25.40.10")

The CATH code is extracted from the FoldSeek target name, which for CATH50
has the format: af_{ACCESSION}_{START}_{END}_{CATHCODE}
""",
    )
    parser.add_argument('-i', '--input', required=True,
                        help='FoldSeek results TSV file')
    parser.add_argument('-c', '--cath-names', required=True,
                        help='CATH names file (cath-names.txt from CATH FTP)')
    parser.add_argument('-o', '--output', default=None,
                        help='Output TSV file (default: stdout)')
    parser.add_argument('--gzip', action='store_true',
                        help='Gzip the output')

    args = parser.parse_args()

    # Parse CATH names
    cath_names = parse_cath_names(args.cath_names)
    print(f"Loaded {len(cath_names)} CATH hierarchy entries", file=sys.stderr)

    # New columns to add
    out_columns = [
        'cath_class',
        'cath_architecture',
        'cath_topology',
        'cath_homologous_superfamily',
        'cath_code',
    ]

    # Open input
    opener_in = gzip.open if args.input.endswith('.gz') else open
    with opener_in(args.input, 'rt') as fin:
        # Read header
        header_line = fin.readline().rstrip('\n')
        header = header_line.split('\t')

        if 'target' not in header:
            print("ERROR: 'target' column not found in input TSV", file=sys.stderr)
            sys.exit(1)

        # Build output header
        new_header = header + out_columns

        # Open output
        if args.output:
            if args.gzip or args.output.endswith('.gz'):
                fout = gzip.open(args.output, 'wt')
            else:
                fout = open(args.output, 'w')
        else:
            fout = sys.stdout

        try:
            fout.write('\t'.join(new_header) + '\n')

            annotated = 0
            total = 0

            for line in fin:
                line = line.rstrip('\n')
                if not line:
                    continue
                fields = line.split('\t')
                total += 1

                # Extract CATH code and build annotation
                target = fields[header.index('target')]
                cath_code = extract_cath_code(target)

                if cath_code:
                    levels = cath_code_to_levels(cath_code)
                    if levels:
                        annotation = []
                        for col in out_columns:
                            if col == 'cath_code':
                                annotation.append(cath_code)
                            else:
                                level_key = col.split('_', 1)[1][0].upper()
                                cath_level = levels.get(level_key, '')
                                annotation.append(cath_names.get(cath_level, ''))
                        annotated += 1
                    else:
                        annotation = [''] * len(out_columns)
                else:
                    annotation = [''] * len(out_columns)

                fout.write('\t'.join(fields + annotation) + '\n')

        finally:
            if args.output and fout is not sys.stdout:
                fout.close()

    print(f"Annotated {annotated}/{total} rows with CATH descriptions", file=sys.stderr)


if __name__ == '__main__':
    main()
