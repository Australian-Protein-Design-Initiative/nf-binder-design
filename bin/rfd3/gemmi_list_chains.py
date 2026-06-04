#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///

"""List unique polymer chain IDs from a structure file via gemmi residues."""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def gemmi_list_chains(structure: Path) -> list[str]:
    """Return sorted unique chain IDs in *structure*.

    Uses ``gemmi residues -c`` with three ``-s`` flags: gemmi allows -s 2x or 3x for
    progressively shorter output; 3x omits residue/atom detail so we only parse chain IDs.
    """
    cmd = ["gemmi", "residues", "-c", "-s", "-s", "-s", str(structure)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        log.error("Command failed: %s\n%s", " ".join(cmd), r.stderr.strip())
        sys.exit(1)
    lines = r.stdout.splitlines()
    if len(lines) < 2:
        return []
    chains: set[str] = set()
    for line in lines[1:]:
        parts = line.split()
        if parts:
            chains.add(parts[0])
    return sorted(chains)


def main() -> int:
    p = argparse.ArgumentParser(description="Print unique chain IDs from a structure file (one per line).")
    p.add_argument("structure", type=Path, help="PDB or mmCIF input")
    args = p.parse_args()
    for ch in gemmi_list_chains(args.structure):
        print(ch)
    return 0


if __name__ == "__main__":
    sys.exit(main())
