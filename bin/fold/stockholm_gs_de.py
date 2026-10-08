#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Pure-Python helper for recovering per-sequence descriptions from a Stockholm
alignment's `#=GS <name> DE <description>` lines.

alphafold.data.parsers.parse_stockholm keeps only the sequence name for each
row's "description" (see af2_msas_to_a3m.py) and silently drops the `DE` free-text
line, which is where jackhmmer's UniProt/UniRef hit headers (`TaxID=`, `RepID=`,
`Tax=`) actually live. This module has no dependency on the `alphafold` package so
it can be unit tested outside the AF2 container.
"""

from __future__ import annotations

import re
from typing import Dict, List

_GS_DE_RE = re.compile(r"^#=GS\s+(\S+)\s+DE\s+(.*?)\s*$")


def parse_gs_de(text: str) -> Dict[str, str]:
    """Map sequence name -> its `#=GS <name> DE <description>` free text.

    A name may have more than one `DE` line (rare); the first one wins, matching
    Stockholm's convention that GS lines are informational annotations, not part
    of the alignment itself.
    """
    gs_de: Dict[str, str] = {}
    for line in text.splitlines():
        m = _GS_DE_RE.match(line)
        if not m:
            continue
        name, desc = m.group(1), m.group(2)
        if name not in gs_de:
            gs_de[name] = desc
    return gs_de


def merge_descriptions(names: List[str], gs_de: Dict[str, str]) -> List[str]:
    """Return `<name> <DE text>` for names with a GS DE line, else `<name>` unchanged."""
    merged = []
    for name in names:
        de = gs_de.get(name)
        merged.append(f"{name} {de}" if de else name)
    return merged
