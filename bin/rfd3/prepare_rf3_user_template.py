#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = []
# ///

"""Prepare a user-supplied RF3 template: select/rename chains to rfd3TargetChain and map template_selection."""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def run_gemmi(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        log.error("Command failed: %s\n%s", " ".join(cmd), r.stderr.strip())
        sys.exit(1)
    return r


def gemmi_list_chains(structure: Path) -> list[str]:
    r = run_gemmi(["gemmi", "residues", "-c", "-s", "-s", "-s", str(structure)])
    lines = r.stdout.splitlines()
    if len(lines) < 2:
        return []
    chains: set[str] = set()
    for line in lines[1:]:
        parts = line.split()
        if parts:
            chains.add(parts[0])
    return sorted(chains)


def structure_format(path: Path) -> str:
    n = path.name.lower()
    if n.endswith(".cif.gz") or n.endswith(".cif"):
        return "cif"
    return "pdb"


def needs_gemmi_normalise(structure: Path) -> bool:
    """Copying would duplicate .gz bytes into a plain .cif/.pdb path; use gemmi instead."""
    return structure.name.lower().endswith(".gz")


def chain_refs_from_token(token: str) -> list[str]:
    """AtomSelection token: chain is first '/'-separated segment, or whole token if no '/'."""
    t = token.strip()
    if not t:
        return []
    if "/" in t:
        return [t.split("/", 1)[0].strip()]
    return [t]


def unique_chain_refs_from_selection(sel: str) -> list[str]:
    seen: list[str] = []
    for part in sel.split(","):
        for ch in chain_refs_from_token(part):
            if ch and ch not in seen:
                seen.append(ch)
    return seen


def map_selection_token(token: str, src: str, dst: str) -> str:
    t = token.strip()
    if not t:
        return t
    if "/" in t:
        chain, rest = t.split("/", 1)
        if chain == src:
            return f"{dst}/{rest}"
        return t
    if t == src:
        return dst
    return t


def map_template_selection(user_sel: Optional[str], src: str, dst: str) -> str:
    if not user_sel or not user_sel.strip():
        return dst
    parts = [p.strip() for p in user_sel.split(",") if p.strip()]
    return ",".join(map_selection_token(p, src, dst) for p in parts)


def resolve_source_chain(
    chains: list[str],
    user_sel: Optional[str],
    target_chain: str,
) -> str:
    if len(chains) == 0:
        log.error("No polymer chains found in template")
        sys.exit(1)
    if len(chains) == 1:
        sole = chains[0]
        if user_sel and user_sel.strip():
            refs = unique_chain_refs_from_selection(user_sel)
            for r in refs:
                if r != sole:
                    log.error(
                        "template_selection references chain %r but the template file has only chain %r",
                        r,
                        sole,
                    )
                    sys.exit(1)
        return sole
    if user_sel and user_sel.strip():
        refs = unique_chain_refs_from_selection(user_sel)
        if len(refs) == 0:
            log.error("Could not parse chain IDs from --template-selection")
            sys.exit(1)
        if len(refs) > 1:
            log.error(
                "Multi-chain template: template_selection must reference exactly one chain "
                "(found %s). Use a single-chain template or one chain in the selection.",
                ",".join(refs),
            )
            sys.exit(1)
        src = refs[0]
        if src not in chains:
            log.error("Chain %r from template_selection not in template chains %s", src, chains)
            sys.exit(1)
        return src
    if target_chain in chains:
        return target_chain
    log.error(
        "Multi-chain template %s: set --rf3_template_selection to the target chain ID in this file "
        "(rfd3 target chain is %r, which is not present).",
        chains,
        target_chain,
    )
    sys.exit(1)


def gemmi_convert(src: Path, dst: Path, extra: Optional[list[str]] = None) -> None:
    cmd: list[str] = ["gemmi", "convert"]
    if extra:
        cmd.extend(extra)
    cmd.extend([str(src), str(dst)])
    run_gemmi(cmd)


def prepare_template(
    structure: Path,
    out_path: Path,
    target_chain: str,
    user_selection: Optional[str],
    selection_out: Path,
) -> None:
    if not structure.exists():
        log.error("Structure not found: %s", structure)
        sys.exit(1)

    out_fmt = structure_format(out_path)
    chains = gemmi_list_chains(structure)
    src = resolve_source_chain(chains, user_selection, target_chain)
    mapped = map_template_selection(user_selection, src, target_chain)
    selection_out.write_text(mapped, encoding="utf-8")
    log.info(
        "RF3 user template: source chain %r -> target %r; template_selection (for RF3 JSON): %r",
        src,
        target_chain,
        mapped,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_pdb = Path("work_rf3_template.pdb")

    if len(chains) == 1:
        if src == target_chain:
            if (
                structure_format(structure) == out_fmt
                and not needs_gemmi_normalise(structure)
            ):
                shutil.copy2(structure, out_path)
            else:
                gemmi_convert(structure, out_path)
        else:
            gemmi_convert(
                structure,
                out_path,
                [f"--rename-chain={src}:{target_chain}"],
            )
        return

    gemmi_convert(structure, tmp_pdb, ["--select", f"//{src}/"])
    if src != target_chain:
        gemmi_convert(
            tmp_pdb,
            out_path,
            [f"--rename-chain={src}:{target_chain}"],
        )
        tmp_pdb.unlink(missing_ok=True)
    else:
        if out_fmt == "cif":
            gemmi_convert(tmp_pdb, out_path)
            tmp_pdb.unlink(missing_ok=True)
        else:
            shutil.move(str(tmp_pdb), str(out_path))


def main() -> int:
    p = argparse.ArgumentParser(
        description="Prepare user RF3 template: align chain IDs to rfd3TargetChain and map template_selection.",
    )
    p.add_argument("--structure", required=True, type=Path, help="Input PDB/mmCIF (template file)")
    p.add_argument("--target-chain", required=True, help="RFD3/RF3 target chain ID after preparation")
    p.add_argument(
        "--template-selection",
        default=None,
        help="User template_selection using chain IDs from the template file (comma-separated)",
    )
    p.add_argument(
        "--template-selection-file",
        type=Path,
        default=None,
        help="File with user template_selection (empty = whole source chain, auto-mapped)",
    )
    p.add_argument("-o", "--output", required=True, type=Path, help="Output template path (e.g. template_rf3.pdb)")
    p.add_argument(
        "--write-mapped-selection",
        required=True,
        type=Path,
        help="Write one-line RF3 template_selection (post-rename chain IDs) for make_rf3_input_spec.py",
    )
    args = p.parse_args()

    user_sel: Optional[str] = args.template_selection
    if args.template_selection_file is not None:
        raw = args.template_selection_file.read_text()
        user_sel = raw if raw else None

    prepare_template(args.structure, args.output, args.target_chain, user_sel, args.write_mapped_selection)
    return 0


if __name__ == "__main__":
    sys.exit(main())
