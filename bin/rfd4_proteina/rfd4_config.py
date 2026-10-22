#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Helper for RFD4-Proteina condition-spec handling in Nextflow.

Modes:
  parse-inputs - Print the resolved input structure path(s) from a spec, one per line.
  stage        - Rewrite input paths in a spec so they point at Nextflow-staged
                 files (basename only), writing JSON.
  generate     - Build a condition spec from --contigs / --hotspot-res, using the
                 upstream default recipe for protein binders.
  infer-chains - Print "target_chain,binder_chain" derived from a contig string.

YAML specs are read when PyYAML is importable (it is, inside the RFD4-Proteina
container); JSON is always readable. Everything this script writes is JSON, which
rfd4-generate accepts wherever it accepts YAML.
"""

import argparse
import json
import logging
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

# Reserved top-level key: shared defaults, not a design (see the condition-spec reference).
GLOBAL_ARGS = "global_args"

# An expandable segment is a length or length range, optionally with a P/R/D chain-type
# suffix: "20-30", "10R", "50-120". A named expandable segment carries a "name:" prefix.
EXPANDABLE_RE = re.compile(r"^\d+(-\d+)?[PRD]?$")
# A preserved segment is a chain id followed by a residue or residue range: "A17-131", "AA52".
PRESERVED_RE = re.compile(r"^(?P<chain>[A-Za-z][A-Za-z0-9_]*?)(?P<start>\d+)(-(?P<end>\d+))?$")


def load_spec(path: str) -> dict:
    text = Path(path).read_text()
    if path.lower().endswith((".yaml", ".yml")):
        try:
            import yaml
        except ImportError:
            raise SystemExit(
                f"Reading {path} needs PyYAML. Convert the spec to JSON, or run this in the "
                "RFD4-Proteina container."
            )
        return yaml.safe_load(text)
    return json.loads(text)


def save_spec(spec: dict, path: str) -> None:
    Path(path).write_text(json.dumps(spec, indent=2) + "\n")
    log.info("Wrote spec to %s", path)


def design_entries(spec: dict):
    """Yield (key, design_dict) for every design, skipping the reserved global_args key."""
    for key, value in spec.items():
        if key == GLOBAL_ARGS:
            continue
        if isinstance(value, dict):
            yield key, value


def _input_list(design: dict) -> List[str]:
    raw = design.get("input")
    if raw is None:
        return []
    return [str(p) for p in raw] if isinstance(raw, list) else [str(raw)]


def parse_inputs(config_path: str, config_dir: Optional[str]) -> List[str]:
    """Resolved input paths, in declaration order, deduplicated.

    `global_args` can carry a shared `input`, so it is read here even though it is not a design.
    """
    spec = load_spec(config_path)
    base = Path(config_dir) if config_dir else Path(config_path).resolve().parent
    seen: Dict[str, None] = {}
    sources = []
    if isinstance(spec.get(GLOBAL_ARGS), dict):
        sources.append(spec[GLOBAL_ARGS])
    sources.extend(design for _, design in design_entries(spec))
    for design in sources:
        for raw in _input_list(design):
            resolved = raw if os.path.isabs(raw) else str((base / raw).resolve())
            seen.setdefault(resolved, None)
    return list(seen)


def stage_spec(config_path: str, output_path: str) -> None:
    """Rewrite every input path to its basename, so a spec works against staged files."""
    spec = load_spec(config_path)
    targets = [spec[GLOBAL_ARGS]] if isinstance(spec.get(GLOBAL_ARGS), dict) else []
    targets.extend(design for _, design in design_entries(spec))
    for design in targets:
        raw = design.get("input")
        if raw is None:
            continue
        if isinstance(raw, list):
            design["input"] = [os.path.basename(str(p)) for p in raw]
        else:
            design["input"] = os.path.basename(str(raw))
    save_spec(spec, output_path)


def split_contig_runs(contig: str) -> List[List[str]]:
    """Split a contig into chain runs on /0 chainbreaks."""
    tokens = [t.strip() for t in contig.split(",")]
    runs: List[List[str]] = [[]]
    for token in tokens:
        if not token:
            continue
        if token == "/0":
            runs.append([])
        else:
            runs[-1].append(token)
    return [run for run in runs if run]


def classify_token(token: str) -> Tuple[str, Optional[str]]:
    """Return (kind, chain_id) where kind is 'preserved' or 'expandable'."""
    if ":" in token or EXPANDABLE_RE.match(token):
        return "expandable", None
    match = PRESERVED_RE.match(token)
    if not match:
        raise SystemExit(f"Could not parse contig token: {token!r}")
    return "preserved", match.group("chain")


def infer_chains(contig: str) -> Tuple[str, str]:
    """Target and binder chain ids for a two-chain binder contig.

    A run containing a preserved segment takes that segment's chain id; a run with none
    is the designed binder, and RFD4-Proteina gives it the first unused letter.
    """
    runs = split_contig_runs(contig)
    if not runs:
        raise SystemExit("Empty contig")

    run_chains = []
    for run in runs:
        chain = None
        for token in run:
            kind, token_chain = classify_token(token)
            if kind == "preserved":
                chain = token_chain
                break
        run_chains.append(chain)

    target_chains = [c for c in run_chains if c]
    if not target_chains:
        raise SystemExit(
            "Contig has no preserved segment, so there is no target chain. A contig made "
            "entirely of expandable segments is rejected by RFD4-Proteina."
        )

    used = set(target_chains)
    binder_chains = []
    for chain in run_chains:
        if chain:
            continue
        fresh = next(c for c in (chr(o) for o in range(ord("A"), ord("Z") + 1)) if c not in used)
        used.add(fresh)
        binder_chains.append(fresh)

    if not binder_chains:
        raise SystemExit(
            "Contig has no designed chain: every run contains a preserved segment. Add an "
            "expandable segment, e.g. \"50-120, /0, A17-131\"."
        )
    if len(binder_chains) > 1 or len(target_chains) > 1:
        log.warning(
            "Contig has %d designed and %d preserved chains; reporting the first of each",
            len(binder_chains),
            len(target_chains),
        )
    return target_chains[0], binder_chains[0]


def target_chain_ids(contig: str) -> List[str]:
    """Every preserved chain id in the contig, in declaration order, deduplicated."""
    chains: Dict[str, None] = {}
    for run in split_contig_runs(contig):
        for token in run:
            kind, chain = classify_token(token)
            if kind == "preserved" and chain:
                chains.setdefault(chain, None)
    return list(chains)


def parse_hotspots(hotspot_res: str) -> List[str]:
    """Turn "A56,A115,A123" into atomworks selections covering each whole residue."""
    selections = []
    for item in hotspot_res.split(","):
        item = item.strip()
        if not item:
            continue
        match = re.match(r"^(?P<chain>[A-Za-z][A-Za-z0-9_]*?)(?P<res>\d+)$", item)
        if not match:
            raise SystemExit(
                f"Could not parse hotspot residue {item!r}; expected a chain id followed by a "
                "residue number, e.g. A56"
            )
        selections.append(f"sel('{match.group('chain')}/*/{match.group('res')}')")
    return selections


def generate_spec(
    design_name: str,
    input_path: str,
    contigs: str,
    hotspot_res: Optional[str],
    output_path: str,
    partial_trajectory_fraction: Optional[float] = None,
) -> None:
    """Build a binder condition spec using the upstream default recipe.

    That recipe is: target backbone geometry as one free-pose coordinate group (C_CRDGRP),
    the target sequence held fixed (C_SEQ), and hotspots where given (C_HOT).
    """
    targets = target_chain_ids(contigs)
    if not targets:
        raise SystemExit("Contig has no preserved segment, so there is no target to condition on")
    target_sel = ",".join(targets) if len(targets) > 1 else targets[0]

    conditions: Dict[str, dict] = {
        "C_CRDGRP": {"1": [{"select": f"sel('{target_sel}/*/*/[N,CA,C,O]')"}]},
        "C_SEQ": {"true": [{"select": f"sel('{target_sel}')"}]},
    }
    if hotspot_res:
        selections = parse_hotspots(hotspot_res)
        if selections:
            conditions["C_HOT"] = {"true": [{"select": s} for s in selections]}

    design = {
        "input": os.path.basename(input_path),
        "contig": contigs,
        "conditions": conditions,
    }
    if partial_trajectory_fraction is not None:
        # Partial diffusion re-noises the input, so its existing conditioning must survive the load.
        design["input_as_annotated_cif"] = True

    save_spec({design_name: design}, output_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="RFD4-Proteina condition-spec helper", allow_abbrev=False
    )
    sub = parser.add_subparsers(dest="mode", required=True)

    p_parse = sub.add_parser("parse-inputs", help="Print resolved input path(s), one per line")
    p_parse.add_argument("config")
    p_parse.add_argument("--config-dir", default=None)

    p_stage = sub.add_parser("stage", help="Rewrite input paths to staged basenames")
    p_stage.add_argument("config")
    p_stage.add_argument("-o", "--output", required=True)

    p_gen = sub.add_parser("generate", help="Build a binder condition spec")
    p_gen.add_argument("--design-name", required=True)
    p_gen.add_argument("--input", required=True)
    p_gen.add_argument("--contigs", required=True)
    p_gen.add_argument("--hotspot-res", default=None)
    p_gen.add_argument("--partial-trajectory-fraction", type=float, default=None)
    p_gen.add_argument("-o", "--output", required=True)

    p_chains = sub.add_parser("infer-chains", help="Print 'target_chain,binder_chain'")
    p_chains.add_argument("--contigs", required=True)

    args = parser.parse_args()

    if args.mode == "parse-inputs":
        for path in parse_inputs(args.config, args.config_dir):
            print(path)
    elif args.mode == "stage":
        stage_spec(args.config, args.output)
    elif args.mode == "generate":
        generate_spec(
            args.design_name,
            args.input,
            args.contigs,
            args.hotspot_res,
            args.output,
            args.partial_trajectory_fraction,
        )
    elif args.mode == "infer-chains":
        target, binder = infer_chains(args.contigs)
        print(f"{target},{binder}")


if __name__ == "__main__":
    main()
