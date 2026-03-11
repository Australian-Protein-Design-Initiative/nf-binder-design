#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Helper script for RFD3 (RFDiffusion3) config handling in Nextflow.

Modes:
  stage      - Rewrite input paths in an existing JSON config so they point
               to Nextflow-staged files (basename only).
  generate   - Create a new RFD3 JSON config from CLI flags, translating
               RFDiffusion v1-style contigs/hotspots to v3 format.
  parse-inputs - Print resolved input file path(s) from config (one per line).
"""

import argparse
import json
import logging
import os
import re
import sys
from pathlib import Path
from typing import Optional

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def load_config(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def save_config(config: dict, path: str) -> None:
    with open(path, "w") as f:
        json.dump(config, f, indent=2)
    log.info("Wrote config to %s", path)


def stage_config(config_path: str, output_path: str) -> None:
    """Rewrite 'input' paths in each spec to use basename only."""
    config = load_config(config_path)

    for key, spec in config.items():
        if isinstance(spec, dict) and "input" in spec:
            original = spec["input"]
            spec["input"] = os.path.basename(original)
            if original != spec["input"]:
                log.info("Rewrote input path: %s -> %s", original, spec["input"])

    save_config(config, output_path)


def convert_v1_contigs_to_v3(contigs: str) -> str:
    """Convert RFDiffusion v1 contig format to v3 format.

    v1: "[A18-132/0 65-120]"  (space-separated, square brackets, /0 attached)
    v3: "A18-132,/0,65-120"   (comma-separated, no brackets, /0 as separate element)
    """
    s = contigs.strip()
    # Strip outer brackets
    s = re.sub(r"^\[+", "", s)
    s = re.sub(r"\]+$", "", s)
    s = s.strip()

    # If already comma-separated (likely v3 format), return as-is
    if "," in s and " " not in s:
        return s

    # Split on spaces (v1 separator), then separate /0 chain breaks
    # from adjacent segments (v1 attaches /0 to preceding element)
    v3_parts = []
    for part in s.split():
        if "/0" in part:
            # Split "A18-132/0" into ["A18-132", "/0"]
            segments = part.split("/0")
            for i, seg in enumerate(segments):
                if seg:
                    v3_parts.append(seg)
                if i < len(segments) - 1:
                    v3_parts.append("/0")
        else:
            v3_parts.append(part)

    return ",".join(v3_parts)


def parse_v1_hotspots(hotspot_str: str) -> dict:
    """Convert v1 hotspot string to v3 select_hotspots dictionary.

    v1: "[A56,A88,A96]" or "A56,A88,A96"
    v3: {"A56": "ALL", "A88": "ALL", "A96": "ALL"}
    """
    s = hotspot_str.strip()
    s = re.sub(r"^\[+", "", s)
    s = re.sub(r"\]+$", "", s)
    s = s.strip()

    if not s:
        return {}

    hotspots = {}
    for residue in s.split(","):
        residue = residue.strip()
        if residue:
            hotspots[residue] = "ALL"

    return hotspots


def generate_config(
    design_name: str,
    input_pdb: str,
    contigs: str,
    hotspot_res: Optional[str] = None,
    partial_t: Optional[float] = None,
    is_non_loopy: Optional[bool] = None,
    translate_v1: bool = True,
) -> dict:
    """Generate an rfd3 JSON config from CLI parameters."""
    spec: dict = {
        "dialect": 2,
        "input": os.path.basename(input_pdb),
    }

    if contigs:
        contig_str = convert_v1_contigs_to_v3(contigs) if translate_v1 else contigs
        spec["contig"] = contig_str

    if hotspot_res:
        hotspots = parse_v1_hotspots(hotspot_res) if translate_v1 else hotspot_res
        if hotspots:
            spec["select_hotspots"] = hotspots
            spec["infer_ori_strategy"] = "hotspots"

    if is_non_loopy is not None:
        spec["is_non_loopy"] = is_non_loopy

    if partial_t is not None:
        spec["partial_t"] = partial_t

    config = {design_name: spec}
    return config


def extract_input_paths(config_path: str, config_dir: Optional[str] = None) -> list[str]:
    """Extract 'input' paths from each spec in an rfd3 JSON config.

    Paths are resolved relative to the config file directory (or config_dir if given).
    Returns unique paths in order of first occurrence.
    """
    path_obj = Path(config_path)
    if config_dir is None:
        config_dir = str(path_obj.parent.resolve())
    else:
        config_dir = str(Path(config_dir).resolve())

    if not path_obj.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    config = load_config(config_path)
    seen: set[str] = set()
    result: list[str] = []
    for key, spec in config.items():
        if isinstance(spec, dict) and "input" in spec:
            raw = spec["input"]
            resolved = Path(raw) if os.path.isabs(raw) else Path(config_dir) / raw
            resolved_str = str(resolved.resolve())
            if resolved_str not in seen:
                seen.add(resolved_str)
                result.append(resolved_str)
    return result


def cmd_parse_inputs(args: argparse.Namespace) -> int:
    paths = extract_input_paths(args.config, args.config_dir)
    for p in paths:
        print(p)
    return 0


def cmd_stage(args: argparse.Namespace) -> int:
    stage_config(args.config, args.output)
    return 0


def cmd_generate(args: argparse.Namespace) -> int:
    is_non_loopy: Optional[bool] = None
    if args.is_non_loopy and args.allow_loopy:
        raise ValueError("--is-non-loopy and --allow-loopy are mutually exclusive")
    if args.is_non_loopy:
        is_non_loopy = True
    elif args.allow_loopy:
        is_non_loopy = False
    config = generate_config(
        design_name=args.design_name,
        input_pdb=args.input_pdb,
        contigs=args.contigs,
        hotspot_res=args.hotspot_res,
        partial_t=args.partial_t,
        is_non_loopy=is_non_loopy,
        translate_v1=not args.no_translate,
    )
    save_config(config, args.output)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="RFDiffusion3 config helper for Nextflow"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- stage subcommand ---
    stage_parser = subparsers.add_parser(
        "stage",
        help="Rewrite input paths in an existing config to use basenames",
    )
    stage_parser.add_argument("config", help="Path to JSON config file")
    stage_parser.add_argument("-o", "--output", required=True, help="Output JSON path")

    # --- parse-inputs subcommand ---
    parse_parser = subparsers.add_parser(
        "parse-inputs",
        help="Print resolved input file path(s) from config (one per line)",
    )
    parse_parser.add_argument("config", help="Path to JSON config file")
    parse_parser.add_argument(
        "--config-dir",
        default=None,
        help="Directory for resolving relative paths (default: config file directory)",
    )

    # --- generate subcommand ---
    gen_parser = subparsers.add_parser(
        "generate",
        help="Generate a new rfd3 JSON config from CLI parameters",
    )
    gen_parser.add_argument("--design-name", default="design", help="Design spec key name")
    gen_parser.add_argument("--input-pdb", required=True, help="Input PDB/CIF file path")
    gen_parser.add_argument("--contigs", required=True, help="Contig string")
    gen_parser.add_argument("--hotspot-res", default=None, help="Hotspot residues")
    gen_parser.add_argument(
        "--partial-t", type=float, default=None,
        help="Partial diffusion noise level (angstroms)",
    )
    gen_parser.add_argument(
        "--is-non-loopy", action="store_true",
        help="Set is_non_loopy to true in the generated config",
    )
    gen_parser.add_argument(
        "--allow-loopy", action="store_true",
        help="Set is_non_loopy to false (allow loopy designs)",
    )
    gen_parser.add_argument(
        "--no-translate", action="store_true",
        help="Don't translate v1 contig/hotspot format to v3",
    )
    gen_parser.add_argument("-o", "--output", required=True, help="Output JSON path")

    args = parser.parse_args()

    if args.command == "stage":
        return cmd_stage(args)
    if args.command == "parse-inputs":
        return cmd_parse_inputs(args)
    if args.command == "generate":
        return cmd_generate(args)

    return 1


if __name__ == "__main__":
    sys.exit(main())
