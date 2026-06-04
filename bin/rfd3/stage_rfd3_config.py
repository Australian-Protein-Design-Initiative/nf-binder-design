#!/usr/bin/env python
# /// script
# requires-python = ">=3.7"
# ///

"""
Helper script for RFD3 (RFDiffusion3) config handling in Nextflow.

Modes:
  stage      - Rewrite input paths in an existing JSON config so they point
               to Nextflow-staged files (basename only). Optionally subsample
               select_hotspots per invocation.
  generate   - Create a new RFD3 JSON config from CLI flags, translating
               RFDiffusion v1-style contigs/hotspots to v3 format.
  parse-inputs - Print resolved input file path(s) from config (one per line).
  infer-binder-chain - Print MPNN binder chain letter (A, B, ...) from contig
               (v3 polymers split on /0; exactly one length-only binder polymer).
  infer-rfd3-chain-pair - Print target,binder chain letters for RFD3/MPNN/RF3
               (two-polymer contigs; optional --binder when MPNN chains are explicit).
"""

import argparse
import json
import logging
import math
import os
import random
import re
import sys
from pathlib import Path
from typing import List, Optional, Set, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def load_config(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def save_config(config: dict, path: str) -> None:
    with open(path, "w") as f:
        json.dump(config, f, indent=2)
    log.info("Wrote config to %s", path)


def subsample_select_hotspots(
    select_hotspots: dict, hotspot_subsample: float, spec_key: str
) -> None:
    """In-place: replace select_hotspots with a random subset of keys (mutates spec)."""
    if not isinstance(select_hotspots, dict) or not select_hotspots:
        return
    keys = [k for k in select_hotspots if k]
    if not keys:
        return
    n_total = len(keys)
    n_keep = max(1, math.ceil(n_total * hotspot_subsample))
    if n_keep >= n_total:
        return
    random.seed()
    chosen = random.sample(keys, n_keep)
    new_map = {k: select_hotspots[k] for k in chosen}
    select_hotspots.clear()
    select_hotspots.update(new_map)
    log.info(
        "Subsampled select_hotspots for %s: kept %d/%d keys: %s",
        spec_key,
        n_keep,
        n_total,
        ",".join(sorted(chosen)),
    )


def stage_config(
    config_path: str, output_path: str, hotspot_subsample: float = 1.0,
) -> None:
    """Rewrite 'input' paths in each spec to use basename only."""
    config = load_config(config_path)

    for key, spec in config.items():
        if isinstance(spec, dict) and "input" in spec:
            original = spec["input"]
            spec["input"] = os.path.basename(original)
            if original != spec["input"]:
                log.info("Rewrote input path: %s -> %s", original, spec["input"])
            if (
                hotspot_subsample < 1.0
                and "select_hotspots" in spec
                and isinstance(spec.get("select_hotspots"), dict)
            ):
                subsample_select_hotspots(spec["select_hotspots"], hotspot_subsample, key)

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


_BINDER_SEG = re.compile(r"^\d+-\d+$")
_TARGET_SEG = re.compile(r"^([A-Za-z]+)(-?\d+)-(-?\d+)$")


def normalise_contig_for_binder_infer(contig: str) -> str:
    """v3 pass-through when comma-separated; otherwise run v1→v3."""
    s = contig.strip()
    if not s:
        return ""
    if "," in s and " " not in s:
        return s
    return convert_v1_contigs_to_v3(s)


def split_polymers_v3(v3_contig: str) -> List[List[str]]:
    tokens = [t.strip() for t in v3_contig.strip().split(",") if t.strip()]
    polymers: List[List[str]] = []
    current: List[str] = []
    for t in tokens:
        if t == "/0":
            polymers.append(current)
            current = []
        else:
            current.append(t)
    polymers.append(current)
    return polymers


def _strip_slash_suffix(part: str) -> str:
    if "/" in part:
        return part.split("/")[0].strip()
    return part.strip()


def is_binder_polymer(segments: List[str]) -> bool:
    if not segments:
        return False
    for seg in segments:
        seg = _strip_slash_suffix(seg)
        if not seg:
            return False
        if _TARGET_SEG.match(seg):
            return False
        if not _BINDER_SEG.match(seg):
            return False
    return True


def infer_binder_chain_from_contig(contig: str) -> str:
    """Single binder polymer (length-only segments); return chr(A + polymer_index)."""
    v3 = normalise_contig_for_binder_infer(contig)
    if not v3:
        raise ValueError("Empty contig after normalisation")

    polymers = [p for p in split_polymers_v3(v3) if p]
    binder_indices = [i for i, p in enumerate(polymers) if is_binder_polymer(p)]
    if len(binder_indices) != 1:
        raise ValueError(
            f"Expected exactly one binder polymer (length-only segments); "
            f"found {len(binder_indices)} in contig {contig!r} (v3={v3!r})"
        )
    return chr(ord("A") + binder_indices[0])


def resolve_rfd3_target_binder_chain_letters(
    contig: str, explicit_binder: Optional[str] = None
) -> Tuple[str, str]:
    """Target and binder chain IDs as RFD3 assigns (first polymer A, second B, …).

    Requires exactly two polymers. Binder index from inference unless ``explicit_binder``
    gives the first MPNN designed chain letter (A or B for two polymers).
    """
    v3 = normalise_contig_for_binder_infer(contig)
    if not v3:
        raise ValueError("Empty contig after normalisation")

    polymers = [p for p in split_polymers_v3(v3) if p]
    if len(polymers) != 2:
        raise ValueError(
            f"Expected exactly two contig polymers for RFD3 target/binder chain IDs; "
            f"got {len(polymers)} in contig {contig!r} (v3={v3!r})"
        )

    if explicit_binder:
        binder = explicit_binder.strip().split(",")[0].strip().upper()
        if len(binder) != 1 or binder not in "AB":
            raise ValueError(
                f"For a two-polymer contig, --binder must be A or B; got {explicit_binder!r}"
            )
        bi = ord(binder) - ord("A")
    else:
        binder = infer_binder_chain_from_contig(contig)
        bi = ord(binder) - ord("A")

    ti = 1 - bi
    target = chr(ord("A") + ti)
    return target, binder


def contig_from_rfd3_config_path(config_path: str) -> str:
    """Read contig from first spec: JSON/YAML object, else whole file as contig string."""
    path = Path(config_path)
    text = path.read_text()
    suf = path.suffix.lower()
    if suf == ".json":
        data = json.loads(text)
        if not isinstance(data, dict) or not data:
            raise ValueError(f"RFD3 config must be a non-empty JSON object: {path}")
        spec = next(iter(data.values()))
        if not isinstance(spec, dict):
            raise ValueError(f"First config entry must be an object: {path}")
        c = spec.get("contig", "")
        if not isinstance(c, str):
            raise ValueError(f"contig must be a string in {path}")
        return c.strip()
    if suf in (".yaml", ".yml"):
        try:
            import yaml  # type: ignore[import-untyped]
        except ImportError as e:
            raise ValueError("YAML config requires PyYAML (pip install pyyaml)") from e
        data = yaml.safe_load(text)
        if not isinstance(data, dict) or not data:
            raise ValueError(f"RFD3 config must be a non-empty mapping: {path}")
        spec = next(iter(data.values()))
        if not isinstance(spec, dict):
            raise ValueError(f"First config entry must be a mapping: {path}")
        c = spec.get("contig", "")
        if not isinstance(c, str):
            raise ValueError(f"contig must be a string in {path}")
        return c.strip()
    return " ".join(text.split()).strip()


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


def extract_input_paths(config_path: str, config_dir: Optional[str] = None) -> List[str]:
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
    seen: Set[str] = set()
    result: List[str] = []
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
    frac = args.hotspot_subsample
    if not (0.0 <= frac <= 1.0):
        log.error("--hotspot-subsample must be between 0.0 and 1.0")
        return 1
    stage_config(args.config, args.output, hotspot_subsample=frac)
    return 0


def cmd_infer_binder_chain(args: argparse.Namespace) -> int:
    try:
        if args.rfd3_config is not None:
            p = Path(args.rfd3_config)
            if not p.is_file():
                log.error("Config not found: %s", p)
                return 1
            contig = contig_from_rfd3_config_path(str(p))
        else:
            contig = args.contig or ""
        chain = infer_binder_chain_from_contig(contig)
    except (ValueError, OSError, json.JSONDecodeError) as e:
        log.error("%s", e)
        return 1
    print(chain)
    return 0


def cmd_infer_rfd3_chain_pair(args: argparse.Namespace) -> int:
    try:
        if args.rfd3_config is not None:
            p = Path(args.rfd3_config)
            if not p.is_file():
                log.error("Config not found: %s", p)
                return 1
            contig = contig_from_rfd3_config_path(str(p))
        else:
            contig = args.contig or ""
        exp = args.binder.strip() if args.binder else None
        target, binder = resolve_rfd3_target_binder_chain_letters(contig, exp)
    except (ValueError, OSError, json.JSONDecodeError) as e:
        log.error("%s", e)
        return 1
    print(f"{target},{binder}")
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
    stage_parser.add_argument(
        "--hotspot-subsample",
        type=float,
        default=1.0,
        metavar="FRACTION",
        help=(
            "Fraction of select_hotspots keys to keep per spec (0.0-1.0); "
            "uses ceil(N*fraction), at least 1. Default 1.0 (no subsampling)."
        ),
    )

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

    infer_parser = subparsers.add_parser(
        "infer-binder-chain",
        help="Print binder chain letter for MPNN --designed_chains from contig or config",
    )
    infer_g = infer_parser.add_mutually_exclusive_group(required=True)
    infer_g.add_argument("--contig", help="RFD3 or v1-style contig string")
    infer_g.add_argument(
        "--rfd3-config",
        type=str,
        dest="rfd3_config",
        help="RFD3 JSON/YAML config path (first entry contig)",
    )

    pair_parser = subparsers.add_parser(
        "infer-rfd3-chain-pair",
        help="Print target,binder chain letters (RFD3 polymer order; two polymers only)",
    )
    pair_g = pair_parser.add_mutually_exclusive_group(required=True)
    pair_g.add_argument("--contig", help="RFD3 or v1-style contig string")
    pair_g.add_argument(
        "--rfd3-config",
        type=str,
        dest="rfd3_config",
        help="RFD3 JSON/YAML config path (first entry contig)",
    )
    pair_parser.add_argument(
        "--binder",
        default=None,
        help="First MPNN designed chain when not using auto (must be A or B)",
    )

    args = parser.parse_args()

    if args.command == "stage":
        return cmd_stage(args)
    if args.command == "parse-inputs":
        return cmd_parse_inputs(args)
    if args.command == "infer-binder-chain":
        return cmd_infer_binder_chain(args)
    if args.command == "infer-rfd3-chain-pair":
        return cmd_infer_rfd3_chain_pair(args)
    if args.command == "generate":
        return cmd_generate(args)

    return 1


if __name__ == "__main__":
    sys.exit(main())
