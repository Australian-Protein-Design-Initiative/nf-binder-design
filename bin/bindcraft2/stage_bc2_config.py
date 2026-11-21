#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Helper script for BindCraft2 (BC2) campaign config handling in Nextflow.

Modes:
  stage        - Rewrite file paths in an existing campaign JSON so they point
                 at Nextflow-staged files (basename only), and optionally
                 subsample each target's hotspots.
  generate     - Create a new campaign JSON from CLI flags, translating
                 BindCraft v1 style target/hotspot options to the BC2 schema.
  parse-inputs - Print resolved input file path(s) from a campaign JSON
                 (one per line).

BC2 resolves `targets[].target_path` and `binder_scaffold` relative to the
campaign JSON's own directory (bindcraft/settings.py read_settings()). Nextflow
stages the referenced files beside the rewritten JSON in the task work
directory, so rewriting them to bare basenames is what makes them resolve.
"""

import argparse
import json
import logging
import math
import os
import random
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Set

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

# The only two settings BC2 treats as campaign-relative file paths.
SCAFFOLD_KEY = "binder_scaffold"
TARGET_PATH_KEY = "target_path"


def load_config(path: str) -> dict:
    with open(path) as f:
        config = json.load(f)
    if not isinstance(config, dict):
        raise ValueError(f"campaign config must be a JSON object: {path}")
    return config


def save_config(config: dict, path: str) -> None:
    # Write to a temporary file and rename over the target. Nextflow stages
    # inputs as symlinks, so a plain open(path, "w") on a name that collides
    # with a staged input would follow the link and silently overwrite the
    # source file outside the work directory. os.replace swaps the symlink
    # itself, and is atomic, so a crash cannot leave a half-written config.
    directory = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(dir=directory, prefix=".stage_bc2_config.", suffix=".json")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(config, f, indent=2)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    log.info("Wrote config to %s", path)


def subsample_hotspots(hotspots: str, fraction: float, seed=None) -> str:
    """Keep ceil(N * fraction) of a comma-separated hotspot list, at least 1.

    Each comma-separated token counts as one item, so a BC2 range such as
    "B12-16" is kept or dropped whole rather than being split. The survivors
    stay in their original order, which keeps the result readable and makes a
    diff against the input obvious.

    Unlike the BindCraft v1 helper, which seeds from system time, the seed is
    passed in: the workflow gives each batch its own campaign_seed, so batches
    get different hotspot subsets but the run as a whole still reproduces.
    """
    items = [h.strip() for h in hotspots.split(",") if h.strip()]
    if not items or fraction >= 1.0:
        return hotspots

    n_keep = max(1, math.ceil(len(items) * fraction))
    rng = random.Random(seed)
    kept_indices = sorted(rng.sample(range(len(items)), n_keep))
    return ",".join(items[i] for i in kept_indices)


def stage_config(
    config_path: str,
    output_path: str,
    hotspot_subsample: float = 1.0,
    subsample_seed: Optional[int] = None,
) -> None:
    """Rewrite campaign-relative file paths to basenames, optionally thinning hotspots."""
    config = load_config(config_path)

    for target in config.get("targets", []) or []:
        if isinstance(target, dict) and target.get(TARGET_PATH_KEY):
            original = str(target[TARGET_PATH_KEY])
            target[TARGET_PATH_KEY] = os.path.basename(original)
            if original != target[TARGET_PATH_KEY]:
                log.info("Rewrote target_path: %s -> %s", original, target[TARGET_PATH_KEY])

    if hotspot_subsample < 1.0:
        for index, target in enumerate(config.get("targets", []) or []):
            if not isinstance(target, dict) or not target.get("hotspots"):
                continue
            original = str(target["hotspots"])
            # Vary the draw per target as well as per batch, or a multi-target
            # campaign would thin every target to the same positions.
            per_target_seed = None if subsample_seed is None else f"{subsample_seed}:{index}"
            thinned = subsample_hotspots(original, hotspot_subsample, per_target_seed)
            target["hotspots"] = thinned
            log.info(
                "Subsampled hotspots for target %d: %s -> %s", index, original, thinned
            )

    if config.get(SCAFFOLD_KEY):
        original = str(config[SCAFFOLD_KEY])
        config[SCAFFOLD_KEY] = os.path.basename(original)
        if original != config[SCAFFOLD_KEY]:
            log.info("Rewrote binder_scaffold: %s -> %s", original, config[SCAFFOLD_KEY])

    save_config(config, output_path)


def extract_input_paths(config_path: str, config_dir: Optional[str] = None) -> List[str]:
    """Resolved file paths referenced by a campaign JSON, in order of first occurrence.

    Returns an empty list for a campaign that only names shipped presets (eg
    {"target": "hPDL1"}), whose paths resolve inside the container instead.
    """
    path_obj = Path(config_path)
    if not path_obj.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    base = Path(config_dir).resolve() if config_dir else path_obj.parent.resolve()
    config = load_config(config_path)

    raw_paths = [
        str(target[TARGET_PATH_KEY])
        for target in (config.get("targets", []) or [])
        if isinstance(target, dict) and target.get(TARGET_PATH_KEY)
    ]
    if config.get(SCAFFOLD_KEY):
        raw_paths.append(str(config[SCAFFOLD_KEY]))

    seen: Set[str] = set()
    result: List[str] = []
    for raw in raw_paths:
        resolved = Path(raw) if os.path.isabs(raw) else base / raw
        resolved_str = str(resolved.resolve())
        if resolved_str not in seen:
            seen.add(resolved_str)
            result.append(resolved_str)
    return result


def parse_binder_length_range(binder_length_range: str) -> List[int]:
    """'60-150' -> [60, 150]; '80' -> [80]; '60,80,100' -> [60, 80, 100].

    BC2 reads binder_lengths as an inclusive range for two values, an exact
    length for one, and a set of discrete choices for more.
    """
    text = (binder_length_range or "").strip()
    if not text:
        raise ValueError("binder_length_range cannot be empty")

    separator = "," if "," in text else "-"
    parts = [p.strip() for p in text.split(separator) if p.strip()]
    if not parts:
        raise ValueError(f"no lengths found in binder_length_range: {binder_length_range!r}")

    try:
        lengths = [int(p) for p in parts]
    except ValueError:
        raise ValueError(
            f"binder_length_range must be integers, eg '60-150' or '80': {binder_length_range!r}"
        )

    if any(v <= 0 for v in lengths):
        raise ValueError(f"binder lengths must be positive: {binder_length_range!r}")
    if separator == "-" and len(lengths) == 2 and lengths[0] > lengths[1]:
        raise ValueError(
            f"binder_length_range start must not exceed end: {binder_length_range!r}"
        )
    return lengths


def generate_config(
    design_name: str,
    input_pdb: str,
    target_chains: Optional[str] = None,
    hotspot_res: Optional[str] = None,
    binder_length_range: Optional[str] = None,
    target_name: Optional[str] = None,
) -> dict:
    """Build a BC2 campaign JSON from BindCraft v1 style flags.

    v1's flat starting_pdb / chains / target_hotspot_residues / lengths become
    BC2's targets list plus binder_lengths.

    project_folder, campaign_seed, max_trajectories and number_of_final_designs
    are deliberately omitted: the workflow supplies them per batch via --set, so
    that a generated config is identical across batches and only the overrides
    differ.
    """
    target: dict = {
        "name": target_name or Path(input_pdb).stem,
        TARGET_PATH_KEY: os.path.basename(input_pdb),
    }
    if target_chains:
        target["chains"] = target_chains
    if hotspot_res:
        target["hotspots"] = hotspot_res

    config: dict = {"campaign_name": design_name, "targets": [target]}
    if binder_length_range:
        config["binder_lengths"] = parse_binder_length_range(binder_length_range)
    return config


def cmd_stage(args: argparse.Namespace) -> int:
    if not 0.0 <= args.hotspot_subsample <= 1.0:
        raise ValueError("--hotspot-subsample must be between 0.0 and 1.0")
    stage_config(
        args.config,
        args.output,
        hotspot_subsample=args.hotspot_subsample,
        subsample_seed=args.subsample_seed,
    )
    return 0


def cmd_parse_inputs(args: argparse.Namespace) -> int:
    for path in extract_input_paths(args.config, args.config_dir):
        print(path)
    return 0


def cmd_generate(args: argparse.Namespace) -> int:
    config = generate_config(
        design_name=args.design_name,
        input_pdb=args.input_pdb,
        target_chains=args.target_chains,
        hotspot_res=args.hotspot_res,
        binder_length_range=args.binder_length_range,
        target_name=args.target_name,
    )
    save_config(config, args.output)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="BindCraft2 campaign config helper for Nextflow")
    subparsers = parser.add_subparsers(dest="command", required=True)

    stage_parser = subparsers.add_parser(
        "stage",
        help="Rewrite campaign file paths to basenames so they resolve beside the staged JSON",
    )
    stage_parser.add_argument("config", help="Path to the campaign JSON")
    stage_parser.add_argument("-o", "--output", required=True, help="Output JSON path")
    stage_parser.add_argument(
        "--hotspot-subsample",
        type=float,
        default=1.0,
        help="Fraction of each target's hotspots to keep, ceil(N*fraction) with a minimum of 1 "
        "(default: 1.0, keep all)",
    )
    stage_parser.add_argument(
        "--subsample-seed",
        type=int,
        default=None,
        help="Seed for the hotspot draw (default: unseeded)",
    )

    parse_parser = subparsers.add_parser(
        "parse-inputs",
        help="Print resolved input file path(s) from a campaign JSON (one per line)",
    )
    parse_parser.add_argument("config", help="Path to the campaign JSON")
    parse_parser.add_argument(
        "--config-dir",
        default=None,
        help="Directory for resolving relative paths (default: the config file's directory)",
    )

    gen_parser = subparsers.add_parser(
        "generate", help="Generate a campaign JSON from CLI parameters"
    )
    gen_parser.add_argument("--design-name", default="bc2_design", help="campaign_name value")
    gen_parser.add_argument("--input-pdb", required=True, help="Target PDB/mmCIF/FASTA file")
    gen_parser.add_argument("--target-name", default=None, help="Target name (default: input stem)")
    gen_parser.add_argument("--target-chains", default=None, help="Target chain(s), eg 'A' or 'A,B'")
    gen_parser.add_argument(
        "--hotspot-res", default=None, help="Hotspot residues, eg 'A54,A56' or 'B12-16'"
    )
    gen_parser.add_argument(
        "--binder-length-range",
        default=None,
        help="'60-150' for a range, '80' for an exact length, '60,80,100' for discrete choices",
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
