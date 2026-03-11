#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
RF3 (RosettaFold3) input generation for Nextflow.

Modes:
  cif   - Produce an RF3 input CIF with optional target chain MSA path
          (_msa_paths_by_chain_id). Legacy; prefer json for new use.
  json  - Generate RF3 input JSON with name and components (seq, optional
          msa_path, chain_id) for rf3 fold inputs=<path>.
"""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def inject_msa_into_cif(structure_cif: Path, target_msa_path: Optional[Path], target_chain: str, out_path: Path) -> None:
    """Write a CIF to out_path that mirrors structure_cif with optional _msa_paths_by_chain_id."""
    with open(structure_cif) as f:
        lines = f.readlines()

    out_lines: List[str] = []
    data_line_idx: Optional[int] = None
    for i, line in enumerate(lines):
        out_lines.append(line)
        if line.startswith("data_") and data_line_idx is None:
            data_line_idx = i
            break

    if data_line_idx is None:
        raise ValueError(f"No data_ block found in {structure_cif}")

    if target_msa_path is not None and target_msa_path.exists():
        msa_abs = str(target_msa_path.resolve())
        out_lines.append("#\n")
        out_lines.append(f"_msa_paths_by_chain_id.{target_chain}   {msa_abs}\n")
        out_lines.append("#\n")

    out_lines.extend(lines[data_line_idx + 1 :])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.writelines(out_lines)
    log.info("Wrote RF3 input CIF to %s", out_path)


def get_chain_sequence_from_cif(structure_cif: Path, chain: str, pdb_to_fasta_script: Path) -> str:
    """Run pdb_to_fasta.py and return the single sequence for the given chain."""
    result = subprocess.run(
        [sys.executable, str(pdb_to_fasta_script), str(structure_cif), "--chains", chain, "-o", "-"],
        capture_output=True,
        text=True,
        check=True,
    )
    seq = ""
    for line in result.stdout.splitlines():
        if line.startswith(">"):
            continue
        seq += line.strip()
    return seq


def generate_rf3_input(
    structure_cif: Path,
    target_chain: str,
    binder_chain: str,
    pdb_to_fasta_script: Path,
    name: str = "rf3_config",
    target_msa_path: Optional[Path] = None,
    template_structure: Optional[Path] = None,
    template_selection: Optional[str] = None,
    basename_paths: bool = False,
) -> dict:
    """Generate RF3 input JSON: name + components with seq/msa_path/chain_id and optional single template path (one per chain)."""
    target_seq = get_chain_sequence_from_cif(structure_cif, target_chain, pdb_to_fasta_script)
    binder_seq = get_chain_sequence_from_cif(structure_cif, binder_chain, pdb_to_fasta_script)
    target_comp: dict = {"seq": target_seq, "chain_id": target_chain}
    if target_msa_path is not None and target_msa_path.exists():
        target_comp["msa_path"] = (
            target_msa_path.name if basename_paths else str(target_msa_path.resolve())
        )
    binder_comp: dict = {"seq": binder_seq, "chain_id": binder_chain}
    components: list[dict] = []

    def make_path(p: Path) -> str:
        if basename_paths:
            return p.name
        return str(p.resolve())

    has_template = template_structure is not None and template_structure.exists()
    if not has_template:
        components.append(target_comp)
    if has_template and target_msa_path is not None and target_msa_path.exists():
        binder_comp["msa_path"] = (
            target_msa_path.name if basename_paths else str(target_msa_path.resolve())
        )
    components.append(binder_comp)
    if has_template:
        components.append({"path": make_path(template_structure)})

    config: dict = {
        "name": name,
        "components": components,
    }
    if has_template and template_selection:
        config["template_selection"] = [template_selection]

    return config


def main() -> int:
    parser = argparse.ArgumentParser(
        description="RF3 (RosettaFold3) input generation: CIF with MSA path or JSON for rf3 fold",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- cif subcommand (legacy) ---
    cif_parser = subparsers.add_parser(
        "cif",
        help="Produce an RF3 input CIF with optional target chain MSA path",
    )
    cif_parser.add_argument("--structure-cif", required=True, type=Path, help="Input structure CIF (target + binder)")
    # Prefer --target-msa for consistency; keep --target-msa-path as a deprecated alias.
    cif_parser.add_argument("--target-msa", default=None, type=Path, help="Path to target chain MSA (.a3m); omit for no MSA")
    cif_parser.add_argument("--target-msa-path", dest="target_msa", default=None, type=Path, help=argparse.SUPPRESS)
    cif_parser.add_argument("--target-chain", default="A", help="Target chain ID for MSA (default: A)")
    cif_parser.add_argument("-o", "--output", required=True, type=Path, help="Output CIF path")

    # --- json subcommand ---
    json_parser = subparsers.add_parser(
        "json",
        help="Generate RF3 input JSON with optional msa_path per component for rf3 fold",
    )
    json_parser.add_argument("--structure-cif", required=True, help="Structure CIF (target + binder)")
    json_parser.add_argument("--target-chain", default="A", help="Target chain ID")
    json_parser.add_argument("--binder-chain", default="B", help="Binder chain ID")
    json_parser.add_argument("--pdb-to-fasta", required=True, help="Path to pdb_to_fasta.py")
    json_parser.add_argument("--name", default="rf3_input", help="Name field in the JSON")
    json_parser.add_argument("--target-msa", default=None, help="Path to target chain MSA (.a3m); optional")
    json_parser.add_argument("--template-structure", default=None, help="Single PDB/CIF to use as RF3 template (e.g. target structure; one per chain)")
    json_parser.add_argument("--template-selection", default=None, help="Optional RF3 template_selection for template chain (e.g. 'A')")
    json_parser.add_argument("--basename-paths", action="store_true", help="Store msa_path and template paths as basenames/relative paths for Nextflow staging")
    json_parser.add_argument("-o", "--output", required=True, help="Output JSON path")

    args = parser.parse_args()

    if args.command == "cif":
        if not args.structure_cif.exists():
            log.error("Structure CIF not found: %s", args.structure_cif)
            return 1
        inject_msa_into_cif(
            args.structure_cif,
            args.target_msa,
            args.target_chain,
            args.output,
        )
        return 0

    if args.command == "json":
        config = generate_rf3_input(
            structure_cif=Path(args.structure_cif),
            target_chain=args.target_chain,
            binder_chain=args.binder_chain,
            pdb_to_fasta_script=Path(args.pdb_to_fasta),
            name=args.name,
            target_msa_path=Path(args.target_msa) if args.target_msa else None,
            template_structure=Path(args.template_structure) if args.template_structure else None,
            template_selection=args.template_selection,
            basename_paths=args.basename_paths,
        )
        with open(args.output, "w") as f:
            # RF3 expects a list of config objects; wrap single config in a list.
            json.dump([config], f, indent=2)
        log.info("Wrote RF3 input JSON list to %s", args.output)
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
