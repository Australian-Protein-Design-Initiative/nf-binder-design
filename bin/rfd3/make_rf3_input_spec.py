#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
RF3 (RosettaFold3) input generation for Nextflow.

Modes:
  json-batch  - One JSON list for multiple CIFs (shared template/MSA): pass
          --structure-cifs and --names in the same order (length must match).
"""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def parse_template_selection_arg(raw: Optional[str]) -> Optional[List[str]]:
    """Split comma-separated RF3 AtomSelection tokens (commas inside a token are not supported)."""
    if raw is None or not str(raw).strip():
        return None
    return [s.strip() for s in str(raw).split(",") if s.strip()]


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
    template_selection: Optional[Sequence[str]] = None,
    basename_paths: bool = False,
    binder_sequence_chain: Optional[str] = None,
) -> dict:
    """Generate RF3 input JSON: name + components with seq/msa_path/chain_id and optional single template path.

    ``template_selection``: list of RF3 AtomSelection tokens (e.g. ``["A"]`` or ``["A/*/1-50", "A/*/80-120"]``).

    Without a template: [target (seq, optional msa_path), binder (seq)].

    With a template: target coordinates come from the template file only — RF3 rejects duplicate
    chain_id across components, so we must not add a separate seq component for the target chain.
    The binder sequence component plus template path; target MSA (if any) is set via top-level
    ``msa_paths`` (merged into chain_info in RF3), mapped to target_chain — never on the binder.

    ``binder_sequence_chain``: chain in the MPNN/RFD3 CIF where the binder sequence is read
    (defaults to ``binder_chain``). JSON ``chain_id`` values use ``target_chain`` and ``binder_chain``
    so RosettaFold3 matches RFDiffusion3 polymer order (template trimmed/renamed to ``target_chain``).
    """
    seq_chain = binder_sequence_chain if binder_sequence_chain else binder_chain
    binder_seq = get_chain_sequence_from_cif(structure_cif, seq_chain, pdb_to_fasta_script)
    binder_comp: dict = {"seq": binder_seq, "chain_id": binder_chain}
    components: list[dict] = []

    def make_path(p: Path) -> str:
        if basename_paths:
            return p.name
        return str(p.resolve())

    def msa_path_str() -> str:
        assert target_msa_path is not None
        return (
            target_msa_path.name if basename_paths else str(target_msa_path.resolve())
        )

    has_template = template_structure is not None and template_structure.exists()
    if has_template:
        components.append(binder_comp)
        components.append({"path": make_path(template_structure)})
    else:
        target_seq = get_chain_sequence_from_cif(structure_cif, target_chain, pdb_to_fasta_script)
        target_comp: dict = {"seq": target_seq, "chain_id": target_chain}
        if target_msa_path is not None and target_msa_path.exists():
            target_comp["msa_path"] = msa_path_str()
        components.append(target_comp)
        components.append(binder_comp)

    config: dict = {
        "name": name,
        "components": components,
    }
    if has_template and template_selection:
        config["template_selection"] = list(template_selection)
    if has_template and target_msa_path is not None and target_msa_path.exists():
        config["msa_paths"] = {target_chain: msa_path_str()}

    return config


def main() -> int:
    parser = argparse.ArgumentParser(
        description="RF3 (RosettaFold3) batch input JSON generation for rf3 fold",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    batch_parser = subparsers.add_parser(
        "json-batch",
        help="Generate RF3 inputs JSON list from multiple structure CIFs (shared options)",
    )
    batch_parser.add_argument(
        "--structure-cifs",
        nargs="+",
        required=True,
        help="Structure CIF paths (space-separated; shell may expand globs before Python)",
    )
    batch_parser.add_argument(
        "--names",
        nargs="+",
        required=True,
        help="RF3 'name' field per CIF, same order and count as --structure-cifs",
    )
    batch_parser.add_argument("--target-chain", default="A", help="Target chain ID")
    batch_parser.add_argument(
        "--binder-chain",
        default="B",
        help="chain_id for binder in RF3 JSON",
    )
    batch_parser.add_argument(
        "--binder-sequence-chain",
        default=None,
        help="Chain in each CIF for binder sequence (default: --binder-chain)",
    )
    batch_parser.add_argument("--pdb-to-fasta", required=True, help="Path to pdb_to_fasta.py")
    batch_parser.add_argument("--target-msa", default=None, help="Shared target MSA (.a3m); optional")
    batch_parser.add_argument("--template-structure", default=None, help="Shared RF3 template path")
    batch_parser.add_argument(
        "--template-selection",
        default=None,
        help="Comma-separated RF3 template_selection AtomSelection tokens",
    )
    batch_parser.add_argument("--basename-paths", action="store_true", help="Basenames for staged paths")
    batch_parser.add_argument("-o", "--output", required=True, help="Output JSON path")

    args = parser.parse_args()

    if args.command == "json-batch":
        cifs = [Path(p) for p in args.structure_cifs]
        names = list(args.names)
        if len(cifs) != len(names):
            log.error("--structure-cifs (%d) and --names (%d) must have the same length", len(cifs), len(names))
            return 1
        for p in cifs:
            if not p.exists():
                log.error("Structure CIF not found: %s", p)
                return 1
        tm = Path(args.target_msa) if args.target_msa else None
        tpl = Path(args.template_structure) if args.template_structure else None
        tsel = parse_template_selection_arg(args.template_selection)
        configs: list[dict] = []
        for cif_path, name in zip(cifs, names):
            configs.append(
                generate_rf3_input(
                    structure_cif=cif_path,
                    target_chain=args.target_chain,
                    binder_chain=args.binder_chain,
                    pdb_to_fasta_script=Path(args.pdb_to_fasta),
                    name=name,
                    target_msa_path=tm,
                    template_structure=tpl,
                    template_selection=tsel,
                    basename_paths=args.basename_paths,
                    binder_sequence_chain=args.binder_sequence_chain,
                )
            )
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            json.dump(configs, f, indent=2)
        log.info("Wrote RF3 batch input JSON (%d examples) to %s", len(configs), out_path)
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
