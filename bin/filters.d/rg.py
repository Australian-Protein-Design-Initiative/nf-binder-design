#!/usr/bin/env python
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "mdanalysis",
#     "gemmi",
#     "pandas",
# ]
# ///

import argparse
import math
import sys
import pandas as pd
import MDAnalysis as mda
import gemmi
from pathlib import Path
import logging


def register_metrics() -> list[str]:
    """Returns a list of metrics that this plugin can calculate."""
    return ["rg"]


def calculate_rg_mmcif(structure_file: str, chain: str = "A") -> float:
    """Calculates radius of gyration for one chain from an mmCIF structure."""
    structure = gemmi.read_structure(structure_file)
    if len(structure) == 0:
        raise ValueError(f"No models found in {structure_file}")

    model = structure[0]
    coords: list[tuple[float, float, float]] = []
    for ch in model:
        if ch.name != chain:
            continue
        for residue in ch:
            for atom in residue:
                pos = atom.pos
                coords.append((pos.x, pos.y, pos.z))

    if not coords:
        raise ValueError(f"Chain {chain} not found in {structure_file}")

    n = float(len(coords))
    cx = sum(p[0] for p in coords) / n
    cy = sum(p[1] for p in coords) / n
    cz = sum(p[2] for p in coords) / n
    rg2 = sum((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 for x, y, z in coords) / n
    return math.sqrt(rg2)


def calculate_rg(structure_file: str, chain: str = "A") -> float:
    """Calculates the radius of gyration for a given chain in a structure file."""
    structure_path = Path(structure_file)
    is_mmcif = structure_path.suffix == ".cif" or structure_path.name.endswith(".cif.gz")

    if is_mmcif:
        return calculate_rg_mmcif(structure_file, chain)

    u = mda.Universe(structure_file)
    selection = u.select_atoms(f"chainID {chain}")
    if not selection:
        raise ValueError(f"Chain {chain} not found in {structure_file}")
    return selection.radius_of_gyration()


def calculate_metrics(
    structure_files: list[str], binder_chains: list[str] = ["A"]
) -> pd.DataFrame:
    """
    Calculates the radius of gyration for a list of structure files.
    """
    # For now, we only handle a single binder chain for simplicity
    if len(binder_chains) > 1:
        logging.warning(
            "This plugin currently only calculates rg for the first specified binder chain."
        )

    chain = binder_chains[0]

    results = []
    for structure_path in structure_files:
        try:
            rg = calculate_rg(str(structure_path), chain)
            results.append({"design_id": Path(structure_path).stem, "rg": rg})
        except Exception as e:
            logging.error(f"Error processing {structure_path}: {e}")
            results.append({"design_id": Path(structure_path).stem, "rg": None})

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.set_index("design_id")
    return df


def main():
    """For standalone execution."""
    parser = argparse.ArgumentParser(
        description="Calculate Radius of Gyration (Rg) for a structure file."
    )
    parser.add_argument("structure_file", help="Path to the structure file.")
    parser.add_argument(
        "--chain", default="A", help="Chain to calculate Rg for. Defaults to 'A'."
    )
    parser.add_argument(
        "--format", default="tsv", choices=["tsv", "json"], help="Output format."
    )

    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, stream=sys.stderr, format="%(levelname)s: %(message)s"
    )

    try:
        rg = calculate_rg(args.structure_file, args.chain)
        design_id = Path(args.structure_file).stem

        if args.format == "json":
            import json

            print(json.dumps({design_id: {"rg": rg}}))
        else:  # tsv
            print(f"design_id\\trg")
            print(f"{design_id}\\t{rg}")

    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
