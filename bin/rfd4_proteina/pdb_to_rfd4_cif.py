#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["gemmi"]
# ///

"""
Convert a PDB to an mmCIF that RFD4-Proteina contigs can address.

RFD4-Proteina resolves a contig against the mmCIF `label_asym_id` and `label_seq_id`, not the
author fields. A plain `gemmi convert` is no good for that on two counts:

  - `label_asym_id` comes out as gemmi's subchain name ("Axp" for chain A's polymer), so a
    contig token like "A18-132" matches nothing.
  - `label_seq_id` is a 1-based index into the entity sequence, so chain A's residue 18 becomes
    label_seq_id 1, and every residue number a user reads off their PDB is shifted.

Both are rewritten here: `label_asym_id` becomes the chain name, and `label_seq_id` becomes the
PDB residue number. A PDB has only author numbering to begin with, so this keeps contigs and
hotspots in the numbering the user already sees, matching how --contigs and --hotspot_res behave
in the other nf-binder-design workflows.

Chains that would collide on `label_asym_id` (a ligand or waters sharing a chain name with a
polymer) keep gemmi's distinct subchain names, since label_asym_id has to stay unique.
"""

import argparse
import logging
import sys
from pathlib import Path

import gemmi

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def convert(input_path: Path, output_path: Path) -> None:
    structure = gemmi.read_structure(str(input_path))
    structure.setup_entities()

    # The first subchain of each chain takes the chain's name; the rest keep gemmi's, so
    # label_asym_id stays unique.
    renames = {}
    claimed = set()
    for model in structure:
        for chain in model:
            first_subchain = next((r.subchain for r in chain), None)
            if first_subchain is not None and chain.name not in claimed:
                renames[first_subchain] = chain.name
                claimed.add(chain.name)
            for residue in chain:
                residue.label_seq = residue.seqid.num

    for model in structure:
        for chain in model:
            for residue in chain:
                residue.subchain = renames.get(residue.subchain, residue.subchain)

    # Entities index subchains by name, so they are remapped rather than re-derived: calling
    # setup_entities() again here would rebuild entities around the new names and write the
    # subchain name into label_entity_id.
    for entity in structure.entities:
        entity.subchains = [renames.get(name, name) for name in entity.subchains]

    # gemmi names entities after a subchain, which writes a non-numeric label_entity_id that
    # atomworks warns about. Number them instead, as mmCIF expects.
    for index, entity in enumerate(structure.entities, start=1):
        entity.name = str(index)

    doc = structure.make_mmcif_document()
    doc.write_file(str(output_path))
    log.info("Wrote %s", output_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a PDB to an mmCIF addressable by RFD4-Proteina contigs",
        allow_abbrev=False,
    )
    parser.add_argument("input", help="Input PDB (or any structure gemmi reads)")
    parser.add_argument("-o", "--output", required=True, help="Output mmCIF path")
    args = parser.parse_args()

    convert(Path(args.input), Path(args.output))


if __name__ == "__main__":
    main()
