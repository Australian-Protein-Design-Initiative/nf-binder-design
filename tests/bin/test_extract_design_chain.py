#!/usr/bin/env python
# /// script
# requires-python = ">=3.10"
# dependencies = ["pytest"]
# ///

"""Tests for bin/foldseek/extract_design_chain.py.

The interesting behaviour is which chain is treated as the binder. The default
smallest-chain heuristic is right for a one-target/one-binder complex, but
BindCraft2 writes every target chain before the binder and the binder may span
several chains, so --last-chains exists to say so explicitly.
"""

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "bin" / "foldseek" / "extract_design_chain.py"
sys.path.insert(0, str(SCRIPT.parent))

import extract_design_chain as edc  # noqa: E402


def _pdb_atom(serial: int, chain: str, resseq: int) -> str:
    return (
        f"ATOM  {serial:>5d}  CA  ALA {chain}{resseq:>4d}    "
        f"{1.0:>8.3f}{2.0:>8.3f}{3.0:>8.3f}{1.0:>6.2f}{0.0:>6.2f}           C\n"
    )


def _write_pdb(path: Path, chain_sizes: dict) -> Path:
    serial = 1
    lines = []
    for chain, n_atoms in chain_sizes.items():
        for residue in range(1, n_atoms + 1):
            lines.append(_pdb_atom(serial, chain, residue))
            serial += 1
    path.write_text("".join(lines) + "END\n")
    return path


def _chains_in(pdb_path: Path) -> list:
    seen = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[21] not in seen:
            seen.append(line[21])
    return seen


def test_default_takes_the_smallest_chain(tmp_path: Path) -> None:
    source = _write_pdb(tmp_path / "complex.pdb", {"A": 50, "B": 20})
    out = tmp_path / "out"
    out.mkdir()

    edc.process_file(str(source), str(out))

    assert _chains_in(out / "complex.pdb") == ["B"]


def test_last_chains_overrides_the_heuristic(tmp_path: Path) -> None:
    """The target is the smaller chain here, so the heuristic would be wrong."""
    source = _write_pdb(tmp_path / "complex.pdb", {"A": 10, "B": 90})
    out = tmp_path / "out"
    out.mkdir()

    edc.process_file(str(source), str(out), last_chains=1)

    assert _chains_in(out / "complex.pdb") == ["B"]


def test_last_chains_keeps_a_multi_chain_binder_together(tmp_path: Path) -> None:
    """A Fab-like binder spans two chains; both must survive."""
    source = _write_pdb(tmp_path / "complex.pdb", {"A": 90, "B": 30, "C": 25})
    out = tmp_path / "out"
    out.mkdir()

    edc.process_file(str(source), str(out), last_chains=2)

    assert _chains_in(out / "complex.pdb") == ["B", "C"]


def test_last_chains_beyond_the_chain_count_keeps_everything(tmp_path: Path) -> None:
    source = _write_pdb(tmp_path / "complex.pdb", {"A": 10, "B": 10})
    out = tmp_path / "out"
    out.mkdir()

    edc.process_file(str(source), str(out), last_chains=5)

    assert _chains_in(out / "complex.pdb") == ["A", "B"]


def test_select_design_chains_uses_file_order_not_sort_order() -> None:
    """Chains are selected by where they appear, not alphabetically."""
    order = ["T", "A", "B"]
    counts = {"T": 100, "A": 20, "B": 20}

    assert edc.select_design_chains(order, counts, last_chains=2) == ["A", "B"]


def test_cli_accepts_last_chains(tmp_path: Path) -> None:
    source = _write_pdb(tmp_path / "complex.pdb", {"A": 10, "B": 90})
    out = tmp_path / "out"

    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(out), "--last-chains", "1", str(source)],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert _chains_in(out / "complex.pdb") == ["B"]


def _write_cif(path: Path, chain_sizes: dict) -> Path:
    header = (
        "data_model\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
    )
    rows = []
    serial = 1
    for chain, n_atoms in chain_sizes.items():
        for residue in range(1, n_atoms + 1):
            rows.append(
                f"ATOM {serial} C CA . ALA {chain} {residue} 1.000 2.000 3.000 1.00 0.00\n"
            )
            serial += 1
    path.write_text(header + "".join(rows) + "#\n")
    return path


def test_mmcif_last_chains_picks_the_trailing_binder(tmp_path: Path) -> None:
    """BindCraft2 writes mmCIF with the target first, and may crop it small."""
    source = _write_cif(tmp_path / "design.cif", {"A": 10, "B": 90})
    out = tmp_path / "out"
    out.mkdir()

    edc.process_file(str(source), str(out), last_chains=1)

    assert _chains_in(out / "design.pdb") == ["B"]


def test_mmcif_default_still_uses_the_heuristic(tmp_path: Path) -> None:
    source = _write_cif(tmp_path / "design.cif", {"A": 90, "B": 10})
    out = tmp_path / "out"
    out.mkdir()

    edc.process_file(str(source), str(out))

    assert _chains_in(out / "design.pdb") == ["B"]
