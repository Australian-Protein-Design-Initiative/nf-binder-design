#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Tests for bin/foldseek/cluster_table.py.

    uv run --with pytest pytest tests/bin/test_cluster_table.py -q

The input format is foldseek easy-cluster's <representative>TAB<member> table, which has
one row per design including the representative's own row.
"""

import csv
import importlib.util
from pathlib import Path
from typing import Dict, List

import pytest

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "foldseek" / "cluster_table.py"
_spec = importlib.util.spec_from_file_location("cluster_table", _MODPATH)
cluster_table = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cluster_table)


def rows_by_id(pairs: List[tuple]) -> Dict[str, Dict[str, str]]:
    return {r["id"]: r for r in cluster_table.build_rows(pairs)}


def test_a_cluster_of_two_and_a_singleton() -> None:
    """The duplicate-design case: foldseek groups the copy with its source."""
    by_id = rows_by_id([("a", "a"), ("a", "b"), ("c", "c")])

    assert len(by_id) == 3
    assert by_id["a"]["foldseek_cluster_size"] == "2"
    assert by_id["b"]["foldseek_cluster_size"] == "2"
    assert by_id["b"]["foldseek_cluster_rep"] == "a"
    assert by_id["c"]["foldseek_cluster_size"] == "1"


def test_only_the_representative_is_flagged() -> None:
    by_id = rows_by_id([("a", "a"), ("a", "b")])
    assert by_id["a"]["foldseek_cluster_is_rep"] == "1"
    assert by_id["b"]["foldseek_cluster_is_rep"] == "0"


def test_clusters_are_numbered_largest_first() -> None:
    """So cluster 1 is the most populated, which is how a diversity summary reads."""
    by_id = rows_by_id([("solo", "solo"), ("big", "big"), ("big", "big2"), ("big", "big3")])
    assert by_id["big"]["foldseek_cluster_id"] == "1"
    assert by_id["solo"]["foldseek_cluster_id"] == "2"


def test_every_design_appears_exactly_once(tmp_path: Path) -> None:
    """A design listed under two representatives would double-count in a join."""
    pairs = [("a", "a"), ("a", "b"), ("c", "c"), ("c", "d")]
    ids = [r["id"] for r in cluster_table.build_rows(pairs)]
    assert sorted(ids) == ["a", "b", "c", "d"]
    assert len(ids) == len(set(ids))


@pytest.mark.parametrize(
    "name,expected",
    [("x.pdb", "x"), ("x.cif", "x"), ("x.cif.gz", "x"), ("x.pdb.gz", "x"), ("x", "x"),
     ("pdl1_b0_pdl1_0.cif", "pdl1_b0_pdl1_0")],
)
def test_extensions_are_stripped_so_ids_join(name: str, expected: str) -> None:
    """The id has to match the design table's `id`, which carries no extension."""
    assert cluster_table.strip_suffix(name) == expected


def test_reads_a_real_foldseek_table(tmp_path: Path) -> None:
    src = tmp_path / "clusters.tsv"
    src.write_text("rep1\trep1\nrep1\tmember1\n\nrep2\trep2\n")
    pairs = cluster_table.read_clusters(src)
    # The blank line is skipped rather than becoming an empty id
    assert pairs == [("rep1", "rep1"), ("rep1", "member1"), ("rep2", "rep2")]


def test_written_table_has_the_documented_columns(tmp_path: Path) -> None:
    src = tmp_path / "clusters.tsv"
    src.write_text("a\ta\na\tb\n")
    out = tmp_path / "table.tsv"
    import sys
    from unittest import mock

    argv = ["cluster_table.py", str(src), "-o", str(out)]
    with mock.patch.object(sys, "argv", argv):
        cluster_table.main()

    with open(out, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == cluster_table.FIELDS
        assert len(list(reader)) == 2
