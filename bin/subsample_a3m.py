#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# ///

"""
ColabFold / AlphaFold2-style MSA row subsampling for an a3m alignment.

Mirrors alphafold.model.tf.data_transforms.sample_msa + crop_extra_msa
(monomer path):

1. Always keep the query (row 0).
2. Shuffle remaining rows without replacement; take the first max_seq as
   selected ("cluster centres").
3. Independently shuffle leftovers and keep up to max_extra_seq as extras.

Uniform over MSA rows (not HHBlits clusters). Bit-exact parity with TF
random_shuffle is not attempted - matching the algorithm is enough for
conformational sampling. Stdlib only so it runs in AF2/Boltz/RF3/Protenix
containers without extra deps.
"""

from __future__ import annotations

import argparse
import logging
import random
import sys
from pathlib import Path
from typing import List, Sequence, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

# (header_without_gt, sequence, 0-based file line of the '>' header)
A3mRecord = Tuple[str, str, int]


def parse_a3m_records(text: str) -> List[A3mRecord]:
    """Return [(header_without_gt, sequence, header_line_idx), ...] in file order."""
    records: List[A3mRecord] = []
    header: str | None = None
    header_line: int = 0
    seq_parts: List[str] = []
    for line_idx, line in enumerate(text.splitlines()):
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts), header_line))
            header = line[1:]
            header_line = line_idx
            seq_parts = []
        elif header is not None:
            seq_parts.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_parts), header_line))
    return records


def format_a3m(records: Sequence[A3mRecord]) -> str:
    lines: List[str] = []
    for header, seq, _line in records:
        lines.append(f">{header}")
        lines.append(seq)
    return "\n".join(lines) + ("\n" if records else "")


def sequence_id(header: str) -> str:
    """First whitespace-delimited token of an a3m/FASTA header (no leading '>')."""
    return header.split()[0] if header.strip() else ""


def write_ids(records: Sequence[A3mRecord], path: Path) -> None:
    """Write '<header_line_0idx>\\t<id>' per kept sequence (original a3m line)."""
    rows: List[str] = []
    for header, _seq, line_idx in records:
        sid = sequence_id(header)
        if sid:
            rows.append(f"{line_idx}\t{sid}\n")
    path.write_text("".join(rows))


def clamp_limits(n: int, max_seq: int, max_extra_seq: int) -> Tuple[int, int]:
    """Clamp like ColabFold when the MSA is shallower than requested."""
    if n < 1:
        raise ValueError("MSA has no sequences")
    max_seq_c = min(max_seq, n)
    remaining = n - max_seq_c
    if remaining <= 0:
        max_extra_c = 0
    else:
        max_extra_c = max(min(remaining, max_extra_seq), 1) if max_extra_seq > 0 else 0
        max_extra_c = min(remaining, max_extra_c)
    return max_seq_c, max_extra_c


def subsample_records(
    records: Sequence[A3mRecord],
    max_seq: int,
    max_extra_seq: int,
    seed: int,
) -> List[A3mRecord]:
    n = len(records)
    max_seq_c, max_extra_c = clamp_limits(n, max_seq, max_extra_seq)
    rng = random.Random(seed)

    if n == 1:
        return [records[0]]

    rest = list(range(1, n))
    rng.shuffle(rest)
    order = [0] + rest
    selected_idx = order[:max_seq_c]
    extra_idx = order[max_seq_c:]
    if max_extra_c > 0 and extra_idx:
        rng.shuffle(extra_idx)
        extra_idx = extra_idx[:max_extra_c]
    else:
        extra_idx = []

    out_idx = list(selected_idx) + list(extra_idx)
    return [records[i] for i in out_idx]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--a3m", required=True, help="Input a3m (query-first)")
    parser.add_argument(
        "--ids-only",
        action="store_true",
        help="Only write sequence IDs from --a3m (no subsample); requires --ids-output",
    )
    parser.add_argument("--max-seq", type=int, help="Selected MSA rows including query")
    parser.add_argument(
        "--max-extra-seq",
        type=int,
        help="Extra MSA rows after selected (second shuffle + crop)",
    )
    parser.add_argument("--seed", type=int, help="RNG seed")
    parser.add_argument(
        "-o",
        "--output",
        help="Output a3m path, or '-' for stdout (required unless --ids-only)",
    )
    parser.add_argument(
        "--ids-output",
        help="Write '<0-based header line>\\t<id>' per sequence (original a3m lines)",
    )
    args = parser.parse_args()

    a3m_path = Path(args.a3m)
    text = a3m_path.read_text()
    records = parse_a3m_records(text)
    if not records:
        log.error("%s contains no sequences", a3m_path)
        return 1

    if args.ids_only:
        if not args.ids_output:
            log.error("--ids-only requires --ids-output")
            return 1
        write_ids(records, Path(args.ids_output))
        log.info("Wrote %d sequence IDs from %s to %s", len(records), a3m_path, args.ids_output)
        return 0

    if args.max_seq is None or args.max_extra_seq is None or args.seed is None or not args.output:
        log.error("--max-seq, --max-extra-seq, --seed and -o/--output are required unless --ids-only")
        return 1
    if args.max_seq < 1:
        log.error("--max-seq must be >= 1")
        return 1
    if args.max_extra_seq < 0:
        log.error("--max-extra-seq must be >= 0")
        return 1

    out_records = subsample_records(records, args.max_seq, args.max_extra_seq, args.seed)
    out_text = format_a3m(out_records)
    log.info(
        "Subsampled %s: %d -> %d sequences (max-seq=%d, max-extra-seq=%d, seed=%d)",
        a3m_path,
        len(records),
        len(out_records),
        args.max_seq,
        args.max_extra_seq,
        args.seed,
    )

    if args.output == "-":
        sys.stdout.write(out_text)
    else:
        Path(args.output).write_text(out_text)

    if args.ids_output:
        write_ids(out_records, Path(args.ids_output))
        log.info("Wrote %d sequence IDs to %s", len(out_records), args.ids_output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
