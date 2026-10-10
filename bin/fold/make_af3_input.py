#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
AlphaFold3 input JSON (alphafold3 dialect, version 2) for fold.nf's ALPHAFOLD3_FOLD.

By default we run AF3 with --run_data_pipeline=false, feeding MSAs computed by the
pipeline's shared MSA stage. In that mode AF3 requires unpairedMsa, pairedMsa and
templates to be non-null for every protein chain (data/featurisation.py), so every
chain gets an unpaired and a paired a3m file (query-only when there is nothing
better) and "templates": [].

--paired-msa-mode empty and --templates-mode search relax that for Germinal
parity: Germinal gives AF3 no pairing at all ("pairedMsa": "") and no templates
key, and lets AF3's data pipeline run so it searches pdb_seqres/mmcif_files for
templates on every chain. --templates-mode search therefore REQUIRES
--run_data_pipeline=true; with the pipeline off AF3 rejects the missing key.

The a3m files are rewritten next to the JSON as chain_<ID>_{unpaired,paired}.a3m and
referenced by basename (*MsaPath resolves relative to the JSON), so the predict task
can subsample chain_A_unpaired.a3m in place without rewriting the JSON:
  - ColabFold '#' header lines and NUL separators are dropped (AF3 requires the
    first record to be the query).
  - Row 0 must equal the chain sequence, else we fail here rather than inside AF3.

Templates: with --templates-dir (bin/fold/match_templates.py output), each chain
listed in --template-chains (default: all) gets its matched templates inline
(`mmcif` plus queryIndices/templateIndices), looked up by sequence md5.

Multimer pairing: AF3 pairs pairedMsa rows across chains by UniProt species mnemonic
parsed from `tr|ACC|NAME_SPECIES` headers - render those with
`bin/fold/msa_taxonomy.py --tool af3` and pass them via --paired-a3m.
"""

import argparse
import hashlib
import json
import logging
import string
import sys
from pathlib import Path
from typing import List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

CHAIN_IDS = string.ascii_uppercase


def parse_fasta_records(fasta_path: Path) -> List[str]:
    """Return one sequence string per FASTA record, in file order (no external deps)."""
    records: List[str] = []
    current: List[str] = []
    for line in fasta_path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                records.append("".join(current))
                current = []
        else:
            current.append(line)
    if current:
        records.append("".join(current))
    return records


def sanitised_name(name: str) -> str:
    """AF3's Input.sanitised_name(): the output directory / file prefix it will use."""
    allowed = set(string.ascii_letters + string.digits + "_-.")
    return "".join(c for c in name.replace(" ", "_") if c in allowed)


def read_a3m_records(text: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    seq: List[str] = []
    for raw in text.replace("\x00", "").splitlines():
        line = raw.rstrip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header = line[1:]
            seq = []
        elif header is not None:
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def clean_a3m(src: Optional[Path], query_seq: str, chain_id: str) -> str:
    """Return AF3-ready a3m text: query first, no '#' lines. Query-only if src is None."""
    if src is None:
        return f">query\n{query_seq}\n"
    records = read_a3m_records(src.read_text())
    if not records:
        log.warning("chain %s: %s has no records - using query-only MSA", chain_id, src.name)
        return f">query\n{query_seq}\n"
    first = records[0][1].replace("-", "")
    if first.upper() != query_seq.upper():
        raise ValueError(
            f"chain {chain_id}: first sequence in {src.name} does not match the FASTA "
            f"sequence (the query must be the first MSA row)"
        )
    return "".join(f">{h}\n{s}\n" for h, s in records)


def _match_per_chain(paths: Optional[List[Path]], n_seq: int, flag: str) -> List[Optional[Path]]:
    if not paths:
        return [None] * n_seq
    if len(paths) != n_seq:
        raise ValueError(
            f"{flag} got {len(paths)} file(s) for {n_seq} chain(s); "
            f"pass exactly one per chain in record order"
        )
    return list(paths)


def chain_templates(templates_dir: Optional[Path], seq: str) -> List[dict]:
    if templates_dir is None or not (templates_dir / "index.json").exists():
        return []
    index = json.loads((templates_dir / "index.json").read_text())
    hits = index.get(hashlib.md5(seq.upper().encode()).hexdigest(), [])
    return [
        {
            "mmcif": (templates_dir / h["cif"]).read_text(),
            "queryIndices": h["query_indices"],
            "templateIndices": h["template_indices"],
        }
        for h in hits
    ]


def make_af3_input(
    fasta_path: Path,
    name: str,
    out_dir: Path,
    seed,
    unpaired_a3m_paths: Optional[List[Path]] = None,
    paired_a3m_paths: Optional[List[Path]] = None,
    templates_dir: Optional[Path] = None,
    template_chains: Optional[List[str]] = None,
    paired_msa_mode: str = "path",
    templates_mode: str = "inline",
) -> dict:
    # `seed` accepts a single int or a list; AF3 runs --num_diffusion_samples
    # structures for each entry of modelSeeds.
    seeds = [seed] if isinstance(seed, int) else list(seed)

    sequences = parse_fasta_records(fasta_path)
    if not sequences:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    if len(sequences) > len(CHAIN_IDS):
        raise ValueError(f"{len(sequences)} chains exceeds the {len(CHAIN_IDS)} supported")

    n = len(sequences)
    unpaired = _match_per_chain(unpaired_a3m_paths, n, "--a3m")
    paired = _match_per_chain(paired_a3m_paths, n, "--paired-a3m")

    entries = []
    for i, seq in enumerate(sequences):
        cid = CHAIN_IDS[i]
        protein = {"id": cid, "sequence": seq}

        unpaired_name = f"chain_{cid}_unpaired.a3m"
        (out_dir / unpaired_name).write_text(clean_a3m(unpaired[i], seq, cid))
        protein["unpairedMsaPath"] = unpaired_name

        if paired_msa_mode == "empty":
            # Germinal parity: no cross-chain pairing at all. "" is an explicitly
            # empty MSA, which is not the same as omitting the key - omitting it
            # would let a running data pipeline go and search for one.
            protein["pairedMsa"] = ""
        else:
            paired_name = f"chain_{cid}_paired.a3m"
            (out_dir / paired_name).write_text(clean_a3m(paired[i], seq, cid))
            protein["pairedMsaPath"] = paired_name

        if templates_mode == "search":
            # Leave the key out entirely so AF3's own data pipeline searches
            # pdb_seqres/mmcif_files for templates. Requires --run_data_pipeline=true,
            # otherwise AF3's featurisation rejects the null.
            pass
        elif templates_mode == "none":
            protein["templates"] = []
        else:
            templates = (chain_templates(templates_dir, seq)
                         if (not template_chains or cid in template_chains) else [])
            if templates:
                log.info("chain %s: %d template(s)", cid, len(templates))
            protein["templates"] = templates

        entries.append({"protein": protein})

    return {
        "name": sanitised_name(name),
        "modelSeeds": seeds,
        "sequences": entries,
        "dialect": "alphafold3",
        "version": 2,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="FASTA file (one record per chain)")
    parser.add_argument("--name", required=True, help="Job name (sanitised as AF3 does; output dir name)")
    parser.add_argument(
        "--a3m",
        nargs="+",
        default=None,
        help="Optional unpaired a3m(s): one per chain in record order (query-only if omitted)",
    )
    parser.add_argument(
        "--paired-a3m",
        dest="paired_a3m",
        nargs="+",
        default=None,
        help="Optional paired a3m(s) with AF3-parseable species headers: one per chain",
    )
    parser.add_argument("--seed", type=int, nargs="+", default=[1],
                        help="modelSeeds entries; several seeds run in one AF3 job (default: 1)")
    parser.add_argument("--templates-dir", dest="templates_dir", default=None, help="Matched-templates directory (bin/fold/match_templates.py)")
    parser.add_argument("--template-chains", dest="template_chains", nargs="+", default=None, help="Chain ID(s) that may be templated (default: all)")
    parser.add_argument(
        "--paired-msa-mode", dest="paired_msa_mode", choices=["path", "empty"], default="path",
        help="path: write a paired a3m per chain and reference it (default). "
             "empty: emit \"pairedMsa\": \"\" so AF3 does no cross-chain pairing, as Germinal does.",
    )
    parser.add_argument(
        "--templates-mode", dest="templates_mode", choices=["inline", "none", "search"], default="inline",
        help="inline: embed templates matched by --templates-dir (default). "
             "none: emit \"templates\": []. "
             "search: omit the key so AF3's own data pipeline searches for templates "
             "(needs --run_data_pipeline=true), as Germinal does.",
    )
    parser.add_argument("-o", "--output", required=True, help="Output JSON path (a3m files are written alongside)")
    args = parser.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    spec = make_af3_input(
        fasta_path=Path(args.fasta),
        name=args.name,
        out_dir=out_path.parent,
        seed=args.seed,
        unpaired_a3m_paths=[Path(p) for p in args.a3m] if args.a3m else None,
        paired_a3m_paths=[Path(p) for p in args.paired_a3m] if args.paired_a3m else None,
        templates_dir=Path(args.templates_dir) if args.templates_dir else None,
        template_chains=args.template_chains,
        paired_msa_mode=args.paired_msa_mode,
        templates_mode=args.templates_mode,
    )
    with open(out_path, "w") as f:
        json.dump(spec, f, indent=2)
    log.info("Wrote AlphaFold3 input JSON (%d chain(s), seeds %s, paired=%s, templates=%s) to %s",
             len(spec["sequences"]), ",".join(str(s) for s in spec["modelSeeds"]),
             args.paired_msa_mode, args.templates_mode, out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
