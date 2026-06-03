#!/usr/bin/env python
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "biopython>=1.75",
# ]
# ///

import argparse
import asyncio
import copy
import logging
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, TextIO, Tuple

from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

logger = logging.getLogger(__name__)

ResidueKey = Tuple[str, int, str]

warnings.filterwarnings("ignore", category=PDBConstructionWarning)
warnings.filterwarnings("ignore", message="WARNING: Unrecognized atom type")
warnings.filterwarnings("ignore", message="WARNING: Negative sasa result!")

# Maximum residue accessible surface area (theoretical) from Tien et al., 2013
# https://doi.org/10.1371/journal.pone.0080635
TIEN_2023_THEORETICAL: dict[str, float] = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLU": 223.0,
    "GLN": 225.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0,
}


@dataclass(frozen=True)
class SiteSpec:
    name: str
    residue_tokens: Tuple[str, ...]


@dataclass
class ResidueDelta:
    key: ResidueKey
    default_label: str
    column_label: str
    resname: str
    delta_angstrom: float
    delta_percent: Optional[float]


@dataclass
class StructureResult:
    pdb_name: str
    target_chains: str
    binder_chains: str
    residue_deltas: Dict[ResidueKey, ResidueDelta]
    total_angstrom: float
    site_angstrom: Dict[str, float]
    site_percent: Dict[str, Optional[float]]


def parse_chain_list(chains_arg: str) -> List[str]:
    return [c.strip() for c in chains_arg.replace(" ", ",").split(",") if c.strip()]


def parse_site_arg(site_str: str) -> SiteSpec:
    if "=" not in site_str:
        raise ValueError(
            f"Invalid --site value '{site_str}': expected NAME=RES1,RES2,..."
        )
    name, residues = site_str.split("=", 1)
    name = name.strip()
    if not name:
        raise ValueError(f"Invalid --site value '{site_str}': empty site name")
    tokens = tuple(r.strip() for r in residues.split(",") if r.strip())
    if not tokens:
        raise ValueError(f"Invalid --site value '{site_str}': no residues listed")
    return SiteSpec(name=name, residue_tokens=tokens)


def structure_chain_ids(structure) -> Set[str]:
    return {chain.get_id() for chain in structure[0]}


def parse_residue_token(token: str, chain_ids: Set[str]) -> ResidueKey:
    for chain_id in sorted(chain_ids, key=len, reverse=True):
        if not token.startswith(chain_id):
            continue
        rest = token[len(chain_id) :]
        try:
            resnum = int(rest)
        except ValueError:
            continue
        return (chain_id, resnum, " ")
    raise ValueError(
        f"Cannot parse residue token '{token}' against chain IDs: {sorted(chain_ids)}"
    )


def default_residue_label(chain_id: str, resnum: int) -> str:
    return f"{chain_id}{resnum}"


def residue_key_sort_key(key: ResidueKey) -> Tuple[str, int, str]:
    return (key[0], key[1], key[2])


def target_chain_ids_in_structure(structure, binder_chains: Set[str]) -> List[str]:
    chain_ids: List[str] = []
    for chain in structure[0]:
        chain_id = chain.get_id()
        if chain_id in binder_chains:
            continue
        if any(is_aa(residue, standard=True) for residue in chain.get_residues()):
            chain_ids.append(chain_id)
    return sorted(chain_ids)


def binder_chain_ids_in_structure(
    structure, binder_chains: Set[str]
) -> List[str]:
    present = structure_chain_ids(structure)
    return sorted(chain_id for chain_id in binder_chains if chain_id in present)


def iter_target_residues(structure, binder_chains: Set[str]) -> Iterable[Tuple[ResidueKey, str]]:
    for chain in structure[0]:
        chain_id = chain.get_id()
        if chain_id in binder_chains:
            continue
        for residue in chain.get_residues():
            if not is_aa(residue, standard=True):
                continue
            res_id = residue.get_id()
            key = (chain_id, res_id[1], res_id[2])
            yield key, residue.get_resname().upper()


def compute_residue_sasa_map(structure, sasa_calculator: ShrakeRupley) -> Dict[ResidueKey, float]:
    sasa_calculator.compute(structure, level="R")
    sasa_by_key: Dict[ResidueKey, float] = {}
    for chain in structure[0]:
        chain_id = chain.get_id()
        for residue in chain.get_residues():
            if not is_aa(residue, standard=True):
                continue
            res_id = residue.get_id()
            key = (chain_id, res_id[1], res_id[2])
            if hasattr(residue, "sasa"):
                sasa_by_key[key] = float(getattr(residue, "sasa"))
            else:
                sasa_by_key[key] = 0.0
    return sasa_by_key


def make_apo_structure(structure, binder_chains: Set[str]):
    apo = copy.deepcopy(structure)
    model = apo[0]
    for chain_id in binder_chains:
        if chain_id in model:
            model.detach_child(chain_id)
    return apo


def delta_percent_for(resname: str, delta_angstrom: float) -> Optional[float]:
    standard_sasa = TIEN_2023_THEORETICAL.get(resname)
    if not standard_sasa or standard_sasa <= 0:
        return None
    return round((delta_angstrom / standard_sasa) * 100, 2)


def site_burial_percent(
    keys: Sequence[ResidueKey],
    residue_deltas: Dict[ResidueKey, ResidueDelta],
) -> Optional[float]:
    """Buried SASA for the site as a fraction of combined max residue SASA."""
    angstrom_sum = sum(residue_deltas[key].delta_angstrom for key in keys)
    theoretical_sum = 0.0
    for key in keys:
        max_sasa = TIEN_2023_THEORETICAL.get(residue_deltas[key].resname)
        if max_sasa and max_sasa > 0:
            theoretical_sum += max_sasa
    if theoretical_sum <= 0:
        return None
    return round((angstrom_sum / theoretical_sum) * 100, 2)


def assign_column_labels(
    keys: Sequence[ResidueKey],
    resnames_by_key: Dict[ResidueKey, str],
    use_residue_names: bool,
) -> Dict[ResidueKey, str]:
    labels: Dict[ResidueKey, str] = {}
    if not use_residue_names:
        for key in keys:
            labels[key] = default_residue_label(key[0], key[1])
        return labels

    proposed: Dict[ResidueKey, str] = {}
    for key in keys:
        aa1 = seq1(resnames_by_key[key], custom_map={"UNK": "X"})
        proposed[key] = f"{aa1}{key[1]}"

    seen: Dict[str, ResidueKey] = {}
    duplicates: Set[str] = set()
    for key, label in proposed.items():
        if label in seen and seen[label] != key:
            duplicates.add(label)
        else:
            seen.setdefault(label, key)

    for key in keys:
        label = proposed[key]
        if label in duplicates:
            fallback = default_residue_label(key[0], key[1])
            logger.warning(
                "Duplicate column label '%s'; using '%s' instead",
                label,
                fallback,
            )
            labels[key] = fallback
        else:
            labels[key] = label
    return labels


def collect_pdb_files(pdb_files_arg: Sequence[str], pdbdir: Optional[str]) -> List[Path]:
    if pdbdir:
        input_dir = Path(pdbdir)
        if not input_dir.is_dir():
            logger.error(
                "Input directory '%s' does not exist or is not a directory.",
                pdbdir,
            )
            sys.exit(1)
        pdb_files = list(input_dir.rglob("*.pdb"))
        if not pdb_files:
            logger.warning("No PDB files found in directory '%s'.", pdbdir)
            sys.exit(0)
        return [p.resolve() for p in pdb_files]

    pdb_files: List[Path] = []
    for pdb_path_str in pdb_files_arg:
        pdb_path = Path(pdb_path_str)
        if pdb_path.is_file() and pdb_path.suffix.lower() == ".pdb":
            pdb_files.append(pdb_path.resolve())
        elif pdb_path.is_file():
            logger.warning(
                "File '%s' is not a PDB file (expected .pdb extension). Skipping.",
                pdb_path_str,
            )
        else:
            logger.warning(
                "File '%s' does not exist or is not a file. Skipping.",
                pdb_path_str,
            )

    if not pdb_files:
        logger.warning("No valid PDB files found in the provided arguments.")
        sys.exit(0)
    return pdb_files


def build_column_labels_from_rows(
    rows: Sequence[StructureResult],
    use_residue_names: bool,
) -> Dict[ResidueKey, str]:
    union_keys: Set[ResidueKey] = set()
    resnames_by_key: Dict[ResidueKey, str] = {}
    for row in rows:
        for key, delta in row.residue_deltas.items():
            union_keys.add(key)
            if key not in resnames_by_key:
                resnames_by_key[key] = delta.resname
    logger.info(
        "Column layout: %d target residue position(s) from %d structure(s).",
        len(union_keys),
        len(rows),
    )
    if not union_keys:
        return {}
    sorted_keys = sorted(union_keys, key=residue_key_sort_key)
    return assign_column_labels(sorted_keys, resnames_by_key, use_residue_names)


def resolve_site_keys(
    sites: Sequence[SiteSpec],
    chain_ids: Set[str],
    binder_chains: Set[str],
    pdb_name: str,
) -> Dict[str, List[ResidueKey]]:
    resolved: Dict[str, List[ResidueKey]] = {}
    for site in sites:
        keys: List[ResidueKey] = []
        for token in site.residue_tokens:
            try:
                key = parse_residue_token(token, chain_ids)
            except ValueError as exc:
                raise ValueError(f"{pdb_name}: {exc}") from exc
            if key[0] in binder_chains:
                raise ValueError(
                    f"{pdb_name}: site '{site.name}' residue '{token}' is on "
                    f"binder chain '{key[0]}'"
                )
            keys.append(key)
        resolved[site.name] = keys
    return resolved


def analyse_structure(
    pdb_file: Path,
    binder_chains: Set[str],
    sites: Sequence[SiteSpec],
    sasa_calculator: ShrakeRupley,
    pdb_parser: PDBParser,
) -> Optional[StructureResult]:
    pdb_name = pdb_file.name
    try:
        structure = pdb_parser.get_structure(pdb_name, str(pdb_file))
    except Exception as exc:
        logger.error("Could not parse PDB file %s: %s", pdb_file, exc)
        return None

    chain_ids = structure_chain_ids(structure)
    missing_binder = binder_chains - chain_ids
    if missing_binder:
        logger.warning(
            "%s: binder chain(s) not found and ignored: %s",
            pdb_name,
            ", ".join(sorted(missing_binder)),
        )

    site_keys_by_name = resolve_site_keys(sites, chain_ids, binder_chains, pdb_name)

    complex_sasa = compute_residue_sasa_map(structure, sasa_calculator)
    apo_structure = make_apo_structure(structure, binder_chains)
    apo_sasa = compute_residue_sasa_map(apo_structure, sasa_calculator)

    residue_deltas: Dict[ResidueKey, ResidueDelta] = {}
    total_angstrom = 0.0

    for key, resname in iter_target_residues(structure, binder_chains):
        default_label = default_residue_label(key[0], key[1])
        apo_value = apo_sasa.get(key, 0.0)
        complex_value = complex_sasa.get(key, 0.0)
        delta_angstrom = round(apo_value - complex_value, 2)
        delta_percent = delta_percent_for(resname, delta_angstrom)
        total_angstrom += delta_angstrom
        residue_deltas[key] = ResidueDelta(
            key=key,
            default_label=default_label,
            column_label=default_label,
            resname=resname,
            delta_angstrom=delta_angstrom,
            delta_percent=delta_percent,
        )

    target_residues_in_structure = {key for key, _ in iter_target_residues(structure, binder_chains)}
    site_angstrom: Dict[str, float] = {}
    site_percent: Dict[str, float] = {}
    for site in sites:
        keys = site_keys_by_name[site.name]
        for key in keys:
            if key not in target_residues_in_structure:
                raise ValueError(
                    f"{pdb_name}: site '{site.name}' residue "
                    f"'{default_residue_label(key[0], key[1])}' not found in structure"
                )
        angstrom_sum = sum(residue_deltas[key].delta_angstrom for key in keys)
        site_angstrom[site.name] = round(angstrom_sum, 2)
        site_percent[site.name] = site_burial_percent(keys, residue_deltas)

    return StructureResult(
        pdb_name=pdb_name,
        target_chains=",".join(target_chain_ids_in_structure(structure, binder_chains)),
        binder_chains=",".join(
            binder_chain_ids_in_structure(structure, binder_chains)
        ),
        residue_deltas=residue_deltas,
        total_angstrom=round(total_angstrom, 2),
        site_angstrom=site_angstrom,
        site_percent=site_percent,
    )


def filter_residue_columns(
    column_labels: Dict[ResidueKey, str],
    rows: Sequence[StructureResult],
    min_change_percent: float,
) -> List[ResidueKey]:
    if min_change_percent <= 0:
        return sorted(column_labels.keys(), key=residue_key_sort_key)

    kept: List[ResidueKey] = []
    dropped = 0
    for key in sorted(column_labels.keys(), key=residue_key_sort_key):
        label = column_labels[key]
        values: List[float] = []
        for row in rows:
            delta = row.residue_deltas.get(key)
            if delta is None or delta.delta_percent is None:
                values.append(float("-inf"))
            else:
                values.append(delta.delta_percent)
        if any(value >= min_change_percent for value in values):
            kept.append(key)
        else:
            dropped += 1
            logger.debug(
                "Dropping columns for %s (all delta_percent < %s)",
                label,
                min_change_percent,
            )

    if min_change_percent > 0:
        logger.info(
            "Dropped %d per-residue column pair(s) below --min-change-percent %s",
            dropped,
            min_change_percent,
        )
    return kept


def format_optional_float(value: Optional[float]) -> str:
    if value is None:
        return ""
    return f"{value:.2f}"


def sort_rows_by_first_site_percent(
    rows: Sequence[StructureResult],
    site_names: Sequence[str],
) -> List[StructureResult]:
    if not site_names:
        return list(rows)
    first_site = site_names[0]

    def sort_key(row: StructureResult) -> float:
        value = row.site_percent.get(first_site)
        if value is None:
            return float("-inf")
        return value

    return sorted(rows, key=sort_key, reverse=True)


def write_output(
    output_path: str,
    rows: Sequence[StructureResult],
    residue_columns: Sequence[ResidueKey],
    column_labels: Dict[ResidueKey, str],
    site_names: Sequence[str],
) -> None:
    header_parts = ["pdb_file", "target_chains", "binder_chains"]
    for key in residue_columns:
        label = column_labels[key]
        header_parts.extend([f"{label}_angstrom", f"{label}_percent"])
    header_parts.append("total_angstrom")
    for site_name in site_names:
        header_parts.extend(
            [f"{site_name}_site_angstrom", f"{site_name}_site_percent"]
        )
    header = "\t".join(header_parts)

    output_file: Optional[TextIO] = None
    if output_path != "-":
        output_file_path = Path(output_path)
        output_file_path.parent.mkdir(parents=True, exist_ok=True)
        output_file = open(output_file_path, "w")
    try:
        if output_path == "-":
            print(header)
        else:
            assert output_file is not None
            output_file.write(header + "\n")

        for row in rows:
            parts = [row.pdb_name, row.target_chains, row.binder_chains]
            for key in residue_columns:
                delta = row.residue_deltas.get(key)
                if delta is None:
                    parts.extend(["", ""])
                else:
                    parts.extend(
                        [
                            format_optional_float(delta.delta_angstrom),
                            format_optional_float(delta.delta_percent),
                        ]
                    )
            parts.append(format_optional_float(row.total_angstrom))
            for site_name in site_names:
                parts.extend(
                    [
                        format_optional_float(row.site_angstrom.get(site_name)),
                        format_optional_float(row.site_percent.get(site_name)),
                    ]
                )
            line = "\t".join(parts)
            if output_path == "-":
                print(line)
            else:
                assert output_file is not None
                output_file.write(line + "\n")
    finally:
        if output_file is not None:
            output_file.close()


def analyse_pdb_worker(
    pdb_path_str: str,
    binder_chains_tuple: Tuple[str, ...],
    sites_tuple: Tuple[SiteSpec, ...],
    probe_radius: float,
    n_points: int,
) -> Optional[StructureResult]:
    pdb_file = Path(pdb_path_str)
    binder_chains = set(binder_chains_tuple)
    pdb_parser = PDBParser(QUIET=True)
    sasa_calculator = ShrakeRupley(
        probe_radius=probe_radius,
        n_points=n_points,
    )
    return analyse_structure(
        pdb_file,
        binder_chains,
        sites_tuple,
        sasa_calculator,
        pdb_parser,
    )


async def analyse_structures_parallel(
    pdb_files: Sequence[Path],
    binder_chains: Set[str],
    sites: Sequence[SiteSpec],
    probe_radius: float,
    n_points: int,
    threads: int,
    sasa_calculator: ShrakeRupley,
    pdb_parser: PDBParser,
) -> List[StructureResult]:
    binder_tuple = tuple(sorted(binder_chains))
    sites_tuple = tuple(sites)

    if threads <= 1 or len(pdb_files) <= 1:
        sequential_rows: List[StructureResult] = []
        for pdb_file in pdb_files:
            logger.info("Processing %s...", pdb_file.name)
            result = analyse_structure(
                pdb_file,
                binder_chains,
                sites,
                sasa_calculator,
                pdb_parser,
            )
            if result is not None:
                sequential_rows.append(result)
        return sequential_rows

    loop = asyncio.get_running_loop()
    with ProcessPoolExecutor(max_workers=threads) as executor:
        for pdb_file in pdb_files:
            logger.info("Processing %s...", pdb_file.name)
        futures = [
            loop.run_in_executor(
                executor,
                analyse_pdb_worker,
                str(pdb_file.resolve()),
                binder_tuple,
                sites_tuple,
                probe_radius,
                n_points,
            )
            for pdb_file in pdb_files
        ]
        results = await asyncio.gather(*futures)

    rows: List[StructureResult] = []
    for result in results:
        if result is not None:
            rows.append(result)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate per-residue delta SASA on target chains when a binder is removed "
            "from a protein complex."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "pdb_files",
        nargs="*",
        type=str,
        help="PDB files to process. Can use shell globs like *.pdb",
    )
    parser.add_argument(
        "--pdbdir",
        type=str,
        help="Directory containing PDB files (scanned recursively). If provided, positional PDB files are ignored.",
    )
    parser.add_argument(
        "--binder-chains",
        required=True,
        type=str,
        help="Comma-separated binder chain ID(s) to remove for the apo target state.",
    )
    parser.add_argument(
        "--site",
        action="append",
        default=[],
        metavar="NAME=RES,...",
        help="Repeatable site subset; sums reported as NAME_site_angstrom and NAME_site_percent.",
    )
    parser.add_argument(
        "--use-residue-names",
        action="store_true",
        help="Use one-letter amino acid residue names in column labels (e.g. K41_percent).",
    )
    parser.add_argument(
        "--min-change-percent",
        type=float,
        default=0.0,
        help=(
            "Drop per-residue columns when every structure's delta_percent is strictly "
            "below this value. Use 0 to keep all per-residue columns."
        ),
    )
    parser.add_argument(
        "--output",
        type=str,
        default="-",
        help="Path to the output TSV file. Use '-' for stdout.",
    )
    parser.add_argument(
        "--probe-radius",
        type=float,
        default=1.4,
        help="Radius of the probe for SASA calculation (in Angstroms).",
    )
    parser.add_argument(
        "--n-points",
        type=int,
        default=100,
        help="Number of points for SASA sphere resolution (higher is more precise but slower).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Worker processes for parallel PDB analysis (default: CPU count).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging output.",
    )
    args = parser.parse_args()

    log_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(
        stream=sys.stderr,
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    binder_chains = set(parse_chain_list(args.binder_chains))
    if not binder_chains:
        logger.error("--binder-chains must list at least one chain ID.")
        sys.exit(1)

    try:
        sites = [parse_site_arg(site_str) for site_str in args.site]
    except ValueError as exc:
        logger.error("%s", exc)
        sys.exit(1)

    pdb_files = collect_pdb_files(args.pdb_files, args.pdbdir)
    logger.info("Found %d PDB file(s) to process.", len(pdb_files))

    pdb_parser = PDBParser(QUIET=True)
    sasa_calculator = ShrakeRupley(
        probe_radius=args.probe_radius,
        n_points=args.n_points,
    )

    threads = args.threads if args.threads is not None else (os.cpu_count() or 1)
    if threads < 1:
        logger.error("--threads must be at least 1.")
        sys.exit(1)

    try:
        rows = asyncio.run(
            analyse_structures_parallel(
                pdb_files,
                binder_chains,
                sites,
                args.probe_radius,
                args.n_points,
                threads,
                sasa_calculator,
                pdb_parser,
            )
        )
    except ValueError as exc:
        logger.error("%s", exc)
        sys.exit(1)

    if not rows:
        logger.error("No structures were processed successfully.")
        sys.exit(1)

    column_labels = build_column_labels_from_rows(rows, args.use_residue_names)
    if not column_labels:
        logger.error("No target residues found in the supplied PDB file(s).")
        sys.exit(1)

    residue_columns = filter_residue_columns(
        column_labels,
        rows,
        args.min_change_percent,
    )
    site_names = [site.name for site in sites]
    rows = sort_rows_by_first_site_percent(rows, site_names)
    write_output(args.output, rows, residue_columns, column_labels, site_names)

    if args.output != "-":
        logger.info("Output written to %s", Path(args.output).resolve())


if __name__ == "__main__":
    main()
