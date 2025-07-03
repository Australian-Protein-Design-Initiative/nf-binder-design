#!/usr/bin/env python
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pandas",
#     "mdanalysis",
# ]
# ///
# NOTE: If running via `uv run`, the dependencies for the filters.d/ plugins must also be here.

import argparse
import sys
import pandas as pd
from pathlib import Path
import shutil
import re
import importlib.util
import logging


def load_plugins(plugins_dir: Path) -> dict:
    """Dynamically loads filter plugins from a directory."""
    plugins = {}
    if not plugins_dir.is_dir():
        logging.warning(f"Plugins directory '{plugins_dir}' not found.")
        return plugins

    for plugin_file in plugins_dir.glob("*.py"):
        try:
            spec = importlib.util.spec_from_file_location(plugin_file.stem, plugin_file)
            if spec and spec.loader:
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                if hasattr(module, "register_metrics"):
                    for metric in module.register_metrics():
                        if metric in plugins:
                            logging.warning(
                                f"Metric '{metric}' from '{plugin_file.name}' conflicts with a previously loaded plugin. Ignoring."
                            )
                        else:
                            plugins[metric] = module
        except Exception as e:
            logging.warning(f"Could not load plugin '{plugin_file.name}': {e}")
    return plugins


def parse_filter_expression(expression: str) -> tuple | None:
    """Parses a filter expression like 'rg<15' into (metric, operator, value)."""
    match = re.match(r"([a-zA-Z0-9_]+)\s*(>=|<=|==|!=|>|<|=)\s*(.+)", expression)
    if match:
        metric, op, value_str = match.groups()
        return metric, op, value_str
    return None


def evaluate_filter(
    df: pd.DataFrame, metric: str, op: str, value_str: str
) -> pd.Series:
    """Evaluates a single filter expression against the DataFrame."""
    if metric not in df.columns:
        raise ValueError(f"Metric '{metric}' not found in data for filtering.")

    col = df[metric]

    # Type coercion logic
    try:
        # Try converting to number first
        value = float(value_str)
        col = pd.to_numeric(col, errors="coerce")
    except ValueError:
        # Fallback to string/boolean
        low_val_str = value_str.lower()
        if low_val_str in ["true", "yes"]:
            value = True
            col = col.astype(bool)
        elif low_val_str in ["false", "no"]:
            value = False
            col = col.astype(bool)
        else:
            # Treat as string
            value = value_str
            col = col.astype(str)

    op_map = {
        ">": lambda a, b: a > b,
        "<": lambda a, b: a < b,
        ">=": lambda a, b: a >= b,
        "<=": lambda a, b: a <= b,
        "==": lambda a, b: a == b,
        "=": lambda a, b: a == b,  # Treat as equality
        "!=": lambda a, b: a != b,
    }
    return op_map[op](col, value)


def main():
    script_dir = Path(__file__).parent.resolve()
    parser = argparse.ArgumentParser(
        description="Filter PDB files based on scores and calculated metrics.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("pdb_files", nargs="*", help="List of PDB files to process.")
    parser.add_argument(
        "--pdb-list", help="Optional path to a file containing one PDB path per line."
    )
    parser.add_argument(
        "--scores",
        action="append",
        dest="scores_files",
        help="Path to a scores file (TSV/CSV). The column used to map to PDB filenames can be specified after a colon (e.g., 'scores.tsv:pdb_name'). Defaults to a 'design_id' or 'pdb' column. Can be specified multiple times.",
    )
    parser.add_argument(
        "--filter",
        action="append",
        dest="filters",
        help='A filter expression (e.g., "rg<30", "plddt>85"). Can be specified multiple times.',
    )
    parser.add_argument(
        "--filter-plugins-dir",
        default=str(script_dir / "filters.d"),
        help="Path to the directory containing filter plugin modules. Defaults to the 'filters.d' directory next to the script.",
    )
    parser.add_argument(
        "--collect-in",
        help="Optional. If provided, copy PDBs to '{path}/accepted' and '{path}/rejected' subdirectories.",
    )
    parser.add_argument(
        "--output",
        default="-",
        help="Path to write the output TSV. Defaults to stdout (-).",
    )
    parser.add_argument(
        "--binder-chains",
        default="A",
        help="Comma-separated list of binder chain IDs. Defaults to 'A'.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, stream=sys.stderr, format="%(levelname)s: %(message)s"
    )

    # --- 1. Input Gathering ---
    pdb_paths = []
    if args.pdb_files:
        pdb_paths.extend(Path(p) for p in args.pdb_files)
    if args.pdb_list:
        with open(args.pdb_list) as f:
            pdb_paths.extend(Path(line.strip()) for line in f if line.strip())

    if not pdb_paths:
        parser.error(
            "No PDB files provided either as positional arguments or via --pdb-list."
        )

    df = pd.DataFrame({"pdb_path": pdb_paths})
    df["design_id"] = df["pdb_path"].apply(lambda p: p.stem)
    df.set_index("design_id", inplace=True)

    logging.info(f"Processing {len(df)} PDB files.")

    # --- 2. Score Loading ---
    if args.scores_files:
        logging.info(f"Loading scores from: {args.scores_files}")
        for score_spec in args.scores_files:
            parts = score_spec.split(":", 1)
            path = parts[0]
            join_col = parts[1] if len(parts) > 1 else "design_id"

            try:
                scores_df = pd.read_csv(path, sep=None, engine="python")

                # Find the correct column to join on, case-insensitive
                id_col_in_scores = None
                for col in scores_df.columns:
                    if col.lower() == join_col.lower():
                        id_col_in_scores = col
                        break

                if not id_col_in_scores:
                    raise ValueError(
                        f"Join column '{join_col}' not found in scores file '{path}'"
                    )

                scores_df.rename(columns={id_col_in_scores: "design_id"}, inplace=True)
                scores_df.set_index("design_id", inplace=True)
                df = df.join(scores_df, how="left")

            except Exception as e:
                parser.error(f"Could not process scores file '{path}': {e}")

    # --- 3. Filter Parsing & 4. Plugin Execution ---
    plugins = load_plugins(Path(args.filter_plugins_dir))
    required_metrics = set()
    parsed_filters = []

    if args.filters:
        logging.info(f"Applying filters: {args.filters}")
        for f_expr in args.filters:
            parsed = parse_filter_expression(f_expr)
            if not parsed:
                parser.error(f"Invalid filter expression syntax: '{f_expr}'")
            metric, _, _ = parsed
            required_metrics.add(metric)
            parsed_filters.append(parsed)

    missing_metrics = required_metrics - set(df.columns)
    metrics_to_calculate = {}  # map of plugin_module -> list_of_metrics

    for metric in missing_metrics:
        plugin_module = plugins.get(metric)
        if not plugin_module:
            parser.error(f"No plugin found to calculate required metric: '{metric}'")
        if plugin_module not in metrics_to_calculate:
            metrics_to_calculate[plugin_module] = []
        metrics_to_calculate[plugin_module].append(metric)

    if metrics_to_calculate:
        logging.info(f"Calculating missing metrics: {list(missing_metrics)}...")
        all_pdb_paths = df["pdb_path"].tolist()
        binder_chains_list = args.binder_chains.split(",")

        for plugin_module, metrics in metrics_to_calculate.items():
            logging.info(
                f"  -> using {Path(plugin_module.__file__).name} for {metrics}"
            )
            calculated_df = plugin_module.calculate_metrics(
                all_pdb_paths, binder_chains=binder_chains_list
            )
            df = df.join(calculated_df, how="left")

    # --- 5. Evaluation ---
    all_pass_columns = []
    if parsed_filters:
        for metric, op, value in parsed_filters:
            pass_col_name = f"{metric}_pass"
            all_pass_columns.append(pass_col_name)
            try:
                df[pass_col_name] = evaluate_filter(df, metric, op, value)
            except ValueError as e:
                parser.error(str(e))

        df["pass"] = df[all_pass_columns].all(axis=1)  # type: ignore
    else:
        # If no filters, everything passes
        df["pass"] = True

    # --- 6. Output Generation ---
    output_handle = open(args.output, "w") if args.output != "-" else sys.stdout
    df.to_csv(output_handle, sep="\t")
    if args.output != "-":
        output_handle.close()
        logging.info(f"Wrote output scores to {args.output}")

    if args.collect_in:
        accepted_dir = Path(args.collect_in) / "accepted"
        rejected_dir = Path(args.collect_in) / "rejected"
        accepted_dir.mkdir(parents=True, exist_ok=True)
        rejected_dir.mkdir(parents=True, exist_ok=True)

        logging.info(f"Collecting PDBs into {args.collect_in}")

        # Resetting index to make 'pdb_path' a regular column for safer iteration
        df_for_copy = df.reset_index()

        for _, row in df_for_copy.iterrows():
            pdb_path = Path(str(row["pdb_path"]))
            dest_dir = accepted_dir if row["pass"] else rejected_dir
            shutil.copy(pdb_path, dest_dir / pdb_path.name)


if __name__ == "__main__":
    main()
