#!/usr/bin/env python3
# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "numpy",
# ]
# ///

"""
Score one AlphaFold2 run for fold.nf: for each structure the run publishes to
fold/predictions/, run ipSAE (bin/ipsae.py) and emit a normalized TSV row (via
bin/parse_fold_confidence.py). Rows for all kept models are concatenated to
stdout (header once), ready for collectFile into af2_fold_scores.tsv.

Which structures are "published" mirrors modules/fold/af2/alphafold2.nf's flat
fold/predictions saveAs rule:
  relax on (default): relaxed_model_*.pdb  (all 5 for keep=all; the best for keep=best)
  --af2_no_relax + keep=all:  unrelaxed_model_*.pdb
  --af2_no_relax + keep=best: ranked_0.pdb (-> ranking_debug order[0] model)
predictions_file = <pred_prefix><cif basename>, matching FoldNaming.af2Prefix.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

BIN = os.path.dirname(os.path.abspath(__file__))


def _suffix_from_model_pdb(name):
    # relaxed_model_1_multimer_v3_pred_0.pdb -> "1_multimer_v3_pred_0"
    base = os.path.basename(name)
    stem = base[: -len(".pdb")]
    for pre in ("relaxed_model_", "unrelaxed_model_"):
        if stem.startswith(pre):
            return stem[len(pre):]
    return None


def _published_structures(run_dir, keep_models, no_relax):
    """Return list of (model_index, cif_basename, pdb_path, suffix)."""
    out = []
    if not no_relax:
        pdbs = sorted(glob.glob(os.path.join(run_dir, "relaxed_model_*.pdb")))
        for p in pdbs:
            suf = _suffix_from_model_pdb(p)
            out.append((suf.split("_")[0], f"relaxed_model_{suf}.cif", p, suf))
    elif keep_models == "best":
        rd = os.path.join(run_dir, "ranking_debug.json")
        with open(rd) as f:
            order = json.load(f)["order"]
        top = order[0]  # e.g. model_1_multimer_v3_pred_0
        suf = top[len("model_"):] if top.startswith("model_") else top
        pdb = os.path.join(run_dir, f"ranked_0.pdb")
        out.append((suf.split("_")[0], "ranked_0.cif", pdb, suf))
    else:
        pdbs = sorted(glob.glob(os.path.join(run_dir, "unrelaxed_model_*.pdb")))
        for p in pdbs:
            suf = _suffix_from_model_pdb(p)
            out.append((suf.split("_")[0], f"unrelaxed_model_{suf}.cif", p, suf))
    return out


def _run_ipsae(pae, pdb, pae_cut, dist_cut):
    """Run ipsae.py in a temp dir; return path to the *_ipsae.tsv (or None)."""
    tmp = tempfile.mkdtemp()
    lp, lpdb = os.path.join(tmp, os.path.basename(pae)), os.path.join(tmp, os.path.basename(pdb))
    shutil.copy(pae, lp)
    shutil.copy(pdb, lpdb)
    r = subprocess.run(
        [sys.executable, os.path.join(BIN, "ipsae.py"), os.path.basename(pae),
         os.path.basename(pdb), str(pae_cut), str(dist_cut)],
        cwd=tmp, capture_output=True, text=True,
    )
    if r.returncode != 0:
        sys.stderr.write(f"ipsae.py failed for {pdb}:\n{r.stderr}\n")
        return None, tmp
    pae_s = str(int(pae_cut)).zfill(2)
    dist_s = str(int(dist_cut)).zfill(2)
    tsv = os.path.join(tmp, f"{os.path.basename(pdb)[:-4]}_{pae_s}_{dist_s}_ipsae.tsv")
    return (tsv if os.path.exists(tsv) else None), tmp


def main():
    p = argparse.ArgumentParser(description="Score one AF2 run into normalized TSV rows")
    p.add_argument("--run-dir", required=True, help="AF2 out/<id> directory")
    p.add_argument("--id", required=True)
    p.add_argument("--pred-prefix", required=True, help="FoldNaming.af2Prefix(meta)")
    p.add_argument("--keep-models", default="all", choices=("all", "best"))
    p.add_argument("--no-relax", action="store_true")
    p.add_argument("--pae-cutoff", type=float, default=10.0)
    p.add_argument("--dist-cutoff", type=float, default=10.0)
    args = p.parse_args()

    structures = _published_structures(args.run_dir, args.keep_models, args.no_relax)
    first = True
    for model_n, cif, pdb, suf in structures:
        pae = os.path.join(args.run_dir, f"pae_model_{suf}.json")
        pkl = os.path.join(args.run_dir, f"result_model_{suf}.pkl")
        ipsae_tsv, tmp = (None, None)
        if os.path.exists(pae) and os.path.exists(pdb):
            ipsae_tsv, tmp = _run_ipsae(pae, pdb, args.pae_cutoff, args.dist_cutoff)
        cmd = [
            sys.executable, os.path.join(BIN, "parse_fold_confidence.py"),
            "--tool", "af2", "--id", args.id, "--model", str(model_n),
            "--original-file", cif, "--predictions-file", f"{args.pred_prefix}{cif}",
            "--pkl", pkl,
        ]
        if ipsae_tsv:
            cmd += ["--ipsae-tsv", ipsae_tsv]
        if not first:
            cmd.append("--no-header")
        out = subprocess.run(cmd, capture_output=True, text=True)
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)
        if out.returncode != 0:
            sys.stderr.write(out.stderr)
            sys.exit(out.returncode)
        sys.stdout.write(out.stdout)
        first = False


if __name__ == "__main__":
    main()
