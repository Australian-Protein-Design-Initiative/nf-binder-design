#!/usr/bin/env python3
"""Call a locally-running OpenFold3 NIM's /predict endpoint and save the
predicted structure(s) plus their confidence scores. Inputs come in via
environment variables (set by openfold3_nim.nf) rather than command-line args,
to avoid multi-layer quoting problems between Nextflow, bash, and this script.
"""
import json
import os
import sys
import urllib.error
import urllib.request
from typing import Dict, List

AA_3_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

SCORE_COLUMNS = [
    "confidence_score",
    "complex_plddt_score",
    "complex_pde_score",
    "ptm_score",
    "iptm_score",
]


def chain_sequences(pdb_path: str) -> Dict[str, str]:
    """One-letter sequence per chain, in the order the chains first appear."""
    sequences: Dict[str, str] = {}
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21]
            # Key on residue number + insertion code so the many ATOM records
            # making up one residue only contribute a single letter.
            residue_key = (chain, line[22:27])
            if residue_key in seen:
                continue
            seen.add(residue_key)
            sequences[chain] = sequences.get(chain, "") + AA_3_1.get(line[17:20].strip(), "X")
    return sequences


def single_sequence_msa(chain: str, sequence: str) -> dict:
    """OpenFold3 NIM has no single-sequence mode - an MSA is mandatory for every
    protein chain - so hand it an alignment holding only the query itself. The
    baseline AF2 initial guess step this replaces also runs without a real MSA,
    so this keeps the two comparable, but it does mean the target chain is
    folded without evolutionary signal.
    """
    return {"main": {"a3m": {"alignment": f">{chain}\n{sequence}", "format": "a3m"}}}


input_pdb_path = os.environ["OF3_INPUT_PDB"]
output_dir = os.environ["OF3_OUTPUT_DIR"]
output_prefix = os.environ["OF3_OUTPUT_PREFIX"]
output_scores_path = os.environ["OF3_OUTPUT_SCORES"]
diffusion_samples = int(os.environ.get("OF3_DIFFUSION_SAMPLES", "1"))

sequences = chain_sequences(input_pdb_path)
if not sequences:
    sys.exit(f"No protein chains found in {input_pdb_path}")

molecules = [
    {
        "type": "protein",
        "id": chain,
        "sequence": sequence,
        "msa": single_sequence_msa(chain, sequence),
    }
    for chain, sequence in sequences.items()
]

payload = {
    "inputs": [
        {
            "input_id": output_prefix,
            "molecules": molecules,
            "diffusion_samples": diffusion_samples,
            "output_format": "pdb",
        }
    ]
}

req = urllib.request.Request(
    "http://localhost:8000/biology/openfold/openfold3/predict",
    data=json.dumps(payload).encode("utf-8"),
    headers={"Content-Type": "application/json"},
    method="POST",
)
try:
    with urllib.request.urlopen(req) as resp:
        result = json.load(resp)
except urllib.error.HTTPError as e:
    # The NIM puts its validation detail in the response body, which urllib
    # otherwise hides behind a bare "HTTP Error 422: Unprocessable Entity".
    sys.stderr.write(f"OpenFold3 NIM returned {e.code}: {e.read().decode(errors='replace')}\n")
    raise

# Ranked best-first by the NIM.
structures = result["outputs"][0]["structures_with_scores"]

rows: List[str] = []
for i, entry in enumerate(structures):
    name = output_prefix if len(structures) == 1 else f"{output_prefix}_s{i}"
    with open(os.path.join(output_dir, f"{name}.pdb"), "w") as f:
        f.write(entry["structure"])
    rows.append("\t".join([f"{name}.pdb"] + [str(entry.get(c, "")) for c in SCORE_COLUMNS]))

with open(output_scores_path, "w") as f:
    f.write("\t".join(["filename"] + SCORE_COLUMNS) + "\n")
    f.write("\n".join(rows) + "\n")

print(f"Chains folded: {', '.join(f'{c}({len(s)})' for c, s in sequences.items())}")
print(f"Wrote {len(structures)} structure(s) for {output_prefix}")
