#!/usr/bin/env python3
"""Call a locally-running ProteinMPNN NIM's /predict endpoint and save the
resulting sequences. Inputs come in via environment variables (set by
dl_binder_design_nim.nf) rather than command-line args, to avoid multi-layer
quoting problems between Nextflow, bash, and this script.
"""
import json
import os
import urllib.request

backbone_pdb_path = os.environ["PMPNN_BACKBONE_PDB"]
design_chain = os.environ["PMPNN_DESIGN_CHAIN"]
sampling_temp = float(os.environ["PMPNN_SAMPLING_TEMP"])
output_fasta_path = os.environ["PMPNN_OUTPUT_FASTA"]

with open(backbone_pdb_path) as f:
    lines = [l.rstrip("\n") for l in f if l.startswith(("ATOM", "TER", "END"))]
pdb_text = "\n".join(lines)

payload = {
    "input_pdb": pdb_text,
    "input_pdb_chains": [design_chain],
    "num_seq_per_target": 1,
    "sampling_temp": [sampling_temp],
    "use_soluble_model": False,
    "ca_only": False,
}

req = urllib.request.Request(
    "http://localhost:8000/biology/ipd/proteinmpnn/predict",
    data=json.dumps(payload).encode("utf-8"),
    headers={"Content-Type": "application/json"},
    method="POST",
)
with urllib.request.urlopen(req) as resp:
    result = json.load(resp)

with open(output_fasta_path, "w") as f:
    f.write(result["mfasta"])

print("Scores:", result.get("scores"))
