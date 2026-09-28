#!/usr/bin/env python3
"""Call a locally-running RFdiffusion NIM's /generate endpoint and save the
resulting PDB. Inputs come in via environment variables (set by
rfdiffusion_nim.nf) rather than command-line args, to avoid multi-layer
quoting problems between Nextflow, bash, and this script.
"""
import json
import os
import urllib.request

input_pdb_path = os.environ["RFD_INPUT_PDB"]
contigs = os.environ["RFD_CONTIGS"]
hotspot_res = os.environ["RFD_HOTSPOT_RES"]
output_pdb_path = os.environ["RFD_OUTPUT_PDB"]

with open(input_pdb_path) as f:
    lines = [l.rstrip("\n") for l in f if l.startswith(("ATOM", "TER", "END"))]
pdb_text = "\n".join(lines)

# hotspot_res arrives as a bracketed string like "[A56]" - the NIM wants a
# plain JSON array of residue IDs, e.g. ["A56"].
hotspot_list = [x.strip() for x in hotspot_res.strip("[]").split(",") if x.strip()]

payload = {
    "input_pdb": pdb_text,
    "contigs": contigs,
    "hotspot_res": hotspot_list,
}

req = urllib.request.Request(
    "http://localhost:8000/biology/ipd/rfdiffusion/generate",
    data=json.dumps(payload).encode("utf-8"),
    headers={"Content-Type": "application/json"},
    method="POST",
)
with urllib.request.urlopen(req) as resp:
    result = json.load(resp)

with open(output_pdb_path, "w") as f:
    f.write(result["output_pdb"])

print("Server-side elapsed ms:", result.get("elapsed_ms"))
