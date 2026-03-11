#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Query the ColabFold MMseqs2 API server for MSA generation.

Implements the same protocol as ColabFold's run_mmseqs2():
https://github.com/sokrypton/ColabFold/blob/main/colabfold/colabfold.py

Reads a FASTA file, submits to the API, polls until complete, downloads the
tar.gz result, extracts and merges uniref + env a3m files, writes one .a3m per
query sequence.
"""

import argparse
import json
import logging
import random
import ssl
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

DEFAULT_HOST = "https://api.colabfold.com"
USER_AGENT = "nf-binder-design/colabfold_remote_msa (Nextflow pipeline)"


def read_fasta_sequences(path: Path) -> List[Tuple[str, str]]:
    """Return list of (header_line, sequence) from a FASTA file."""
    pairs: List[Tuple[str, str]] = []
    current_header: str | None = None
    current_seq: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    pairs.append((current_header, "".join(current_seq).replace(" ", "")))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            pairs.append((current_header, "".join(current_seq).replace(" ", "")))
    return pairs


def _urlopen(req: urllib.request.Request, timeout: int = 30) -> urllib.request.addinfourl:
    return urllib.request.urlopen(req, timeout=timeout, context=ssl.create_default_context())


def submit(seqs: List[str], mode: str, host_url: str, headers: dict, N: int = 101) -> dict:
    query = ""
    for i, seq in enumerate(seqs):
        query += f">{N + i}\n{seq}\n"
    data = urllib.parse.urlencode({"q": query, "mode": mode}).encode()
    req = urllib.request.Request(
        f"{host_url}/ticket/msa",
        data=data,
        headers={**headers, "Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    for _ in range(5):
        try:
            with _urlopen(req, timeout=30) as res:
                text = res.read().decode()
            return json.loads(text)
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
            log.warning("Error submitting: %s. Retrying...", e)
            time.sleep(5 + random.randint(0, 5))
        except json.JSONDecodeError:
            raise RuntimeError(f"Server did not return JSON: {text[:200]}")
    raise RuntimeError("Failed to submit to ColabFold MSA server after retries")


def poll_status(ticket_id: str, host_url: str, headers: dict) -> dict:
    req = urllib.request.Request(f"{host_url}/ticket/{ticket_id}", headers=headers)
    for _ in range(5):
        try:
            with _urlopen(req, timeout=30) as res:
                return json.loads(res.read().decode())
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
            log.warning("Error fetching status: %s. Retrying...", e)
            time.sleep(5 + random.randint(0, 5))
        except json.JSONDecodeError:
            return {"status": "ERROR"}
    return {"status": "ERROR"}


def download_result(ticket_id: str, out_path: Path, host_url: str, headers: dict) -> None:
    req = urllib.request.Request(f"{host_url}/result/download/{ticket_id}", headers=headers)
    for _ in range(5):
        try:
            with _urlopen(req, timeout=120) as res:
                with open(out_path, "wb") as f:
                    f.write(res.read())
            return
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
            log.warning("Error downloading result: %s. Retrying...", e)
            time.sleep(5 + random.randint(0, 5))
    raise RuntimeError("Failed to download result after retries")


def run_remote_msa(
    seqs: List[str],
    prefix: str,
    host_url: str = DEFAULT_HOST,
    use_env: bool = True,
    use_filter: bool = True,
) -> List[str]:
    """
    Submit sequences to ColabFold MSA API, wait for completion, extract and merge a3m.
    Returns one merged a3m string per input sequence (same order as seqs).
    """
    mode = "env" if (use_filter and use_env) else "all" if use_filter else "env-nofilter" if use_env else "nofilter"
    path = f"{prefix}_{mode}"
    Path(path).mkdir(parents=True, exist_ok=True)
    tar_gz_file = Path(f"{path}/out.tar.gz")

    headers = {"User-Agent": USER_AGENT}
    N = 101
    seqs_unique: List[str] = []
    for s in seqs:
        if s not in seqs_unique:
            seqs_unique.append(s)
    Ms = [N + seqs_unique.index(s) for s in seqs]

    if not tar_gz_file.is_file():
        out = submit(seqs_unique, mode, host_url, headers, N)
        while out.get("status") in ("UNKNOWN", "RATELIMIT"):
            sleep_t = 5 + random.randint(0, 5)
            log.info("Sleeping %ds (status: %s)", sleep_t, out.get("status"))
            time.sleep(sleep_t)
            out = submit(seqs_unique, mode, host_url, headers, N)
        if out.get("status") == "ERROR":
            raise RuntimeError(
                "MMseqs2 API returned ERROR. Check input is valid protein sequence."
            )
        if out.get("status") == "MAINTENANCE":
            raise RuntimeError("MMseqs2 API is under maintenance. Try again later.")
        ticket_id = out["id"]
        while out.get("status") in ("UNKNOWN", "RUNNING", "PENDING"):
            t = 5 + random.randint(0, 5)
            log.info("Waiting for job (status: %s)...", out.get("status"))
            time.sleep(t)
            out = poll_status(ticket_id, host_url, headers)
        if out.get("status") == "ERROR":
            raise RuntimeError("MMseqs2 API job failed.")
        if out.get("status") != "COMPLETE":
            raise RuntimeError(f"Unexpected status: {out.get('status')}")
        download_result(ticket_id, tar_gz_file, host_url, headers)

    with tarfile.open(tar_gz_file) as tar:
        tar.extractall(path=path)
    path_p = Path(path)
    uniref = next(path_p.rglob("uniref.a3m"), None)
    env_a3m = next(path_p.rglob("bfd.mgnify30.metaeuk30.smag30.a3m"), None) if use_env else None
    a3m_files = [str(uniref)] if uniref else []
    if env_a3m:
        a3m_files.append(str(env_a3m))
    if not a3m_files:
        raise FileNotFoundError(f"No .a3m files found in {path} after extracting {tar_gz_file}")

    a3m_lines: Dict[int, List[str]] = {}
    for a3m_file in a3m_files:
        ap = Path(a3m_file)
        if not ap.exists():
            continue
        update_M = True
        M = None
        with open(ap) as f:
            for line in f:
                if "\x00" in line:
                    line = line.replace("\x00", "")
                    update_M = True
                if line.startswith(">") and update_M:
                    try:
                        M = int(line[1:].rstrip().split()[0])
                    except (ValueError, IndexError):
                        pass
                    update_M = False
                if M is not None:
                    if M not in a3m_lines:
                        a3m_lines[M] = []
                    a3m_lines[M].append(line)

    return ["".join(a3m_lines.get(n, [])) for n in Ms]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Query ColabFold MMseqs2 API for MSA and write .a3m file(s)",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="Input FASTA file (one or more sequences)",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output path: file (single seq) or directory (multiple seqs)",
    )
    parser.add_argument(
        "--host",
        default=DEFAULT_HOST,
        help="ColabFold API base URL [default: %(default)s]",
    )
    parser.add_argument(
        "--no-env",
        action="store_true",
        help="Disable environmental database",
    )
    parser.add_argument(
        "--no-filter",
        action="store_true",
        help="Disable MSA filtering",
    )
    args = parser.parse_args()

    if not args.fasta.exists():
        log.error("FASTA not found: %s", args.fasta)
        return 1

    pairs = read_fasta_sequences(args.fasta)
    if not pairs:
        log.error("No sequences found in %s", args.fasta)
        return 1
    seqs = [s for _, s in pairs]
    prefix = "colabfold_remote"

    a3m_strings = run_remote_msa(
        seqs,
        prefix,
        host_url=args.host,
        use_env=not args.no_env,
        use_filter=not args.no_filter,
    )

    # Decide base output path. If a directory is given, name the .a3m after the input FASTA
    # so multiple targets can coexist without clashing (e.g. PDL1_A.fasta -> PDL1_A.a3m).
    if args.output.suffix == ".a3m":
        out_path = args.output
    elif args.output.suffix:
        out_path = args.output
    else:
        base = args.fasta.stem
        out_path = Path(args.output) / f"{base}.a3m"
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if len(a3m_strings) == 1:
        with open(out_path, "w") as f:
            f.write(a3m_strings[0])
        log.info("Wrote %s", out_path)
    else:
        for i, content in enumerate(a3m_strings):
            p = out_path.parent / f"{out_path.stem}_{i}{out_path.suffix}"
            with open(p, "w") as f:
                f.write(content)
            log.info("Wrote %s", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
