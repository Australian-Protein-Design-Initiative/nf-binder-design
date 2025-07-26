#!/usr/bin/env python3

import subprocess
import sys
import argparse
import logging
import os
import random
import time
import re
from typing import List, Tuple, Optional


def get_available_gpus() -> List[int]:
    """Get list of all available GPU indices."""
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=index", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            check=True,
        )
        return [
            int(gpu_id.strip())
            for gpu_id in result.stdout.strip().split("\n")
            if gpu_id.strip()
        ]
    except subprocess.CalledProcessError:
        logging.error("Failed to query available GPUs")
        sys.exit(1)


def get_gpu_memory_usage(gpu_id: int) -> Tuple[int, int]:
    """Get memory usage for a specific GPU.

    Returns:
        Tuple of (memory_used_mb, memory_total_mb)
    """
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "-i",
                str(gpu_id),
                "--query-gpu=memory.used,memory.total",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        memory_line = result.stdout.strip()
        memory_used, memory_total = [int(x.strip()) for x in memory_line.split(",")]
        return memory_used, memory_total
    except subprocess.CalledProcessError:
        logging.error(f"Failed to query memory for GPU {gpu_id}")
        sys.exit(1)
    except ValueError:
        logging.error(
            f"Failed to parse memory output for GPU {gpu_id}: {result.stdout}"
        )
        sys.exit(1)


def parse_gpu_candidates(gpu_candidates: str) -> List[int]:
    """Parse GPU candidates string into list of GPU IDs."""
    if not gpu_candidates or gpu_candidates.lower() == "all":
        return get_available_gpus()

    try:
        return [
            int(gpu_id.strip())
            for gpu_id in gpu_candidates.split(",")
            if gpu_id.strip()
        ]
    except ValueError:
        logging.error(f"Invalid GPU candidates format: {gpu_candidates}")
        sys.exit(1)


def get_gpu_processes(gpu_id: int) -> List[int]:
    """Get list of process IDs running on a specific GPU."""
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "-i",
                str(gpu_id),
                "--query-compute-apps=pid",
                "--format=csv,noheader",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        pids = []
        for line in result.stdout.strip().split("\n"):
            line = line.strip()
            if line and line != "[Not Supported]":
                try:
                    pids.append(int(line))
                except ValueError:
                    continue
        return pids
    except subprocess.CalledProcessError:
        logging.warning(f"Failed to query processes for GPU {gpu_id}")
        return []


def get_process_cmdline(pid: int) -> str:
    """Get command line for a specific process ID."""
    # Try using ps aux command first - this should work even in containers
    try:
        logging.info(f"Using ps aux to get cmdline for PID {pid}")
        result = subprocess.run(
            ["ps", "aux"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0 and result.stdout.strip():
            # Parse ps aux output to find the specific PID
            for line in result.stdout.strip().split("\n"):
                if line.strip():
                    # Split by whitespace and get PID from second column
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            line_pid = int(parts[1])
                            if line_pid == pid:
                                # Reconstruct command line from remaining parts
                                cmdline = " ".join(
                                    parts[10:]
                                )  # Command starts at column 11
                                logging.info(
                                    f"PS aux successful for PID {pid}: {cmdline}"
                                )
                                return cmdline
                        except ValueError:
                            continue
            logging.info(f"PID {pid} not found in ps aux output")
        else:
            logging.info(
                f"PS aux command failed: returncode={result.returncode}, stderr='{result.stderr.strip()}'"
            )
    except subprocess.TimeoutExpired:
        logging.info(f"PS aux command timeout")
    except Exception as e:
        logging.info(f"PS aux error: {e}")

    # Fallback: try the /proc approach
    try:
        cmdline_path = f"/proc/{pid}/cmdline"
        logging.info(f"Trying /proc fallback for PID {pid}: {cmdline_path}")

        if not os.path.exists(cmdline_path):
            logging.info(f"Cmdline path does not exist: {cmdline_path}")
        else:
            with open(cmdline_path, "rb") as f:
                cmdline_bytes = f.read()
                logging.info(f"Read {len(cmdline_bytes)} bytes from {cmdline_path}")

                if len(cmdline_bytes) > 0:
                    # Replace null bytes with spaces for readable command line
                    cmdline = (
                        cmdline_bytes.replace(b"\x00", b" ")
                        .decode("utf-8", errors="ignore")
                        .strip()
                    )
                    logging.info(f"/proc fallback successful for PID {pid}: {cmdline}")
                    return cmdline
                else:
                    logging.info(
                        f"Empty cmdline file for PID {pid} - process might have exited"
                    )

    except PermissionError as e:
        logging.info(f"Permission denied reading cmdline for PID {pid}: {e}")
    except (OSError, IOError) as e:
        logging.info(f"OS/IO error reading cmdline for PID {pid}: {e}")
    except Exception as e:
        logging.info(f"Unexpected error reading cmdline for PID {pid}: {e}")

    logging.info(f"Failed to get cmdline for PID {pid} using all methods")
    return ""


def is_gpu_busy_with_process(gpu_id: int, process_pattern: str) -> bool:
    """Check if GPU is running any process matching the specified pattern in command line."""
    if not process_pattern:
        return False

    try:
        # Compile regex pattern
        pattern = re.compile(process_pattern)
    except re.error as e:
        logging.error(f"Invalid regex pattern '{process_pattern}': {e}")
        sys.exit(1)

    pids = get_gpu_processes(gpu_id)
    if pids:
        logging.info(f"Found {len(pids)} process(es) on GPU {gpu_id}: {pids}")

        for pid in pids:
            cmdline = get_process_cmdline(pid)
            logging.info(f"  Process {pid}: '{cmdline}'")

            if cmdline:  # We got a command line
                if pattern.search(cmdline):
                    logging.info(
                        f"MATCH FOUND: Process {pid} - {cmdline[:80]} found on GPU {gpu_id} - excluding this GPU"
                    )
                    return True
                else:
                    logging.info(f"  Process {pid} does NOT match exclude pattern")
            else:
                logging.info(
                    f"  WARNING: Process {pid} has empty command line - possibly due to container PID namespace isolation"
                )

        logging.info(f"No matching processes found on GPU {gpu_id}")
    else:
        logging.info(f"No processes running on GPU {gpu_id}")
    return False


def find_least_utilized_gpu(
    gpu_candidates: List[int], exclude_pattern: Optional[str] = None
) -> Optional[int]:
    """Find the GPU with the least memory utilization."""
    if not gpu_candidates:
        logging.error("No GPU candidates provided")
        sys.exit(1)

    best_gpu = None
    min_memory_used = float("inf")
    available_gpus = []

    for gpu_id in gpu_candidates:
        try:
            # Check if GPU is busy with the specified process
            if exclude_pattern and is_gpu_busy_with_process(gpu_id, exclude_pattern):
                continue

            memory_used, memory_total = get_gpu_memory_usage(gpu_id)
            logging.info(f"GPU {gpu_id}: {memory_used} MB / {memory_total} MB used")
            available_gpus.append(gpu_id)

            if memory_used < min_memory_used:
                min_memory_used = memory_used
                best_gpu = gpu_id

        except Exception as e:
            logging.warning(f"Skipping GPU {gpu_id} due to error: {e}")
            continue

    if best_gpu is None:
        if exclude_pattern and len(available_gpus) == 0:
            logging.error(
                f"No usable GPU found - all candidates are busy with {exclude_pattern}"
            )
        else:
            logging.error("No usable GPU found")
        sys.exit(1)

    logging.info(f"Selected GPU {best_gpu} with {min_memory_used} MB used")
    return best_gpu


def main():
    parser = argparse.ArgumentParser(
        description="Find a GPU not running one of our tasks, and with least VRAM utilization"
    )
    parser.add_argument(
        "gpu_candidates",
        nargs="?",
        default="all",
        help="Comma-separated list of GPU candidates or 'all' (default: all)",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )
    parser.add_argument(
        "--exclude",
        "-e",
        help="Exclude GPUs running processes matching this regex pattern in command line",
    )
    parser.add_argument(
        "--random-wait",
        type=float,
        metavar="SECONDS",
        help="Wait a random time (0 to SECONDS) before checking GPUs to avoid race conditions",
    )

    args = parser.parse_args()

    # Set up logging
    log_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=log_level, format="%(message)s", stream=sys.stderr)

    # Random wait to avoid race conditions
    if args.random_wait:
        wait_time = random.uniform(0, args.random_wait) + random.uniform(0, 1)
        logging.info(f"Waiting for {wait_time:.2f} seconds")
        time.sleep(wait_time)

    # Parse GPU candidates
    gpu_candidates = parse_gpu_candidates(args.gpu_candidates)
    logging.info(f"GPU candidates: {gpu_candidates}")

    # Find optimal GPU
    optimal_gpu = find_least_utilized_gpu(gpu_candidates, args.exclude)

    # Output only the GPU number (for easy shell capture)
    print(optimal_gpu)


if __name__ == "__main__":
    main()
