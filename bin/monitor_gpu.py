#!/usr/bin/env python3
# /// script
# dependencies = []
# ///

import subprocess
import sys
import argparse
import logging
import signal
import time
import csv
import os
from typing import Optional


def query_gpu_metrics() -> Optional[list]:
    """Query GPU metrics using nvidia-smi.
    
    Returns:
        List of dictionaries with GPU metrics, or None on error.
        Each dict contains: timestamp, name, pci_bus_id, temperature_gpu,
        utilization_gpu, utilization_memory, memory_used, memory_total, memory_free
    """
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=timestamp,name,pci.bus_id,temperature.gpu,utilization.gpu,utilization.memory,memory.used,memory.total,memory.free",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        )
        
        gpus = []
        for line in result.stdout.strip().split("\n"):
            if not line.strip():
                continue
            
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 9:
                logging.warning(f"Skipping malformed line: {line}")
                continue
            
            gpus.append({
                "timestamp": parts[0],
                "name": parts[1],
                "pci_bus_id": parts[2],
                "temperature_gpu": parts[3],
                "utilization_gpu": parts[4],
                "utilization_memory": parts[5],
                "memory_used": parts[6],
                "memory_total": parts[7],
                "memory_free": parts[8],
            })
        
        return gpus
    
    except subprocess.CalledProcessError as e:
        logging.error(f"nvidia-smi failed: {e.stderr}")
        return None
    except subprocess.TimeoutExpired:
        logging.error("nvidia-smi timed out")
        return None
    except Exception as e:
        logging.error(f"Unexpected error querying GPU: {e}")
        return None


def get_job_id() -> str:
    """Get job ID from environment variables (SLURM_JOB_ID or PBS_JOBID)."""
    return os.environ.get("SLURM_JOB_ID") or os.environ.get("PBS_JOBID") or ""


def write_csv_header(writer: csv.DictWriter, process_name: str, task_hash: str, task_index: str, job_id: str) -> None:
    """Write CSV header with metadata columns."""
    fieldnames = [
        "timestamp",
        "gpu_name",
        "pci_bus_id",
        "temperature_gpu",
        "utilization_gpu",
        "utilization_memory",
        "memory_used",
        "memory_total",
        "memory_free",
        "process_name",
        "task_hash",
        "task_index",
        "hpc_job_id",
    ]
    writer.writeheader()


def write_gpu_row(writer: csv.DictWriter, gpu_data: dict, process_name: str, task_hash: str, task_index: str, job_id: str) -> None:
    """Write a single GPU metrics row with metadata."""
    writer.writerow({
        "timestamp": gpu_data["timestamp"],
        "gpu_name": gpu_data["name"],
        "pci_bus_id": gpu_data["pci_bus_id"],
        "temperature_gpu": gpu_data["temperature_gpu"],
        "utilization_gpu": gpu_data["utilization_gpu"],
        "utilization_memory": gpu_data["utilization_memory"],
        "memory_used": gpu_data["memory_used"],
        "memory_total": gpu_data["memory_total"],
        "memory_free": gpu_data["memory_free"],
        "process_name": process_name,
        "task_hash": task_hash,
        "task_index": task_index,
        "hpc_job_id": job_id,
    })


class GPUMonitor:
    """GPU monitoring loop that writes metrics to CSV."""
    
    def __init__(self, output_file: str, process_name: str, task_hash: str, task_index: str, job_id: str, interval: float = 1.0):
        self.output_file = output_file
        self.process_name = process_name
        self.task_hash = task_hash
        self.task_index = task_index
        self.job_id = job_id
        self.interval = interval
        self.running = True
        self.csv_file = None
        self.writer = None
        
        signal.signal(signal.SIGTERM, self._signal_handler)
        signal.signal(signal.SIGINT, self._signal_handler)
    
    def _signal_handler(self, signum, frame):
        """Handle shutdown signals gracefully."""
        logging.info(f"Received signal {signum}, shutting down...")
        self.running = False
    
    def run(self) -> None:
        """Run the monitoring loop."""
        try:
            self.csv_file = open(self.output_file, "w", newline="")
            self.writer = csv.DictWriter(
                self.csv_file,
                fieldnames=[
                    "timestamp",
                    "gpu_name",
                    "pci_bus_id",
                    "temperature_gpu",
                    "utilization_gpu",
                    "utilization_memory",
                    "memory_used",
                    "memory_total",
                    "memory_free",
                    "process_name",
                    "task_hash",
                    "task_index",
                    "hpc_job_id",
                ],
                lineterminator="\n",
            )
            write_csv_header(self.writer, self.process_name, self.task_hash, self.task_index, self.job_id)
            
            while self.running:
                gpus = query_gpu_metrics()
                if gpus:
                    for gpu in gpus:
                        write_gpu_row(self.writer, gpu, self.process_name, self.task_hash, self.task_index, self.job_id)
                    self.csv_file.flush()
                else:
                    logging.warning("Failed to query GPU metrics, continuing...")
                
                time.sleep(self.interval)
        
        except KeyboardInterrupt:
            logging.info("Interrupted by user")
        except Exception as e:
            logging.error(f"Error in monitoring loop: {e}")
            sys.exit(1)
        finally:
            if self.csv_file:
                self.csv_file.close()


def main():
    parser = argparse.ArgumentParser(
        description="Monitor GPU utilisation metrics and write to CSV with task metadata"
    )
    parser.add_argument(
        "--process-name",
        required=True,
        help="Nextflow process name (e.g., BINDCRAFT)",
    )
    parser.add_argument(
        "--task-hash",
        required=True,
        help="Nextflow task hash",
    )
    parser.add_argument(
        "--task-index",
        required=True,
        help="Nextflow task index",
    )
    parser.add_argument(
        "--job-id",
        default=None,
        help="Job ID (overrides SLURM_JOB_ID/PBS_JOBID environment variables)",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=1.0,
        help="Sampling interval in seconds (default: 1.0)",
    )
    parser.add_argument(
        "--output",
        default="gpu_stats.csv",
        help="Output CSV file path (default: gpu_stats.csv)",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    
    args = parser.parse_args()
    
    log_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=log_level, format="%(message)s", stream=sys.stderr)
    
    job_id = args.job_id if args.job_id is not None else get_job_id()
    
    monitor = GPUMonitor(
        output_file=args.output,
        process_name=args.process_name,
        task_hash=args.task_hash,
        task_index=args.task_index,
        job_id=job_id,
        interval=args.interval,
    )
    
    monitor.run()


if __name__ == "__main__":
    main()
