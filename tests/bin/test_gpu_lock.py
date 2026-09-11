#!/usr/bin/env python
# /// script
# requires-python = ">=3.10"
# dependencies = ["pytest"]
# ///

"""Tests for bin/gpu_lock.sh.

The property that matters is mutual exclusion under a genuine race, so these
tests launch real concurrent processes rather than asserting on the script's
text. `nvidia-smi` is stubbed so the tests run anywhere, including CI without
GPUs.
"""

import os
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

GPU_LOCK = Path(__file__).resolve().parents[2] / "bin" / "gpu_lock.sh"


@pytest.fixture
def env(tmp_path):
    """PATH with a stub nvidia-smi reporting two GPUs."""
    binp = tmp_path / "bin"
    binp.mkdir()
    smi = binp / "nvidia-smi"
    smi.write_text("#!/bin/bash\necho 0\necho 1\n")
    smi.chmod(0o755)
    e = dict(os.environ)
    e["PATH"] = f"{binp}:{e['PATH']}"
    return e


def _worker_script(tmp_path, candidates="0,1", slots=1, timeout=30, hold=0.6):
    """A task that claims a GPU, records its occupancy interval, then exits."""
    w = tmp_path / "worker.sh"
    w.write_text(textwrap.dedent(f"""\
        #!/bin/bash
        source {GPU_LOCK}
        nfbd_acquire_gpu "{candidates}" "{tmp_path}/locks" {slots} {timeout} >/dev/null 2>&1 \
            || {{ echo "FAIL" >> {tmp_path}/events; exit 1; }}
        echo "$(date +%s.%N) START $CUDA_VISIBLE_DEVICES" >> {tmp_path}/events
        sleep {hold}
        echo "$(date +%s.%N) END $CUDA_VISIBLE_DEVICES" >> {tmp_path}/events
    """))
    w.chmod(0o755)
    return w


def _intervals(tmp_path):
    """Parse events into {gpu: [(start, end), ...]}."""
    per_gpu = {}
    open_start = {}
    for line in (tmp_path / "events").read_text().splitlines():
        if line.strip() == "FAIL":
            pytest.fail("a worker failed to acquire a GPU")
        ts, kind, gpu = line.split()
        if kind == "START":
            open_start.setdefault(gpu, []).append(float(ts))
        else:
            start = open_start[gpu].pop()
            per_gpu.setdefault(gpu, []).append((start, float(ts)))
    return per_gpu


def _max_overlap(intervals):
    """Peak number of simultaneously open intervals."""
    events = [(s, 1) for s, _ in intervals] + [(e, -1) for _, e in intervals]
    events.sort()
    cur = peak = 0
    for _, delta in events:
        cur += delta
        peak = max(peak, cur)
    return peak


@pytest.mark.skipif(shutil.which("flock") is None, reason="flock not available")
def test_never_double_books_a_gpu(tmp_path, env):
    """16 tasks racing for 2 GPUs must never put two on one card."""
    worker = _worker_script(tmp_path)
    procs = [subprocess.Popen([str(worker)], env=env) for _ in range(16)]
    for p in procs:
        assert p.wait(timeout=120) == 0

    per_gpu = _intervals(tmp_path)
    assert set(per_gpu) == {"0", "1"}, "both GPUs should be used"
    for gpu, ivs in per_gpu.items():
        assert _max_overlap(ivs) == 1, f"GPU {gpu} was double-booked"
    assert sum(len(v) for v in per_gpu.values()) == 16


@pytest.mark.skipif(shutil.which("flock") is None, reason="flock not available")
def test_slots_per_device_allows_sharing(tmp_path, env):
    """slots_per_device=2 permits exactly two concurrent tasks per GPU."""
    worker = _worker_script(tmp_path, slots=2)
    procs = [subprocess.Popen([str(worker)], env=env) for _ in range(8)]
    for p in procs:
        assert p.wait(timeout=120) == 0

    per_gpu = _intervals(tmp_path)
    peak = max(_max_overlap(ivs) for ivs in per_gpu.values())
    assert peak == 2, f"expected 2 concurrent tasks per GPU, saw {peak}"


@pytest.mark.skipif(shutil.which("flock") is None, reason="flock not available")
def test_times_out_when_no_slot_frees(tmp_path, env):
    """A task waiting on a permanently held single GPU fails rather than hanging."""
    lock_dir = tmp_path / "locks"
    lock_dir.mkdir()
    # Hold GPU 0's only slot for longer than the waiter's timeout.
    holder = subprocess.Popen(
        ["bash", "-c", f'exec 9>"{lock_dir}/gpu0.slot0.lock"; flock 9; sleep 30'],
        env=env,
    )
    try:
        result = subprocess.run(
            ["bash", "-c",
             f'source {GPU_LOCK}; nfbd_acquire_gpu "0" "{lock_dir}" 1 3'],
            env=env, capture_output=True, text=True, timeout=60,
        )
        assert result.returncode == 1
        assert "no GPU slot free" in result.stderr
    finally:
        holder.kill()
        holder.wait()


@pytest.mark.skipif(shutil.which("flock") is None, reason="flock not available")
def test_lock_released_when_task_is_killed(tmp_path, env):
    """A SIGKILLed task must not leave its GPU claimed."""
    lock_dir = tmp_path / "locks"
    lock_dir.mkdir()
    victim = subprocess.Popen(
        ["bash", "-c",
         f'source {GPU_LOCK}; nfbd_acquire_gpu "0" "{lock_dir}" 1 10 >/dev/null; sleep 60'],
        env=env,
    )
    # Wait for the claim to actually be taken before killing.
    deadline = __import__("time").time() + 20
    while __import__("time").time() < deadline:
        probe = subprocess.run(
            ["bash", "-c",
             f'exec 9>"{lock_dir}/gpu0.slot0.lock"; flock -n 9 && echo FREE || echo HELD'],
            env=env, capture_output=True, text=True)
        if "HELD" in probe.stdout:
            break
        __import__("time").sleep(0.2)
    else:
        pytest.fail("victim never acquired the lock")

    victim.kill()
    victim.wait()

    result = subprocess.run(
        ["bash", "-c", f'source {GPU_LOCK}; nfbd_acquire_gpu "0" "{lock_dir}" 1 5'],
        env=env, capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, "lock was not released when the holder was killed"


@pytest.mark.skipif(shutil.which("flock") is None, reason="flock not available")
def test_fills_every_gpu_before_doubling_up(tmp_path, env):
    """With slots>1, an idle GPU must be preferred over a second slot on a busy one.

    This is the failure the whole script exists to prevent, and it is invisible
    to the exclusivity tests: ordering the slot list GPU-major keeps mutual
    exclusion perfectly intact while still stacking two tasks on one card and
    leaving the other idle.

    Placement rotates the GPU order per task, so a single two-task trial splits
    correctly about half the time even when the ordering is wrong. The test
    therefore repeats the trial: two concurrent tasks must land on different
    GPUs *every* time. Under the GPU-major bug each trial is a coin flip, so
    ten trials fail with probability ~1 - 2^-10.
    """
    for trial in range(10):
        run = tmp_path / f"t{trial}"
        run.mkdir()
        worker = _worker_script(run, slots=4, hold=0.8)
        procs = [subprocess.Popen([str(worker)], env=env) for _ in range(2)]
        for p in procs:
            assert p.wait(timeout=120) == 0
        per_gpu = _intervals(run)
        counts = {g: len(v) for g, v in per_gpu.items()}
        assert set(per_gpu) == {"0", "1"}, (
            f"trial {trial}: two concurrent tasks stacked on one GPU instead of "
            f"using both, got {counts}"
        )


@pytest.mark.skipif(shutil.which("flock") is None, reason="flock not available")
def test_second_slot_only_used_once_every_gpu_is_busy(tmp_path, env):
    """A third task may double up, but only after both cards are occupied."""
    worker = _worker_script(tmp_path, slots=4, hold=1.2)
    procs = [subprocess.Popen([str(worker)], env=env) for _ in range(3)]
    for p in procs:
        assert p.wait(timeout=120) == 0

    per_gpu = _intervals(tmp_path)
    counts = {g: len(v) for g, v in per_gpu.items()}
    assert set(counts) == {"0", "1"}
    assert sorted(counts.values()) == [1, 2], (
        f"three tasks over two GPUs should place 1 and 2, got {counts}"
    )
