#!/usr/bin/env python3
from __future__ import annotations
import os, sys, shutil, subprocess, tempfile
from pathlib import Path
from contextlib import contextmanager
from typing import Iterable, Optional

def say(msg: str) -> None:
    print(f"[humanfilt] {msg}", file=sys.stderr)

def available_threads() -> int:
    for k in ("HUMANFILT_THREADS", "THREADS", "CPU_COUNT", "SLURM_CPUS_PER_TASK", "NSLOTS", "OMP_NUM_THREADS"):
        v = os.environ.get(k)
        if v and v.isdigit():
            return max(1, int(v))
    return max(1, os.cpu_count() or 1)

def ensure_tools(names: Iterable[str]) -> None:
    missing = [n for n in names if shutil.which(n) is None]
    if missing:
        raise SystemExit("Missing required tools: " + ", ".join(missing))

@contextmanager
def temp_workdir(prefix: str = "hf_"):
    d = tempfile.mkdtemp(prefix=prefix)
    try:
        yield d
    finally:
        try:
            shutil.rmtree(d, ignore_errors=True)
        except Exception:
            pass

def _count_with_seqkit(path: Path, threads: Optional[int] = None) -> Optional[int]:
    try:
        cmd = ["seqkit", "stats", "-T"]
        if threads and threads > 0:
            cmd += ["-j", str(threads)]
        cmd += [str(path)]
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        lines = [ln for ln in res.stdout.strip().splitlines() if ln.strip()]
        if len(lines) < 2:
            return None
        row = lines[1].split("\t")
        if len(row) >= 4:
            return int(row[3])
    except Exception:
        return None
    return None

def _count_with_wc(path: Path, threads: Optional[int] = None) -> Optional[int]:
    try:
        if str(path).endswith(".gz"):
            pigz = shutil.which("pigz")
            if pigz:
                cmd = [pigz]
                if threads and threads > 0:
                    cmd += ["-p", str(threads)]
                cmd += ["-cd", str(path)]
            else:
                gz = shutil.which("gzip")
                cmd = [gz, "-cd", str(path)] if gz else ["zcat", str(path)]
        else:
            cmd = ["cat", str(path)]
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
        out, _ = p2.communicate(); p1.wait()
        if p2.returncode == 0 and out.strip().isdigit():
            return int(out.strip()) // 4
    except Exception:
        return None
    return None

def count_reads(pathlike, threads: Optional[int] = None) -> int:
    p = Path(pathlike)
    n = _count_with_seqkit(p, threads)
    if n is not None:
        return n
    n = _count_with_wc(p, threads)
    if n is not None:
        return n
    say(f"WARNING: could not count reads for {p}; returning 0")
    return 0
