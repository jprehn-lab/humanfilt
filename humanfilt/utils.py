from pathlib import Path
import gzip, shutil, subprocess, sys

def nreads_fastq(path: Path) -> int:
    p = Path(path)
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt", errors="ignore") as fh:
        lines = 0
        for _ in fh:
            lines += 1
    return lines // 4

def ensure_tools(required=(), anyof=()):
    missing = [c for c in required if shutil.which(c) is None]
    if missing:
        sys.exit(f"Missing required tools: {', '.join(missing)}")
    if anyof and not any(shutil.which(c) for c in anyof):
        sys.exit(f"Need at least one of: {', '.join(anyof)}")

def call(cmd, **kw):
    cmd = [str(x) for x in cmd]
    return subprocess.check_call(cmd, **kw)
