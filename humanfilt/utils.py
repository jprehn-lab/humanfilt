from pathlib import Path
import gzip, shutil, subprocess, os

def have(exe: str) -> bool:
    return shutil.which(exe) is not None

def ensure_tools(required=(), anyof=()):
    missing = [c for c in required if not have(c)]
    if missing:
        raise SystemExit(f"[humanfilt] Missing tools: {', '.join(missing)}")
    if anyof:
        if not any(have(c) for c in anyof):
            raise SystemExit(f"[humanfilt] Need at least one of: {', '.join(anyof)}")

def call(cmd, **kw):
    return subprocess.check_call([str(x) for x in cmd], **kw)

def nreads_fastq(path: Path) -> int:
    """Fast read count: prefers seqkit, falls back to counting lines/4 (supports .gz)."""
    path = Path(path)
    if have("seqkit"):
        try:
            out = subprocess.check_output(["seqkit","stats","-T",str(path)], text=True)
            line = [l for l in out.splitlines() if l.strip()][-1]
            return int(line.split('\t')[3])  # num_seqs
        except Exception:
            pass
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="ignore") as fh:
        lines = sum(1 for _ in fh)
    return lines // 4

def user_cache_dir() -> Path:
    xdg = os.environ.get("XDG_DATA_HOME")
    return Path(xdg)/"humanfilt" if xdg else Path.home()/".local"/"share"/"humanfilt"

def detect_threads(default_if_unknown: int = 8) -> int:
    for v in ("HUMANFILT_THREADS","THREADS","OMP_NUM_THREADS","CPU_COUNT"):
        val = os.environ.get(v)
        if val and str(val).isdigit():
            return max(1, int(val))
    try:
        c = os.cpu_count() or default_if_unknown
        return max(1, int(c))
    except Exception:
        return default_if_unknown

