#!/usr/bin/env python3
"""
humanfilt.downloader
Downloads reference bundles from Zenodo, verifies (if .sha256 sidecars exist),
extracts them into a writable cache, and writes config.json with discovered paths.
"""

import os
import sys
import json
import hashlib
import subprocess
from pathlib import Path

# ---------------------------
# Defaults & environment
# ---------------------------

# Hard-default Zenodo record; user may override via env at runtime
os.environ.setdefault("HUMANFILT_ZENODO_RECORD", "17020482")
RECORD_ID = os.environ["HUMANFILT_ZENODO_RECORD"]

# Bundles expected in the Zenodo record (filenames must match exactly)
BUNDLES = [
    {"name": "grch38.tar.zst",                         "tag": "grch38"},
    {"name": "t2t.tar.zst",                            "tag": "t2t"},
    {"name": "hprc.tar.zst",                           "tag": "hprc"},
    {"name": "gencode.tar.zst",                        "tag": "gencode"},
    {"name": "univec.tar.zst",                         "tag": "univec"},
    {"name": "humanfilt-kraken2-human-202409.tar.zst", "tag": "kraken2"},
]

# ---------------------------
# Small helpers
# ---------------------------

def _print(msg: str):
    print(f"[humanfilt] {msg}", file=sys.stderr)

def _user_cache_dir() -> Path:
    xdg = os.environ.get("XDG_DATA_HOME")
    return Path(xdg) / "humanfilt" if xdg else Path.home() / ".local" / "share" / "humanfilt"

def _ensure_writable(p: Path) -> bool:
    try:
        p.mkdir(parents=True, exist_ok=True)
        t = p / ".hf_write_test"
        t.write_text("ok")
        t.unlink()
        return True
    except Exception:
        return False

def _curl(url: str, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    _print(f"downloading {url} -> {dest}")
    # NOTE: requires 'curl' (declared in conda recipe)
    subprocess.check_call(["curl", "-L", "--fail", "--progress-bar", "-o", str(tmp), url])
    tmp.replace(dest)

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def _try_verify_with_sidecar(archive: Path):
    """
    If <file>.sha256 exists in the Zenodo record, verify the archive.
    If not present (404), skip verification silently.
    """
    sidecar_url = f"https://zenodo.org/records/{RECORD_ID}/files/{archive.name}.sha256?download=1"
    sidecar = archive.with_suffix(archive.suffix + ".sha256.txt")
    try:
        _curl(sidecar_url, sidecar)
        want = Path(sidecar).read_text().strip().split()[0]
        got = _sha256(archive)
        if want != got:
            raise RuntimeError(f"Checksum mismatch for {archive.name}: got {got}, expected {want}")
        _print(f"sha256 OK for {archive.name}")
    except subprocess.CalledProcessError:
        _print(f"no sidecar sha256 for {archive.name}; skipping verify")
    finally:
        sidecar.unlink(missing_ok=True)

def _extract(archive: Path, into: Path):
    """
    Extracts tar archives, supporting .tar.zst (zstd), .tar.gz and .tar.
    Requires 'tar' (and 'zstd' for .zst) — both declared in the conda recipe.
    """
    into.mkdir(parents=True, exist_ok=True)
    name = archive.name.lower()
    _print(f"extracting {archive.name} -> {into}")
    if name.endswith(".tar.zst") or name.endswith(".tzst") or name.endswith(".zst"):
        subprocess.check_call(["tar", "-I", "zstd", "-xf", str(archive), "-C", str(into)])
    elif name.endswith(".tar.gz") or name.endswith(".tgz"):
        subprocess.check_call(["tar", "-xzf", str(archive), "-C", str(into)])
    elif name.endswith(".tar"):
        subprocess.check_call(["tar", "-xf", str(archive), "-C", str(into)])
    else:
        raise RuntimeError(f"Unsupported archive format: {archive}")

def _scan_cfg(base: Path) -> dict:
    """
    Walks the extracted tree and deduces paths needed by the pipelines.
    Returns a config dict written to base/config.json
    """
    cfg = {"data_dir": str(base.resolve())}

    # GRCh38 fasta (and fai typically nearby)
    grch_fa = next((p for p in base.rglob("GCF_*GRCh38*.f*a*")), None) \
           or next((p for p in base.rglob("*grch38*.f*a*")), None)
    cfg["grch38_fa"] = str(grch_fa) if grch_fa else ""

    # T2T bowtie2 prefix (.1.bt2 file without the suffix)
    t2t_bt2_1 = next((p for p in base.rglob("*.1.bt2")), None)
    cfg["t2t_bt2_prefix"] = str(t2t_bt2_1)[:-5] if t2t_bt2_1 else ""

    # HPRC fasta & minimap2 index
    hprc_fa = next((p for p in base.rglob("*hprc*.*fa*")), None)
    hprc_mmi = next((p for p in base.rglob("*hprc*.mmi")), None)
    cfg["hprc_fa"] = str(hprc_fa) if hprc_fa else ""
    cfg["hprc_mmi"] = str(hprc_mmi) if hprc_mmi else ""

    # GENCODE transcripts fasta & mmi
    gen_fa = next((p for p in base.rglob("*gencode*transcripts*.f*a*")), None)
    gen_mmi = next((p for p in base.rglob("*gencode*transcripts*.mmi")), None)
    cfg["gencode_tx_fa"] = str(gen_fa) if gen_fa else ""
    cfg["gencode_tx_mmi"] = str(gen_mmi) if gen_mmi else ""

    # UniVec core file
    univec = next((p for p in base.rglob("UniVec_Core*")), None)
    cfg["univec_ref"] = str(univec) if univec else ""

    # Kraken2 DB directory (folder with hash.k2d/opts.k2d/taxo.k2d)
    kdir = ""
    kroot = base / "kraken2"
    if kroot.exists():
        for d in kroot.iterdir():
            if d.is_dir() and all((d / f).exists() for f in ("hash.k2d", "opts.k2d", "taxo.k2d")):
                kdir = str(d)
                break
    cfg["kraken2_db"] = kdir

    return cfg

# ---------------------------
# Public API
# ---------------------------

def ensure_data_ready(data_dir=None, force: bool = False) -> dict:
    """
    Ensure all Zenodo bundles are present & extracted, then write config.json.
    Returns the config dict.
    """
    requested = Path(data_dir).expanduser() if data_dir else _user_cache_dir()
    base = requested if _ensure_writable(requested) else _user_cache_dir()
    if not _ensure_writable(base):
        raise SystemExit(f"Cannot write to '{requested}' or fallback '{base}'")

    bundles_dir = base / "bundles"
    bundles_dir.mkdir(parents=True, exist_ok=True)

    for spec in BUNDLES:
        url = f"https://zenodo.org/records/{RECORD_ID}/files/{spec['name']}?download=1"
        dst = bundles_dir / spec["name"]
        if force and dst.exists():
            _print(f"removing (force) {dst}")
            dst.unlink()
        if not dst.exists():
            _curl(url, dst)
            _try_verify_with_sidecar(dst)
        _extract(dst, base)

    cfg = _scan_cfg(base)
    cfg_path = base / "config.json"
    cfg_path.write_text(json.dumps(cfg, indent=2))
    _print(f"config -> {cfg_path}")
    return cfg

# ---------------------------
# CLI (for manual testing)
# ---------------------------

def _parse_cli(argv):
    import argparse
    ap = argparse.ArgumentParser(description="Download/extract HumanFilt references from Zenodo.")
    ap.add_argument("--data-dir", default=None, help="Cache directory (default: ~/.local/share/humanfilt)")
    ap.add_argument("--force", action="store_true", help="Re-download archives if present")
    return ap.parse_args(argv)

if __name__ == "__main__":
    args = _parse_cli(sys.argv[1:])
    ensure_data_ready(args.data_dir, force=args.force)

