#!/usr/bin/env python3
"""
Download/extract HumanFilt reference bundles from Zenodo, auto-discover paths,
and validate the final config. Robust to layouts like t2t/bowtie2/*.bt2.
"""

from __future__ import annotations
import os
import sys
import json
import hashlib
import subprocess
from pathlib import Path
from typing import Optional, Dict

# ---------------------------
# Defaults & environment
# ---------------------------

# Baked-in record; user may override at runtime
os.environ.setdefault("HUMANFILT_ZENODO_RECORD", "17020482")
RECORD_ID = os.environ["HUMANFILT_ZENODO_RECORD"]

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

def _say(msg: str) -> None:
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

def _curl(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    _say(f"downloading {url} -> {dest}")
    subprocess.check_call(["curl", "-L", "--fail", "--progress-bar", "-o", str(tmp), url])
    tmp.replace(dest)

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def _try_verify_with_sidecar(archive: Path) -> None:
    """Try <file>.sha256 sidecar from Zenodo; if not present, skip."""
    sidecar_url = f"https://zenodo.org/records/{RECORD_ID}/files/{archive.name}.sha256?download=1"
    sidecar = archive.with_suffix(archive.suffix + ".sha256.txt")
    try:
        _curl(sidecar_url, sidecar)
        want = Path(sidecar).read_text().strip().split()[0]
        got = _sha256(archive)
        if want != got:
            raise RuntimeError(f"Checksum mismatch for {archive.name}: got {got}, expected {want}")
        _say(f"sha256 OK for {archive.name}")
    except subprocess.CalledProcessError:
        _say(f"no sidecar sha256 for {archive.name}; skipping verify")
    finally:
        sidecar.unlink(missing_ok=True)

def _extract(archive: Path, into: Path) -> None:
    """Extract .tar.zst/.tar.gz/.tar using tar (and zstd)."""
    into.mkdir(parents=True, exist_ok=True)
    name = archive.name.lower()
    _say(f"extracting {archive.name} -> {into}")
    if name.endswith(".tar.zst") or name.endswith(".tzst") or name.endswith(".zst"):
        subprocess.check_call(["tar", "-I", "zstd", "-xf", str(archive), "-C", str(into)])
    elif name.endswith(".tar.gz") or name.endswith(".tgz"):
        subprocess.check_call(["tar", "-xzf", str(archive), "-C", str(into)])
    elif name.endswith(".tar"):
        subprocess.check_call(["tar", "-xf", str(archive), "-C", str(into)])
    else:
        raise RuntimeError(f"Unsupported archive format: {archive}")

# ---------------------------
# Presence checks per bundle
# ---------------------------

def _bundle_present(tag: str, base: Path) -> bool:
    """Return True if the extracted content for this tag appears present."""
    if tag == "grch38":
        d = base / "grch38"
        if not d.exists():
            return False
        fa = next(d.glob("*.fna"), None)
        have_idx = bool(list(d.glob("*.bwt"))) and bool(list(d.glob("*.amb")))
        return fa is not None and have_idx

    if tag == "t2t":
        d = base / "t2t"
        if not d.exists():
            return False
        # recursive check for any *.1.bt2 (works with t2t/bowtie2/*)
        return next(d.rglob("*.1.bt2"), None) is not None

    if tag == "hprc":
        d = base / "hprc"
        return (d / "hprc-v1.1-mc-grch38.gfa.fa").exists() and (d / "hprc-v1.1-mc-grch38.mmi").exists()

    if tag == "gencode":
        d = base / "gencode"
        return (d / "gencode.v44.transcripts.fa").exists() and (d / "gencode_v44.mmi").exists()

    if tag == "univec":
        d = base / "univec"
        return (d / "UniVec_Core").exists() or bool(list(d.glob("*.nsq")))

    if tag == "kraken2":
        d = base / "kraken2" / "human_202409"
        return all((d / f).exists() for f in ("hash.k2d", "opts.k2d", "taxo.k2d"))

    return False

# ---------------------------
# Discovery
# ---------------------------

def _derive_bwa_prefix_for_fasta(fa: Path) -> Optional[Path]:
    """
    Find BWA index prefix for a given FASTA.
    We check:
      - <fa>.* (prefix includes .fna)
      - <fa without extension>.*
      - bwa/<fa>.* and bwa/<stem>.*
    """
    d = fa.parent
    stem_no_ext = fa.with_suffix("")  # drop .fna
    candidates = [
        fa,
        stem_no_ext,
        d / "bwa" / fa.name,
        d / "bwa" / stem_no_ext.name,
    ]
    need = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    for pref in candidates:
        if all((Path(str(pref) + suf)).exists() for suf in need):
            return Path(pref)
    any_bwt = next(d.glob("*.bwt"), None)
    if any_bwt:
        return Path(str(any_bwt)[:-4])
    return None

def _scan_cfg(base: Path) -> Dict[str, str]:
    """
    Build config dict with discovered paths/prefixes.
    Matches your on-disk layout.
    """
    cfg: Dict[str, str] = {"data_dir": str(base.resolve())}

    # GRCh38
    gr_dir = base / "grch38"
    gr_fa = next(gr_dir.glob("*.fna"), None) if gr_dir.exists() else None
    cfg["grch38_fa"] = str(gr_fa) if gr_fa else ""
    gr_pref = _derive_bwa_prefix_for_fasta(gr_fa) if gr_fa else None
    cfg["grch38_bwa_prefix"] = str(gr_pref) if gr_pref else ""

    # T2T Bowtie2 prefix — search recursively and strip the exact ".1.bt2" suffix
    t2_dir = base / "t2t"
    t2_bt2_1 = next(t2_dir.rglob("*.1.bt2"), None) if t2_dir.exists() else None
    if t2_bt2_1:
        s = str(t2_bt2_1)
        cfg["t2t_bt2_prefix"] = s.removesuffix(".1.bt2")
    else:
        cfg["t2t_bt2_prefix"] = ""

    # HPRC
    hp_dir = base / "hprc"
    hp_fa = hp_dir / "hprc-v1.1-mc-grch38.gfa.fa"
    hp_mmi = hp_dir / "hprc-v1.1-mc-grch38.mmi"
    cfg["hprc_fa"] = str(hp_fa) if hp_fa.exists() else ""
    cfg["hprc_mmi"] = str(hp_mmi) if hp_mmi.exists() else ""

    # GENCODE
    ge_dir = base / "gencode"
    ge_fa = ge_dir / "gencode.v44.transcripts.fa"
    ge_mmi = ge_dir / "gencode_v44.mmi"
    cfg["gencode_tx_fa"] = str(ge_fa) if ge_fa.exists() else ""
    cfg["gencode_tx_mmi"] = str(ge_mmi) if ge_mmi.exists() else ""

    # UniVec
    uv_dir = base / "univec"
    uv = uv_dir / "UniVec_Core"
    cfg["univec_ref"] = str(uv) if uv.exists() else str(next(uv_dir.glob("*.nsq"), Path("")))

    # Kraken2
    k_dir = base / "kraken2" / "human_202409"
    if all((k_dir / f).exists() for f in ("hash.k2d", "opts.k2d", "taxo.k2d")):
        cfg["kraken2_db"] = str(k_dir)
    else:
        cfg["kraken2_db"] = ""

    return cfg

# Public: non-downloading config load/scan
def load_or_scan_config(data_dir: Optional[str] = None) -> Dict[str, str]:
    """
    Load config.json if present in the chosen data dir; otherwise scan the
    directory contents to derive paths. Does NOT attempt any downloads.
    """
    requested = Path(data_dir).expanduser() if data_dir else _user_cache_dir()
    base = requested
    cfg_path = base / "config.json"
    if cfg_path.exists():
        try:
            return json.loads(cfg_path.read_text())
        except Exception:
            # fall back to scanning
            pass
    return _scan_cfg(base)

# ---------------------------
# Validation
# ---------------------------

def _must_bt2(prefix: str) -> None:
    need = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
    miss = [prefix + s for s in need if not Path(prefix + s).exists()]
    if miss:
        raise SystemExit("[humanfilt] Missing bowtie2 index files for T2T:\n  " + "\n  ".join(miss))
    # Sanity: confirm the index is readable by bowtie2-inspect
    try:
        subprocess.check_call(["bowtie2-inspect", "-n", prefix], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        sizes = {s: (Path(prefix + s).stat().st_size if Path(prefix + s).exists() else 0) for s in need}
        raise SystemExit(
            "[humanfilt] Corrupted or unreadable bowtie2 index for T2T. "
            f"Prefix: {prefix}\n  Sizes: " + ", ".join([f"{k}:{v}" for k,v in sizes.items()]) +
            "\n  Try re-running 'humanfilt setup --force' or rebuilding the index."
        )

def _must_bwa(prefix: str) -> None:
    need = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    miss = [prefix + s for s in need if not Path(prefix + s).exists()]
    if miss:
        raise SystemExit("[humanfilt] Missing BWA index files for GRCh38:\n  " + "\n  ".join(miss))

def _must_exist(path: str, label: str) -> None:
    if not path or not Path(path).exists():
        raise SystemExit(f"[humanfilt] Missing or not found: {label} -> {path}")

def validate_cfg_or_die(cfg: Dict[str, str]) -> None:
    # For prefixes, do NOT require the bare prefix path to exist; check index pieces instead.
    _must_exist(cfg.get("grch38_fa", ""), "grch38_fa")
    _must_exist(cfg.get("hprc_fa", ""), "hprc_fa")
    _must_exist(cfg.get("univec_ref", ""), "univec_ref")
    _must_exist(cfg.get("kraken2_db", ""), "kraken2_db")
    _must_bwa(cfg.get("grch38_bwa_prefix", ""))
    _must_bt2(cfg.get("t2t_bt2_prefix", ""))

# ---------------------------
# Public API
# ---------------------------

def ensure_data_ready(data_dir: Optional[str] = None, force: bool = False) -> Dict[str, str]:
    """
    Ensure bundles are present/extracted, write config.json, return config.
    Skips download/extract if we detect content is already present.
    """
    requested = Path(data_dir).expanduser() if data_dir else _user_cache_dir()
    base = requested if _ensure_writable(requested) else _user_cache_dir()
    if not _ensure_writable(base):
        raise SystemExit(f"Cannot write to '{requested}' or fallback '{base}'")

    bundles_dir = base / "bundles"
    bundles_dir.mkdir(parents=True, exist_ok=True)

    for spec in BUNDLES:
        name = spec["name"]
        tag = spec["tag"]
        archive = bundles_dir / name

        if not force and _bundle_present(tag, base):
            _say(f"already present: {name}")
        else:
            if not archive.exists() or force:
                url = f"https://zenodo.org/records/{RECORD_ID}/files/{name}?download=1"
                _curl(url, archive)
                _try_verify_with_sidecar(archive)
            _extract(archive, base)
            if not _bundle_present(tag, base):
                raise SystemExit(f"After extracting {name}, expected files for '{tag}' were not found.")

    cfg = _scan_cfg(base)
    cfg_path = base / "config.json"
    cfg_path.write_text(json.dumps(cfg, indent=2))
    _say(f"config -> {cfg_path}")
    return cfg

# ---------------------------
# CLI (manual testing)
# ---------------------------

def _parse_cli(argv):
    import argparse
    ap = argparse.ArgumentParser(description="Download/extract HumanFilt references from Zenodo.")
    ap.add_argument("--data-dir", default=None, help="Cache directory (default: ~/.local/share/humanfilt)")
    ap.add_argument("--force", action="store_true", help="Re-download archives if present")
    ap.add_argument("--validate", action="store_true", help="Validate config after (or without) extraction")
    return ap.parse_args(argv)

if __name__ == "__main__":
    args = _parse_cli(sys.argv[1:])
    cfg = ensure_data_ready(args.data_dir, force=args.force)
    if args.validate:
        validate_cfg_or_die(cfg)
