import os, json, subprocess, hashlib, sys
from pathlib import Path
from .utils import user_cache_dir
import os, json, subprocess, hashlib, sys
from pathlib import Path
from .utils import user_cache_dir

# Hard-default Zenodo record; user can still override via env

os.environ.setdefault("HUMANFILT_ZENODO_RECORD", "17020482")
RECORD_ID = os.environ("HUMANFILT_ZENODO_RECORD")

# Archives present in that record
BUNDLES = [
    {"name":"grch38.tar.zst",                         "tag":"grch38"},
    {"name":"t2t.tar.zst",                            "tag":"t2t"},
    {"name":"hprc.tar.zst",                           "tag":"hprc"},
    {"name":"gencode.tar.zst",                        "tag":"gencode"},
    {"name":"univec.tar.zst",                         "tag":"univec"},
    {"name":"humanfilt-kraken2-human-202409.tar.zst", "tag":"kraken2"},
]

def _curl(url: str, dest: Path):
    print(f"[humanfilt] downloading {url} -> {dest}", file=sys.stderr)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    subprocess.check_call(["curl","-L","--fail","--progress-bar","-o",str(tmp), url])
    tmp.rename(dest)

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

def _try_verify_with_sidecar(archive: Path):
    """If file.sha256 exists in Zenodo, verify; otherwise skip verification."""
    url = f"https://zenodo.org/records/{RECORD_ID}/files/{archive.name}.sha256?download=1"
    side = archive.with_suffix(archive.suffix + ".sha256.txt")
    try:
        _curl(url, side)
        want = Path(side).read_text().strip().split()[0]
        got  = _sha256(archive)
        if want != got:
            raise RuntimeError(f"Checksum mismatch for {archive.name}: {got} != {want}")
        print(f"[humanfilt] sha256 OK for {archive.name}", file=sys.stderr)
    except subprocess.CalledProcessError:
        print(f"[humanfilt] no sidecar sha256 for {archive.name}; skipping verify", file=sys.stderr)
    finally:
        side.unlink(missing_ok=True)

def _extract(archive: Path, into: Path):
    print(f"[humanfilt] extracting {archive.name}", file=sys.stderr)
    subprocess.check_call(["tar","-I","zstd","-xf",str(archive),"-C",str(into)])

def _ensure_writable(p: Path) -> bool:
    try:
        p.mkdir(parents=True, exist_ok=True)
        t = p/".hf_write_test"
        t.write_text("ok"); t.unlink()
        return True
    except Exception:
        return False

def _scan_cfg(base: Path) -> dict:
    """Infer final paths after extraction."""
    cfg = {"data_dir": str(base.resolve())}

    # GRCh38 fasta
    grch_fa = next((p for p in base.rglob("GCF_*GRCh38*.f*a*")), None) \
           or next((p for p in base.rglob("*grch38*.f*a*")), None)
    cfg["grch38_fa"] = str(grch_fa) if grch_fa else ""

    # T2T bowtie2 prefix
    t2t_bt2_1 = next((p for p in base.rglob("*.1.bt2")), None)
    cfg["t2t_bt2_prefix"] = str(t2t_bt2_1)[:-5] if t2t_bt2_1 else ""

    # HPRC fasta + mmi
    hprc_fa  = next((p for p in base.rglob("*hprc*.*fa*")), None)
    hprc_mmi = next((p for p in base.rglob("*hprc*.mmi")), None)
    cfg["hprc_fa"]  = str(hprc_fa) if hprc_fa else ""
    cfg["hprc_mmi"] = str(hprc_mmi) if hprc_mmi else ""

    # GENCODE transcripts + mmi
    gen_fa  = next((p for p in base.rglob("*gencode*transcripts*.f*a*")), None)
    gen_mmi = next((p for p in base.rglob("*gencode*transcripts*.mmi")), None)
    cfg["gencode_tx_fa"]  = str(gen_fa) if gen_fa else ""
    cfg["gencode_tx_mmi"] = str(gen_mmi) if gen_mmi else ""

    # UniVec_Core
    univec = next((p for p in base.rglob("UniVec_Core*")), None)
    cfg["univec_ref"] = str(univec) if univec else ""

    # Kraken2 DB directory
    kdir = ""
    kroot = base/"kraken2"
    if kroot.exists():
        for d in kroot.iterdir():
            if d.is_dir() and all((d/f).exists() for f in ("hash.k2d","opts.k2d","taxo.k2d")):
                kdir = str(d); break
    cfg["kraken2_db"] = kdir
    return cfg

def ensure_data_ready(data_dir=None, force=False):
    requested = Path(data_dir).expanduser() if data_dir else user_cache_dir()
    data = requested if _ensure_writable(requested) else user_cache_dir()
    if not _ensure_writable(data):
        raise SystemExit(f"[humanfilt] cannot write to '{requested}' or fallback '{data}'")

    bundles_dir = data/"bundles"
    bundles_dir.mkdir(parents=True, exist_ok=True)

    for spec in BUNDLES:
        url = f"https://zenodo.org/records/{RECORD_ID}/files/{spec['name']}?download=1"
        dst = bundles_dir/spec["name"]
        if force and dst.exists():
            dst.unlink()
        if not dst.exists():
            _curl(url, dst)
            _try_verify_with_sidecar(dst)
        _extract(dst, data)

    cfg = _scan_cfg(data)
    (data/"config.json").write_text(json.dumps(cfg, indent=2))
    print(f"[humanfilt] config -> {(data/'config.json')}", file=sys.stderr)
    return cfg

