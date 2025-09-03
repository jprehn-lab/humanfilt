import os, json, subprocess, hashlib, shutil
from pathlib import Path

# Default: your published Zenodo record with all bundles
DEFAULT_RECORD = os.environ.get("HUMANFILT_ZENODO_RECORD", "17020482")
RECORDS = {
    "refs":   os.environ.get("HUMANFILT_ZENODO_RECORD_REFS", DEFAULT_RECORD),
    "kraken": os.environ.get("HUMANFILT_ZENODO_RECORD_KRAKEN", DEFAULT_RECORD),
}

DEFAULT_DATA = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share")) / "humanfilt"
CFG_NAME = "config.json"

# Files present in record 17020482 (leave sha256 "" to skip verify)
FILES = [
    {"record":"refs",   "name":"grch38.tar.zst",                         "dest":"bundles/grch38.tar.zst",                         "sha256":""},
    {"record":"refs",   "name":"t2t.tar.zst",                            "dest":"bundles/t2t.tar.zst",                            "sha256":""},
    {"record":"refs",   "name":"hprc.tar.zst",                           "dest":"bundles/hprc.tar.zst",                           "sha256":""},
    {"record":"refs",   "name":"gencode.tar.zst",                        "dest":"bundles/gencode.tar.zst",                        "sha256":""},
    {"record":"refs",   "name":"univec.tar.zst",                         "dest":"bundles/univec.tar.zst",                         "sha256":""},
    {"record":"kraken", "name":"humanfilt-kraken2-human-202409.tar.zst", "dest":"bundles/humanfilt-kraken2-human-202409.tar.zst", "sha256":""},
]

def _curl(url: str, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    subprocess.check_call(["curl","-L","--fail","-o",str(tmp), url])
    tmp.rename(dest)

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""): h.update(chunk)
    return h.hexdigest()

def _extract_tar(archive: Path, into: Path):
    into.mkdir(parents=True, exist_ok=True)
    s = str(archive)
    if s.endswith(".tar.zst"):
        subprocess.check_call(["tar","-I","zstd","-xf",s,"-C",str(into)])
    elif s.endswith(".tar.gz") or s.endswith(".tgz"):
        subprocess.check_call(["tar","-zxf",s,"-C",str(into)])
    elif s.endswith(".tar"):
        subprocess.check_call(["tar","-xf",s,"-C",str(into)])
    else:
        raise RuntimeError(f"Unsupported archive type: {archive}")

def _ensure_tree_after_extract(data: Path):
    # normalize: data/refs/{grch38,t2t,hprc,gencode,univec}
    for comp in ["grch38","t2t","hprc","gencode","univec"]:
        top = data/comp
        target = data/"refs"/comp
        if top.exists() and not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(top), str(target))
    # ensure kraken2/human_*
    k2 = data/"kraken2"
    if not k2.exists():
        for p in data.glob("human_*"):
            k2.mkdir(parents=True, exist_ok=True)
            shutil.move(str(p), str(k2/p.name))
            break

def _find_under(base: Path, globs):
    for g in globs:
        for p in base.glob(g):
            return p
    return None

def _bt2_prefix_from_file(bt2_file: Path) -> str:
    suf = [".rev.2.bt2",".rev.1.bt2",".4.bt2",".3.bt2",".2.bt2",".1.bt2"]
    name = bt2_file.name
    for s in suf:
        if name.endswith(s):
            return str(bt2_file.with_name(name[:-len(s)]))
    return str(bt2_file)

def ensure_data_ready(data_dir=None, force=False):
    data = Path(data_dir) if data_dir else DEFAULT_DATA
    data.mkdir(parents=True, exist_ok=True)

    cfg_path = data/CFG_NAME
    if cfg_path.exists() and not force:
        return json.loads(cfg_path.read_text())

    # download + extract
    for spec in FILES:
        rec_id = RECORDS[spec["record"]]
        url = f"https://zenodo.org/records/{rec_id}/files/{spec['name']}?download=1"
        dest = data/spec["dest"]
        if force and dest.exists():
            dest.unlink()
        if not dest.exists():
            _curl(url, dest)
        if spec.get("sha256"):
            got = _sha256(dest)
            if got != spec["sha256"]:
                dest.unlink(missing_ok=True)
                raise RuntimeError(f"Checksum mismatch for {dest}")
        _extract_tar(dest, data)

    _ensure_tree_after_extract(data)

    # resolve config
    grch = data/"refs"/"grch38"
    grch_fa = _find_under(grch, ["*.fa","*.fna","*.fa.gz","*.fna.gz"])
    amb = _find_under(grch/"bwa", ["*.amb"]) or _find_under(grch, ["*.amb"])
    grch_bwa_prefix = str(amb.with_suffix("")) if amb else (str(grch_fa) if grch_fa else "")

    t2t = data/"refs"/"t2t"
    t2t_fa = _find_under(t2t, ["*.fa","*.fna","*.fa.gz","*.fna.gz"])
    bt2_any = _find_under(t2t/"bowtie2", ["*.bt2"]) or _find_under(t2t, ["*.bt2"])
    t2t_bt2_prefix = _bt2_prefix_from_file(bt2_any) if bt2_any else ""

    hprc = data/"refs"/"hprc"
    hprc_fa = _find_under(hprc, ["*.fa","*.fna","*.fa.gz","*.fna.gz"])
    hprc_mmi = _find_under(hprc, ["*.mmi"])

    gencode = data/"refs"/"gencode"
    gencode_fa = _find_under(gencode, ["*transcripts*.fa","*transcripts*.fa.gz","*.fa","*.fa.gz"])
    gencode_mmi = _find_under(gencode, ["*.mmi"])

    univec = data/"refs"/"univec"
    uni = _find_under(univec, ["UniVec_Core*","*.fa","*.fa.gz"])

    kroot = data/"kraken2"
    kdb = None
    cand = list(kroot.glob("human_*"))
    if cand: kdb = cand[0]

    cfg = {
        "data_dir":          str(data.resolve()),
        "univec_ref":        str(uni) if uni else "",
        "kraken2_db":        str(kdb) if kdb else "",
        "grch38_fa":         str(grch_fa) if grch_fa else "",
        "grch38_bwa_prefix": grch_bwa_prefix,
        "t2t_fa":            str(t2t_fa) if t2t_fa else "",
        "t2t_bt2_prefix":    t2t_bt2_prefix,
        "hprc_fa":           str(hprc_fa) if hprc_fa else "",
        "hprc_mmi":          str(hprc_mmi) if hprc_mmi else "",
        "gencode_tx_fa":     str(gencode_fa) if gencode_fa else "",
        "gencode_tx_mmi":    str(gencode_mmi) if gencode_mmi else "",
    }
    (data/CFG_NAME).write_text(json.dumps(cfg, indent=2))
    return cfg
