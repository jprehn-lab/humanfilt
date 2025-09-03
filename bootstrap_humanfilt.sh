#!/usr/bin/env bash
set -euo pipefail

say(){ printf "• %s\n" "$*"; }

mkdir -p humanfilt recipe tests

########################################
# Python package
########################################
say "writing humanfilt/__init__.py"
cat > humanfilt/__init__.py <<'PYEOF'
__all__ = []
PYEOF

say "writing humanfilt/utils.py"
cat > humanfilt/utils.py <<'PYEOF'
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
PYEOF

say "writing humanfilt/downloader.py"
cat > humanfilt/downloader.py <<'PYEOF'
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
        subprocess.check_ call(["tar","-I","zstd","-xf",s,"-C",str(into)])
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
PYEOF

say "writing humanfilt/pipeline_wgs.py"
cat > humanfilt/pipeline_wgs.py <<'PYEOF'
from pathlib import Path
import subprocess, shutil, tempfile
from .utils import nreads_fastq, ensure_tools

def _trim_galore_paired(out, r1, r2, threads, q, l):
    subprocess.check_call([
        "trim_galore","--paired","--cores",str(threads),
        "--quality", str(q), "--length", str(l),
        "--output_dir", str(out),
        str(r1), str(r2)
    ])
    # canonical names from Trim Galore
    base = r1.name.replace("_R1.fastq.gz","").replace("_R1.fastq","")
    tr_r1 = out/f"{base}_R1_val_1.fq.gz"
    tr_r2 = out/f"{base}_R2_val_2.fq.gz"
    if not tr_r1.exists(): tr_r1 = tr_r1.with_suffix("")  # .fq
    if not tr_r2.exists(): tr_r2 = tr_r2.with_suffix("")
    return tr_r1, tr_r2, base

def _trim_galore_single(out, r, threads, q, l):
    subprocess.check_call([
        "trim_galore","--cores",str(threads),
        "--quality", str(q), "--length", str(l),
        "--output_dir", str(out),
        str(r)
    ])
    base = r.name.replace(".fastq.gz","").replace(".fastq","")
    tr = out/f"{base}_trimmed.fq.gz"
    if not tr.exists(): tr = tr.with_suffix("")  # .fq
    if not tr.exists():
        cand = list(out.glob(f"{base}*trimmed.fq*"))
        if not cand: raise SystemExit(f"No trimmed file for {r}")
        tr = cand[0]
    return tr, base

def run_wgs(input_dir, output_dir, report_csv, threads, cfg, *,
            layout: str = "paired", trim_quality: int = 20, trim_length: int = 20):
    """
    layout: 'paired' or 'single'
    """
    ensure_tools(
        required=["bwa","bowtie2","minimap2","samtools","kraken2","trim_galore","bbduk.sh"],
        anyof=["fastuniq","dedupe.sh"]
    )

    inp = Path(input_dir); out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    rep = Path(report_csv)

    if layout == "paired":
        if not rep.exists():
            rep.write_text(
                "Sample,Initial_R1,Initial_R2,Trimmed_R1,Trimmed_R2,VectorClean_R1,VectorClean_R2,"
                "Dedup_R1,Dedup_R2,Kraken2_R1,Kraken2_R2,BWA_R1,BWA_R2,T2T_R1,T2T_R2,HPRC_R1,HPRC_R2\n"
            )
        pairs = {}
        for r1 in sorted(list(inp.glob("*_R1.fastq.gz")) + list(inp.glob("*_R1.fastq"))):
            r2 = Path(str(r1).replace("_R1.fastq.gz","_R2.fastq.gz").replace("_R1.fastq","_R2.fastq"))
            if r2.exists():
                sample = r1.name.replace("_R1.fastq.gz","").replace("_R1.fastq","")
                pairs[sample] = (r1, r2)

        for sample, (r1, r2) in pairs.items():
            print(f"[WGS-PE] {sample}")
            init_r1, init_r2 = nreads_fastq(r1), nreads_fastq(r2)

            tr_r1, tr_r2, base = _trim_galore_paired(out, r1, r2, threads, trim_quality, trim_length)
            trc_r1, trc_r2 = nreads_fastq(tr_r1), nreads_fastq(tr_r2)

            vec_r1 = out/f"{base}_vec_R1.fq";  vec_r2 = out/f"{base}_vec_R2.fq"
            bb_stats = out/f"{base}_bbduk_stats.txt"
            subprocess.check_call([
                "bbduk.sh", f"in1={tr_r1}", f"in2={tr_r2}",
                f"out1={vec_r1}", f"out2={vec_r2}",
                f"ref={cfg['univec_ref']}", "k=27", "hdist=1",
                f"stats={bb_stats}", f"threads={threads}"
            ])
            vcc_r1, vcc_r2 = nreads_fastq(vec_r1), nreads_fastq(vec_r2)

            dd_r1 = out/f"{base}_dedup_R1.fq"; dd_r2 = out/f"{base}_dedup_R2.fq"
            if shutil.which("fastuniq"):
                import tempfile
                with tempfile.TemporaryDirectory() as td:
                    pairlist = Path(td)/"pairs.txt"
                    pairlist.write_text(f"{vec_r1}\n{vec_r2}\n")
                    subprocess.check_call(["fastuniq","-i",str(pairlist),"-t","q","-o",str(dd_r1),"-p",str(dd_r2)])
            else:
                subprocess.check_call(["dedupe.sh", f"in1={vec_r1}", f"in2={vec_r2}", f"out1={dd_r1}", f"out2={dd_r2}",
                                       "ac=f", f"threads={threads}"])
            ddc_r1, ddc_r2 = nreads_fastq(dd_r1), nreads_fastq(dd_r2)

            kr_un_r1 = out/f"{base}_kraken2_R_1.fq"
            kr_un_r2 = out/f"{base}_kraken2_R_2.fq"
            kr_report = out/f"{base}_kraken2_report.txt"
            subprocess.check_call([
                "kraken2","--db",cfg["kraken2_db"],"--threads",str(threads),
                "--report",str(kr_report),"--paired",
                "--unclassified-out", str(out/f"{base}_kraken2_R#.fq"),
                str(dd_r1), str(dd_r2)
            ])
            krc_r1, krc_r2 = nreads_fastq(kr_un_r1), nreads_fastq(kr_un_r2)

            bwa_un_r1 = out/f"{base}_un_grch38_bwa_R1.fq"
            bwa_un_r2 = out/f"{base}_un_grch38_bwa_R2.fq"
            p1 = subprocess.Popen(["bwa","mem","-t",str(threads), cfg["grch38_bwa_prefix"], str(kr_un_r1), str(kr_un_r2)],
                                  stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["samtools","view","-b","-f","12","-F","256"], stdin=p1.stdout, stdout=subprocess.PIPE)
            subprocess.check_call(["samtools","fastq","-1",str(bwa_un_r1),"-2",str(bwa_un_r2),
                                   "-0","/dev/null","-s","/dev/null","-n"], stdin=p2.stdout)
            bwc_r1, bwc_r2 = nreads_fastq(bwa_un_r1), nreads_fastq(bwa_un_r2)

            subprocess.check_call([
                "bowtie2","-p",str(threads), "-x", cfg["t2t_bt2_prefix"],
                "-1",str(bwa_un_r1),"-2",str(bwa_un_r2),
                "--very-sensitive","--un-conc",str(out/f"{base}_un_t2t_R%.fq"), "-S","/dev/null"
            ])
            t2t_un_r1 = out/f"{base}_un_t2t_R1.fq"
            t2t_un_r2 = out/f"{base}_un_t2t_R2.fq"
            t2c_r1, t2c_r2 = nreads_fastq(t2t_un_r1), nreads_fastq(t2t_un_r2)

            hprc_un_r1 = out/f"{base}_un_hprc_R1.fq"
            hprc_un_r2 = out/f"{base}_un_hprc_R2.fq"
            p3 = subprocess.Popen(["minimap2","-t",str(threads),"-ax","sr",(cfg.get("hprc_mmi") or cfg["hprc_fa"]),
                                   str(t2t_un_r1), str(t2t_un_r2)], stdout=subprocess.PIPE)
            p4 = subprocess.Popen(["samtools","view","-b","-f","12","-F","256"], stdin=p3.stdout, stdout=subprocess.PIPE)
            subprocess.check_call(["samtools","fastq","-1",str(hprc_un_r1),"-2",str(hprc_un_r2),
                                   "-0","/dev/null","-s","/dev/null","-n"], stdin=p4.stdout)
            hpc_r1, hpc_r2 = nreads_fastq(hprc_un_r1), nreads_fastq(hprc_un_r2)

            with open(rep, "a") as rfh:
                rfh.write(",".join(map(str, [
                  sample, init_r1, init_r2, trc_r1, trc_r2, vcc_r1, vcc_r2, ddc_r1, ddc_r2,
                  krc_r1, krc_r2, bwc_r1, bwc_r2, t2c_r1, t2c_r2, hpc_r1, hpc_r2
                ])) + "\n")

            for f in [kr_un_r1, kr_un_r2, t2t_un_r1, t2t_un_r2, bwa_un_r1, bwa_un_r2, bb_stats, kr_report,
                      vec_r1, vec_r2, dd_r1, dd_r2]:
                Path(f).unlink(missing_ok=True)

    else:
        if not rep.exists():
            rep.write_text("Sample,Initial,Trimmed,VectorClean,Dedup,Kraken2,BWA,T2T,HPRC\n")
        singles = [f for f in sorted(list(Path(input_dir).glob("*.fastq*"))) if "_R2" not in f.name]
        for r in singles:
            base = r.name.replace(".fastq.gz","").replace(".fastq","").replace("_R1","")
            print(f"[WGS-SE] {base}")
            init = nreads_fastq(r)

            tr, base = _trim_galore_single(Path(output_dir), r, threads, trim_quality, trim_length)
            trc = nreads_fastq(tr)

            vec = Path(output_dir)/f"{base}_vec.fq"
            bb_stats = Path(output_dir)/f"{base}_bbduk_stats.txt"
            subprocess.check_call([
                "bbduk.sh", f"in={tr}", f"out={vec}",
                f"ref={cfg['univec_ref']}", "k=27", "hdist=1",
                f"stats={bb_stats}", f"threads={threads}"
            ])
            vcc = nreads_fastq(vec)

            dd = Path(output_dir)/f"{base}_dedup.fq"
            subprocess.check_call(["dedupe.sh", f"in={vec}", f"out={dd}", "ac=f", f"threads={threads}"])
            ddc = nreads_fastq(dd)

            kr_un = Path(output_dir)/f"{base}_kraken2.fastq"
            kr_rep = Path(output_dir)/f"{base}_kraken2_report.txt"
            subprocess.check_call([
                "kraken2","--db",cfg["kraken2_db"],"--threads",str(threads),
                "--report",str(kr_rep),"--unclassified-out",str(kr_un), str(dd)
            ])
            krc = nreads_fastq(kr_un)

            bwa_un = Path(output_dir)/f"{base}_un_grch38_bwa.fq"
            p1 = subprocess.Popen(["bwa","mem","-t",str(threads), cfg["grch38_bwa_prefix"], str(kr_un)],
                                  stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["samtools","view","-b","-F","256"], stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(["samtools","view","-b","-f","4"], stdin=p2.stdout, stdout=subprocess.PIPE)
            with open(bwa_un, "wb") as fout:
                subprocess.check_call(["samtools","fastq","-"], stdin=p3.stdout, stdout=fout)
            bwc = nreads_fastq(bwa_un)

            t2t_un = Path(output_dir)/f"{base}_un_t2t.fq"
            subprocess.check_call([
                "bowtie2","-p",str(threads), "-x", cfg["t2t_bt2_prefix"],
                "-U", str(bwa_un),
                "--very-sensitive","--un", str(t2t_un), "-S","/dev/null"
            ])
            t2c = nreads_fastq(t2t_un)

            hprc_un = Path(output_dir)/f"{base}_un_hprc.fq"
            p4 = subprocess.Popen(["minimap2","-t",str(threads),"-ax","sr",(cfg.get("hprc_mmi") or cfg["hprc_fa"]),
                                   str(t2t_un)], stdout=subprocess.PIPE)
            p5 = subprocess.Popen(["samtools","view","-b","-F","256"], stdin=p4.stdout, stdout=subprocess.PIPE)
            p6 = subprocess.Popen(["samtools","view","-b","-f","4"], stdin=p5.stdout, stdout=subprocess.PIPE)
            with open(hprc_un, "wb") as fout:
                subprocess.check_call(["samtools","fastq","-"], stdin=p6.stdout, stdout=fout)
            hpc = nreads_fastq(hprc_un)

            with open(rep, "a") as rfh:
                rfh.write(f"{base},{init},{trc},{vcc},{ddc},{krc},{bwc},{t2c},{hpc}\n")

            for f in [kr_un, kr_rep, vec, dd, bb_stats, bwa_un, t2t_un]:
                Path(f).unlink(missing_ok=True)
PYEOF

say "writing humanfilt/pipeline_rnaseq.py"
cat > humanfilt/pipeline_rnaseq.py <<'PYEOF'
from pathlib import Path
import re, subprocess
from .utils import nreads_fastq, ensure_tools

def _sample_from_name(p: Path) -> str:
    base = p.name
    base = re.sub(r"\.gcap[^.]*", "", base)
    base = base.replace(".R1.fastq.gz","").replace(".R1.fastq","").replace(".fastq.gz","").replace(".fastq","")
    return base

def run_rnaseq(input_dir, output_dir, report_csv, threads, cfg, *,
               rna_preset="illumina", pattern=None, trim_quality: int = 20, trim_length: int = 20):
    """RNA-seq is single-end in this pipeline."""
    ensure_tools(required=["minimap2","samtools","kraken2","trim_galore","fastqc","bbduk.sh","dedupe.sh"])
    inp = Path(input_dir); out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    report = Path(report_csv)
    if not report.exists():
        report.write_text("File,Initial_fastq,Trimmed_fastq,VectorClean,Dedup,Kraken2,Gencode,GRCh38,T2T,HPRC\n")
    patt = pattern or "*.R1.fastq.gz"
    files = sorted(inp.glob(patt)) or sorted(inp.glob(patt.replace(".gz","")))
    if not files:
        files = [f for f in sorted(list(inp.glob("*.fastq*"))) if "_R2" not in f.name]
    mm_preset = "map-ont" if rna_preset == "ont" else "sr"

    for r1 in files:
        sample = _sample_from_name(r1)
        print(f"[RNA] {sample} ({rna_preset})")
        init = nreads_fastq(r1)

        subprocess.check_call([
            "trim_galore","--fastqc","--cores",str(threads),
            "--quality", str(trim_quality), "--length", str(trim_length),
            "--output_dir", str(out), str(r1)
        ])
        tr = out/f"{sample}_trimmed.fq.gz"
        if not tr.exists(): tr = out/f"{sample}_trimmed.fq"
        if not tr.exists():
            cand = list(out.glob(f"{sample}*trimmed.fq*"))
            if not cand:
                print(f"[{sample}] No trimmed file; skipping."); continue
            tr = cand[0]
        trc = nreads_fastq(tr)

        vec = out/f"{sample}_vec.fq"
        bb_stats = out/f"{sample}_bbduk_stats.txt"
        subprocess.check_call([
            "bbduk.sh", f"in={tr}", f"out={vec}", f"ref={cfg['univec_ref']}",
            "k=27", "hdist=1", f"stats={bb_stats}", f"threads={threads}"
        ])
        vcc = nreads_fastq(vec)

        dd = out/f"{sample}_dedup.fq"
        subprocess.check_call(["dedupe.sh", f"in={vec}", f"out={dd}", "ac=f", f"threads={threads}"])
        ddc = nreads_fastq(dd)

        kr_un = out/f"{sample}_kraken2.fastq"
        kr_rep = out/f"{sample}_kraken2_report.txt"
        subprocess.check_call([
            "kraken2","--db",cfg["kraken2_db"],"--threads",str(threads),
            "--report",str(kr_rep),"--unclassified-out",str(kr_un), str(dd)
        ])
        krc = nreads_fastq(kr_un)
        current = kr_un

        def _filter(tag, ref_or_mmi, curr):
            un = out/f"{sample}_un_{tag}.fastq"
            p1 = subprocess.Popen(["minimap2","-t",str(threads),"-ax",mm_preset, ref_or_mmi, str(curr)], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["samtools","view","-b","-F","256"], stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(["samtools","view","-b","-f","4"], stdin=p2.stdout, stdout=subprocess.PIPE)
            with open(un, "wb") as fout:
                subprocess.check_call(["samtools","fastq","-"], stdin=p3.stdout, stdout=fout)
            return un, nreads_fastq(un)

        gencode_ref = cfg.get("gencode_tx_mmi") or cfg["gencode_tx_fa"]
        current, c_gencode = _filter("gencode", gencode_ref, current)
        grch38_ref = cfg.get("grch38_mmi") or cfg["grch38_fa"]
        current, c_grch38 = _filter("grch38", grch38_ref, current)
        t2t_ref = cfg.get("t2t_mmi") or cfg["t2t_fa"]
        current, c_t2t = _filter("t2t", t2t_ref, current)
        hprc_ref = cfg.get("hprc_mmi") or cfg["hprc_fa"]
        current, c_hprc = _filter("hprc", hprc_ref, current)

        with open(report, "a") as r:
            r.write(f"{sample},{init},{trc},{vcc},{ddc},{krc},{c_gencode},{c_grch38},{c_t2t},{c_hprc}\n")

        for f in [kr_un, kr_rep, vec, dd, bb_stats]:
            Path(f).unlink(missing_ok=True)
PYEOF

say "writing humanfilt/cli.py"
cat > humanfilt/cli.py <<'PYEOF'
import argparse, json
from pathlib import Path
from .downloader import ensure_data_ready
from .pipeline_wgs import run_wgs
from .pipeline_rnaseq import run_rnaseq

def _parser():
    p = argparse.ArgumentParser(
        prog="humanfilt",
        description="Human read filtering (WGS paired/single & RNA-seq single) with first-run Zenodo downloads",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("setup", help="Download/refresh references & DBs from Zenodo")
    s.add_argument("--data-dir", default=None, help="Cache directory (default: ~/.local/share/humanfilt)")
    s.add_argument("--force", action="store_true", help="Re-download even if present")

    r = sub.add_parser("run", help="Run the pipeline")
    r.add_argument("--mode", choices=["wgs","rnaseq"], required=True)
    r.add_argument("--input", required=True)
    r.add_argument("--output", required=True)
    r.add_argument("--report", required=True)
    r.add_argument("--threads", type=int, default=8)
    r.add_argument("--data-dir", default=None)
    r.add_argument("--kraken2-db", default=None, help="Override Kraken2 DB path")

    # Trim Galore controls (requested)
    r.add_argument("--trim-quality", type=int, default=20,
                   help="Trim Galore --quality <INT>. Default: 20.")
    r.add_argument("--trim-length", type=int, default=20,
                   help="Trim Galore --length <INT>. Default: 20 (0 disables read discard).")

    # WGS layout (paired or single). RNA-seq is always single.
    r.add_argument("--wgs-layout", choices=["paired","single"], default="paired",
                   help="Only used when --mode wgs (default: paired).")

    # RNA-seq specifics
    r.add_argument("--rna-preset", choices=["illumina","ont"], default="illumina",
                   help="Minimap2 preset for RNA stage (illumina=sr, ont=map-ont)")
    r.add_argument("--pattern", default=None,
                   help="RNA input pattern (default: *.R1.fastq.gz; fallback: any *.fastq* not matching _R2).")
    return p

def main():
    args = _parser().parse_args()

    if args.cmd == "setup":
        cfg = ensure_data_ready(args.data_dir, force=args.force)
        print(json.dumps(cfg, indent=2))
        return

    cfg = ensure_data_ready(args.data_dir, force=False)
    if args.kraken2_db:
        cfg["kraken2_db"] = args.kraken2_db

    Path(args.output).mkdir(parents=True, exist_ok=True)

    if args.mode == "wgs":
        run_wgs(args.input, args.output, args.report, args.threads, cfg,
                layout=args.wgs_layout, trim_quality=args.trim_quality, trim_length=args.trim_length)
    else:
        run_rnaseq(args.input, args.output, args.report, args.threads, cfg,
                   rna_preset=args.rna_preset, pattern=args.pattern,
                   trim_quality=args.trim_quality, trim_length=args.trim_length)

if __name__ == "__main__":
    main()
PYEOF

########################################
# Build system & recipe
########################################
say "writing pyproject.toml"
cat > pyproject.toml <<'TOMLEOF'
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "humanfilt"
version = "0.6.0"
description = "WGS paired/single & RNA-seq single filtering with first-run Zenodo download (GRCh38/T2T/HPRC/GENCODE/UniVec + Kraken2)."
readme = "README.md"
requires-python = ">=3.9"
license = {text = "MIT"}
authors = [{name="Maria S. Frolova"}]

[project.scripts]
humanfilt = "humanfilt.cli:main"

[tool.setuptools.packages.find]
where = ["."]
include = ["humanfilt*"]
TOMLEOF

say "writing recipe/meta.yaml"
cat > recipe/meta.yaml <<'YAMEOF'
package:
  name: humanfilt
  version: "0.6.0"

source:
  path: ..

build:
  noarch: python
  script: bash recipe/build.sh
  entry_points:
    - humanfilt = humanfilt.cli:main

requirements:
  host:
    - python >=3.9
    - pip
  run:
    - python >=3.9
    - bwa
    - bowtie2
    - minimap2
    - samtools
    - kraken2
    - trim-galore
    - fastqc
    - bbmap
    - fastuniq
    - pigz
    - curl
    - coreutils
    - tar
    - zstd

test:
  commands:
    - humanfilt --help

about:
  home: https://github.com/jprehn-lab/humanfilt
  license: MIT
  summary: "WGS paired/single & RNA-seq single decontamination with first-run Zenodo download (record 17020482)."
YAMEOF

say "writing recipe/build.sh"
cat > recipe/build.sh <<'BASHEOF'
#!/usr/bin/env bash
set -eux
$PYTHON -m pip install --no-deps --no-build-isolation .
BASHEOF
chmod +x recipe/build.sh

########################################
# Tests, README, gitignore, license
########################################
say "writing tests/test_cli.sh"
cat > tests/test_cli.sh <<'BASHEOF'
#!/usr/bin/env bash
set -eux
humanfilt --help
BASHEOF
chmod +x tests/test_cli.sh

say "writing README.md"
cat > README.md <<'MDEOF'
# humanfilt

Conda-installable CLI for human read decontamination.

- Tools: bwa, bowtie2, minimap2, samtools, kraken2, trim-galore, fastqc, bbmap, fastuniq
- First-run download from **Zenodo record `17020482`** (split tar.zst bundles)
- Vector removal (UniVec via **BBDuk**) and **dedup** after Trim Galore
- **WGS (paired)**: Trim → UniVec → Dedup → Kraken2 → BWA(GRCh38) → Bowtie2(T2T) → minimap2(HPRC)
- **WGS (single)**: Trim → UniVec → Dedup → Kraken2 → BWA(GRCh38) → Bowtie2(T2T) → minimap2(HPRC)
- **RNA-seq (single)**: Trim → UniVec → Dedup → Kraken2 → minimap2(GENCODE → GRCh38 → T2T → HPRC)

## Trim Galore controls
- `--trim-quality <INT>` (default **20**)  
  Trim low-quality ends (BWA algorithm).
- `--trim-length <INT>` (default **20**, **0** disables)  
  Discard reads shorter than INT after trimming (quality/adapter).

## Quickstart
```bash
pip install -e .

export HUMANFILT_ZENODO_RECORD=17020482
humanfilt setup                 # downloads + extracts into ~/.local/share/humanfilt by default

# WGS paired-end
humanfilt run --mode wgs --wgs-layout paired \
  --input /path/in --output /path/out --report /path/out/report.csv --threads 16 \
  --trim-quality 20 --trim-length 20

# WGS single-end
humanfilt run --mode wgs --wgs-layout single \
  --input /path/in --output /path/out --report /path/out/report.csv --threads 16

# RNA-seq single-end
humanfilt run --mode rnaseq \
  --input /path/in --output /path/out --report /path/out/report.csv --threads 16 \
  --rna-preset illumina --trim-quality 20 --trim-length 20

