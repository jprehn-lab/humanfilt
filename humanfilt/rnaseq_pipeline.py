#!/usr/bin/env python3
from __future__ import annotations

import os
import csv
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from .utils import say as _say, ensure_tools, available_threads


def _count_reads_fastq(path: Path, threads: int) -> int:
    try:
        # prefer seqkit
        cmd = ["seqkit", "stats", "-T", "-j", str(threads), str(path)]
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        lines = [ln for ln in res.stdout.strip().splitlines() if ln.strip()]
        if len(lines) >= 2:
            return int(lines[1].split("\t")[3])
    except Exception:
        pass
    try:
        # fallback: gzip -cd | awk
        if str(path).endswith(".gz"):
            p1 = subprocess.Popen(["gzip", "-cd", str(path)], stdout=subprocess.PIPE)
        else:
            p1 = subprocess.Popen(["cat", str(path)], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["awk", 'NR%4==1{c++} END{print c+0}'], stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
        out, _ = p2.communicate()
        p1.wait()
        return int(out.strip()) if out and out.strip().isdigit() else 0
    except Exception:
        return 0


def _find_trimmed_single(tmp: Path, base: str) -> Path:
    cand = [tmp / f"{base}_trimmed.fq.gz", tmp / f"{base}_trimmed.fq"]
    for c in cand:
        if c.exists():
            return c
    hits = list(tmp.glob("*.f*q*"))
    if not hits:
        raise FileNotFoundError("Trim Galore output not found")
    return sorted(hits, key=lambda p: p.stat().st_mtime, reverse=True)[0]


def run_rnaseq(input_dir: str, output_dir: str, report_csv: str, threads: Optional[int], cfg: dict,
               trim_quality: int = 20, trim_length: int = 20, save_bams: bool = False):
    """
    RNA-seq (single-end) pipeline (HISAT/Bowtie2/T2T/HPRC), modeled on rna_hisat.sh.
    Writes only final HPRC-unmapped FASTQs to output_dir; intermediates/logs in tmp are removed.
    Requires env overrides for HISAT2 and Bowtie2 GENCODE indices, and T2T MMI if not present in config.
    """
    t = threads or available_threads()
    ensure_tools(["trim_galore", "bbduk.sh", "kraken2", "hisat2", "bowtie2", "minimap2", "samtools", "gzip"])

    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)

    # Resolve references
    kraken_db = cfg.get("kraken2_db", "")
    univec = cfg.get("univec_ref", "")
    hisat2_prefix = os.environ.get("HUMANFILT_HISAT2_PREFIX") or os.environ.get("HISAT2_PREFIX", "")
    bt2_gencode_prefix = os.environ.get("HUMANFILT_BT2_GENCODE_PREFIX") or os.environ.get("BT2_GENCODE_PREFIX", "")
    t2t_mmi = os.environ.get("HUMANFILT_T2T_MMI") or os.environ.get("T2T_MMI", "")
    hprc_mmi = cfg.get("hprc_mmi", "")

    missing = []
    if not kraken_db or not Path(kraken_db).exists():
        missing.append(f"kraken2_db -> {kraken_db}")
    if not univec or not Path(univec).exists():
        missing.append(f"univec_ref -> {univec}")
    if not hisat2_prefix:
        missing.append("HISAT2_PREFIX (env) not set")
    if not bt2_gencode_prefix:
        missing.append("BT2_GENCODE_PREFIX (env) not set")
    if not t2t_mmi or not Path(t2t_mmi).exists():
        missing.append(f"T2T_MMI (env) -> {t2t_mmi}")
    if not hprc_mmi or not Path(hprc_mmi).exists():
        missing.append(f"hprc_mmi -> {hprc_mmi}")
    if missing:
        _say("Missing or not found:\n  - " + "\n  - ".join(missing))
        raise SystemExit(1)

    reads = sorted(list(Path(input_dir).glob("*.fastq")) +
                   list(Path(input_dir).glob("*.fastq.gz")) +
                   list(Path(input_dir).glob("*.fq")) +
                   list(Path(input_dir).glob("*.fq.gz")))
    if not reads:
        _say(f"No FASTQ files found in {input_dir}")
        return

    with open(report_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File","Initial","TrimGalore","UniVec","Kraken2_unclassified","HISAT2_unmapped","Bowtie2_unmapped","T2T_unmapped","HPRC_unmapped"])

        for inp in reads:
            base = inp.name
            base = base.split(".fastq")[0].split(".fq")[0]
            _say(f"==== RNA :: {base} ====")
            with tempfile.TemporaryDirectory(prefix=f"hf_rna_{base}_") as tmpdir:
                tmp = Path(tmpdir)

                init_n = _count_reads_fastq(inp, t)
                current = str(inp)

                # Trim Galore
                _say("  [1/7] Trim Galore…")
                trim_outdir = tmp / "trim"
                trim_outdir.mkdir(parents=True, exist_ok=True)
                with open(tmp / f"{base}.trimgalore.log", "w") as lg:
                    subprocess.check_call(["trim_galore", "--cores", str(t), "--quality", str(trim_quality),
                                           "--length", str(trim_length), "--output_dir", str(trim_outdir), current], stdout=lg, stderr=lg)
                trimmed = _find_trimmed_single(trim_outdir, base)
                tg_out = tmp / f"{base}.trimmed.fq.gz"
                if str(trimmed).endswith(".gz"):
                    subprocess.check_call(["cp", "-f", str(trimmed), str(tg_out)])
                else:
                    with open(tg_out, "wb") as fhgz:
                        subprocess.check_call(["gzip", "-c", str(trimmed)], stdout=fhgz)
                tg_n = _count_reads_fastq(tg_out, t)
                current = str(tg_out)

                # BBDuk UniVec
                _say("  [2/7] BBduk UniVec…")
                univec_out = tmp / f"{base}.univec.fq.gz"
                with open(tmp / f"{base}.univec.bbduk.txt", "w") as lg:
                    subprocess.check_call(["bbduk.sh", f"in={current}", f"out={univec_out}", f"ref={univec}",
                                           "k=27", "hdist=1", f"threads={t}", f"stats={tmp}/{base}.univec.bbduk.txt", "overwrite=t"],
                                          stdout=lg, stderr=lg)
                univec_n = _count_reads_fastq(univec_out, t)
                current = str(univec_out)

                # Kraken2 keep UNCLASSIFIED
                _say("  [3/7] Kraken2 (keep UNCLASSIFIED)…")
                kr_tmp = tmp / f"{base}.kr_unclassified.fq"
                with open(tmp / f"{base}.kraken2.report.txt", "w") as rep, open(tmp / f"{base}.kraken2.classification.txt", "w") as cls:
                    subprocess.check_call(["kraken2", "--db", kraken_db, "--threads", str(t), "--report", str(rep.name),
                                           "--unclassified-out", str(kr_tmp), current], stdout=cls, stderr=subprocess.DEVNULL)
                kr_un = tmp / f"{base}.kr_unclassified.fq.gz"
                with open(kr_un, "wb") as fhgz:
                    if kr_tmp.exists() and kr_tmp.stat().st_size > 0:
                        subprocess.check_call(["gzip", "-c", str(kr_tmp)], stdout=fhgz)
                    else:
                        fhgz.write(b"")
                kr_n = _count_reads_fastq(kr_un, t)
                current = str(kr_un)

                # HISAT2 keep UNMAPPED
                _say("  [4/7] HISAT2 (keep UNMAPPED)…")
                hisat_un = tmp / f"{base}.hisat_unmapped.fq.gz"
                with open(tmp / f"{base}.hisat2.log", "w") as lg:
                    subprocess.check_call(["hisat2", "-p", str(t), "-x", hisat2_prefix, "-U", current,
                                           "--no-unal", "--un-gz", str(hisat_un), "-S", "/dev/null"], stderr=lg, stdout=lg)
                hisat_n = _count_reads_fastq(hisat_un, t)
                current = str(hisat_un)

                # Bowtie2 vs GENCODE keep UNMAPPED
                _say("  [5/7] Bowtie2 (GENCODE; keep UNMAPPED)…")
                bt2_un = tmp / f"{base}.bt2_unmapped.fq.gz"
                with open(tmp / f"{base}.bowtie2.log", "w") as lg:
                    subprocess.check_call(["bowtie2", "-p", str(t), "-x", bt2_gencode_prefix, "-U", current,
                                           "--very-sensitive", "--un-gz", str(bt2_un), "-S", "/dev/null"], stderr=lg, stdout=lg)
                bt2_n = _count_reads_fastq(bt2_un, t)
                current = str(bt2_un)

                # minimap2 vs T2T keep UNMAPPED
                _say("  [6/7] minimap2 vs T2T (keep UNMAPPED)…")
                t2t_un = tmp / f"{base}.t2t_unmapped.fq.gz"
                with open(tmp / f"{base}.t2t.mm2.log", "w") as lg, open(t2t_un, "wb") as fhgz:
                    p1 = subprocess.Popen(["minimap2", "-t", str(t), "-ax", "sr", t2t_mmi, current], stdout=subprocess.PIPE, stderr=lg, text=False)
                    p2 = subprocess.Popen(["samtools", "view", "-@", str(t), "-b", "-f", "4"], stdin=p1.stdout, stdout=subprocess.PIPE)
                    p3 = subprocess.Popen(["samtools", "fastq", "-@", str(t), "-"], stdin=p2.stdout, stdout=subprocess.PIPE)
                    subprocess.check_call(["gzip"], stdin=p3.stdout, stdout=fhgz)
                    p1.wait(); p2.wait(); p3.wait()
                t2t_n = _count_reads_fastq(t2t_un, t)
                current = str(t2t_un)

                # minimap2 vs HPRC keep UNMAPPED (final output in output_dir)
                _say("  [7/7] minimap2 vs HPRC (keep UNMAPPED)…")
                hprc_un = Path(output_dir) / f"{base}.hprc_unmapped.fq.gz"
                with open(tmp / f"{base}.hprc.mm2.log", "w") as lg, open(hprc_un, "wb") as fhgz:
                    p1 = subprocess.Popen(["minimap2", "-t", str(t), "-ax", "sr", hprc_mmi, current], stdout=subprocess.PIPE, stderr=lg, text=False)
                    p2 = subprocess.Popen(["samtools", "view", "-@", str(t), "-b", "-f", "4"], stdin=p1.stdout, stdout=subprocess.PIPE)
                    p3 = subprocess.Popen(["samtools", "fastq", "-@", str(t), "-"], stdin=p2.stdout, stdout=subprocess.PIPE)
                    subprocess.check_call(["gzip"], stdin=p3.stdout, stdout=fhgz)
                    p1.wait(); p2.wait(); p3.wait()
                hprc_n = _count_reads_fastq(hprc_un, t)

                w.writerow([base, init_n, tg_n, univec_n, kr_n, hisat_n, bt2_n, t2t_n, hprc_n])

