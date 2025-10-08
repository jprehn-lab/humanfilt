#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple

from humanfilt.utils import (
    say as _say,
    ensure_tools,
    count_reads,
    available_threads,
)


def _pair_for_r1(p: Path) -> Path:
    return Path(re.sub(r"(_R)1(\.f(?:ast)?q(?:\.gz)?)$", r"\g<1>2\g<2>", str(p), flags=re.IGNORECASE))


def _find_inputs_paired(folder: Path) -> List[Tuple[Path, Path, str]]:
    pats = ["*_R1.fq", "*_R1.fq.gz", "*_R1.fastq", "*_R1.fastq.gz"]
    r1s: List[Path] = []
    for pat in pats:
        r1s.extend(folder.glob(pat))
    seen = set()
    items: List[Tuple[Path, Path, str]] = []
    for r1 in sorted(r1s):
        r2 = _pair_for_r1(r1)
        if not r2.exists():
            continue
        sample = re.sub(r"(_R)1(\.f(?:ast)?q(?:\.gz)?)$", "", r1.name, flags=re.IGNORECASE)
        key = (str(r1), str(r2))
        if key in seen:
            continue
        seen.add(key)
        items.append((r1, r2, sample))
    return items



def _safe_unlink(p: Path) -> None:
    try:
        p.unlink()
    except Exception:
        pass


def _trim_galore_pair(r1: Path, r2: Path, sample: str, tmp: Path, t: int,
                      trim_quality: int, trim_length: int, lg) -> Tuple[Path, Path]:
    cmd = [
        "trim_galore", "--paired",
        "--cores", str(t),
        "--quality", str(trim_quality),
        "--length",  str(trim_length),
        "--output_dir", str(tmp),
        str(r1), str(r2),
    ]
    subprocess.check_call(cmd, stdout=lg, stderr=lg)

    # Derive outputs by scanning the per-sample temp directory for Trim Galore outputs
    # This is robust regardless of input basenames (original or Kraken-filtered).
    c1 = sorted(list(tmp.glob("*_val_1.fq.gz")) + list(tmp.glob("*_val_1.fq")), key=lambda p: p.stat().st_mtime, reverse=True)
    c2 = sorted(list(tmp.glob("*_val_2.fq.gz")) + list(tmp.glob("*_val_2.fq")), key=lambda p: p.stat().st_mtime, reverse=True)
    if not c1 or not c2:
        raise FileNotFoundError(f"Trim Galore outputs not found in {tmp}")
    return c1[0], c2[0]



def _kraken2_paired(inp_r1: Path, inp_r2: Path, sample: str, tmp: Path, t: int, db: str, lg) -> Tuple[Path, Path]:
    out_pat = tmp / f"{sample}.kr_R#.fq"
    rep = tmp / f"{sample}.kraken_report.txt"
    cmd = [
        "kraken2", "--db", db, "--threads", str(t), "--report", str(rep),
        "--paired", "--unclassified-out", str(out_pat),
        str(inp_r1), str(inp_r2),
    ]
    subprocess.check_call(cmd, stdout=lg, stderr=lg)
    cand_r1 = [tmp / f"{sample}.kr_R1.fq", tmp / f"{sample}.kr_R1.fastq", tmp / f"{sample}.kr_R_1.fq", tmp / f"{sample}.kr_R_1.fastq"]
    cand_r2 = [tmp / f"{sample}.kr_R2.fq", tmp / f"{sample}.kr_R2.fastq", tmp / f"{sample}.kr_R_2.fq", tmp / f"{sample}.kr_R_2.fastq"]
    r1 = next((p for p in cand_r1 if p.exists()), None)
    r2 = next((p for p in cand_r2 if p.exists()), None)
    if not (r1 and r2):
        raise FileNotFoundError(f"Kraken2 did not create expected unclassified outputs; looked for: {cand_r1 + cand_r2}")
    return r1, r2



def _fastuniq_pair(dd_in_r1: Path, dd_in_r2: Path, sample: str, tmp: Path, lg) -> Tuple[Path, Path]:
    fu_r1, fu_r2 = dd_in_r1, dd_in_r2
    cleanup = []
    if str(dd_in_r1).endswith(".gz"):
        fu_r1 = tmp / f"{sample}.R1.fq"
        fu_r2 = tmp / f"{sample}.R2.fq"
        subprocess.check_call(["gunzip", "-c", str(dd_in_r1)], stdout=open(fu_r1, "wb"), stderr=lg)
        subprocess.check_call(["gunzip", "-c", str(dd_in_r2)], stdout=open(fu_r2, "wb"), stderr=lg)
        cleanup.extend([fu_r1, fu_r2])
    lst = tmp / "pairs.txt"
    lst.write_text(f"{fu_r1}\n{fu_r2}\n")
    out1 = tmp / f"{sample}.dedup_R1.fq"
    out2 = tmp / f"{sample}.dedup_R2.fq"
    try:
        subprocess.check_call(["fastuniq", "-i", str(lst), "-t", "q", "-o", str(out1), "-p", str(out2)], stdout=lg, stderr=lg)
    finally:
        for p in cleanup:
            _safe_unlink(p)
        _safe_unlink(lst)
    return out1, out2



def _bbduk_univec(r1: Path, r2: Path, sample: str, tmp: Path, t: int, univec: str, lg) -> Tuple[Path, Path]:
    out1 = tmp / f"{sample}.vec_R1.fq"
    out2 = tmp / f"{sample}.vec_R2.fq"
    stats = tmp / f"{sample}.bbduk_stats.txt"
    subprocess.check_call([
        "bbduk.sh",
        f"in1={r1}", f"in2={r2}",
        f"out1={out1}", f"out2={out2}",
        f"ref={univec}", "k=27", "hdist=1",
        f"stats={stats}", f"threads={t}",
    ], stdout=lg, stderr=lg)
    return out1, out2



def _bwa_unmapped_pair(prefix_or_fa: str, r1: Path, r2: Path,
                       out1: Path, out2: Path, t: int, lg):
    p1 = subprocess.Popen(["bwa", "mem", "-t", str(t), prefix_or_fa, str(r1), str(r2)], stdout=subprocess.PIPE, stderr=lg, text=False)
    p2 = subprocess.Popen(["samtools", "view", "-@", str(t), "-b", "-f", "12", "-F", "256"], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=lg, text=False)
    try:
        subprocess.check_call(["samtools", "fastq", "-@", str(t), "-1", str(out1), "-2", str(out2), "-0", "/dev/null", "-s", "/dev/null", "-n"], stdin=p2.stdout, stdout=lg, stderr=lg)
    finally:
        try:
            if p1.stdout:
                p1.stdout.close()
        except Exception:
            pass
        try:
            if p2.stdout:
                p2.stdout.close()
        except Exception:
            pass
        rc2 = p2.wait(); rc1 = p1.wait()
    if rc2 != 0:
        raise subprocess.CalledProcessError(rc2, "samtools view")
    if rc1 != 0:
        raise subprocess.CalledProcessError(rc1, "bwa mem")


def _bowtie2_t2t_unmapped(bt2_prefix: str, r1: Path, r2: Path, sample: str, tmp: Path, t: int, lg) -> Tuple[Path, Path]:
    out_prefix = tmp / f"{sample}.un_t2t"
    cmd = ["bowtie2", "-p", str(t), "-x", bt2_prefix, "-1", str(r1), "-2", str(r2), "--very-sensitive", "--un-conc", str(out_prefix), "-S", "/dev/null"]
    subprocess.check_call(cmd, stdout=lg, stderr=lg)
    cand1 = Path(str(out_prefix) + ".1"); cand2 = Path(str(out_prefix) + ".2")
    if not (cand1.exists() and cand2.exists()):
        suffix = "".join(out_prefix.suffixes)
        stem = out_prefix.stem
        alt1 = out_prefix.with_name(stem + ".1" + suffix)
        alt2 = out_prefix.with_name(stem + ".2" + suffix)
        if alt1.exists() and alt2.exists():
            cand1, cand2 = alt1, alt2
        else:
            raise FileNotFoundError(f"Bowtie2 --un-conc outputs not found: {cand1}, {cand2} or {alt1}, {alt2}.")
    return cand1, cand2


def _minimap2_hprc_unmapped(hprc_fa: str, r1: Path, r2: Path, out1: Path, out2: Path, t: int, lg):
    p1 = subprocess.Popen(["minimap2", "-t", str(t), "-ax", "sr", hprc_fa, str(r1), str(r2)], stdout=subprocess.PIPE, stderr=lg, text=False)
    p2 = subprocess.Popen(["samtools", "view", "-@", str(t), "-b", "-f", "12", "-F", "256"], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=lg, text=False)
    try:
        subprocess.check_call(["samtools", "fastq", "-@", str(t), "-1", str(out1), "-2", str(out2), "-0", "/dev/null", "-s", "/dev/null", "-n"], stdin=p2.stdout, stdout=lg, stderr=lg)
    finally:
        try:
            if p1.stdout:
                p1.stdout.close()
        except Exception:
            pass
        try:
            if p2.stdout:
                p2.stdout.close()
        except Exception:
            pass
        rc2 = p2.wait(); rc1 = p1.wait()
    if rc2 != 0:
        raise subprocess.CalledProcessError(rc2, "samtools view")
    if rc1 != 0:
        raise subprocess.CalledProcessError(rc1, "minimap2")


def run_wgs(input_dir: str, output_dir: str, report_csv: str, threads: int,
            cfg: dict,
            trim_quality: int = 20, trim_length: int = 20,
            keep_temp: bool = False):
    t = threads or available_threads()
    ensure_tools(["bwa", "samtools", "bowtie2", "minimap2", "trim_galore", "fastuniq", "bbduk.sh", "dedupe.sh", "kraken2", "seqkit"])

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    gr_fa = cfg.get("grch38_fa", "")
    bwa_prefix = cfg.get("grch38_bwa_prefix") or gr_fa
    t2t_bt2 = cfg.get("t2t_bt2_prefix", "")
    hprc_fa = cfg.get("hprc_fa", "")
    kr_db   = cfg.get("kraken2_db", "")
    univec  = cfg.get("univec_ref", "")

    missing = []
    if not gr_fa or not Path(gr_fa).exists():
        missing.append(f"grch38_fa -> {gr_fa}")
    if not t2t_bt2 or not Path(t2t_bt2 + ".1.bt2").exists():
        missing.append(f"t2t_bt2_prefix (prefix) -> {t2t_bt2}")
    if not hprc_fa or not Path(hprc_fa).exists():
        missing.append(f"hprc_fa -> {hprc_fa}")
    if not kr_db or not Path(kr_db).exists():
        missing.append(f"kraken2_db -> {kr_db}")
    if not univec or not Path(univec).exists():
        missing.append(f"univec_ref -> {univec}")
    if missing:
        _say("Missing or not found:\n  - " + "\n  - ".join(missing))
        raise SystemExit(1)

    inp = Path(input_dir)
    items = _find_inputs_paired(inp)
    if not items:
        _say(f"No paired FASTQ files found in {input_dir}")
        return

    with open(report_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Sample","Initial_R1","Initial_R2","Kraken_R1","Kraken_R2","TrimGalore_R1","TrimGalore_R2","Dedup_R1","Dedup_R2","BBDuk_R1","BBDuk_R2","BWA_GRCh38_R1","BWA_GRCh38_R2","Bowtie2_T2T_R1","Bowtie2_T2T_R2","Minimap2_HPRC_R1","Minimap2_HPRC_R2"])
        fh.flush(); os.fsync(fh.fileno())

        for rec in items:
            r1, r2, sample = rec

            _say(f"--- WGS paired :: {sample} ---")
            if keep_temp:
                tmp = out_path / f"{sample}__tmp"
                tmp.mkdir(parents=True, exist_ok=True)
                _tmp_cm = None
            else:
                _tmp_cm = tempfile.TemporaryDirectory(prefix=f"hf_{sample}_", dir=str(out_path))
                tmp = Path(_tmp_cm.name)
            log = tmp / f"{sample}.log"
            with open(log, "a") as lg:
                if True:
                    init_r1 = count_reads(r1, t)
                    init_r2 = count_reads(r2, t)

                    # Kraken2 first to drop clear human reads early
                    kr1, kr2 = _kraken2_paired(r1, r2, sample, tmp, t, kr_db, lg)
                    kr1_n = count_reads(kr1, t); kr2_n = count_reads(kr2, t)

                    # Trim after Kraken2 on the remaining reads
                    tg1, tg2 = _trim_galore_pair(kr1, kr2, sample, tmp, t, trim_quality, trim_length, lg)
                    tg1_n = count_reads(tg1, t); tg2_n = count_reads(tg2, t)
                    _safe_unlink(kr1); _safe_unlink(kr2)

                    dd1, dd2 = _fastuniq_pair(tg1, tg2, sample, tmp, lg)
                    dd1_n = count_reads(dd1, t); dd2_n = count_reads(dd2, t)
                    _safe_unlink(tg1); _safe_unlink(tg2)

                    vec1, vec2 = _bbduk_univec(dd1, dd2, sample, tmp, t, univec, lg)
                    vec1_n = count_reads(vec1, t); vec2_n = count_reads(vec2, t)
                    _safe_unlink(dd1); _safe_unlink(dd2)

                    un1 = tmp / f"{sample}.un_grch38_R1.fq"
                    un2 = tmp / f"{sample}.un_grch38_R2.fq"
                    try:
                        _bwa_unmapped_pair(bwa_prefix, vec1, vec2, un1, un2, t, lg)
                    except subprocess.CalledProcessError:
                        _say("[WGS] ERROR: BWA step failed; aborting sample.")
                        if not keep_temp:
                            try:
                                _tmp_cm.cleanup()  # type: ignore
                            except Exception:
                                pass
                        continue
                    un1_n = count_reads(un1, t); un2_n = count_reads(un2, t)
                    _safe_unlink(vec1); _safe_unlink(vec2)

                    try:
                        t2t1, t2t2 = _bowtie2_t2t_unmapped(t2t_bt2, un1, un2, sample, tmp, t, lg)
                    except Exception as e:
                        _say(f"[WGS] WARNING: Bowtie2 step failed for {sample} ({e}); passing through.")
                        t2t1, t2t2 = un1, un2
                    t2t1_n = count_reads(t2t1, t); t2t2_n = count_reads(t2t2, t)
                    if t2t1 is not un1:
                        _safe_unlink(un1); _safe_unlink(un2)

                    out1 = out_path / f"{sample}_un_hprc_R1.fq"
                    out2 = out_path / f"{sample}_un_hprc_R2.fq"
                    try:
                        _minimap2_hprc_unmapped(hprc_fa, t2t1, t2t2, out1, out2, t, lg)
                    except subprocess.CalledProcessError:
                        _say(f"[WGS] ERROR: Minimap2 step failed for {sample}; aborting sample.")
                        continue
                    out1_n = count_reads(out1, t); out2_n = count_reads(out2, t)
                    if t2t1 is not un1:
                        _safe_unlink(t2t1); _safe_unlink(t2t2)

                    w.writerow([sample, init_r1, init_r2, kr1_n, kr2_n, tg1_n, tg2_n,
                                dd1_n, dd2_n, vec1_n, vec2_n,
                                un1_n, un2_n, t2t1_n, t2t2_n, out1_n, out2_n])
                    fh.flush(); os.fsync(fh.fileno())

                

            if _tmp_cm is not None:
                _tmp_cm.cleanup()
    _say(f"WGS done. Summary: {report_csv}")
