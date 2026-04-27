#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import re
import subprocess
import tempfile
from contextlib import suppress
from pathlib import Path
from typing import List, Optional, Tuple

from humanfilt.utils import (
    say as _say,
    ensure_tools,
    count_reads,
    available_threads,
)


_FASTQ_EXT_RE = re.compile(r"(?i)\.f(?:ast)?q(?:\.gz)?$")
_R_READ_RE = re.compile(r"(?i)^(?P<sample>.+?)(?P<token>[_\.-]r(?P<read>[12])(?:_001)?)$")
_NUM_READ_RE = re.compile(r"(?i)^(?P<sample>.+?)(?P<token>[_\.-](?P<read>[12]))$")


def _parse_paired_fastq_name(p: Path) -> Optional[Tuple[str, int, str]]:
    m = _FASTQ_EXT_RE.search(p.name)
    if not m:
        return None
    stem = p.name[:m.start()]
    for rx, family in ((_R_READ_RE, "r"), (_NUM_READ_RE, "n")):
        mm = rx.match(stem)
        if mm:
            sample = mm.group("sample").rstrip("._-")
            if sample:
                return sample, int(mm.group("read")), family
    return None


def _find_inputs_paired(folder: Path) -> List[Tuple[Path, Path, str]]:
    files = sorted(
        p for p in folder.iterdir()
        if p.is_file() and _FASTQ_EXT_RE.search(p.name)
    )
    parsed = []
    for p in files:
        info = _parse_paired_fastq_name(p)
        if info is not None:
            sample, read, family = info
            parsed.append((sample, read, family, p))

    grouped = {}
    for sample, read, family, path in parsed:
        grp = grouped.setdefault(sample, {"families": set(), 1: [], 2: []})
        grp["families"].add(family)
        grp[read].append(path)

    issues: List[str] = []
    items: List[Tuple[Path, Path, str]] = []
    for sample in sorted(grouped):
        grp = grouped[sample]
        families = grp["families"]
        r1s = sorted(grp[1])
        r2s = sorted(grp[2])
        if len(families) > 1:
            issues.append(
                f"{sample}: mixed naming conventions detected ({', '.join(sorted(families))}); files: "
                + ", ".join(p.name for p in sorted(r1s + r2s))
            )
            continue
        if len(r1s) != 1 or len(r2s) != 1:
            parts = []
            if len(r1s) != 1:
                parts.append(f"R1 candidates={len(r1s)}")
            if len(r2s) != 1:
                parts.append(f"R2 candidates={len(r2s)}")
            issues.append(
                f"{sample}: could not infer a unique pair ({', '.join(parts)}); files: "
                + ", ".join(p.name for p in sorted(r1s + r2s))
            )
            continue
        items.append((r1s[0], r2s[0], sample))

    if issues:
        raise ValueError(
            "Input pairing failed.\n"
            "Supported naming includes *_R1/*_R2, *_R1_001/*_R2_001, *_1/*_2, and *.1/*.2.\n  - "
            + "\n  - ".join(issues)
        )
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
        "--dont_gzip",
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
    out1 = tmp / f"{sample}.trim_R1.fq"
    out2 = tmp / f"{sample}.trim_R2.fq"
    if out1.exists():
        out1.unlink()
    if out2.exists():
        out2.unlink()
    c1[0].rename(out1)
    c2[0].rename(out2)
    return out1, out2



def _kraken2_paired(inp_r1: Path, inp_r2: Path, sample: str, tmp: Path, t: int, db: str, lg) -> Tuple[Path, Path]:
    out_pat = tmp / f"{sample}.kr_R#.fq"
    rep = tmp / f"{sample}.kraken_report.txt"
    cmd = [
        "kraken2", "--db", db, "--threads", str(t), "--report", str(rep),
        "--paired", "--unclassified-out", str(out_pat),
        str(inp_r1), str(inp_r2),
    ]
    subprocess.check_call(cmd, stdout=lg, stderr=lg)
    std_r1 = tmp / f"{sample}.kr_R1.fq"
    std_r2 = tmp / f"{sample}.kr_R2.fq"
    legacy_pairs = [
        (tmp / f"{sample}.kr_R1.fastq", std_r1),
        (tmp / f"{sample}.kr_R_1.fq", std_r1),
        (tmp / f"{sample}.kr_R_1.fastq", std_r1),
        (tmp / f"{sample}.kr_R2.fastq", std_r2),
        (tmp / f"{sample}.kr_R_2.fq", std_r2),
        (tmp / f"{sample}.kr_R_2.fastq", std_r2),
    ]

    if std_r1.exists() and std_r2.exists():
        r1, r2 = std_r1, std_r2
    else:
        # Normalize any alternative names Kraken created
        for src, dst in legacy_pairs:
            if src.exists():
                if dst.exists():
                    dst.unlink()
                src.rename(dst)

        r1 = std_r1 if std_r1.exists() else None
        r2 = std_r2 if std_r2.exists() else None

    if not (r1 and r2):
        # Kraken2 omits the --unclassified-out files when every read is classified.
        # Emit empty placeholders so downstream steps still run.
        std_r1.touch(exist_ok=True)
        std_r2.touch(exist_ok=True)
        _say(f"[Kraken2] no unclassified reads for {sample}; created empty FASTQs")
        return std_r1, std_r2

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
    fq1 = tmp / f"{sample}.un_t2t_R1.fq"
    fq2 = tmp / f"{sample}.un_t2t_R2.fq"
    if fq1.exists():
        fq1.unlink()
    if fq2.exists():
        fq2.unlink()
    cand1.rename(fq1)
    cand2.rename(fq2)
    return fq1, fq2


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
            save_intermediates: bool = False,
            save_logs: bool = False):
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
    try:
        items = _find_inputs_paired(inp)
    except ValueError as e:
        _say(str(e))
        raise SystemExit(1)
    if not items:
        _say(f"No paired FASTQ files found in {input_dir}")
        return

    tmp_root = out_path / ("intermediates" if save_intermediates else "tmp")
    tmp_root.mkdir(parents=True, exist_ok=True)

    def _stage_unlink(*paths: Path) -> None:
        if save_intermediates:
            return
        for p in paths:
            _safe_unlink(p)

    with open(report_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Sample","Initial_R1","Initial_R2","Kraken_R1","Kraken_R2","TrimGalore_R1","TrimGalore_R2","Dedup_R1","Dedup_R2","BBDuk_R1","BBDuk_R2","BWA_GRCh38_R1","BWA_GRCh38_R2","Bowtie2_T2T_R1","Bowtie2_T2T_R2","Minimap2_HPRC_R1","Minimap2_HPRC_R2"])
        fh.flush(); os.fsync(fh.fileno())

        for rec in items:
            r1, r2, sample = rec

            _say(f"--- WGS paired :: {sample} ---")
            if save_intermediates:
                tmp = tmp_root / f"{sample}__tmp"
                tmp.mkdir(parents=True, exist_ok=True)
                _tmp_cm = None
            else:
                _tmp_cm = tempfile.TemporaryDirectory(prefix=f"hf_{sample}_", dir=str(tmp_root))
                tmp = Path(_tmp_cm.name)
            log_handle = None
            if save_logs:
                log = tmp / f"{sample}.log"
                log_handle = open(log, "a")
            lg = log_handle
            try:
                init_r1 = count_reads(r1, t)
                init_r2 = count_reads(r2, t)

                # Kraken2 first to drop clear human reads early
                kr1, kr2 = _kraken2_paired(r1, r2, sample, tmp, t, kr_db, lg)
                kr1_n = count_reads(kr1, t); kr2_n = count_reads(kr2, t)

                # Trim after Kraken2 on the remaining reads
                tg1, tg2 = _trim_galore_pair(kr1, kr2, sample, tmp, t, trim_quality, trim_length, lg)
                tg1_n = count_reads(tg1, t); tg2_n = count_reads(tg2, t)
                _stage_unlink(kr1, kr2)

                dd1, dd2 = _fastuniq_pair(tg1, tg2, sample, tmp, lg)
                dd1_n = count_reads(dd1, t); dd2_n = count_reads(dd2, t)
                _stage_unlink(tg1, tg2)

                vec1, vec2 = _bbduk_univec(dd1, dd2, sample, tmp, t, univec, lg)
                vec1_n = count_reads(vec1, t); vec2_n = count_reads(vec2, t)
                _stage_unlink(dd1, dd2)

                un1 = tmp / f"{sample}.un_grch38_R1.fq"
                un2 = tmp / f"{sample}.un_grch38_R2.fq"
                try:
                    _bwa_unmapped_pair(bwa_prefix, vec1, vec2, un1, un2, t, lg)
                except subprocess.CalledProcessError:
                    _say("[WGS] ERROR: BWA step failed; aborting sample.")
                    if _tmp_cm is not None:
                        with suppress(Exception):
                            _tmp_cm.cleanup()
                    continue
                un1_n = count_reads(un1, t); un2_n = count_reads(un2, t)
                _stage_unlink(vec1, vec2)

                try:
                    t2t1, t2t2 = _bowtie2_t2t_unmapped(t2t_bt2, un1, un2, sample, tmp, t, lg)
                except Exception as e:
                    _say(f"[WGS] WARNING: Bowtie2 step failed for {sample} ({e}); passing through.")
                    t2t1, t2t2 = un1, un2
                t2t1_n = count_reads(t2t1, t); t2t2_n = count_reads(t2t2, t)
                if t2t1 is not un1:
                    _stage_unlink(un1, un2)

                out1 = out_path / f"{sample}_un_hprc_R1.fq"
                out2 = out_path / f"{sample}_un_hprc_R2.fq"
                try:
                    _minimap2_hprc_unmapped(hprc_fa, t2t1, t2t2, out1, out2, t, lg)
                except subprocess.CalledProcessError:
                    _say(f"[WGS] ERROR: Minimap2 step failed for {sample}; aborting sample.")
                    continue
                out1_n = count_reads(out1, t); out2_n = count_reads(out2, t)
                if t2t1 is not un1:
                    _stage_unlink(t2t1, t2t2)

                w.writerow([sample, init_r1, init_r2, kr1_n, kr2_n, tg1_n, tg2_n,
                            dd1_n, dd2_n, vec1_n, vec2_n,
                            un1_n, un2_n, t2t1_n, t2t2_n, out1_n, out2_n])
                fh.flush(); os.fsync(fh.fileno())
            finally:
                if log_handle:
                    log_handle.close()

            if _tmp_cm is not None:
                _tmp_cm.cleanup()

    if not save_intermediates:
        with suppress(Exception):
            tmp_root.rmdir()
    _say(f"WGS done. Summary: {report_csv}")
