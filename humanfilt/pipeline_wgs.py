#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple

from .utils import (
    say as _say,
    ensure_tools,
    count_reads,
    available_threads,
)

# -------------------------------
# Small local helpers
# -------------------------------

def _ls_dir(p: Path) -> str:
    """Return a simple 'ls' listing of a directory for debugging."""
    try:
        lines = []
        for q in sorted(p.iterdir()):
            try:
                sz = q.stat().st_size
                lines.append(f"{q.name}\t{sz}")
            except Exception:
                lines.append(q.name)
        return "\n".join(lines)
    except Exception as e:
        return f"<ls failed: {e}>"

def _pick_first_existing(candidates):
    """Return the first existing path from a list, else raise FileNotFoundError."""
    for p in candidates:
        if Path(p).exists():
            return Path(p)
    raise FileNotFoundError("None of the candidate files exist.")

_R1_PAT = re.compile(r"(_R)1(\.f(?:ast)?q(?:\.gz)?)$", re.IGNORECASE)

def _pair_for_r1(p: Path) -> Path:
    """Turn ..._R1.fastq(.gz) into ..._R2.fastq(.gz) robustly."""
    return Path(_R1_PAT.sub(r"\g<1>2\g<2>", str(p)))

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

def _find_inputs_single(folder: Path) -> List[Tuple[Path, str]]:
    pats = ["*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"]
    files: List[Path] = []
    for pat in pats:
        files.extend(folder.glob(pat))
    items: List[Tuple[Path, str]] = []
    for r in sorted(files):
        nm = re.sub(r"(\.f(?:ast)?q(?:\.gz)?)$", "", r.name, flags=re.IGNORECASE)
        items.append((r, nm))
    return items

# -------------------------------
# Steps
# -------------------------------

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

    tg_r1 = tmp / f"{sample}_R1_val_1.fq.gz"
    tg_r2 = tmp / f"{sample}_R2_val_2.fq.gz"
    if not tg_r1.exists():
        tg_r1 = tmp / f"{sample}_R1_val_1.fq"
    if not tg_r2.exists():
        tg_r2 = tmp / f"{sample}_R2_val_2.fq"
    return tg_r1, tg_r2

def _trim_galore_single(r: Path, sample: str, tmp: Path, t: int,
                        trim_quality: int, trim_length: int, lg) -> Path:
    cmd = [
        "trim_galore",
        "--cores", str(t),
        "--quality", str(trim_quality),
        "--length",  str(trim_length),
        "--output_dir", str(tmp),
        str(r),
    ]
    subprocess.check_call(cmd, stdout=lg, stderr=lg)

    cand = [
        tmp / f"{sample}_trimmed.fq.gz",
        tmp / f"{sample}_trimmed.fq",
    ]
    for c in cand:
        if c.exists():
            return c
    # Fallback to the newest FASTQ-ish file
    return max(list(tmp.glob("*.f*q*")), key=lambda p: p.stat().st_mtime)

def _fastuniq(dd_in_r1: Path, dd_in_r2: Path, sample: str, tmp: Path, lg) -> Tuple[Path, Path]:
    """FastUniq requires plain FASTQ; if gz, stream-decompress to tmp first."""
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
        subprocess.check_call(
            ["fastuniq", "-i", str(lst), "-t", "q", "-o", str(out1), "-p", str(out2)],
            stdout=lg, stderr=lg
        )
    finally:
        for p in cleanup:
            try: p.unlink()
            except FileNotFoundError: pass
        try: lst.unlink()
        except FileNotFoundError: pass
    return out1, out2

def _bbduk_univec(r1: Path, r2: Path, sample: str, tmp: Path, t: int, univec: str, lg) -> Tuple[Path, Path]:
    out1 = tmp / f"{sample}.vec_R1.fq"
    out2 = tmp / f"{sample}.vec_R2.fq"
    stats = tmp / f"{sample}.bbduk_stats.txt"
    subprocess.check_call(
        [
            "bbduk.sh",
            f"in1={r1}", f"in2={r2}",
            f"out1={out1}", f"out2={out2}",
            f"ref={univec}", "k=27", "hdist=1",
            f"stats={stats}", f"threads={t}",
        ],
        stdout=lg, stderr=lg
    )
    return out1, out2

def _kraken2_paired(dd_r1: Path, dd_r2: Path, sample: str, tmp: Path, t: int, db: str, lg) -> Tuple[Path, Path]:
    out_pat = tmp / f"{sample}.kr_R#.fq"
    rep = tmp / f"{sample}.kraken_report.txt"
    cmd = [
        "kraken2", "--db", db, "--threads", str(t), "--report", str(rep),
        "--paired", "--unclassified-out", str(out_pat),
        str(dd_r1), str(dd_r2),
    ]
    subprocess.check_call(cmd, stdout=lg, stderr=lg)

    cand_r1 = [
        tmp / f"{sample}.kr_R1.fq",
        tmp / f"{sample}.kr_R1.fastq",
        tmp / f"{sample}.kr_R_1.fq",
        tmp / f"{sample}.kr_R_1.fastq",
    ]
    cand_r2 = [
        tmp / f"{sample}.kr_R2.fq",
        tmp / f"{sample}.kr_R2.fastq",
        tmp / f"{sample}.kr_R_2.fq",
        tmp / f"{sample}.kr_R_2.fastq",
    ]
    try:
        r1 = _pick_first_existing(cand_r1)
        r2 = _pick_first_existing(cand_r2)
        return r1, r2
    except FileNotFoundError:
        _say(
            "[WGS] WARNING: Kraken2 did not create expected unclassified files.\n"
            f"  looked for: {cand_r1 + cand_r2}\n"
            f"  temp dir listing:\n{_ls_dir(tmp)}\n"
            "Proceeding without Kraken2 filtering (pass-through)."
        )
        return dd_r1, dd_r2

def _bwa_unmapped_pair(prefix_or_fa: str, r1: Path, r2: Path,
                       out1: Path, out2: Path, t: int, lg):
    """
    bwa mem -> both-unmapped -> fastq
    `prefix_or_fa` can be the fasta path or the BWA index prefix.
    """
    p1 = subprocess.Popen(
        ["bwa", "mem", "-t", str(t), prefix_or_fa, str(r1), str(r2)],
        stdout=subprocess.PIPE, stderr=lg, text=False
    )
    p2 = subprocess.Popen(
        ["samtools", "view", "-@", str(t), "-b", "-f", "12", "-F", "256"],
        stdin=p1.stdout, stdout=subprocess.PIPE, stderr=lg, text=False
    )
    rc1 = rc2 = 0
    try:
        subprocess.check_call(
            ["samtools", "fastq", "-@", str(t),
             "-1", str(out1), "-2", str(out2),
             "-0", "/dev/null", "-s", "/dev/null", "-n"],
            stdin=p2.stdout, stdout=lg, stderr=lg
        )
    finally:
        try:
            if p1.stdout: p1.stdout.close()
        except Exception:
            pass
        try:
            if p2.stdout: p2.stdout.close()
        except Exception:
            pass
        rc2 = p2.wait()
        rc1 = p1.wait()
    if rc2 != 0:
        raise subprocess.CalledProcessError(rc2, "samtools view")
    if rc1 != 0:
        raise subprocess.CalledProcessError(rc1, "bwa mem")

def _bowtie2_t2t_unmapped(bt2_prefix: str, r1: Path, r2: Path,
                          sample: str, tmp: Path, t: int, lg) -> Tuple[Path, Path]:
    """
    Run bowtie2 vs T2T and collect both-unmapped pairs.
    Note: bowtie2 --un-conc writes two files by appending .1 and .2 to the path provided.
    We avoid printf-style patterns here for maximum compatibility.
    """
    out_prefix = tmp / f"{sample}.un_t2t"
    cmd = [
        "bowtie2", "-p", str(t), "-x", bt2_prefix,
        "-1", str(r1), "-2", str(r2),
        "--very-sensitive",
        "--un-conc", str(out_prefix),
        "-S", "/dev/null",
    ]
    subprocess.check_call(cmd, stdout=lg, stderr=lg)

    cand1 = Path(str(out_prefix) + ".1")
    cand2 = Path(str(out_prefix) + ".2")
    if not (cand1.exists() and cand2.exists()):
        suffix = "".join(out_prefix.suffixes)
        stem = out_prefix.stem
        alt1 = out_prefix.with_name(stem + ".1" + suffix)
        alt2 = out_prefix.with_name(stem + ".2" + suffix)
        if alt1.exists() and alt2.exists():
            cand1, cand2 = alt1, alt2
        else:
            raise FileNotFoundError(
                f"Bowtie2 --un-conc outputs not found: {cand1}, {cand2} or {alt1}, {alt2}. tmp listing:\n{_ls_dir(tmp)}"
            )
    return cand1, cand2

def _minimap2_hprc_unmapped(hprc_fa: str, r1: Path, r2: Path,
                            out1: Path, out2: Path, t: int, lg):
    p1 = subprocess.Popen(
        ["minimap2", "-t", str(t), "-ax", "sr", hprc_fa, str(r1), str(r2)],
        stdout=subprocess.PIPE, stderr=lg, text=False
    )
    p2 = subprocess.Popen(
        ["samtools", "view", "-@", str(t), "-b", "-f", "12", "-F", "256"],
        stdin=p1.stdout, stdout=subprocess.PIPE, stderr=lg, text=False
    )
    rc1 = rc2 = 0
    try:
        subprocess.check_call(
            ["samtools", "fastq", "-@", str(t),
             "-1", str(out1), "-2", str(out2),
             "-0", "/dev/null", "-s", "/dev/null", "-n"],
            stdin=p2.stdout, stdout=lg, stderr=lg
        )
    finally:
        try:
            if p1.stdout: p1.stdout.close()
        except Exception:
            pass
        try:
            if p2.stdout: p2.stdout.close()
        except Exception:
            pass
        rc2 = p2.wait()
        rc1 = p1.wait()
    if rc2 != 0:
        raise subprocess.CalledProcessError(rc2, "samtools view")
    if rc1 != 0:
        raise subprocess.CalledProcessError(rc1, "minimap2")

# -------------------------------
# Public entry
# -------------------------------

def run_wgs(input_dir: str, output_dir: str, report_csv: str, threads: int,
            cfg: dict, layout: str = "paired",
            trim_quality: int = 20, trim_length: int = 20,
            keep_temp: bool = False, save_bams: bool = False):
    """
    WGS pipeline.
    - layout: 'paired' or 'single'.
    - Only the final HPRC-unmapped reads are written to output_dir.
    - Intermediates live in a temp dir and are removed automatically.
    """
    t = threads or available_threads()

    ensure_tools([
        "bwa", "samtools", "bowtie2", "minimap2",
        "trim_galore", "fastuniq", "bbduk.sh", "dedupe.sh", "kraken2", "seqkit"
    ])

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Refs from cfg
    gr_fa = cfg.get("grch38_fa", "")
    bwa_prefix = cfg.get("grch38_bwa_prefix") or gr_fa
    t2t_bt2 = cfg.get("t2t_bt2_prefix", "")
    hprc_fa = cfg.get("hprc_fa", "")
    kr_db   = cfg.get("kraken2_db", "")
    univec  = cfg.get("univec_ref", "")

    missing = []
    if not gr_fa or not Path(gr_fa).exists():
        missing.append(f"grch38_fa -> {gr_fa}")
    # Require bowtie2 T2T index (prefix)
    if not t2t_bt2 or not Path(t2t_bt2 + ".1.bt2").exists():
        missing.append(f"t2t_bt2_prefix (prefix) -> {t2t_bt2}")
    if not hprc_fa or not Path(hprc_fa).exists():
        missing.append(f"hprc_fa -> {hprc_fa}")
    if not kr_db or not Path(kr_db).exists():
        missing.append(f"kraken2_db -> {kr_db}")
    if not univec or not Path(univec).exists():
        missing.append(f"univec_ref -> {univec}")
    # Ensure BWA index exists for the given prefix (or fasta path used as prefix)
    bwa_idx_exts = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    if not bwa_prefix or not all(Path(bwa_prefix + ext).exists() for ext in bwa_idx_exts):
        missing.append(
            f"grch38_bwa_prefix index files not found for prefix '{bwa_prefix}' "
            f"(need {', '.join(bwa_idx_exts)})"
        )

    if missing:
        _say("Missing or not found:\n  - " + "\n  - ".join(missing))
        raise SystemExit(1)

    inp = Path(input_dir)
    if layout == "paired":
        items = _find_inputs_paired(inp)
        if not items:
            _say(f"No paired R1 files found in {input_dir}")
            return
    else:
        items = _find_inputs_single(inp)
        if not items:
            _say(f"No FASTQ files found in {input_dir}")
            return

    with open(report_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        if layout == "paired":
            w.writerow([
                "Sample","Initial_R1","Initial_R2","TrimGalore_R1","TrimGalore_R2",
                "Dedup_R1","Dedup_R2","BBDuk_R1","BBDuk_R2","Kraken2_R1","Kraken2_R2",
                "BWA_GRCh38_R1","BWA_GRCh38_R2","Bowtie2_T2T_R1","Bowtie2_T2T_R2",
                "Minimap2_HPRC_R1","Minimap2_HPRC_R2"
            ])
        else:
            w.writerow([
                "Sample","Initial","TrimGalore","Dedup","BBDuk","Kraken2",
                "BWA_GRCh38","Minimap2_HPRC"
            ])

        for rec in items:
            if layout == "paired":
                r1, r2, sample = rec
            else:
                r, sample = rec

            _say(f"--- WGS {layout} :: {sample} ---")
            if keep_temp:
                tmp = out_path / f"{sample}__tmp"
                tmp.mkdir(parents=True, exist_ok=True)
                _tmp_cm = None
            else:
                _tmp_cm = tempfile.TemporaryDirectory(prefix=f"hf_{sample}_")
                tmp = Path(_tmp_cm.name)
            log = tmp / f"{sample}.log"
            with open(log, "a") as lg:
                if layout == "paired":
                    init_r1 = count_reads(r1, t)
                    init_r2 = count_reads(r2, t)

                    tg1, tg2 = _trim_galore_pair(r1, r2, sample, tmp, t, trim_quality, trim_length, lg)
                    tg1_n = count_reads(tg1, t); tg2_n = count_reads(tg2, t)

                    dd1, dd2 = _fastuniq(tg1, tg2, sample, tmp, lg)
                    dd1_n = count_reads(dd1, t); dd2_n = count_reads(dd2, t)

                    vec1, vec2 = _bbduk_univec(dd1, dd2, sample, tmp, t, univec, lg)
                    vec1_n = count_reads(vec1, t); vec2_n = count_reads(vec2, t)

                    kr1, kr2 = _kraken2_paired(vec1, vec2, sample, tmp, t, kr_db, lg)
                    kr1_n = count_reads(kr1, t); kr2_n = count_reads(kr2, t)

                    _say(f"[WGS] BWA GRCh38 (both-unmapped)… :: {sample}")
                    un1 = tmp / f"{sample}.un_grch38_R1.fq"
                    un2 = tmp / f"{sample}.un_grch38_R2.fq"
                    try:
                        _bwa_unmapped_pair(bwa_prefix, kr1, kr2, un1, un2, t, lg)
                    except subprocess.CalledProcessError:
                        _say("[WGS] ERROR: BWA step failed; aborting sample.")
                        continue
                    un1_n = count_reads(un1, t); un2_n = count_reads(un2, t)

                    _say(f"[WGS] Bowtie2 T2T (both-unmapped)… :: {sample}")
                    try:
                        t2t1, t2t2 = _bowtie2_t2t_unmapped(t2t_bt2, un1, un2, sample, tmp, t, lg)
                    except Exception as e:
                        _say(f"[WGS] WARNING: Bowtie2 step failed for {sample} ({e}); passing through.")
                        t2t1, t2t2 = un1, un2
                    t2t1_n = count_reads(t2t1, t); t2t2_n = count_reads(t2t2, t)

                    _say(f"[WGS] Minimap2 HPRC (both-unmapped)… :: {sample}")
                    out1 = out_path / f"{sample}_un_hprc_R1.fq"
                    out2 = out_path / f"{sample}_un_hprc_R2.fq"
                    try:
                        _minimap2_hprc_unmapped(hprc_fa, t2t1, t2t2, out1, out2, t, lg)
                    except subprocess.CalledProcessError:
                        _say(f"[WGS] ERROR: Minimap2 step failed for {sample}; aborting sample.")
                        continue
                    out1_n = count_reads(out1, t); out2_n = count_reads(out2, t)

                    w.writerow([sample, init_r1, init_r2, tg1_n, tg2_n,
                                dd1_n, dd2_n, vec1_n, vec2_n, kr1_n, kr2_n,
                                un1_n, un2_n, t2t1_n, t2t2_n, out1_n, out2_n])

                else:
                    # Single-end full pipeline mirroring paired steps
                    init = count_reads(r, t)

                    tg = _trim_galore_single(r, sample, tmp, t, trim_quality, trim_length, lg)
                    tg_n = count_reads(tg, t)

                    # Deduplicate single-end using bbmap's dedupe.sh
                    dd = _dedupe_single(tg, sample, tmp, t, lg)
                    dd_n = count_reads(dd, t)

                    vec = _bbduk_univec_single(dd, sample, tmp, t, univec, lg)
                    vec_n = count_reads(vec, t)

                    kr = _kraken2_single(vec, sample, tmp, t, kr_db, lg)
                    kr_n = count_reads(kr, t)

                    _say(f"[WGS] BWA GRCh38 (unmapped)… :: {sample}")
                    un_gr = tmp / f"{sample}.un_grch38_bwa.fq"
                    bwa_bam = (out_path / f"{sample}.bwa.bam") if save_bams else (tmp / f"{sample}.bwa.bam")
                    try:
                        _bwa_unmapped_single(bwa_prefix, kr, un_gr, t, lg, bam_out=bwa_bam)
                    except subprocess.CalledProcessError:
                        _say(f"[WGS] ERROR: BWA step failed for {sample}; aborting sample.")
                        continue
                    un_gr_n = _count_unmapped_from_bam(bwa_bam, t, lg)
                    # Save stage FASTQ to output directory as well
                    try:
                        shutil.copyfile(un_gr, out_path / f"{sample}_un_grch38.fq")
                    except Exception:
                        pass

                    _say(f"[WGS] Minimap2 HPRC (unmapped)… :: {sample}")
                    out = out_path / f"{sample}_un_hprc.fq"
                    mm2_hprc_bam = (out_path / f"{sample}.mm2_hprc.bam") if save_bams else (tmp / f"{sample}.mm2_hprc.bam")
                    try:
                        _minimap2_unmapped_single(hprc_fa, un_gr, out, t, "sr", lg, bam_out=mm2_hprc_bam)
                    except subprocess.CalledProcessError:
                        _say(f"[WGS] ERROR: Minimap2 (HPRC) failed for {sample}; aborting sample.")
                        continue
                    out_n = _count_unmapped_from_bam(mm2_hprc_bam, t, lg)

                    w.writerow([sample, init, tg_n, dd_n, vec_n, kr_n, un_gr_n, out_n])

            if _tmp_cm is not None:
                _tmp_cm.cleanup()
    _say(f"WGS done. Summary: {report_csv}")

# -------------------------------
# Single-end helpers
# -------------------------------

def _dedupe_single(r: Path, sample: str, tmp: Path, t: int, lg) -> Path:
    out = tmp / f"{sample}.dedup.fq"
    subprocess.check_call([
        "dedupe.sh", f"in={r}", f"out={out}", f"threads={t}", "overwrite=t"
    ], stdout=lg, stderr=lg)
    return out

def _bbduk_univec_single(r: Path, sample: str, tmp: Path, t: int, univec: str, lg) -> Path:
    out = tmp / f"{sample}.vec.fq"
    stats = tmp / f"{sample}.bbduk_stats.txt"
    subprocess.check_call(
        [
            "bbduk.sh",
            f"in={r}", f"out={out}",
            f"ref={univec}", "k=27", "hdist=1",
            f"stats={stats}", f"threads={t}",
        ],
        stdout=lg, stderr=lg
    )
    return out

def _kraken2_single(r: Path, sample: str, tmp: Path, t: int, db: str, lg) -> Path:
    out = tmp / f"{sample}.kr.fq"
    rep = tmp / f"{sample}.kraken_report.txt"
    cmd = [
        "kraken2", "--db", db, "--threads", str(t), "--report", str(rep),
        "--unclassified-out", str(out), str(r)
    ]
    subprocess.check_call(cmd, stdout=lg, stderr=lg)
    return out

def _count_unmapped_from_bam(bam: Path, t: int, lg) -> int:
    """Count unmapped reads (-f 4) in a BAM using samtools view -c -f 4."""
    try:
        out = subprocess.check_output(["samtools", "view", "-@", str(t), "-c", "-f", "4", str(bam)], text=True, stderr=lg)
        s = out.strip()
        return int(s) if s.isdigit() else 0
    except Exception:
        return 0

def _bwa_unmapped_single(prefix_or_fa: str, r: Path, out: Path, t: int, lg, bam_out: Path | None = None):
    """
    Keep ALL unmapped reads (flag 0x4). Do NOT exclude QC-fail/duplicates/etc.
    bwa mem | [optional BAM save] | samtools view -f 4 | samtools fastq
    """
    bwa_cmd = ["bwa", "mem", "-t", str(t), prefix_or_fa, str(r)]

    if bam_out is not None:
        # bwa -> BAM (with header)
        p1 = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=lg, text=False)
        p2 = subprocess.Popen(
            ["samtools", "view", "-h", "-@", str(t), "-b", "-", "-o", str(bam_out)],
            stdin=p1.stdout, stderr=lg, text=False
        )
        try:
            if p1.stdout: p1.stdout.close()
        except Exception:
            pass
        rc2 = p2.wait(); rc1 = p1.wait()
        if rc2 != 0:
            raise subprocess.CalledProcessError(rc2, "samtools view (SAM->BAM)")
        if rc1 != 0:
            raise subprocess.CalledProcessError(rc1, "bwa mem")

        # BAM -> unmapped-only (ALL) -> FASTQ
        p3 = subprocess.Popen(
            ["samtools", "view", "-@", str(t), "-b", "-f", "4", str(bam_out)],
            stdout=subprocess.PIPE, stderr=lg, text=False
        )
        try:
            subprocess.check_call(
                ["samtools", "fastq", "-@", str(t), "-n", "-", "-o", str(out)],
                stdin=p3.stdout, stdout=lg, stderr=lg
            )
        finally:
            try:
                if p3.stdout: p3.stdout.close()
            except Exception:
                pass
            rc3 = p3.wait()
        if rc3 != 0:
            raise subprocess.CalledProcessError(rc3, "samtools view (filter)")

    else:
        # Streamed: bwa -> view -f 4 -> fastq
        p1 = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=lg, text=False)
        p2 = subprocess.Popen(
            ["samtools", "view", "-@", str(t), "-b", "-f", "4", "-"],
            stdin=p1.stdout, stdout=subprocess.PIPE, stderr=lg, text=False
        )
        try:
            subprocess.check_call(
                ["samtools", "fastq", "-@", str(t), "-n", "-", "-o", str(out)],
                stdin=p2.stdout, stdout=lg, stderr=lg
            )
        finally:
            try:
                if p1.stdout: p1.stdout.close()
            except Exception:
                pass
            try:
                if p2.stdout: p2.stdout.close()
            except Exception:
                pass
            rc2 = p2.wait(); rc1 = p1.wait()
        if rc2 != 0:
            raise subprocess.CalledProcessError(rc2, "samtools view (filter)")
        if rc1 != 0:
            raise subprocess.CalledProcessError(rc1, "bwa mem")

def _bowtie2_t2t_unmapped_single(bt2_prefix: str, r: Path, sample: str, tmp: Path, t: int, lg, bam_out: Path | None = None) -> Path:
    """
    Run Bowtie2 vs T2T:
    - If bam_out is provided: first write a full BAM for provenance.
    - Then robustly collect unmapped reads via --un to FASTQ.
    """
    out = tmp / f"{sample}.un_t2t.fq"

    # Save BAM if requested
    if bam_out is not None:
        p1 = subprocess.Popen(["bowtie2", "-p", str(t), "-x", bt2_prefix, "-U", str(r), "--very-sensitive", "-S", "-"],
                              stdout=subprocess.PIPE, stderr=lg, text=False)
        p2 = subprocess.Popen(["samtools", "view", "-h", "-@", str(t), "-b", "-", "-o", str(bam_out)],
                              stdin=p1.stdout, stderr=lg, text=False)
        try:
            if p1.stdout:
                p1.stdout.close()
        except Exception:
            pass
        rc2 = p2.wait(); rc1 = p1.wait()
        if rc2 != 0:
            raise subprocess.CalledProcessError(rc2, "samtools view (SAM->BAM)")
        if rc1 != 0:
            raise subprocess.CalledProcessError(rc1, "bowtie2")

    # Robust unmapped FASTQ via --un
    cmd = [
        "bowtie2", "-p", str(t), "-x", bt2_prefix,
        "-U", str(r), "--very-sensitive", "--un", str(out),
        "-S", "/dev/null",
    ]
    subprocess.check_call(cmd, stdout=lg, stderr=lg)
    if not out.exists():
        raise FileNotFoundError(f"Bowtie2 --un output not found: {out}")
    return out

def _minimap2_unmapped_single(hprc_fa: str, r: Path, out: Path, t: int, preset: str, lg, bam_out: Path | None = None):
    """
    Keep ALL unmapped (flag 0x4) from minimap2 output. No -F mask.
    """
    p1 = subprocess.Popen(["minimap2", "-t", str(t), "-ax", "sr", hprc_fa, str(r)],
                          stdout=subprocess.PIPE, stderr=lg, text=False)
    if bam_out is not None:
        # Save BAM
        p_save = subprocess.Popen(["samtools", "view", "-@", str(t), "-b", "-", "-o", str(bam_out)],
                                  stdin=p1.stdout, stderr=lg, text=False)
        try:
            if p1.stdout: p1.stdout.close()
        except Exception:
            pass
        rc_save = p_save.wait(); rc_mm2 = p1.wait()
        if rc_save != 0:
            raise subprocess.CalledProcessError(rc_save, "samtools view (SAM->BAM)")
        if rc_mm2 != 0:
            raise subprocess.CalledProcessError(rc_mm2, "minimap2")

        # Count (optional)
        try:
            cnt = subprocess.check_output(["samtools", "view", "-@", str(t), "-c", "-f", "4", str(bam_out)], text=True).strip()
            try:
                lg.write(f"minimap2_unmapped_count\t{cnt}\n"); lg.flush()
            except Exception:
                pass
        except Exception:
            pass

        # Filter unmapped -> FASTQ
        p2 = subprocess.Popen(["samtools", "view", "-@", str(t), "-b", "-f", "4", str(bam_out)],
                              stdout=subprocess.PIPE, stderr=lg, text=False)
        try:
            subprocess.check_call(["samtools", "fastq", "-@", str(t), "-n", "-", "-o", str(out)],
                                  stdin=p2.stdout, stdout=lg, stderr=lg)
        finally:
            try:
                if p2.stdout: p2.stdout.close()
            except Exception:
                pass
            rc2 = p2.wait()
        if rc2 != 0:
            raise subprocess.CalledProcessError(rc2, "samtools view (filter)")
    else:
        # Stream: minimap2 | view -f 4 | fastq
        p2 = subprocess.Popen(["samtools", "view", "-@", str(t), "-b", "-f", "4", "-"],
                              stdin=p1.stdout, stdout=subprocess.PIPE, stderr=lg, text=False)
        rc1 = rc2 = 0
        try:
            subprocess.check_call(["samtools", "fastq", "-@", str(t), "-n", "-", "-o", str(out)],
                                  stdin=p2.stdout, stdout=lg, stderr=lg)
        finally:
            try:
                if p1.stdout: p1.stdout.close()
            except Exception:
                pass
            try:
                if p2.stdout: p2.stdout.close()
            except Exception:
                pass
            rc2 = p2.wait(); rc1 = p1.wait()
        if rc2 != 0:
            raise subprocess.CalledProcessError(rc2, "samtools view")
        if rc1 != 0:
            raise subprocess.CalledProcessError(rc1, "minimap2")
