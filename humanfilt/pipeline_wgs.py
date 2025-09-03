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
