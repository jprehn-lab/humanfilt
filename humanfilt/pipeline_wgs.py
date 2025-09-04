from pathlib import Path
import subprocess, tempfile, shutil
from .utils import ensure_tools, nreads_fastq

def run_wgs(input_dir, output_dir, report_csv, threads, cfg,
            layout="paired", trim_q=20, trim_len=20, kraken2_override=None):
    """WGS pipeline.
    layout: 'paired' or 'single'
    """
    ensure_tools(
        required=["trim_galore","bbduk.sh","kraken2","bwa","bowtie2","minimap2","samtools","fastqc","seqkit"],
        anyof=["fastuniq","dedupe.sh"]
    )

    inp = Path(input_dir); out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    report = Path(report_csv)
    if layout == "paired":
        header = ("Sample,Initial_R1,Initial_R2,Trim_R1,Trim_R2,Dedup_R1,Dedup_R2,"
                  "BBDuk_R1,BBDuk_R2,Kraken2_R1,Kraken2_R2,BWA_R1,BWA_R2,"
                  "T2T_R1,T2T_R2,HPRC_R1,HPRC_R2\n")
    else:
        header = "Sample,Initial,Trim,Dedup,BBDuk,Kraken2,GRCh38,T2T,HPRC\n"
    if not report.exists():
        report.write_text(header)

    def count(p): return nreads_fastq(Path(p))

    if layout == "paired":
        # discover R1/R2 pairs
        pairs = {}
        for r1 in list(inp.glob("*_R1.fq*")) + list(inp.glob("*_R1.fastq*")):
            r2 = Path(str(r1).replace("_R1.","_R2."))
            if r2.exists():
                sample = (r1.name
                          .replace("_R1.fastq.gz","")
                          .replace("_R1.fastq","")
                          .replace("_R1.fq.gz","")
                          .replace("_R1.fq",""))
                pairs[sample] = (r1, r2)

        for sample, (r1, r2) in sorted(pairs.items()):
            print(f"[WGS] {sample} (paired)")
            init_r1, init_r2 = count(r1), count(r2)

            # 1) Trim Galore
            subprocess.check_call(["trim_galore","--paired","--cores",str(threads),
                                   "--quality",str(trim_q),"--length",str(trim_len),
                                   "--output_dir",str(out), str(r1), str(r2)])
            tr_r1 = out/f"{sample}_R1_val_1.fq.gz"; tr_r2 = out/f"{sample}_R2_val_2.fq.gz"
            if not tr_r1.exists(): tr_r1 = Path(str(tr_r1).removesuffix(".gz"))
            if not tr_r2.exists(): tr_r2 = Path(str(tr_r2).removesuffix(".gz"))
            trc_r1, trc_r2 = count(tr_r1), count(tr_r2)

            # 2) Dedup (fastuniq preferred, else dedupe.sh)
            dd_r1 = out/f"{sample}_dedup_R1.fq"; dd_r2 = out/f"{sample}_dedup_R2.fq"
            if shutil.which("fastuniq"):
                # fastuniq needs uncompressed inputs
                with tempfile.TemporaryDirectory() as td:
                    r1u = Path(td)/"R1.fq"; r2u = Path(td)/"R2.fq"
                    subprocess.check_call(["bash","-lc",f"pigz -dc {tr_r1} > {r1u}"])
                    subprocess.check_call(["bash","-lc",f"pigz -dc {tr_r2} > {r2u}"])
                    lst = Path(td)/"pairs.txt"; lst.write_text(f"{r1u}\n{r2u}\n")
                    subprocess.check_call(["fastuniq","-i",str(lst),"-t","q","-o",str(dd_r1),"-p",str(dd_r2)])
            else:
                subprocess.check_call(["dedupe.sh", f"in1={tr_r1}", f"in2={tr_r2}",
                                       f"out1={dd_r1}", f"out2={dd_r2}", "ac=f", f"threads={threads}"])
            ddc_r1, ddc_r2 = count(dd_r1), count(dd_r2)

            # 3) BBDuk (UniVec)
            vec_r1 = out/f"{sample}_bbduk_R1.fq"; vec_r2 = out/f"{sample}_bbduk_R2.fq"
            subprocess.check_call(["bbduk.sh", f"in1={dd_r1}", f"in2={dd_r2}",
                                   f"out1={vec_r1}", f"out2={vec_r2}",
                                   f"ref={cfg['univec_ref']}", "k=27","hdist=1", f"threads={threads}"])
            vcc_r1, vcc_r2 = count(vec_r1), count(vec_r2)

            # 4) Kraken2 (keep unclassified)
            krdb = kraken2_override or cfg["kraken2_db"]
            kr_un_r1 = out/f"{sample}_kraken2_R_1.fq"
            kr_un_r2 = out/f"{sample}_kraken2_R_2.fq"
            subprocess.check_call(["kraken2","--db",krdb,"--threads",str(threads),
                                   "--paired","--unclassified-out", str(out/f"{sample}_kraken2_R#.fq"),
                                   str(vec_r1), str(vec_r2)])
            krc_r1, krc_r2 = count(kr_un_r1), count(kr_un_r2)

            # 5) BWA vs GRCh38 → both-unmapped pairs
            bwa_un_r1 = out/f"{sample}_un_grch38_bwa_R1.fq"
            bwa_un_r2 = out/f"{sample}_un_grch38_bwa_R2.fq"
            p1 = subprocess.Popen(["bwa","mem","-t",str(threads), cfg["grch38_fa"], str(kr_un_r1), str(kr_un_r2)],
                                  stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["samtools","view","-@",str(threads),"-b","-f","12","-F","256"],
                                  stdin=p1.stdout, stdout=subprocess.PIPE)
            subprocess.check_call(["samtools","fastq","-@",str(threads),
                                   "-1",str(bwa_un_r1),"-2",str(bwa_un_r2),
                                   "-0","/dev/null","-s","/dev/null","-n"], stdin=p2.stdout)
            bwc_r1, bwc_r2 = count(bwa_un_r1), count(bwa_un_r2)

            # 6) Bowtie2 vs T2T → --un-conc retains both-unmapped
            subprocess.check_call(["bowtie2","-p",str(threads),"-x",cfg["t2t_bt2_prefix"],
                                   "-1",str(bwa_un_r1),"-2",str(bwa_un_r2),
                                   "--very-sensitive","--un-conc",str(out/f"{sample}_un_t2t_R%.fq"), "-S","/dev/null"])
            t2t_un_r1 = out/f"{sample}_un_t2t_R1.fq"; t2t_un_r2 = out/f"{sample}_un_t2t_R2.fq"
            t2c_r1, t2c_r2 = count(t2t_un_r1), count(t2t_un_r2)

            # 7) minimap2 vs HPRC → both-unmapped pairs
            hprc_un_r1 = out/f"{sample}_un_hprc_R1.fq"
            hprc_un_r2 = out/f"{sample}_un_hprc_R2.fq"
            ref_hprc = cfg.get("hprc_mmi") or cfg["hprc_fa"]
            p3 = subprocess.Popen(["minimap2","-t",str(threads),"-ax","sr", ref_hprc, str(t2t_un_r1), str(t2t_un_r2)],
                                  stdout=subprocess.PIPE)
            p4 = subprocess.Popen(["samtools","view","-@",str(threads),"-b","-f","12","-F","256"],
                                  stdin=p3.stdout, stdout=subprocess.PIPE)
            subprocess.check_call(["samtools","fastq","-@",str(threads),
                                   "-1",str(hprc_un_r1),"-2",str(hprc_un_r2),
                                   "-0","/dev/null","-s","/dev/null","-n"], stdin=p4.stdout)
            hpc_r1, hpc_r2 = count(hprc_un_r1), count(hprc_un_r2)

            with open(report, "a") as r:
                r.write(",".join(map(str, [sample, init_r1, init_r2, trc_r1, trc_r2, ddc_r1, ddc_r2,
                                           vcc_r1, vcc_r2, krc_r1, krc_r2, bwc_r1, bwc_r2,
                                           t2c_r1, t2c_r2, hpc_r1, hpc_r2])) + "\n")

    else:
        # SINGLE-END
        files = sorted(list(inp.glob("*.fastq*")) + list(inp.glob("*.fq*")))
        for r1 in files:
            sample = r1.name.rsplit(".",1)[0]
            print(f"[WGS] {sample} (single)")
            init = count(r1)

            subprocess.check_call(["trim_galore","--cores",str(threads),
                                   "--quality",str(trim_q),"--length",str(trim_len),
                                   "--output_dir",str(out), str(r1)])
            tr = next((p for p in [out/f"{sample}_trimmed.fq.gz", out/f"{sample}_trimmed.fq"] if p.exists()), None)
            if not tr:
                cand = list(out.glob(f"{sample}*trim*fq*"))
                if not cand:
                    print(f"[humanfilt] no trimmed file for {sample}; skip")
                    continue
                tr = cand[0]
            trc = count(tr)

            dd = out/f"{sample}_dedup.fq"
            subprocess.check_call(["dedupe.sh", f"in={tr}", f"out={dd}", "ac=f", f"threads={threads}"])
            ddc = count(dd)

            vec = out/f"{sample}_bbduk.fq"
            subprocess.check_call(["bbduk.sh", f"in={dd}", f"out={vec}",
                                   f"ref={cfg['univec_ref']}", "k=27","hdist=1", f"threads={threads}"])
            vcc = count(vec)

            krdb = kraken2_override or cfg["kraken2_db"]
            kr_un = out/f"{sample}_kraken2.fastq"
            subprocess.check_call(["kraken2","--db",krdb,"--threads",str(threads),
                                   "--unclassified-out",str(kr_un), str(vec)])
            krc = count(kr_un)
            current = kr_un

            # minimap2 chain: GRCh38 -> T2T -> HPRC
            def _filter(tag, ref):
                un = out/f"{sample}_un_{tag}.fastq"
                p1 = subprocess.Popen(["minimap2","-t",str(threads),"-ax","sr", ref, str(current)], stdout=subprocess.PIPE)
                p2 = subprocess.Popen(["samtools","view","-@",str(threads),"-b","-F","256"], stdin=p1.stdout, stdout=subprocess.PIPE)
                p3 = subprocess.Popen(["samtools","view","-@",str(threads),"-b","-f","4"], stdin=p2.stdout, stdout=subprocess.PIPE)
                subprocess.check_call(["samtools","fastq","-0","/dev/null","-s","/dev/null","-n","-"],
                                      stdin=p3.stdout, stdout=open(un,"wb"))
                return un

            current = _filter("grch38", cfg["grch38_fa"]); c_grch = count(current)
            t2t_fa = next((p for p in Path(cfg["data_dir"]).rglob("T2T*.*fa*")), Path(cfg["grch38_fa"]))
            current = _filter("t2t", str(t2t_fa)); c_t2t = count(current)
            ref_hprc = cfg.get("hprc_mmi") or cfg["hprc_fa"]
            current = _filter("hprc", ref_hprc); c_hprc = count(current)

            with open(report, "a") as r:
                r.write(f"{sample},{init},{trc},{ddc},{vcc},{krc},{c_grch},{c_t2t},{c_hprc}\n")

