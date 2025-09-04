from pathlib import Path
import subprocess
from .utils import ensure_tools, nreads_fastq

def run_rnaseq(input_dir, output_dir, report_csv, threads, cfg,
               rna_preset="illumina", trim_q=20, trim_len=20,
               pattern=None, kraken2_override=None):
    """RNA-seq is always single-end. Presets: 'illumina' -> sr, 'ont' -> map-ont."""
    ensure_tools(required=["trim_galore","bbduk.sh","kraken2","minimap2","samtools","fastqc","dedupe.sh","seqkit"])

    mm = "map-ont" if rna_preset == "ont" else "sr"
    inp = Path(input_dir); out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    report = Path(report_csv)
    if not report.exists():
        report.write_text("Sample,Initial,Trim,Dedup,BBDuk,Kraken2,GENCODE,GRCh38,T2T,HPRC\n")

    patt = pattern or "*.R1.fastq.gz"
    files = sorted(list(inp.glob(patt)) or list(inp.glob(patt.replace(".gz",""))))
    if not files:
        print(f"[humanfilt] no files matched {patt} in {inp}")
        return

    def count(p): return nreads_fastq(Path(p))

    for r1 in files:
        sample = r1.name.split(".R1.fastq")[0]
        print(f"[RNA] {sample}")
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

        # minimap2 filters (GENCODE -> GRCh38 -> T2T -> HPRC)
        def _filter(tag, ref, preset):
            un = out/f"{sample}_un_{tag}.fastq"
            p1 = subprocess.Popen(["minimap2","-t",str(threads),"-ax",preset, ref, str(kr_un)], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["samtools","view","-@",str(threads),"-b","-F","256"], stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(["samtools","view","-@",str(threads),"-b","-f","4"], stdin=p2.stdout, stdout=subprocess.PIPE)
            subprocess.check_call(["samtools","fastq","-0","/dev/null","-s","/dev/null","-n","-"],
                                  stdin=p3.stdout, stdout=open(un,"wb"))
            return un

        gencode_ref = cfg.get("gencode_tx_mmi") or cfg["gencode_tx_fa"]
        current = _filter("gencode", gencode_ref, mm); c_gencode = count(current)
        current = _filter("grch38", cfg["grch38_fa"], mm); c_grch = count(current)
        t2t_fa = next((p for p in Path(cfg["data_dir"]).rglob("T2T*.*fa*")), Path(cfg["grch38_fa"]))
        current = _filter("t2t", str(t2t_fa), mm); c_t2t = count(current)
        ref_hprc = cfg.get("hprc_mmi") or cfg["hprc_fa"]
        current = _filter("hprc", ref_hprc, mm); c_hprc = count(current)

        with open(report, "a") as r:
            r.write(f"{sample},{init},{trc},{ddc},{vcc},{krc},{c_gencode},{c_grch},{c_t2t},{c_hprc}\n")

