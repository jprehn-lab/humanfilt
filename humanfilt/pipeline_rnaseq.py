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
