#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
from pathlib import Path

# Default Zenodo record baked in; users may still override via env
os.environ.setdefault("HUMANFILT_ZENODO_RECORD", "17020482")

from .downloader import ensure_data_ready, validate_cfg_or_die, load_or_scan_config
from .pipeline_wgs import run_wgs
try:
    # Optional; we may hide RNA in some builds
    from .pipeline_rnaseq import run_rnaseq  # type: ignore
except Exception:
    run_rnaseq = None

def build_parser():
    p = argparse.ArgumentParser(
        prog="humanfilt",
        description="Human read filtering (WGS & RNA-seq) with first-run Zenodo downloads",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("setup", help="Download/refresh references & DBs from Zenodo")
    s.add_argument("--data-dir", default=None, help="Cache dir (default: ~/.local/share/humanfilt)")
    s.add_argument("--force", action="store_true", help="Re-download archives even if present")

    enable_rna = bool(os.environ.get("HUMANFILT_ENABLE_RNA"))
    r = sub.add_parser("run", help="Run the pipeline (auto-download refs on first run)")
    r.add_argument("--mode", choices=["wgs", "rnaseq"] if enable_rna else ["wgs"], required=True)
    r.add_argument("--input", required=True)
    r.add_argument("--output", required=True)
    r.add_argument("--report", required=True)
    r.add_argument("--threads", type=int, default=None)
    r.add_argument("--trim-quality", type=int, default=20)
    r.add_argument("--trim-length", type=int, default=20)
    r.add_argument("--data-dir", default=None)
    r.add_argument("--kraken2-db", default=None)
    # (no RNA-specific CLI args needed)
    r.add_argument("--no-auto-setup", action="store_true")
    r.add_argument("--keep-temp", action="store_true", help="Keep per-sample temp dirs for debugging")
    r.add_argument("--save-bams", action="store_true", help="Save alignment BAMs for single-end (BWA/mm2)")

    if False:
        # Lightweight bash-based RNA test (no UniVec/BBDuk) per user request
        rb = sub.add_parser("rna-test", help="Run bash RNA test pipeline (trim->kraken2->BWA->HPRC)")
        rb.add_argument("--input", required=True)
        rb.add_argument("--output", required=True)
        rb.add_argument("--report", required=True)
        rb.add_argument("--threads", type=int, default=None)
        rb.add_argument("--data-dir", default=None)
        rb.add_argument("--no-auto-setup", action="store_true")
        rb.add_argument("--save-bams", action="store_true", help="Also save full BAMs for BWA and HPRC")
        rb.add_argument("--grch38-prefix", dest="grch38_prefix", default=None, help="Override GRCh38 BWA index prefix")
        rb.add_argument("--grch38-fasta", dest="grch38_fasta", default=None, help="Override GRCh38 FASTA path")
        rb.add_argument("--hprc-fasta", dest="hprc_fasta", default=None, help="Override HPRC FASTA path")
        rb.add_argument("--kraken2-db", dest="kraken2_db", default=None, help="Override Kraken2 DB path")

        # NEXT variant wrapper using external sequence
        rn = sub.add_parser("rna-test-next", help="Run external NEXT RNA pipeline wrapper")
        rn.add_argument("--input", required=True)
        rn.add_argument("--output", required=True)
        rn.add_argument("--report", required=True)
        rn.add_argument("--threads", type=int, default=None)
        rn.add_argument("--data-dir", default=None)
        rn.add_argument("--no-auto-setup", action="store_true")
        rn.add_argument("--save-bams", action="store_true")
        rn.add_argument("--script", default="/mnt/dat5/colon_cancer/cleaning_fastq_next.sh")
        rn.add_argument("--grch38-prefix", dest="grch38_prefix", default=None)
        rn.add_argument("--grch38-fasta", dest="grch38_fasta", default=None)
        rn.add_argument("--hprc-fasta", dest="hprc_fasta", default=None)
        rn.add_argument("--kraken2-db", dest="kraken2_db", default=None)
    return p

def main():
    args = build_parser().parse_args()

    if args.cmd == "setup":
        cfg = ensure_data_ready(args.data_dir, force=args.force)
        validate_cfg_or_die(cfg)
        print("Setup complete.", file=sys.stderr)
        return

    # run
    if args.no_auto_setup:
        # Do not auto-download; only use existing cache/config if present
        cfg = load_or_scan_config(args.data_dir)
    else:
        # Ensure refs are present, downloading if missing
        cfg = ensure_data_ready(args.data_dir, force=False)

    if getattr(args, "kraken2_db", None):
        cfg["kraken2_db"] = args.kraken2_db

    # Validate only for WGS; external RNA script manages its own refs
    if getattr(args, "mode", None) == "wgs":
        validate_cfg_or_die(cfg)

    threads = args.threads
    # (RNA test branches removed)
    elif args.mode == "wgs":
        run_wgs(
            args.input, args.output, args.report, threads, cfg,
            trim_quality=args.trim_quality, trim_length=args.trim_length,
            keep_temp=args.keep_temp, save_bams=args.save_bams,
        )
    else:
        # Hidden by default; enable with HUMANFILT_ENABLE_RNA=1
        if not os.environ.get("HUMANFILT_ENABLE_RNA"):
            print("RNA mode is disabled in this build (coming soon). Set HUMANFILT_ENABLE_RNA=1 to enable.", file=sys.stderr)
            sys.exit(2)
        from .rnaseq_pipeline import run_rnaseq
        run_rnaseq(
            args.input, args.output, args.report, threads, cfg,
            trim_quality=args.trim_quality, trim_length=args.trim_length,
            save_bams=args.save_bams,
        )

if __name__ == "__main__":
    main()
