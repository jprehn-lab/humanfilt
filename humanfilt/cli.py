#!/usr/bin/env python3
import os
import sys
import argparse

# Reuse the downloader/config from the installed humanfilt package
from humanfilt.downloader import ensure_data_ready, validate_cfg_or_die, load_or_scan_config
from pathlib import Path as _Path
from .pipeline_wgs import run_wgs


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="humanfilt",
        description="Paired-end WGS human read filtering for Illumina inputs.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("setup", help="Download and prepare reference bundles; records chosen data dir")
    s.add_argument("--data-dir", default=None, help="Destination directory for references (default: user cache)")
    s.add_argument("--threads", type=int, default=None)
    s.add_argument("--force", action="store_true")

    r = sub.add_parser("run", help="Run the paired-end WGS pipeline")
    r.add_argument("--mode", choices=["wgs", "rna-seq"], required=True)
    r.add_argument("--input", required=True)
    r.add_argument("--output", required=True)
    r.add_argument("--report", required=True)
    r.add_argument("--threads", type=int, default=None)
    r.add_argument("--trim-quality", type=int, default=20)
    r.add_argument("--trim-length", type=int, default=20)
    r.add_argument("--data-dir", default=None)
    r.add_argument("--kraken2-db", dest="kraken2_db", default=None)
    r.add_argument("--no-auto-setup", action="store_true")
    r.add_argument("--keep-temp", action="store_true", help="Keep per-sample temp dirs for debugging")
    return p

    # (Download subcommand is defined below in main if needed)


def main() -> int:
    args = build_parser().parse_args()

    if args.cmd == "setup":
        thr = args.threads
        cfg = ensure_data_ready(args.data_dir, threads=thr, force=args.force)
        validate_cfg_or_die(cfg)
        base = cfg.get("data_dir", "") or (os.environ.get("XDG_DATA_HOME") and str(_Path(os.environ.get("XDG_DATA_HOME")) / "humanfilt") or str(_Path.home() / ".local" / "share" / "humanfilt"))
        print(base)
        return 0

    # Prepare refs/config for run
    if args.no_auto_setup:
        cfg = load_or_scan_config(args.data_dir)
    else:
        cfg = ensure_data_ready(args.data_dir, threads=args.threads, force=False)

    if getattr(args, "kraken2_db", None):
        cfg["kraken2_db"] = args.kraken2_db

    if args.mode == "rna-seq":
        print("--mode rna-seq is recognized but not implemented yet.", file=sys.stderr)
        return 2

    validate_cfg_or_die(cfg)
    run_wgs(
        args.input, args.output, args.report, args.threads, cfg,
        trim_quality=args.trim_quality, trim_length=args.trim_length,
        keep_temp=args.keep_temp,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
