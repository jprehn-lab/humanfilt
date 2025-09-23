#!/usr/bin/env python3
import os
import sys
import argparse
from pathlib import Path  # optional
# import subprocess  # only if you actually use it

# Default Zenodo record baked in; users may still override via env
os.environ.setdefault("HUMANFILT_ZENODO_RECORD", "17020482")

from .downloader import ensure_data_ready, validate_cfg_or_die, load_or_scan_config
from .pipeline_wgs import run_wgs
try:
    # Optional; we may hide RNA in some builds
    from .pipeline_rnaseq import run_rnaseq  # type: ignore
except Exception:
    run_rnaseq = None


def build_parser() -> argparse.ArgumentParser:
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
    r.add_argument("--kraken2-db", dest="kraken2_db", default=None)
    r.add_argument("--no-auto-setup", action="store_true")
    r.add_argument("--keep-temp", action="store_true", help="Keep per-sample temp dirs for debugging")
    r.add_argument("--save-bams", action="store_true", help="Save alignment BAMs for single-end (BWA/mm2)")
    return p


def main() -> int:
    args = build_parser().parse_args()

    if args.cmd == "setup":
        cfg = ensure_data_ready(args.data_dir, force=args.force)
        validate_cfg_or_die(cfg)
        print("Setup complete.", file=sys.stderr)
        return 0

    # Prepare refs/config for 'run'
    if args.no_auto_setup:
        cfg = load_or_scan_config(args.data_dir)
    else:
        cfg = ensure_data_ready(args.data_dir, force=False)

    if getattr(args, "kraken2_db", None):
        cfg["kraken2_db"] = args.kraken2_db

    threads = args.threads

    if args.mode == "wgs":
        validate_cfg_or_die(cfg)  # only WGS needs this here
        run_wgs(
            args.input, args.output, args.report, threads, cfg,
            trim_quality=args.trim_quality, trim_length=args.trim_length,
            keep_temp=args.keep_temp, save_bams=args.save_bams,
        )
        return 0

    elif args.mode == "rnaseq":
        if run_rnaseq is None:
            print("RNA mode is disabled in this build (coming soon).", file=sys.stderr)
            return 2
        # Call your RNA-seq runner (adjust signature if different)
        return int(run_rnaseq(
            args.input, args.output, args.report, threads, cfg
        ) or 0)

    else:
        print("Unknown mode", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())

