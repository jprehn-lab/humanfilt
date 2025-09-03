import argparse, json
from pathlib import Path
from .downloader import ensure_data_ready
from .pipeline_wgs import run_wgs
from .pipeline_rnaseq import run_rnaseq

def _parser():
    p = argparse.ArgumentParser(
        prog="humanfilt",
        description="Human read filtering (WGS paired/single & RNA-seq single) with first-run Zenodo downloads",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("setup", help="Download/refresh references & DBs from Zenodo")
    s.add_argument("--data-dir", default=None, help="Cache directory (default: ~/.local/share/humanfilt)")
    s.add_argument("--force", action="store_true", help="Re-download even if present")

    r = sub.add_parser("run", help="Run the pipeline")
    r.add_argument("--mode", choices=["wgs","rnaseq"], required=True)
    r.add_argument("--input", required=True)
    r.add_argument("--output", required=True)
    r.add_argument("--report", required=True)
    r.add_argument("--threads", type=int, default=8)
    r.add_argument("--data-dir", default=None)
    r.add_argument("--kraken2-db", default=None, help="Override Kraken2 DB path")

    # Trim Galore controls (requested)
    r.add_argument("--trim-quality", type=int, default=20,
                   help="Trim Galore --quality <INT>. Default: 20.")
    r.add_argument("--trim-length", type=int, default=20,
                   help="Trim Galore --length <INT>. Default: 20 (0 disables read discard).")

    # WGS layout (paired or single). RNA-seq is always single.
    r.add_argument("--wgs-layout", choices=["paired","single"], default="paired",
                   help="Only used when --mode wgs (default: paired).")

    # RNA-seq specifics
    r.add_argument("--rna-preset", choices=["illumina","ont"], default="illumina",
                   help="Minimap2 preset for RNA stage (illumina=sr, ont=map-ont)")
    r.add_argument("--pattern", default=None,
                   help="RNA input pattern (default: *.R1.fastq.gz; fallback: any *.fastq* not matching _R2).")
    return p

def main():
    args = _parser().parse_args()

    if args.cmd == "setup":
        cfg = ensure_data_ready(args.data_dir, force=args.force)
        print(json.dumps(cfg, indent=2))
        return

    cfg = ensure_data_ready(args.data_dir, force=False)
    if args.kraken2_db:
        cfg["kraken2_db"] = args.kraken2_db

    Path(args.output).mkdir(parents=True, exist_ok=True)

    if args.mode == "wgs":
        run_wgs(args.input, args.output, args.report, args.threads, cfg,
                layout=args.wgs_layout, trim_quality=args.trim_quality, trim_length=args.trim_length)
    else:
        run_rnaseq(args.input, args.output, args.report, args.threads, cfg,
                   rna_preset=args.rna_preset, pattern=args.pattern,
                   trim_quality=args.trim_quality, trim_length=args.trim_length)

if __name__ == "__main__":
    main()
