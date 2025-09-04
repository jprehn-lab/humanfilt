import argparse
from .downloader import ensure_data_ready
from .pipeline_wgs import run_wgs
from .pipeline_rnaseq import run_rnaseq
from .utils import detect_threads
import os


os.environ.setdefault("HUMANFILT_ZENODO_RECORD", "17020482")


def _parser():
    p = argparse.ArgumentParser(
        prog="humanfilt",
        description="Human read filtering (WGS & RNA-seq) with first-run Zenodo downloads",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("setup", help="Download/refresh references & DBs from Zenodo")
    s.add_argument("--data-dir", default=None, help="Cache directory (default: ~/.local/share/humanfilt)")
    s.add_argument("--force", action="store_true", help="Re-download even if present")

    r = sub.add_parser("run", help="Run the pipeline")
    r.add_argument("--mode", choices=["wgs","rnaseq"], required=True)
    r.add_argument("--input", required=True, help="Input folder with FASTQs")
    r.add_argument("--output", required=True, help="Output folder")
    r.add_argument("--report", required=True, help="CSV report path")
    r.add_argument("--threads", type=int, default=0, help="Override threads (default: auto from env/CPU)")
    r.add_argument("--trim-quality", type=int, default=20, help="Trim Galore: --quality INT (default 20)")
    r.add_argument("--trim-length", type=int, default=20, help="Trim Galore: --length  INT (default 20)")
    r.add_argument("--data-dir", default=None, help="Override cache folder")
    r.add_argument("--kraken2-db", default=None, help="Override Kraken2 DB path")
    r.add_argument("--wgs-layout", choices=["paired","single"], default="paired", help="WGS layout")
    r.add_argument("--rna-preset", choices=["illumina","ont"], default="illumina",
                   help="Minimap2 preset for RNA filters (illumina=sr, ont=map-ont)")
    r.add_argument("--pattern", default=None, help="RNA input pattern (default: *.R1.fastq.gz)")
    return p

def main():
    args = _parser().parse_args()
    if args.cmd == "setup":
        ensure_data_ready(args.data_dir, force=args.force)
        return

    threads = args.threads if args.threads and args.threads > 0 else detect_threads()
    cfg = ensure_data_ready(args.data_dir, force=False)

    if args.mode == "wgs":
        run_wgs(args.input, args.output, args.report, threads, cfg,
                layout=args.wgs_layout, trim_q=args.trim_quality, trim_len=args.trim_length,
                kraken2_override=args.kraken2_db)
    else:
        run_rnaseq(args.input, args.output, args.report, threads, cfg,
                   rna_preset=args.rna_preset, trim_q=args.trim_quality, trim_len=args.trim_length,
                   pattern=args.pattern, kraken2_override=args.kraken2_db)

if __name__ == "__main__":
    main()

