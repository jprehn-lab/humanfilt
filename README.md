# humanfilt

Conda CLI for WGS human‐read decontamination.

- Tool dependencies are installed by Conda.
- References/DBs are auto‑downloaded on first run via `humanfilt setup` from Zenodo record `17020482`.
- Threads: auto‑detected from env (`HUMANFILT_THREADS`, `THREADS`, `OMP_NUM_THREADS`, `CPU_COUNT`) or CPU count; override with `--threads`.

## Quick start

1) Create environment (dev example)
```bash
conda create -n hf -c conda-forge -c bioconda \
  python=3.10 bwa bowtie2 minimap2 samtools kraken2 trim-galore \
  bbmap fastuniq seqkit curl zstd tar pigz -y
conda activate hf
pip install -e .
```

2) Download references (first time)
```bash
humanfilt setup
```

3) Run WGS (paired‑end)
```bash
humanfilt run \
  --mode wgs \
  --input  /path/to/paired_fastqs \
  --output /path/to/out \
  --report /path/to/report.csv \
  --threads 16
```

Outputs:
- Final HPRC‑unmapped FASTQs per sample in `--output`.
- CSV summary with per‑step counts at `--report`.

Recommended workflow
- Install once via Conda, then run `humanfilt setup` to populate your per‑user cache.
- After that, use `humanfilt run` normally — no extra flags needed. The tool reuses the cached references and does not re-download.

## WGS pipeline (paired)

Per sample:
1. Trim Galore (`--trim-quality`, `--trim-length`)
2. De‑duplicate (FastUniq)
3. UniVec screen (BBDuk)
4. Kraken2 (keep UNCLASSIFIED)
5. BWA vs GRCh38 (keep both‑unmapped pairs)
6. Bowtie2 vs T2T (keep both‑unmapped pairs)
7. Minimap2 vs HPRC (keep both‑unmapped pairs) → final FASTQs

## Command options

Common
- `--threads N`              number of threads (auto if omitted)
- `--no-auto-setup`          don’t download; use existing cache
- `--keep-temp`              keep per‑sample temp dirs/logs (debug)

Trim Galore
- `--trim-quality N`         Phred cutoff (default: 20)
- `--trim-length N`          minimum read length (default: 20)

Kraken2
- `--kraken2-db PATH`        override database path for this run

## References & cache

On first run, `humanfilt setup` downloads and indexes:
- GRCh38 (BWA),
- T2T (Bowtie2),
- HPRC (minimap2 .mmi),
- UniVec (BBDuk),
- Kraken2 (human DB).

Cache location: `~/.local/share/humanfilt` by default. Override with `--data-dir` or set `HUMANFILT_ZENODO_RECORD` to point to an alternate Zenodo record ID.

## Examples

Run with custom trim settings
```bash
humanfilt run --mode wgs \
  --input /data/wgs/paired \
  --output /data/wgs/out \
  --report /data/wgs/report.csv \
  --threads 24 \
  --trim-quality 25 --trim-length 30
```

Skip downloads (use pre‑populated cache)
```bash
humanfilt run --mode wgs \
  --input /data/wgs/paired \
  --output /data/wgs/out \
  --report /data/wgs/report.csv \
  --threads 24 \
  --no-auto-setup
```
