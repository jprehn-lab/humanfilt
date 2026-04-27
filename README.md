# HumanFilt

HumanFilt removes human reads from paired-end whole-genome sequencing (WGS) FASTQ files and writes cleaned FASTQs plus a per-sample summary report.

- Install with Conda.
- Run `humanfilt setup` once to download reference files.
- Run `humanfilt run` on a folder of paired FASTQs.
- License: MIT.

## Quick start

1. Create and activate the environment

```bash
conda create -n hf -y -c conda-forge -c bioconda humanfilt=1.1.0
conda activate hf
```

2. Download references (first time only)

```bash
humanfilt setup
```

To store references in a custom location:

```bash
humanfilt setup --data-dir /path/to/humanfilt-data
```

3. Run HumanFilt on paired-end WGS FASTQs

```bash
humanfilt run \
  --mode wgs \
  --input /path/to/paired_fastqs \
  --output /path/to/out \
  --report /path/to/report.csv \
  --threads 16
```

**By default**: All intermediates **except final output** are deleted after successful completion.

**To keep everything**:
```bash
humanfilt run ... --save-intermediates
```

**File naming**: Each intermediate reflects the filtering step completed (`kraken2`, `trimmed`, `grch38`, `t2t`, `final`)


## What you get

After a successful run, HumanFilt writes:

- Final non-human FASTQ files for each sample in `--output`
- A CSV report with per-step read counts in `--report`
- Temporary working files, which are removed automatically unless you keep them with `--save-intermediates`



## Common examples

Run with custom trimming:

```bash
humanfilt run \
  --mode wgs \
  --input /data/wgs/paired \
  --output /data/wgs/out \
  --report /data/wgs/report.csv \
  --threads 16 \
  --trim-quality 25 \
  --trim-length 30
```

Keep intermediate files and logs:

```bash
humanfilt run \
  --mode wgs \
  --input /data/wgs/paired \
  --output /data/wgs/out \
  --report /data/wgs/report.csv \
  --threads 24 \
  --save-intermediates \
  --save-logs
```

## Input files

`humanfilt run --mode wgs` expects paired-end FASTQ files in a single input directory.

Supported read-1 naming patterns include:

- `*_R1.fq`
- `*_R1.fq.gz`
- `*_R1.fastq`
- `*_R1.fastq.gz`
- `*_R1_001.fq`
- `*_R1_001.fq.gz`
- `*_R1_001.fastq`
- `*_R1_001.fastq.gz`
- `*_1.fq`
- `*_1.fq.gz`
- `*_1.fastq`
- `*_1.fastq.gz`
- `*.1.fq`
- `*.1.fq.gz`
- `*.1.fastq`
- `*.1.fastq.gz`

Each sample must have one matching read-2 file with the same sample stem, for example `_R2`, `_R2_001`, `_2`, or `.2`.

Examples:

- `sampleA_R1.fastq.gz` + `sampleA_R2.fastq.gz`
- `sampleA_R1_001.fastq.gz` + `sampleA_R2_001.fastq.gz`
- `sampleA_1.fastq.gz` + `sampleA_2.fastq.gz`
- `sampleA.1.fastq.gz` + `sampleA.2.fastq.gz`
- `patient7_R1.fq` + `patient7_R2.fq`


## References and cache

`humanfilt setup` downloads and indexes the references used by the WGS pipeline, including GRCh38, T2T, HPRC, GENCODE transcripts, UniVec, and a Kraken2 human database.

Default cache location:

- `$XDG_DATA_HOME/humanfilt` if `XDG_DATA_HOME` is set
- Otherwise `~/.local/share/humanfilt`

You can override this with:

```bash
humanfilt setup --data-dir PATH
```

or

```bash
humanfilt run --data-dir PATH
```

After `humanfilt setup --data-dir PATH`, HumanFilt stores a small pointer so later runs can reuse that location automatically.

## Main options

### Setup

- `--data-dir PATH` – reference/cache location
- `--threads N` – number of threads to use
- `--force` – re-download and rebuild references

### Run

- `--mode {wgs,rna-seq}` – pipeline mode (`rna-seq` currently not implemented)
- `--data-dir PATH` – use a specific reference/cache directory
- `--threads N` – number of threads
- `--no-auto-setup` – do not download references automatically
- `--save-intermediates` / `--keep-temp` – keep intermediate FASTQs
- `--save-logs` – keep per-sample logs

### Trimming

- `--trim-quality N` – Phred quality cutoff (default: 20)
- `--trim-length N` – minimum read length (default: 20)

### Kraken2

- `--kraken2-db PATH` – use a custom Kraken2 database for this run

## WGS pipeline

For each sample, HumanFilt runs:

1. Kraken2 (keep unclassified reads)
2. Trim Galore
3. FastUniq deduplication
4. UniVec screening with BBDuk
5. BWA against GRCh38
6. Bowtie2 against T2T
7. Minimap2 against HPRC

The final output is the paired-end reads that remain unmapped after all filtering steps.

## Notes

- HumanFilt auto-detects available CPU threads if `--threads` is not set.
- `rna-seq` is recognized as a mode but is not implemented yet.

## License

MIT License. See `LICENSE`.
