# humanfilt

Conda-installable CLI for human read decontamination.

- Tools: bwa, bowtie2, minimap2, samtools, kraken2, trim-galore, fastqc, bbmap, fastuniq
- First-run download from **Zenodo record `17020482`** (split tar.zst bundles)
- Vector removal (UniVec via **BBDuk**) and **dedup** after Trim Galore
- **WGS (paired)**: Trim → UniVec → Dedup → Kraken2 → BWA(GRCh38) → Bowtie2(T2T) → minimap2(HPRC)
- **WGS (single)**: Trim → UniVec → Dedup → Kraken2 → BWA(GRCh38) → Bowtie2(T2T) → minimap2(HPRC)
- **RNA-seq (single)**: Trim → UniVec → Dedup → Kraken2 → minimap2(GENCODE → GRCh38 → T2T → HPRC)

## Trim Galore controls
- `--trim-quality <INT>` (default **20**)  
  Trim low-quality ends (BWA algorithm).
- `--trim-length <INT>` (default **20**, **0** disables)  
  Discard reads shorter than INT after trimming (quality/adapter).

## Quickstart
```bash
pip install -e .

export HUMANFILT_ZENODO_RECORD=17020482
humanfilt setup                 # downloads + extracts into ~/.local/share/humanfilt by default

# WGS paired-end
humanfilt run --mode wgs --wgs-layout paired \
  --input /path/in --output /path/out --report /path/out/report.csv --threads 16 \
  --trim-quality 20 --trim-length 20

# WGS single-end
humanfilt run --mode wgs --wgs-layout single \
  --input /path/in --output /path/out --report /path/out/report.csv --threads 16

# RNA-seq single-end
humanfilt run --mode rnaseq \
  --input /path/in --output /path/out --report /path/out/report.csv --threads 16 \
  --rna-preset illumina --trim-quality 20 --trim-length 20

