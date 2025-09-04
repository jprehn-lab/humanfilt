# humanfilt

Conda/Bioconda CLI for WGS & RNA-seq decontamination.

- Tool dependencies are installed by Conda.
- **References/DBs download on first run** via `humanfilt setup` from Zenodo record `17020482`.
- Threads: auto-detected from env (`HUMANFILT_THREADS`, `THREADS`, `OMP_NUM_THREADS`, `CPU_COUNT`) or CPU count; override with `--threads`.

## Install (dev)
```bash
conda create -n hf -c conda-forge -c bioconda \
  python=3.10 bwa bowtie2 minimap2 samtools kraken2 trim-galore fastqc \
  bbmap fastuniq seqkit curl zstd tar pigz -y
conda activate hf
pip install -e .

