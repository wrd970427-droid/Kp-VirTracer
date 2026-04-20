# Kp-VirTracer

Kp-VirTracer is an R package + command-line workflow for analyzing
**Klebsiella pneumoniae** genome assemblies and tracing virulence-related
transmission signals.

It uses a dual-layer architecture:

- Conda environment layer: external tools (`kleborate`, `mob_recon`, `blastn`,
  `makeblastdb`, `fastANI`, `prodigal`)
- R package layer: workflow orchestration, parsing, summarization, and CLI entry

## Project Overview

- Input: a directory containing `*.fa`, `*.fasta`, or `*.fna` files
- Minimum sample count: `>= 2`
- Output: standardized per-stage directories and summary tables

## Installation

### 1) Create conda environment

```bash
git clone <repo_url>
cd Kp-VirTracer
bash install_env.sh
conda activate kpvirtracer
```

### 2) Install the R package

```bash
Rscript install_package.R
```

Optional GitHub install:

```r
remotes::install_github("YOUR_GITHUB_USERNAME/Kp-VirTracer")
```

## Prerequisite: Build Virulence BLAST DB (Required)

Before running Kp-VirTracer, you must create a nucleotide BLAST database
from your virulence FASTA file.

```bash
makeblastdb -in /path/to/virulence_genes.fasta -dbtype nucl -out /path/to/db/virulence_db
```

Important:

- Runtime argument must be the BLAST DB prefix (`/path/to/db/virulence_db`)
- Do not pass FASTA path to `--virulence-db`
- Kp-VirTracer does not build this DB at runtime

## Run

### Rscript entry (recommended)

```bash
Rscript run_kpvirtracer.R \
  --input /path/to/fasta_dir \
  --output /path/to/output_dir \
  --virulence-db /path/to/db/virulence_db \
  --threads 8
```

### R interactive

```r
library(KpVirTracer)
run_kp_virtracer(
  input = "/path/to/fasta_dir",
  output = "/path/to/output_dir",
  virulence_db = "/path/to/db/virulence_db",
  threads = 8
)
```

### Package CLI entry

```bash
KpVirTracer \
  --input /path/to/fasta_dir \
  --output /path/to/output_dir \
  --virulence-db /path/to/db/virulence_db \
  --threads 8
```

## External Dependencies

The R package does not auto-install external tools. Install them in the conda
environment first:

- `kleborate`
- `mob_recon` (`mob_suite`)
- `blastn`
- `makeblastdb`
- `fastANI`
- `prodigal`

## Output Structure

```text
output/
|- logs/
|- temp/
|- sample_manifest.tsv
|- 01_kleborate/
|- 02_mobsuite/
|- 03_blast_virulence/
|- 04_ani/
|- 05_annotation/
|- 06_hgt/
`- summary/
```

Key result files:

- `summary/sample_summary.tsv`
- `summary/virulence_hits.tsv`
- `04_ani/ani_results.tsv`
- `06_hgt/hgt_events.tsv`
- `summary/final_summary.tsv`

## Environment Check

```bash
bash scripts/check_env.sh
```

## Smoke Test

```bash
bash scripts/smoke_test.sh
```
