# Kp-VirTracer

Kp-VirTracer is an R package + command-line workflow for tracing virulence gene
transmission signals in **Klebsiella pneumoniae** genome assemblies.

## What This Tool Does

- Input: a folder containing `*.fa`, `*.fasta`, or `*.fna` assemblies
- Minimum sample count: `>= 2`
- Core analyses: Kleborate, MOB-suite, virulence BLAST, ANI, summary reporting
- Output: structured stage-by-stage results and final summary tables

## Architecture

- Environment layer (Conda/Mamba): external bioinformatics tools
- R package layer: workflow orchestration, parsing, integration, and CLI

---

## Installation

### 1) Clone repository

```bash
git clone https://github.com/wrd970427-droid/Kp-VirTracer.git
cd Kp-VirTracer
```

### 2) Create environment (Recommended: `mamba`)

This project is designed to create the runtime environment from `environment.yml`.

```bash
mamba env create -f environment.yml
mamba activate kpvirtracer
```

If `mamba` is unavailable, use:

```bash
conda env create -f environment.yml
conda activate kpvirtracer
```

### 3) Install the R package locally

```bash
Rscript install_package.R
```

---

## Required Prerequisite: Build Virulence BLAST DB

Before running Kp-VirTracer, you must build a nucleotide BLAST database from
your virulence FASTA file:

```bash
makeblastdb -in /path/to/virulence_genes.fasta -dbtype nucl -out /path/to/db/virulence_db
```

Important rules:

- `--virulence-db` must be the BLAST DB prefix from `makeblastdb -out`
- Do not pass FASTA directly to `--virulence-db`
- The pipeline does not create this database at runtime

---

## Quick Start

### Command line (recommended)

```bash
Rscript run_kpvirtracer.R \
  --input /path/to/fasta_dir \
  --output /path/to/output_dir \
  --virulence-db /path/to/db/virulence_db \
  --threads 8
```

### Package CLI entry

```bash
KpVirTracer \
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

---

## Output Layout

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

---

## Environment Check

```bash
bash scripts/check_env.sh
```

## Smoke Test

```bash
bash scripts/smoke_test.sh
```
