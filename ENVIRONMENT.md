# Environment version pins (backup)

Snapshot date: **2026-07-14**  
Validated with the current Kp-VirTracer pipeline on linux-64.

## Which file to use

| File | Purpose | When to use |
|------|---------|-------------|
| [`environment.yml`](environment.yml) | Flexible ranges | Quick / exploratory installs |
| [`environment.pinned.yml`](environment.pinned.yml) | Direct deps pinned | Recommended rebuild of `kpvirtracer` |
| [`environment.lock.yml`](environment.lock.yml) | Full freeze (all transitive deps) | Disaster recovery / bit-for-bit recreate |
| [`envs/copla.pinned.yml`](envs/copla.pinned.yml) | COPLA direct deps pinned | Rebuild optional `copla` env |
| [`envs/copla.lock.yml`](envs/copla.lock.yml) | COPLA full freeze | Exact COPLA recreate |

## Recreate commands

```bash
# Main analysis environment (preferred pinned recipe)
mamba env create -f environment.pinned.yml
mamba activate kpvirtracer
Rscript install_package.R

# Or full freeze
mamba env create -f environment.lock.yml

# Optional COPLA companion env
mamba env create -f envs/copla.pinned.yml
# or: mamba env create -f envs/copla.lock.yml
```

## Critical runtime pins (kpvirtracer)

| Component | Version |
|-----------|---------|
| Python | 3.11.15 |
| R | 4.5.3 |
| Kleborate | 3.2.4 |
| MOB-suite | 3.1.9 |
| BLAST+ | 2.17.0 |
| FastANI | 1.34 |
| Prodigal | 2.6.3 |
| Biostrings | 2.78.0 |

## Critical runtime pins (copla)

| Component | Version |
|-----------|---------|
| Python | **3.8.19** (do not upgrade) |
| graph-tool | 2.37 |
| BLAST+ | 2.9.0 |
| Prodigal | 2.6.3 |
| plasmidfinder | 2.1.6 |
| NumPy | 1.24.4 |
| pandas | 2.0.3 |
| Biopython | 1.83 |

## Policy

Do **not** upgrade Python or core packages inside these environments without an explicit decision and re-validation of Kleborate, MOB-suite, BLAST, ANI, Prodigal, and COPLA PTU outputs.
