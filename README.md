# Kp-VirTracer

Kp-VirTracer is a reproducible R-based workflow for tracing virulence signals and relatedness in *Klebsiella pneumoniae* genome assemblies, with explicit chromosome/plasmid-aware interpretation.

## Why This Project

Hypervirulence in *K. pneumoniae* can be associated with plasmid-borne, chromosome-borne, or mixed virulence determinants. Kp-VirTracer provides a practical pipeline to:

- screen virulence loci by BLAST
- separate chromosome and plasmid evidence using MOB-suite contig typing
- infer plasmid PTU using COPLA (when available)
- classify each sample into `p-hvKp`, `c-hvKp`, `pc-hvKp`, or `nKp`
- run ANI only on virulent strains using chromosome-only sequences

## 研究背景与意义（中文）

在近缘肺炎克雷伯菌（尤其是高风险谱系）中，近期HGT识别面临一个典型难题：

- 宿主间进化距离很近，近期转移片段可提供的系统发育信息位点有限
- 在这种场景下，`gene tree` 与 `species tree` 的“一致”并不总能证明“纯垂直继承”
- 传统“树冲突法”并非错误，而是对“近缘、近期、模块化、质粒介导传播”的检测功效会下降

因此，本项目采用的是“互补证据框架”，而不是“替代系统发育”：

1. 用系统发育信息作为宿主背景约束（host background constraint）。
2. 用区室与载体证据补足盲区：染色体/质粒分区、局部共线性、PTU/replicon、oriT/relaxase/T4SS、IS与融合线索。
3. 用机制一致性支持传播推断：不仅回答“像不像HGT”，也回答“为什么能够发生传播”。

这套框架的实际价值在于：

- 对近缘菌中近期毒力模块传播更敏感
- 对“non-conjugative 但可被 helper 动员”的真实场景更有解释力
- 输出可直接用于流调汇报、风险分层和后续实验验证设计

简要定位：

- Kp-VirTracer is a **phylogeny-constrained, mechanism-aware** workflow for recent plasmid-borne virulence dissemination inference.

## 理论边界与互补分析框架（详细版）

### 核心判断（适用于论文与汇报）

- 传统系统发育方法在近缘肺炎克雷伯菌中并不是“错误”，而是对近期HGT存在“系统性漏检”风险。
- 在“同一物种、近缘谱系、近期传播、质粒模块化载荷”的场景下，`gene tree` 与 `species tree` 的一致性不能直接当作“纯垂直继承”的正证据。
- 最稳妥的研究定位是：以系统发育作为宿主背景约束，再引入区室、局部上下文、质粒骨架与动员机制证据做互补推断。

### 为什么会漏检：可识别性而非单纯算力问题

- 供体和受体过近时，转移片段可提供的变异位点太少，难以重建稳定基因树。
- 近期转移后的模块若继续在受体分支中垂直传播，最终树形可能仍与 `species tree` 相容。
- 因此“缺乏冲突”常常表示“分辨率不足”，而不一定表示“没有HGT”。

### 近缘Kp中的机制事实（与本项目直接相关）

- pLVPK-like 等毒力质粒虽然常被标为 non-conjugative / non-mobilizable，但在真实生态中可被 helper 质粒动员。
- IS 介导重排、融合质粒、耐药-毒力共整合是高频机制线索。
- 所以需要同时追踪：`区室位置 + 骨架归属 + 动员元件 + 局部共线性`，而不只看单基因树冲突。

### 方法学评价矩阵（从“能否识别近期近缘HGT”角度）

| 方法类别 | 对近期近缘HGT敏感性 | 特异性 | 主要数据要求 | 主要局限 | 在本项目中的定位 |
|---|---:|---:|---|---|---|
| species-tree vs gene-tree 冲突 / reconciliation | 低-中 | 中-高 | 高质量核心比对与正交基因树 | 受 gene tree 噪声影响大 | 作为背景约束，不单独裁决 |
| presence-absence / profile | 中 | 中 | 泛基因组矩阵 | 对机制解释能力有限 | 大队列初筛 |
| 区室感知 + 局部共线性 | 高 | 中-高 | 可靠注释与邻域信息 | 阈值敏感，需对照 | 近期事件识别主轴 |
| PTU / replicon / backbone 分型 | 中-高 | 中 | 质粒序列或高可信 bins | 新骨架可能无PTU | 追踪载体流动主轴 |
| oriT/relaxase/T4CP/T4SS + IS/融合注释 | 对“能否传播”极高 | 中 | 完整或近完整质粒序列 | 需结合上下文解释 | 机制解释主轴 |
| hybrid long-read plasmidomics | 真值层 | 最高 | 长读长+短读长 | 成本较高 | 关键节点确认层 |

### 证据链建议（可直接用于结果章节）

- `同PTU/同backbone`
- `模块边界存在IS痕迹或融合线索`
- `oriT/relaxase/T4SS或helper条件可满足`
- `跨宿主背景但局部上下文高度相似`
- `宿主树不能支持单次祖先继承`

当以上证据共同成立时，应优先表述为：

- `consistent with recent horizontal dissemination`
- `recent-HGT-like`（而非“仅凭单证据即排除垂直遗传”）

### 推荐可复现工作流（面向大队列）

1. 样本质控与去冗余：QC、污染检查、MLST/Kleborate、同ST内近重复折叠。  
2. 区室划分：MOB-suite 进行 plasmid binning 与 mobility 预测（用于初筛）。  
3. 质粒骨架分型：COPLA/PTU + replicon/backbone 聚类。  
4. 宿主背景系统发育：核心SNP + 去重组（如 Gubbins）构建背景树。  
5. 模块级HGT初筛：以毒力模块为单位做 identity/coverage/跨背景筛选。  
6. 局部上下文复核：共线性与邻域保守性分析，做阈值敏感性。  
7. 传播机制注释：oriT/relaxase/T4CP/T4SS、IS、fusion/cointegration。  
8. 真值确认：关键节点做 hybrid long-read 组装验证。  

### 审稿友好型验证策略（建议在补充材料中实现）

- 阈值敏感性：例如 SI cutoff 6/7/8、coverage/identity 多阈值复跑。
- 空模型比较：随机置换宿主标签或邻域，评估是否显著高于随机。
- 负对照模块：选择 housekeeping 或稳定染色体标记，验证框架不会“过判HGT”。
- 关键机制验证：至少对代表性“non-mobilizable × helper”组合做实验验证。

### 写作表达模板（可直接复用）

中文：

1. 本研究并不否定系统发育在宿主背景重建中的价值，而是将其作为约束层，并以区室/骨架/机制证据补足近期传播盲区。  
2. 在近缘背景中，gene tree 与 species tree 的一致性不应直接解释为垂直继承的正证据。  
3. 我们将相关事件表述为“recent-HGT-like / consistent with recent horizontal dissemination”，而非单证据绝对判定。  

English:

1. We treat host phylogeny as a background constraint and complement it with compartment-aware and mechanism-aware evidence.  
2. In closely related strains, tree concordance should not be interpreted as definitive evidence of strict vertical inheritance.  
3. We report events as consistent with recent horizontal dissemination rather than claiming exclusion of vertical inheritance from a single signal.  

### 开放问题与局限

- 无完整质粒装配时，fusion/cointegration/shared-backbone 结论属于高可信推断，而非最终证明。
- 局部共线性破坏不专属于HGT，也可能来自重排/插入/缺失，需联合区室与机制证据解释。
- oriT/relaxase/T4SS 的存在表示“潜在可传播性”，不等于“已发生传播”。
- 公共数据库存在时间、地区和克隆组成偏倚，结论应以趋势与框架表达为主。

## Workflow Overview

```mermaid
flowchart TD
    A["Input assemblies (*.fa/*.fasta/*.fna)"] --> B["Kleborate"]
    A --> C["MOB-suite (mob_recon)"]
    A --> D["Virulence BLAST"]
    C --> E["Contig location map (chromosome/plasmid)"]
    D --> F["Virulence hits with location labels"]
    E --> F
    F --> G["Virulence type classification: p-hvKp / c-hvKp / pc-hvKp / nKp"]
    G --> H["Select virulent samples only"]
    C --> I["Chromosome contig IDs"]
    H --> J["Build chromosome-only FASTA per virulent sample"]
    I --> J
    J --> K["fastANI (virulent chromosomes only)"]
    F --> L["Summary reports"]
    G --> L
    K --> L
```

## Key Features

- Chromosome/plasmid-aware virulence interpretation.
- Explicit virulence type output:
  - `p-hvKp`: virulence genes on plasmid only
  - `c-hvKp`: virulence genes on chromosome only
  - `pc-hvKp`: virulence genes on both
  - `nKp`: no virulence hit
- ANI analysis constrained to biologically relevant comparisons:
  - virulent strains only
  - chromosome sequences only
- PTU assignment integrates COPLA output into `plasmid_manifest.tsv` (optional, non-blocking).
- CLI parameter overrides for both BLAST and ANI while keeping sensible defaults from config.

## Installation

### 1) Clone

```bash
git clone https://github.com/wrd970427-droid/Kp-VirTracer.git
cd Kp-VirTracer
```

### 2) Create Conda/Mamba Environment

```bash
mamba env create -f environment.yml
mamba activate kpvirtracer
```

If `mamba` is unavailable:

```bash
conda env create -f environment.yml
conda activate kpvirtracer
```

### 3) Install R Package

```bash
Rscript install_package.R
```

Note: if you modify package source code, reinstall before running CLI entry scripts.

### 4) Optional: Install COPLA In a Dedicated Conda Environment

COPLA is not installed through this repository `environment.yml`. Install it separately following the official project style:

```bash
PROJECT_ROOT_DIRECTORY=~/COPLA
git clone https://github.com/santirdnd/COPLA ${PROJECT_ROOT_DIRECTORY}
cd ${PROJECT_ROOT_DIRECTORY}

conda env create -f copla.environment.yml -n copla
```

If needed by your COPLA setup, install MacSyFinder in its dedicated environment:

```bash
conda env create -f macsyfinder.environment.yml -n macsyfinder
```

Download COPLA databases:

```bash
cd ${PROJECT_ROOT_DIRECTORY}
bin/download_Copla_databases.sh
head databases/Copla_RS84/CoplaDB.fofn
```

Recommended post-install check:

```bash
bin/post_install_test.sh
```

Then configure Kp-VirTracer to call this COPLA environment (see config section below).

## Required Prerequisite: Build Virulence BLAST DB

```bash
makeblastdb \
  -in /path/to/virulence_genes.fasta \
  -dbtype nucl \
  -out /path/to/db/virulence_db
```

Rules:

- pass BLAST DB prefix to `--virulence-db`
- do not pass raw FASTA to `--virulence-db`

## COPLA Configuration In Kp-VirTracer

Edit your config file (or `inst/extdata/config.yaml.example`) and set:

```yaml
copla:
  enabled: true
  conda_bin: "conda"
  conda_env: "copla"
  python_bin: "python3"
  script_path: "/abs/path/to/COPLA/bin/copla.py"
  pickle_path: "/abs/path/to/COPLA/databases/Copla_RS84/RS84f_sHSBM.pickle"
  fofn_path: "/abs/path/to/COPLA/databases/Copla_RS84/CoplaDB.fofn"
  topology: "linear"
```

### How Users Can Locate Their Own COPLA Paths

Run the following commands on the target server:

```bash
# 1) Locate COPLA script
find / -type f -name "copla.py" 2>/dev/null

# 2) Locate COPLA model pickle
find / -type f -name "RS84f_sHSBM.pickle" 2>/dev/null

# 3) Locate COPLA database fofn
find / -type f -name "CoplaDB.fofn" 2>/dev/null

# 4) Locate conda executable
which conda

# 5) Confirm COPLA environment exists
conda env list | grep -E '^copla\\s'
```

Use the discovered absolute paths in the `copla:` config section.

Runtime behavior:

- if COPLA is fully configured, PTU is inferred and written to `02_mobsuite/plasmid_manifest.tsv`
- if COPLA is missing or misconfigured, pipeline does not stop; a warning is logged and PTU fields remain `NA`

Quick validation before running pipeline:

```bash
conda run -n copla python3 /abs/path/to/COPLA/bin/copla.py --help
```

## Quick Start

```bash
Rscript run_kpvirtracer.R \
  --input /path/to/fasta_dir \
  --output /path/to/output_dir \
  --virulence-db /path/to/db/virulence_db \
  --threads 8
```

### Optional BLAST Overrides

```bash
Rscript run_kpvirtracer.R \
  --input /path/to/fasta_dir \
  --output /path/to/output_dir \
  --virulence-db /path/to/db/virulence_db \
  --blast-task blastn \
  --blast-evalue 1e-10 \
  --blast-min-identity 90 \
  --blast-min-coverage 80 \
  --threads 8
```

### Optional ANI Overrides

```bash
Rscript run_kpvirtracer.R \
  --input /path/to/fasta_dir \
  --output /path/to/output_dir \
  --virulence-db /path/to/db/virulence_db \
  --ani-relatedness 99 \
  --ani-min-fraction 0.2 \
  --ani-frag-len 3000 \
  --ani-kmer 16 \
  --threads 8
```

## Output Structure

```text
output/
|- logs/
|- temp/
|- sample_manifest.tsv
|- 01_kleborate/
|- 02_mobsuite/
|- 03_blast_virulence/
|  |- logs/
|  |- raw/
|  |- filtered/
|  |- by_location/
|  `- summary/
|- 04_ani/
|- 05_annotation/
|- 06_hgt/
`- summary/
```

Key files:

- `02_mobsuite/plasmid_manifest.tsv`
- `summary/sample_summary.tsv`
- `summary/virulence_type_summary.tsv`
- `summary/virulence_hits.tsv`
- `03_blast_virulence/summary/sample_virulence_profile.tsv`
- `03_blast_virulence/summary/virulence_type_calls.tsv`
- `04_ani/ani_results.tsv`
- `05_annotation/annotation_manifest.tsv`
- `05_annotation/annotation_summary.tsv`
- `06_hgt/hgt_events.tsv`
- `summary/final_summary.tsv`

## Result Tables And Column Dictionary

### `02_mobsuite/plasmid_manifest.tsv`

- `Sample_ID`: sample identifier derived from input FASTA filename.
- `MOB_Status`: execution status of `mob_recon` (`ok`, `failed`, `not_run`).
- `plasmid_id`: plasmid/contig identifier reported by MOB-suite.
- `replicon_type`: replicon typing result.
- `relaxase`: relaxase typing result.
- `mobility`: predicted mobility class from MOB-suite.
- `predicted_host_range`: host range estimate from MOB-suite.
- `PTU`: PTU assignment from COPLA (if available); `NA` when unassigned/unavailable.
- `PTU_Score`: COPLA confidence score; higher indicates more stable assignment.
- `PTU_Host_Range`: COPLA host range annotation (if reported).
- `PTU_Notes`: COPLA notes/warnings (for example, unassigned cluster size notes).

Interpretation:

- `PTU=NA` with notes like "could not be assigned" usually means the plasmid does not map confidently to a named PTU.
- `MOB_Status=ok` but `PTU=NA` is possible and biologically plausible.

### `03_blast_virulence/summary/virulence_hits.tsv`

- `Sample_ID`: sample identifier.
- `Gene`: virulence gene matched in BLAST database.
- `Location`: inferred location (`chromosome` or `plasmid`) using MOB contig typing.
- `Contig_ID`: matched query contig identifier.
- `Start`, `End`, `Strand`: hit coordinates and orientation on query contig.
- `Identity`: BLAST percent identity.
- `Coverage`: estimated coverage percentage for the hit.
- `Evalue`: BLAST e-value.
- `Bitscore`: BLAST bitscore.
- `Plasmid_ID`: reserved field (currently may be `NA`).
- `Source_FASTA`: original assembly FASTA path.

Interpretation:

- High-confidence hits should satisfy your configured identity/coverage/evalue thresholds.
- `Location` is central for downstream `p-hvKp/c-hvKp/pc-hvKp` classification.

### `03_blast_virulence/summary/sample_virulence_profile.tsv`

- `Sample_ID`: sample identifier.
- `Virulence_Genes_All`: all detected virulence genes.
- `Chromosomal_Genes`: subset detected on chromosome contigs.
- `Plasmid_Genes`: subset detected on plasmid contigs.
- `n_chrom_hits`: count of chromosome-assigned hits.
- `n_plasmid_hits`: count of plasmid-assigned hits.
- `Virulence_Location`: `chromosome`, `plasmid`, `both`, or `absent`.
- `Virulence_Type`: final per-sample type (`c-hvKp`, `p-hvKp`, `pc-hvKp`, `nKp`).

### `03_blast_virulence/summary/virulence_type_calls.tsv`

- `Sample_ID`: sample identifier.
- `Virulence_Location`: location category used for typing.
- `Virulence_Type`: final virulence type.
- `n_chrom_hits`: chromosome hit count.
- `n_plasmid_hits`: plasmid hit count.

### `summary/sample_summary.tsv`

- `Sample_ID`: sample identifier.
- `Virulence_Type`: one of `p-hvKp/c-hvKp/pc-hvKp/nKp`.
- `hvKp_Status`: compact status (`hvKp` vs `nKp`).
- `Virulence_Genes`, `Chromosomal_Genes`, `Plasmid_Genes`: gene lists used for typing.
- `Virulence_Location`: `chromosome/plasmid/both/absent`.
- `Plasmid_Count`: number of plasmid rows in plasmid manifest.
- `Replicon_Types`: merged replicon annotations.
- `PTU_Types`: merged PTU annotations.
- `Kleborate_Virulence_Score`: Kleborate virulence score when available.

### `summary/virulence_type_summary.tsv`

- `Virulence_Type`: virulence class.
- `n_samples`: number of samples in each class.

### `04_ani/ani_results.tsv`

- `Sample_A`, `Sample_B`: compared virulent samples.
- `ANI`: average nucleotide identity.
- `FragmentsMapped`: mapped fragment count reported by FastANI.
- `FragmentsTotal`: total fragments considered.
- `Relatedness`: `Related`, `Unrelated`, or `Unknown` (based on ANI threshold and output availability).

Interpretation:

- ANI runs only for virulent samples and only on chromosome-only FASTA.
- `Relatedness=Related` means ANI met/exceeded `--ani-relatedness`.

### `05_annotation/annotation_manifest.tsv`

- `Sample_ID`: sample identifier.
- `Contig_ID`: contig identifier.
- `molecule_type`: `chromosome` or `plasmid` (from MOB contig report).
- `n_virulence_hits`: count of virulence hits on this contig.
- `Virulence_Genes`: comma-separated virulence genes detected on this contig.
- `Prodigal_Status`: status of prodigal run for that sample (`ok`, `failed`, `not_run`).

### `05_annotation/annotation_summary.tsv`

- `Sample_ID`: sample identifier.
- `n_contigs`: number of contigs included in annotation manifest.
- `n_virulence_contigs`: number of contigs with at least one virulence hit.
- `Prodigal_Status`: sample-level prodigal status.

### `06_hgt/hgt_events.tsv`

- `Sample_A`, `Sample_B`: sample pair.
- `Gene`: shared virulence gene considered for event inference.
- `Location_A`, `Location_B`: gene location in each sample.
- `Plasmid_A`, `Plasmid_B`: plasmid flag placeholders for plasmid-associated entries.
- `HGT_Type`: event heuristic label (`location_switch`, `shared_plasmid_gene`, `shared_chromosomal_gene`).
- `Confirmed`: heuristic confirmation flag using ANI/location consistency.
- `pident`, `qcov`: aggregated identity/coverage evidence from virulence hits.
- `Synteny_Index`: placeholder/basic synteny score (current implementation baseline).
- `Evidence_Files`: primary files used for inference.

### `summary/final_summary.tsv`

- `n_samples`: total number of analyzed samples.
- `n_virulence`: number of samples with virulence type not equal to `nKp`.
- `n_ani_pairs`: number of ANI pairwise comparisons performed.
- `n_hgt_events`: number of inferred HGT rows.

## Virulence Typing Rules

- BLAST hits are assigned to chromosome/plasmid by MOB-suite `contig_report.txt`.
- BLAST stage writes layered outputs:
  - `raw/`: raw `outfmt 6` results per sample
  - `filtered/`: threshold-filtered results per sample
  - `by_location/`: chromosome/plasmid split results per sample
  - `summary/`: merged hit table and sample-level virulence calls
- Classification per sample:
  - `p-hvKp`: plasmid-only virulence hits
  - `c-hvKp`: chromosome-only virulence hits
  - `pc-hvKp`: both chromosome and plasmid virulence hits
  - `nKp`: no virulence hit

## ANI Rules

- ANI is run only for virulent samples (`p-hvKp`, `c-hvKp`, `pc-hvKp`).
- ANI input is reconstructed chromosome-only FASTA from MOB-suite chromosome contig IDs.
- Relatedness is reported by threshold (`--ani-relatedness`, default from config).

## PTU Calling (COPLA)

- `mob_recon` is used for plasmid reconstruction and plasmid/chromosome labeling.
- PTU is called by COPLA from extracted plasmid FASTA per sample via:
  - `conda run -n <copla_env> python3 <copla.py> <query_fasta> <RS84f_sHSBM.pickle> <CoplaDB.fofn> <output_dir>`
- COPLA is optional:
  - if available and configured, PTU/score are filled into `02_mobsuite/plasmid_manifest.tsv`
  - if unavailable/misconfigured, workflow continues and PTU fields remain `NA`

## Environment Validation

```bash
bash scripts/check_env.sh
```

## Smoke Test

```bash
bash scripts/smoke_test.sh
```

## Methodological Context

Kp-VirTracer integrates ideas and tools commonly used in published microbial genomics pipelines, including:

- Kleborate
- MOB-suite
- BLAST+
- FastANI

For manuscript writing or formal citation sections, please cite the corresponding original tool publications used in your analysis environment/version.
