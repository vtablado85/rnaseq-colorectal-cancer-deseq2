# rnaseq-colorectal-cancer-deseq2

> **Status:** work in progress

Differential expression analysis of colorectal cancer RNA-seq data with DESeq2, using a paired design (same patient: normal colon, primary tumour and liver metastasis).

## Biological question

Which genes and pathways change from normal colon tissue to primary colorectal tumour, and from primary tumour to liver metastasis?

## Data

- **Dataset:** [GSE50760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760) (GEO)
- 54 samples from 18 colorectal cancer patients: normal colon, primary tumour and liver metastasis for each patient
- Illumina HiSeq 2000, paired-end RNA-seq
- Raw data are **not** stored in this repository; `scripts/01_download.R` downloads them into `data/raw/`

## Planned workflow

1. Download counts and sample metadata (GEOquery)
2. Quality control and exploratory analysis (PCA, sample distances)
3. Differential expression with DESeq2, design `~ patient + tissue`
4. Functional enrichment (clusterProfiler)
5. Report with Quarto

## Repository structure

```
data/raw/        downloaded data (ignored by git)
data/processed/  intermediate objects (ignored by git)
scripts/         numbered analysis scripts
results/figures/ plots
results/tables/  DE results and enrichment tables
```

## How to reproduce

_To be completed._

## Results

_To be completed._
