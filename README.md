# EasyRNA-Seq Integrated Platform v2.0

A single Shiny application covering the full RNA-seq pipeline from FASTQ to figures.

## Features

| Stage | Tools / Modules |
|---|---|
| Upstream QC | fastp |
| Alignment | STAR, Salmon |
| Quantification | featureCounts (Docker-based) |
| Filtering & Normalization | TMM, CPM |
| Dimensionality Reduction | PCA, t-SNE, UMAP |
| Differential Expression | edgeR |
| Gene Set Enrichment | fgsea + MSigDB |
| Pathway Analysis | clusterProfiler (GO, KEGG) |
| Time-series Analysis | maSigPro |
| Gene ID Conversion | org.Hs.eg.db / biomaRt |

## Prerequisites

- [Docker Desktop](https://www.docker.com/products/docker-desktop/) — required for upstream steps (QC, alignment, quantification)
- R 4.2 or later
- RStudio

## Quick Start

1. Open `app.R` in RStudio.
2. Click **Run App**.

> Docker Desktop must be running before launching any upstream steps.

## R Package Installation

```r
install.packages(c(
  "shiny", "shinydashboard", "DT", "plotly",
  "ggplot2", "pheatmap", "Rtsne", "umap"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "edgeR", "fgsea", "clusterProfiler",
  "org.Hs.eg.db", "maSigPro"
))
```

## Input Format

**Upstream pipeline**

- Paired-end FASTQ files (`.fastq.gz`)

**Downstream (count matrix + metadata)**

- Count matrix: tab-separated file; rows = genes, columns = samples; first column = gene IDs
- Metadata: tab-separated file with a `sample` column and a `group` column

## Test Data

`cd4_cd14.tsv` — Human CD4 T-cell vs CD14 monocyte count matrix (~26,000 genes).
Use this file to explore all downstream modules without running the upstream pipeline.

## Repository Structure

```
toll_rnaseq/
├── app.R              # Main integrated Shiny application
├── R/                 # Module source files
├── scripts/           # Helper scripts
├── cd4_cd14.tsv       # Test count matrix
├── install_pkg.R      # Package installation helper
└── toll_rnaseq.Rproj  # RStudio project file
```

## License

MIT License — see LICENSE for details.
