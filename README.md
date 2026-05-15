# EasyRNA-Seq

A Shiny-based GUI platform for RNA-seq analysis designed for wet-lab scientists with no command-line experience. Covers the full pipeline from raw FASTQ files to publication-ready figures.

---

## Features

| Module | Description |
|--------|-------------|
| **Data Upload** | Load a count matrix (TSV/CSV) and sample metadata in a point-and-click interface |
| **Filtering** | Remove low-count genes with interactive threshold controls |
| **Normalization** | TMM and CPM normalization via edgeR |
| **Dimension Reduction** | PCA, t-SNE, and UMAP plots (interactive, powered by plotly) |
| **Differential Expression** | DEG analysis using edgeR; volcano plots, MA plots, and result tables |
| **GSEA** | Gene Set Enrichment Analysis with fgsea and MSigDB gene sets |
| **GO/KEGG Enrichment** | Over-representation analysis with clusterProfiler; dot plots and bar charts |
| **Time-series Analysis** | Temporal DEG detection and cluster profiling using maSigPro |

---

## Prerequisites

### Option A — Docker (recommended)

- [Docker Desktop](https://docs.docker.com/desktop/) installed and running
- No R installation required

### Option B — Local RStudio

- R 4.2 or later
- RStudio
- Required packages (see [Local Setup](#option-b--local-rstudio-1))

---

## Quick Start

### Option A — Docker

```bash
cd after_count
docker-compose up --build
```

Open your browser at **http://localhost:3838**

The first build downloads and compiles all R packages and may take 15–30 minutes. Subsequent starts are instant.

To stop: press `Ctrl+C` in the terminal.

### Option B — Local RStudio

1. Open `after_count/app.R` in RStudio (downstream analysis only), or open the root `app.R` for the integrated platform.
2. Install dependencies:

```r
install.packages(c(
  "shiny", "plotly", "DT", "shinycssloaders",
  "pheatmap", "ggplot2", "tibble", "writexl",
  "dplyr", "tidyr", "Rtsne", "umap"
))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "edgeR", "limma", "fgsea", "clusterProfiler",
  "org.Hs.eg.db", "org.Mm.eg.db", "maSigPro"
))
```

3. Click **Run App**.

---

## Input Format

### Count Matrix

- Format: tab-separated (`.tsv`) or comma-separated (`.csv`)
- Rows: genes (row names = gene IDs, e.g. Ensembl or gene symbols)
- Columns: samples
- Values: raw integer counts (do **not** pre-normalize)

Example:

```
gene_id        Sample1  Sample2  Sample3  Sample4
ENSG00000001   120      98       210      187
ENSG00000002   0        2        1        0
...
```

### Metadata

- Format: CSV
- One row per sample; must include at least a sample ID column and a group/condition column
- Column names must match the count matrix column names

Example:

```
sample,condition
Sample1,CD4
Sample2,CD4
Sample3,CD14
Sample4,CD14
```

---

## Test Data

A demo count matrix is included for immediate use:

- **File**: `cd4_cd14.tsv`
- **Content**: Human CD4 T cells vs CD14 monocytes, ~26,000 genes

Load this file in the Data Upload module to explore all analysis features without your own data.

---

## Repository Structure

```
toll_rnaseq/
├── app.R                   # Integrated platform (upstream + downstream)
├── after_count/            # Downstream statistical analysis app
│   ├── app.R
│   ├── Dockerfile
│   ├── docker-compose.yml
│   └── R/                  # Analysis modules
├── before_count/           # Upstream preprocessing app (fastp, STAR/Salmon, featureCounts)
│   └── EasyRNASeq_Preprocessor/
├── cd4_cd14.tsv            # Test data
└── install_pkg.R           # Package installer script
```

---

## License

License TBD.
