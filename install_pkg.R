# install_pkg.R

# 1. CRANパッケージのインストール
cran_packages <- c(
  "shiny", "shinycssloaders", "shinyjs", "DT", "plotly", 
  "dplyr", "purrr", "data.table", "tibble", "tidyr", 
  "ggplot2", "pheatmap", "RColorBrewer", "igraph", "Rtsne", "umap", 
  "msigdbr", "ggrepel", "rhandsontable", "writexl", "matrixStats",
  "remotes", "readr", "readxl", "testit", "data.tree", "limSolve",
  "e1071", "rlang", "stringr", "magrittr"
)

new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages)

# 2. Bioconductorパッケージのインストール
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

options(repos = BiocManager::repositories(version = "3.20"))

bioc_packages <- c(
  "edgeR", "DESeq2", "fgsea", "clusterProfiler", 
  "enrichplot", "AnnotationDbi", "maSigPro",
  "BiocParallel", "preprocessCore", "Biobase", "biomaRt", "sva",
  "quantiseqr", "HDF5Array", "GSVA", "genefilter",
  # 生物種ごとのアノテーションDB (必要に応じて追加・コメントアウト)
  "org.Hs.eg.db",  # ヒト
  "org.Mm.eg.db"   # マウス
)

new_bioc_packages <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(new_bioc_packages)) BiocManager::install(new_bioc_packages, version = "3.20", update = FALSE, ask = FALSE)

# 3. GitHubパッケージのインストール (Deconvolution等)
github_packages <- c(
  "GfellerLab/EPIC",
  "grst/MCPcounter",
  "grst/mMCPcounter",
  "omnideconv/ComICS",
  "dviraran/xCell",
  "cansysbio/ConsensusTME",
  "omnideconv/immunedeconv"
)

github_pkg_names <- c(
  "EPIC",
  "MCPcounter",
  "mMCPcounter",
  "ComICS",
  "xCell",
  "ConsensusTME",
  "immunedeconv"
)

for (i in seq_along(github_packages)) {
  if (!requireNamespace(github_pkg_names[[i]], quietly = TRUE)) {
    remotes::install_github(github_packages[[i]], dependencies = FALSE, upgrade = "never")
  }
}

if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  stop("immunedeconv のインストールに失敗しました。直前のエラーを確認してください。", call. = FALSE)
}

message("すべてのパッケージのインストールチェックが完了しました！")
