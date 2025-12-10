# install_pkg.R

# 1. CRANパッケージのインストール
cran_packages <- c(
  "shiny", "shinycssloaders", "shinyjs", "DT", "plotly", 
  "dplyr", "purrr", "data.table", "tibble", "tidyr", 
  "ggplot2", "pheatmap", "RColorBrewer", "igraph", "Rtsne", "umap", 
  "msigdbr", "ggrepel", "rhandsontable", "writexl", "matrixStats"
)

new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages)

# 2. Bioconductorパッケージのインストール
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_packages <- c(
  "edgeR", "DESeq2", "fgsea", "clusterProfiler", 
  "enrichplot", "AnnotationDbi", "maSigPro",
  # 生物種ごとのアノテーションDB (必要に応じて追加・コメントアウト)
  "org.Hs.eg.db",  # ヒト
  "org.Mm.eg.db"   # マウス
)

new_bioc_packages <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(new_bioc_packages)) BiocManager::install(new_bioc_packages, update = FALSE)

message("すべてのパッケージのインストールチェックが完了しました！")
