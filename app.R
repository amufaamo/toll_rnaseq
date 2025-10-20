# app.R (モジュール統合修正版 - セッション管理をタブ1に移動)

# --- 1. アプリケーション設定 ---
options(shiny.maxRequestSize = 100*1024^2) # 例: 100MB

# --- 2. 必要なパッケージをロード ---
library(shiny)
library(shinycssloaders)
library(shinyjs)
library(DT)
library(plotly)
library(dplyr)
library(purrr)
library(data.table)
library(tibble)
library(edgeR)
library(ggplot2)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(Rtsne)
library(umap)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(ggrepel) # 遺伝子ラベルの重なり回避に必要
library(maSigPro) # 時系列解析用に追加
# ★★★ 生物種に合わせて OrgDb パッケージをロード ★★★
# library(org.Hs.eg.db) # ヒトの場合
library(org.Mm.eg.db) # マウスの場合 (例として残しています)

# --- 3. ヘルパー関数 ---
# (もしあれば)

# --- 4. モジュールファイルを読み込み ---
required_modules <- c("R/module_data_upload_metadata.R", "R/module_filtering.R",
                      "R/module_processing.R", "R/module_dimension_reduction.R",
                      "R/module_deg_analysis.R", "R/module_gsea.R",
                      "R/module_go_enrichment_integrated.R",
                      "R/module_timeseries_analysis.R"
)
missing_files <- required_modules[!file.exists(required_modules)]
if(length(missing_files) > 0) {
  stop("以下のモジュールファイルが見つかりません:\n", paste(missing_files, collapse="\n"))
}
# モジュールファイルを読み込む
source("R/module_data_upload_metadata.R")
source("R/module_filtering.R")
source("R/module_processing.R")
source("R/module_dimension_reduction.R")
source("R/module_deg_analysis.R")
source("R/module_gsea.R")
source("R/module_go_enrichment_integrated.R")
source("R/module_timeseries_analysis.R")

# --- 5. UI 定義 ---
ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel("RNA-seq データ解析パイプライン"),
  tabsetPanel(
    id = "mainTabs",
    tabPanel("1. データ入力とセッション管理", dataUploadMetadataUI("dataTab")),
    tabPanel("2. フィルタリング", filteringUI("filterTab")),
    tabPanel("3. データ前処理", processingUI("procTab")),
    tabPanel("4. サンプル間類似性", dimensionReductionUI("dimRedTab")),
    tabPanel("5. DEG解析", degAnalysisUI("degTab")),
    tabPanel("6. GSEA", gseaUI("gseaTab")),
    tabPanel("7. GO Enrichment", goEnrichmentIntegratedUI("go_module")),
    tabPanel("8. 時系列解析", timeseriesAnalysisUI("timeseriesTab"))
  )
)

# --- 6. Server 定義 ---
server <- function(input, output, session) {
  
  # --- 6.1 共有リアクティブ値の初期化 ---
  rv <- reactiveValues(
    merged_data = NULL,
    sample_metadata = NULL,
    file_info = NULL,
    gene_lengths = NULL,
    filtered_keep = NULL,
    background_genes_original = NULL,
    deg_results = NULL,
    selected_species = NULL,
    current_gene_id_type = "ENTREZID"
  )
  # --- 6.2 各モジュールサーバー関数の呼び出し ---
  
  # 1. データ入力とセッション管理
  dataUploadMetadataServer("dataTab", rv)
  
  # 2. フィルタリング
  filteringServer("filterTab", rv)
  
  # --- filtered_keep が更新されたら、背景遺伝子リスト(Entrez ID)を更新 ---
  observe({
    req(rv$merged_data, rv$filtered_keep)
    
    gene_ids_from_merged_data <- rv$merged_data$Geneid
    
    if (!is.null(gene_ids_from_merged_data) &&
        length(gene_ids_from_merged_data) == length(rv$filtered_keep)) {
      
      rv$background_genes_original <- gene_ids_from_merged_data[rv$filtered_keep]
      message("[App] Updated rv$background_genes_original with EntrezIDs based on filtering. Count: ", length(rv$background_genes_original))
      
    } else {
      warning("Length mismatch between rv$merged_data$Geneid and rv$filtered_keep, or gene IDs from merged_data is NULL. Cannot update background_genes_original correctly.")
      rv$background_genes_original <- NULL
    }
  })
  
  
  # 3. データ前処理
  processingServer("procTab", rv)
  
  # 4. サンプル間類似性
  dimensionReductionServer("dimRedTab", rv)
  
  # 5. DEG解析
  degAnalysisServer("degTab", rv)
  
  # 6. GSEA
  gseaServer("gseaTab", rv)
  
  # 7. GOエンリッチメント解析
  goEnrichmentIntegratedServer(
    id = "go_module",
    deg_results_reactive = reactive(rv$deg_results),
    background_genes_reactive = reactive(rv$background_genes_original),
    selected_species_code_reactive = reactive(rv$selected_species)
  )
  
  # 8. 時系列解析
  timeseriesAnalysisServer("timeseriesTab", rv)
  
}

# --- 7. アプリケーションの実行 ---
shinyApp(ui = ui, server = server)

# --- End of app.R ---