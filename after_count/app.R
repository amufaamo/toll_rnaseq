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
library(rhandsontable) # ★ これを追加！
# [Phase 1] command log が使用
library(jsonlite)
library(yaml)

# ★★★ 生物種に合わせて OrgDb パッケージをロード ★★★
library(org.Mm.eg.db) # マウスの場合 (例として残しています)

# --- 3. ヘルパー関数 ---
# (もしあれば)

# --- 4. モジュールファイルを読み込み ---
required_modules <- c("R/module_data_upload_metadata_new.R", "R/module_filtering.R",
                      "R/module_processing.R", "R/module_dimension_reduction.R",
                      "R/module_deg_analysis.R", "R/module_gsea.R",
                      "R/module_go_enrichment_integrated.R",
                      "R/module_timeseries_analysis.R",
                      "R/module_gene_barplot_swap.R",
                      "R/module_figure_enrichment.R"
)
missing_files <- required_modules[!file.exists(required_modules)]
if(length(missing_files) > 0) {
  stop("以下のモジュールファイルが見つかりません:\n", paste(missing_files, collapse="\n"))
}
# モジュールファイルを読み込む
source("R/module_data_upload_metadata_new.R")
source("R/module_filtering.R")
source("R/module_processing.R")
source("R/module_dimension_reduction.R")
source("R/module_deg_analysis.R")
source("R/module_gsea.R")
source("R/module_go_enrichment_integrated.R")
source("R/module_timeseries_analysis.R")
source("R/module_gene_barplot_swap.R")
source("R/module_figure_enrichment.R")
# [Phase 1] command layer + command log tab (additive only)
source("R/command_core.R")
source("R/module_command_log.R")
# GTFアノテーション パース + 共通ID表示変換ヘルパー
source("R/gtf_utils.R")

# --- 4.5 名前空間衝突の解消 ---
# AnnotationDbi / S4Vectors / clusterProfiler 等が dplyr の動詞を上書きするため、
# 全モジュール読み込み後に dplyr 版を明示的に最優先で再束縛する。
# (例: dplyr::select が AnnotationDbi::select に隠れ "継承メソッドが見付かりません"、
#      dplyr::rename が隠れ rename(logFC = log2FoldChange) で "object not found" になる問題)
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
slice <- dplyr::slice
count <- dplyr::count
first <- dplyr::first

# --- 5. UI 定義 ---
ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel("RNA-seq データ解析パイプライン v1.0.0"),
  tabsetPanel(
    id = "mainTabs",
    tabPanel("1. データ入力とセッション管理", dataUploadMetadataUI("dataTab")),
    tabPanel("2. フィルタリング", filteringUI("filterTab")),
    tabPanel("3. データ前処理", processingUI("procTab")),
    tabPanel("4. サンプル間類似性", dimensionReductionUI("dimRedTab")),
    tabPanel("5. DEG解析", degAnalysisUI("degTab")),
    tabPanel("6. GSEA", gseaUI("gseaTab")),
    tabPanel("7. GO Enrichment", goEnrichmentIntegratedUI("go_module")),
    tabPanel("8. 時系列解析", timeseriesAnalysisUI("timeseriesTab")),
    tabPanel("9. 鍵遺伝子 Barplot & スワップ検証", geneBarplotSwapUI("swapTab")),
    tabPanel("10. 作図 & 非モデル生物エンリッチメント", figureEnrichmentUI("figEnrichTab")),
    # [Phase 1] new tab
    tabPanel("11. コマンドログ", commandLogUI("cmdLogTab"))
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
    gene_annotation = NULL,   # GTFアップロード由来: data.frame(Geneid, gene_name, biotype, length)
    filtered_keep = NULL,
    background_genes_original = NULL,
    deg_results = NULL,
    selected_species = NULL,
    current_gene_id_type = "ENTREZID",
    # [Phase 1] command log state (additive)
    command_log = tibble::tibble(
      timestamp   = as.POSIXct(character()),
      session_id  = character(),
      op          = character(),
      status      = character(),
      message     = character(),
      params_json = character(),
      result_json = character()
    ),
    session_id = format(Sys.time(), "%Y%m%d_%H%M%S")
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
    selected_species_code_reactive = reactive(rv$selected_species),
    gene_annotation_reactive = reactive(rv$gene_annotation)
  )
  
  # 8. 時系列解析
  timeseriesAnalysisServer("timeseriesTab", rv)

  # 9. 鍵遺伝子 Barplot & スワップ検証
  geneBarplotSwapServer("swapTab", rv)

  # 10. 作図 & 非モデル生物エンリッチメント
  figureEnrichmentServer("figEnrichTab")

  # 11. [Phase 1] コマンドログ
  commandLogServer("cmdLogTab", rv)

}

# --- 7. アプリケーションの実行 ---
shinyApp(ui = ui, server = server)

# --- End of app.R ---
