# R/module_deg_analysis.R
# (フィルタリング指標としてFDR/P-valueを選択可能にし、プロットも連動させる修正版)

library(shiny)
library(edgeR)
library(DESeq2) # DESeq2パッケージをロード
library(dplyr)
library(tidyr)
library(DT)
library(plotly)
library(ggplot2)
library(shinycssloaders)
library(AnnotationDbi) # ID変換に必要
library(writexl) # Excelダウンロードに必要
library(tibble) # rownames_to_columnのため
library(pheatmap) # ヒートマップ用

# --- Species to OrgDb package mapping (他のモジュールと共通) ---
orgdb_species_map <- c(
  "Homo_sapiens" = "org.Hs.eg.db",
  "Mus_musculus" = "org.Mm.eg.db",
  "Rattus_norvegicus" = "org.Rn.eg.db",
  "Drosophila_melanogaster" = "org.Dm.eg.db",
  "Caenorhabditis_elegans" = "org.Ce.eg.db",
  "Danio_rerio" = "org.Dr.eg.db",
  "Saccharomyces_cerevisiae" = "org.Sc.sgd.db",
  "Arabidopsis_thaliana" = "org.At.tair.db",
  "Bos_taurus" = "org.Bt.eg.db",
  "Gallus_gallus" = "org.Gg.eg.db",
  "Canis_familiaris" = "org.Cf.eg.db",
  "Macaca_mulatta" = "org.Mmu.eg.db",
  "Pan_troglodytes" = "org.Pt.eg.db",
  "Sus_scrofa" = "org.Ss.eg.db",
  "Xenopus_laevis" = "org.Xl.eg.db",
  "Anopheles_gambiae" = "org.Ag.eg.db",
  "Escherichia_coli_K12" = "org.EcK12.eg.db",
  "Escherichia_coli_Sakai" = "org.EcSakai.eg.db",
  "Plasmodium_falciparum" = "org.Pf.plasmo.db",
  "Myxococcus_xanthus_DK122" = "org.Mxanthus.db"
)

# (Helper functionは変更なし)
detect_gene_id_type <- function(ids) {
  ids_clean <- na.omit(ids[ids != "" & !is.na(ids)])
  if (length(ids_clean) == 0) {
    return("UNKNOWN")
  }
  n_total <- length(ids_clean)
  n_sample <- min(n_total, 1000)
  ids_sample <- sample(ids_clean, n_sample)
  if (mean(grepl("^ENS[A-Z0-9]*[FPTG]\\d{10,}(\\.\\d+)?$", ids_sample, ignore.case = TRUE)) > 0.8) {
    return("ENSEMBL")
  }
  if (mean(grepl("^[NX][CMPWTZ]_[0-9]+(\\.\\d+)?$", ids_sample, ignore.case = TRUE)) > 0.8) {
    return("REFSEQ")
  }
  if (mean(grepl("^[0-9]+$", ids_sample)) > 0.9) {
    return("ENTREZID")
  }
  if (mean(grepl("^([A-Za-z][A-Za-z0-9-]*[A-Za-z0-9]|[A-Za-z])$", ids_sample)) > 0.7 &&
    mean(grepl("^ENS", ids_sample, ignore.case = TRUE)) < 0.2 &&
    mean(grepl("^[NX][M_]", ids_sample, ignore.case = TRUE)) < 0.2 &&
    mean(grepl("^[0-9]+$", ids_sample)) < 0.2) {
    return("SYMBOL")
  }
  return("UNKNOWN")
}


degAnalysisUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 4,
      accordion(
        open = c("group_settings", "thresholds"),

        accordion_panel(
          title = "比較グループ設定",
          value = "group_settings",
          icon  = icon("users"),
          radioButtons(ns("analysis_type"), "解析タイプ:",
            choices = c("標準比較 (Pairwise / Interaction)" = "std", "多群比較 (LRT: ANOVA-like)" = "lrt"),
            selected = "std", inline = TRUE
          ),
          tags$div(class = "alert alert-warning", style = "padding: 8px; font-size: 0.85rem;",
            icon("lightbulb"), HTML(" <b>標準比較:</b> 2群比較・相互作用。<b>多群比較 (LRT):</b> 3群以上のANOVA的検定。")
          ),
          uiOutput(ns("degGroupSelectionUI"))
        ),

        accordion_panel(
          title = "解析手法・バッチ補正",
          value = "method_settings",
          icon  = icon("cogs"),
          radioButtons(ns("deg_method"), "統計アルゴリズム:",
            choices = c("edgeR (推奨)" = "edgeR", "DESeq2 (頑健性高)" = "DESeq2"),
            selected = "edgeR", inline = TRUE
          ),
          helpText(icon("info-circle"), " edgeRは高速。DESeq2はサンプル数が多い場合により頑健。"),
          hr(),
          checkboxInput(ns("use_batch"), "バッチ補正を行う (Include Batch in model)", value = FALSE),
          conditionalPanel(
            condition = paste0("input['", ns("use_batch"), "'] == true"),
            selectInput(ns("batch_col"), "バッチ項の列 (Batch column):", choices = NULL)
          )
        ),

        accordion_panel(
          title = "フィルタリング閾値",
          value = "thresholds",
          icon  = icon("filter"),
          helpText(icon("question-circle"), " 以下の閾値はテーブルとプロットに即時反映されます。"),
          radioButtons(ns("sig_metric"), "有意差の指標:",
            choices = c("FDR (adjusted P-value)" = "FDR", "P-value" = "PValue"),
            selected = "FDR", inline = TRUE
          ),
          conditionalPanel(
            condition = paste0("input['", ns("sig_metric"), "'] == 'FDR'"),
            numericInput(ns("degFDR"), "FDR 閾値:", value = 0.05, min = 0, max = 1, step = 0.01)
          ),
          conditionalPanel(
            condition = paste0("input['", ns("sig_metric"), "'] == 'PValue'"),
            numericInput(ns("degPValue"), "P-value 閾値:", value = 0.05, min = 0, max = 1, step = 0.01)
          ),
          numericInput(ns("degLogFC"), "Log2 Fold Change 閾値 (|LogFC| >):", value = 1, min = 0, step = 0.1)
        ),

        accordion_panel(
          title = "実行・ダウンロード",
          value = "run_download",
          icon  = icon("play"),
          actionButton(ns("runDEG"), "DEG解析実行", icon = icon("play"), class = "btn-primary w-100"),
          hr(),
          conditionalPanel(
            condition = paste0("input['", ns("analysis_type"), "'] == 'std'"),
            h6("ボルケーノプロット: ハイライト遺伝子 (任意)"),
            uiOutput(ns("highlightGenesUI")),
            hr()
          ),
          selectInput(ns("deg_id_display_type"), "結果テーブルの遺伝子IDタイプ:",
            choices = c("Gene Symbol" = "SYMBOL", "Entrez ID (内部ID)" = "ENTREZID"),
            selected = "SYMBOL"
          ),
          p(class = "text-muted small", "全遺伝子リスト（フィルタリングなし）をダウンロード:"),
          downloadButton(ns("downloadExcelResults"), "Excel (.xlsx)", icon = icon("file-excel"), class = "btn-sm w-100 mb-1"),
          br(),
          downloadButton(ns("downloadCsvResults"), "CSV (.csv)", icon = icon("file-csv"), class = "btn-sm w-100")
        )
      )
    ),
    mainPanel(
      width = 8,
      h4("解析サマリー"),
      uiOutput(ns("summary_boxes_ui")),
      withSpinner(verbatimTextOutput(ns("degSummary")), type = 6),
      hr(),
      conditionalPanel(
        condition = paste0("input['", ns("analysis_type"), "'] == 'std'"),
        h4("MAプロット (旧MDプロット)"),
        helpText("赤点・青点が有意に発現変動している遺伝子を示します。"),
        withSpinner(plotOutput(ns("degMDPlot")), type = 6),
        hr(),
        h4("ボルケーノプロット"),
        helpText("有意な発現変動遺伝子をLogFCと選択した指標でプロットします。"),
        withSpinner(plotlyOutput(ns("degVolcanoPlot")), type = 6),
        hr()
      ),
      h4("ヒートマップ (Top変動遺伝子)"),
      helpText("FDR/P-valueでソートされた上位遺伝子の発現パターンを表示します。"),
      numericInput(ns("heatmapTopN"), "表示する上位遺伝子数:", value = 50, min = 1, step = 1),
      numericInput(ns("heatmapFontSize"), "文字の大きさを変更:", value = 10, min = 1, step = 1),
      br(),
      downloadButton(ns("downloadHeatmapPDF"), "ヒートマップをPDFで保存 (.pdf)", icon = icon("file-pdf")),
      br(), br(),
      withSpinner(plotOutput(ns("degHeatmap"), height = "600px"), type = 6),
      hr(),
      h4("発現パターンのクラスタリング (k-means)"),
      helpText("有意変動遺伝子(DEGs)を類似したパターンのグループに分類し、群間のトレンドを可視化します。(多群比較で特に有効です)"),
      checkboxInput(ns("show_elbow"), "エルボー法を表示 (最適なクラスター数 k の探索)", value = FALSE),
      numericInput(ns("kmeans_k"), "クラスター数 (k):", value = 4, min = 2, max = 20, step = 1),
      actionButton(ns("send_to_go"), "現在のクラスターをGO解析へ送信", icon = icon("paper-plane")),
      downloadButton(ns("downloadKmeansCsv"), "結果をダウンロード (.csv)", icon = icon("file-csv")),
      br(), br(),
      withSpinner(plotOutput(ns("kmeansPlot"), height = "600px"), type = 6),
      hr(),
      h4("差次発現遺伝子テーブル (Top Tags)"),
      helpText("選択した指標とLogFCの閾値でフィルタリングされた結果が表示されます。"),
      withSpinner(DTOutput(ns("degResultTable")), type = 6),
      h5("Up-regulated Genes (表示IDタイプ適用)"),
      withSpinner(verbatimTextOutput(ns("upGenesList")), type = 6),
      h5("Down-regulated Genes (表示IDタイプ適用)"),
      withSpinner(verbatimTextOutput(ns("downGenesList")), type = 6)
    )
  )
}

degAnalysisServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- 名前空間衝突の解消 (ローカル束縛) ---
    # Shiny の loadSupport が R/ を自動再ソースし、AnnotationDbi/S4Vectors 等の
    # library() が app.R のグローバル再束縛後に dplyr の動詞を再び上書きするため、
    # モジュール実行環境(クロージャ)内で dplyr 版をローカルに束縛して確実に最優先化する。
    select <- dplyr::select
    rename <- dplyr::rename
    filter <- dplyr::filter
    mutate <- dplyr::mutate
    arrange <- dplyr::arrange
    slice  <- dplyr::slice

    detected_deg_input_id_type <- reactiveVal("ENTREZID")
    # 表示IDタイプの選択肢: GTFアップロード時は rv$gtf_id_choices を、無ければ既定3択を使う
    deg_id_choices <- function() {
      ch <- rv$gtf_id_choices
      if (is.null(ch)) ch <- gtf_display_id_choices(NULL)
      ch
    }
    observeEvent(rv$merged_data, {
      req(rv$merged_data)
      detected_deg_input_id_type("ENTREZID")
      ch <- deg_id_choices()
      updateSelectInput(session, "deg_id_display_type", choices = ch$choices, selected = ch$selected)
    })
    # GTFが (カウントデータの後に) アップロード/変更されたら表示IDタイプの選択肢も追従させる
    observeEvent(rv$gtf_id_choices, {
      ch <- deg_id_choices()
      updateSelectInput(session, "deg_id_display_type", choices = ch$choices, selected = ch$selected)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)
    # 共通ヘルパー annotate_display_ids に委譲 (GTFアノテーション優先 -> OrgDb -> 生ID)
    convert_entrez_ids_for_display <- function(entrez_ids, target_display_type, selected_species_code) {
      if (target_display_type == "ENTREZID" || is.null(target_display_type) || !nzchar(target_display_type)) {
        return(as.character(entrez_ids))
      }
      # 生物種未設定 (例: 外部DE結果アップロードで species 未指定) のときは
      # ID変換せず元のGeneidをそのまま表示する (req で停止させない)。
      if (is.null(selected_species_code) || !nzchar(selected_species_code)) {
        return(as.character(entrez_ids))
      }
      annotate_display_ids(entrez_ids, target_display_type, selected_species_code,
                           gene_annotation = rv$gene_annotation, orgdb_map = orgdb_species_map)
    }

    output$degGroupSelectionUI <- renderUI({
      req(rv$sample_metadata, input$analysis_type)
      active_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      factor_cols <- setdiff(colnames(active_metadata), c("id", "current_name", "active", "time"))

      # バッチ選択用の列を更新
      updateSelectInput(session, "batch_col", choices = factor_cols, selected = if (length(factor_cols) > 1) factor_cols[2] else factor_cols[1])

      if (input$analysis_type == "std") {
        tagList(
          selectizeInput(ns("std_group_cols"), "解析に使用するグループ列 (最大2列。1列=ペアワイズ, 2列=相互作用):",
            choices = factor_cols, selected = factor_cols[1], multiple = TRUE, options = list(maxItems = 2)
          ),
          uiOutput(ns("dynamic_std_ui"))
        )
      } else if (input$analysis_type == "lrt") {
        tagList(
          selectInput(ns("lrt_group_col"), "解析に使用するグループ列:", choices = factor_cols, selected = factor_cols[1]),
          uiOutput(ns("dynamic_lrt_ui"))
        )
      }
    })

    output$dynamic_std_ui <- renderUI({
      req(input$std_group_cols, rv$sample_metadata)
      active_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      cols <- input$std_group_cols

      if (length(cols) == 1) {
        c_groups <- unique(active_metadata[[cols[1]]])
        selected_target <- c_groups[1]
        selected_ref <- if (length(c_groups) > 1) c_groups[2] else c_groups[1]
        tagList(
          helpText("ペアワイズ比較: 選択されたグループで指定の2群間を比較します。"),
          selectInput(ns("target_group"), "Target (ターゲット):", choices = c_groups, selected = selected_target),
          selectInput(ns("reference_group"), "Reference (基準):", choices = c_groups, selected = selected_ref)
        )
      } else if (length(cols) == 2) {
        col1 <- cols[1]
        col2 <- cols[2]
        f1_levels <- unique(active_metadata[[col1]])
        f2_levels <- unique(active_metadata[[col2]])
        tagList(
          helpText("相互作用モデル (Interaction): 2要素の組み合わせによる特異的な変動を検出します。"),
          selectInput(ns("f1_ref"), paste0("要因1 (", col1, ") のベースライン*:"), choices = f1_levels, selected = f1_levels[1]),
          selectInput(ns("f1_target"), paste0("要因1 (", col1, ") のターゲット:"), choices = f1_levels, selected = if (length(f1_levels) > 1) f1_levels[2] else f1_levels[1]),
          selectInput(ns("f2_ref"), paste0("要因2 (", col2, ") のベースライン*:"), choices = f2_levels, selected = f2_levels[1]),
          selectInput(ns("f2_target"), paste0("要因2 (", col2, ") のターゲット:"), choices = f2_levels, selected = if (length(f2_levels) > 1) f2_levels[2] else f2_levels[1])
        )
      }
    })

    output$dynamic_lrt_ui <- renderUI({
      req(input$lrt_group_col, rv$sample_metadata)
      active_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      c_groups <- unique(active_metadata[[input$lrt_group_col]])
      tagList(
        checkboxGroupInput(ns("lrt_groups"), "解析に含めるレベル (3群以上推奨):",
          choices = c_groups, selected = c_groups, inline = TRUE
        ),
        helpText("選択したレベル間で発現が異なる遺伝子を検出します (ANOVA的)。")
      )
    })

    output$highlightGenesUI <- renderUI({
      all_entrez_ids_in_deg <- NULL
      if (!is.null(rv$deg_results) && !is.null(rv$deg_results$top_tags) && "Geneid" %in% colnames(rv$deg_results$top_tags$table)) {
        all_entrez_ids_in_deg <- unique(as.character(rv$deg_results$top_tags$table$Geneid))
      }
      display_type <- input$deg_id_display_type
      id_type_label <- switch(display_type,
        "SYMBOL" = "Gene Symbol",
        "ENTREZID" = "Entrez ID",
        "GENENAME" = "Gene Name",
        display_type
      )
      highlight_label <- paste0("ハイライトする遺伝子を選択 (", id_type_label, "で検索):")
      choices_for_ui <- NULL
      if (!is.null(all_entrez_ids_in_deg)) {
        choices_for_ui <- convert_entrez_ids_for_display(all_entrez_ids_in_deg, display_type, rv$selected_species)
        choices_for_ui <- sort(unique(na.omit(choices_for_ui[!grepl("\\(変換不可\\)$", choices_for_ui)])))
      }
      if (is.null(choices_for_ui)) {
        return(tags$p(em("DEG解析を実行すると、ここで遺伝子を選択できます。")))
      }
      selectizeInput(ns("highlight_genes_select"),
        label = highlight_label, choices = choices_for_ui,
        selected = NULL, multiple = TRUE, options = list(placeholder = "入力または選択...", plugins = list("remove_button"))
      )
    })

    observeEvent(input$runDEG, {
      req(rv$merged_data, rv$sample_metadata, input$deg_method, input$analysis_type)
      is_lrt <- input$analysis_type == "lrt"
      is_interaction <- input$analysis_type == "std" && length(input$std_group_cols) == 2
      is_pairwise <- input$analysis_type == "std" && length(input$std_group_cols) == 1

      if (is_pairwise) req(input$target_group, input$reference_group)
      if (is_lrt) req(input$lrt_groups, input$lrt_group_col)
      if (is_interaction) req(input$std_group_cols, input$f1_ref, input$f2_ref)

      sig_metric_val <- if (input$sig_metric == "FDR") input$degFDR else input$degPValue
      logfc_val <- input$degLogFC
      shiny::validate(
        shiny::need(is.numeric(sig_metric_val) && sig_metric_val >= 0 && sig_metric_val <= 1, "エラー: 有意差の閾値は0から1の間の数値で入力してください。"),
        shiny::need(is.numeric(logfc_val) && logfc_val >= 0, "エラー: LogFC閾値は0以上の数値で入力してください。")
      )

      active_samples_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      shiny::validate(shiny::need(nrow(active_samples_metadata) > 0, "解析対象のサンプルが選択されていません。"))
      active_sample_names <- active_samples_metadata$current_name

      counts_df_iso_entrez <- isolate(rv$merged_data)
      counts_df_iso_entrez_active <- counts_df_iso_entrez[, c("Geneid", intersect(colnames(counts_df_iso_entrez), active_sample_names)), drop = FALSE]
      samples_metadata_iso <- active_samples_metadata

      keep_vector_initial_filter <- isolate(rv$filtered_keep)
      shiny::validate(shiny::need(is.data.frame(counts_df_iso_entrez_active) && "Geneid" %in% colnames(counts_df_iso_entrez_active), "カウントデータが不正です。"))

      count_matrix_full_entrez <- as.matrix(counts_df_iso_entrez_active[, -1])
      rownames(count_matrix_full_entrez) <- counts_df_iso_entrez_active$Geneid
      shiny::validate(shiny::need(length(keep_vector_initial_filter) == nrow(isolate(rv$merged_data)), "フィルタリング情報が古いため、フィルタリングタブを再確認してください。"))

      count_matrix_ui_filtered_entrez <- count_matrix_full_entrez[keep_vector_initial_filter, , drop = FALSE]
      shiny::validate(shiny::need(nrow(count_matrix_ui_filtered_entrez) > 0, "フィルタリングの結果、遺伝子が残りませんでした。"))

      if (is_interaction) {
        valid_samples <- colnames(count_matrix_ui_filtered_entrez)
        count_matrix_for_deg <- count_matrix_ui_filtered_entrez[, valid_samples, drop = FALSE]
        metadata_for_deg <- samples_metadata_iso %>% filter(current_name %in% valid_samples)

        col1 <- input$std_group_cols[1]
        col2 <- input$std_group_cols[2]
        f1 <- factor(metadata_for_deg[[col1]])
        f1 <- droplevels(f1)
        if (input$f1_ref %in% levels(f1)) f1 <- relevel(f1, ref = input$f1_ref)

        f2 <- factor(metadata_for_deg[[col2]])
        f2 <- droplevels(f2)
        if (input$f2_ref %in% levels(f2)) f2 <- relevel(f2, ref = input$f2_ref)

        group_factor <- NULL
        comparison_label <- paste("Interaction:", col1, "*", col2)
      } else if (is_lrt) {
        col_lrt <- input$lrt_group_col
        selected_lrt_groups <- input$lrt_groups
        shiny::validate(shiny::need(length(selected_lrt_groups) >= 2, "LRT解析には2群以上を選択してください。"))
        samples_to_keep <- samples_metadata_iso$current_name[samples_metadata_iso[[col_lrt]] %in% selected_lrt_groups]
        valid_samples <- intersect(colnames(count_matrix_ui_filtered_entrez), samples_to_keep)
        shiny::validate(shiny::need(length(valid_samples) >= 2, "選択されたグループにサンプルが不足しています。"))
        count_matrix_for_deg <- count_matrix_ui_filtered_entrez[, valid_samples, drop = FALSE]
        metadata_for_deg <- samples_metadata_iso %>% filter(current_name %in% valid_samples)
        group_factor <- factor(metadata_for_deg[[col_lrt]])

        target_group_name <- NULL
        reference_group_name <- NULL
        comparison_label <- paste(selected_lrt_groups, collapse = " / ")
      } else if (is_pairwise) {
        col1 <- input$std_group_cols[1]
        target_group_name <- input$target_group
        reference_group_name <- input$reference_group
        samples_to_keep <- samples_metadata_iso$current_name[samples_metadata_iso[[col1]] %in% c(target_group_name, reference_group_name)]
        valid_samples <- intersect(colnames(count_matrix_ui_filtered_entrez), samples_to_keep)
        shiny::validate(shiny::need(length(valid_samples) >= 2 && length(unique(samples_metadata_iso[[col1]][samples_metadata_iso$current_name %in% valid_samples])) == 2, "各比較グループに最低1サンプル必要です。"))
        count_matrix_for_deg <- count_matrix_ui_filtered_entrez[, valid_samples, drop = FALSE]
        metadata_for_deg <- samples_metadata_iso %>% filter(current_name %in% valid_samples)
        group_factor <- factor(metadata_for_deg[[col1]]) %>% relevel(ref = reference_group_name)
        comparison_label <- paste0(target_group_name, " vs ", reference_group_name)
      }

      deg_results_list <- tryCatch(
        {
          # ★★★ 解析タイプ・手法による分岐 ★★★
          if (is_interaction) {
            # === 相互作用モデル (強制的にDESeq2を使用) ===
            colData <- data.frame(row.names = metadata_for_deg$current_name, factor1 = f1, factor2 = f2)
            if (input$use_batch && !is.null(input$batch_col)) {
              colData$batch <- factor(metadata_for_deg[[input$batch_col]])
              design_formula <- ~ batch + factor1 + factor2 + factor1:factor2
            } else {
              design_formula <- ~ factor1 + factor2 + factor1:factor2
            }
            dds <- DESeqDataSetFromMatrix(countData = round(count_matrix_for_deg), colData = colData, design = design_formula)
            keep <- rowSums(counts(dds)) >= 10
            dds <- dds[keep, ]
            shiny::validate(shiny::need(nrow(dds) > 0, "フィルタリングの結果、遺伝子が0になりました。"))
            dds <- DESeq(dds)
            res_names <- resultsNames(dds)

            # ユーザーが指定したターゲットレベルに基づくInteraction項の名前を構築
            # 例: "factor1one.factor2plus" または "factor1_one_vs_ctrl.factor2_plus_vs_wo" など
            # DESeq2内部では通常 "factor1[Target].factor2[Target]" の形式になります。
            expected_name <- paste0("factor1", input$f1_target, ".factor2", input$f2_target)

            # もし期待する名前が見つからなければ、部分一致で探す
            if (expected_name %in% res_names) {
              int_name <- expected_name
            } else {
              int_name <- res_names[grepl(paste0("factor1", input$f1_target, ".*factor2", input$f2_target, "|factor2", input$f2_target, ".*factor1", input$f1_target), res_names)]
              if (length(int_name) > 1) int_name <- int_name[1]
              if (length(int_name) == 0) {
                # それでも見つからなければフォールバック
                int_name <- res_names[grepl("factor1.*factor2|factor2.*factor1", res_names)]
                if (length(int_name) > 1) int_name <- int_name[1]
                if (length(int_name) == 0) int_name <- res_names[length(res_names)]
              }
            }

            res <- results(dds, name = int_name)

            res_df <- as.data.frame(res) %>%
              rownames_to_column("Geneid") %>%
              rename(logFC = log2FoldChange, PValue = pvalue, FDR = padj)
            res_df$FDR[is.na(res_df$FDR)] <- 1

            sig_col_name <- input$sig_metric
            up_genes <- res_df %>%
              filter(!is.na(!!sym(sig_col_name)), !!sym(sig_col_name) < sig_metric_val, !is.na(logFC), logFC > logfc_val) %>%
              pull(Geneid) %>%
              unique()
            down_genes <- res_df %>%
              filter(!is.na(!!sym(sig_col_name)), !!sym(sig_col_name) < sig_metric_val, !is.na(logFC), logFC < -logfc_val) %>%
              pull(Geneid) %>%
              unique()

            comp_lbl <- c(paste("Interaction", int_name), "Baseline")

            list(
              normalized_counts = log2(counts(dds, normalized = TRUE) + 1), analysis_method = "DESeq2", analysis_type = "interaction",
              deseq_res = res, top_tags = list(table = res_df), comparison = comp_lbl,
              significant_up_genes = up_genes, significant_down_genes = down_genes, significant_genes = c(up_genes, down_genes),
              background_genes_original = rownames(dds), metric_at_run = sig_col_name,
              sig_threshold_at_run = sig_metric_val, logfc_threshold_at_run = logfc_val, status = "analysis_completed"
            )
          } else if (is_lrt) {
            # === LRT (多群比較) ===
            if (input$deg_method == "edgeR") {
              dge <- DGEList(counts = count_matrix_for_deg, group = group_factor, genes = data.frame(Geneid = rownames(count_matrix_for_deg)))
              if (input$use_batch && !is.null(input$batch_col)) {
                batch <- factor(metadata_for_deg[[input$batch_col]])
                design <- model.matrix(~ batch + group_factor)
                coef_indices <- (nlevels(batch) + 1):ncol(design)
              } else {
                design <- model.matrix(~ group_factor)
                coef_indices <- 2:ncol(design)
              }
              dge <- dge[filterByExpr(dge, design = design), , keep.lib.sizes = FALSE]
              shiny::validate(shiny::need(nrow(dge) > 0, "filterByExprの結果、遺伝子が0になりました。"))
              dge <- calcNormFactors(dge)
              dge <- estimateDisp(dge, design, robust = TRUE)
              fit <- glmQLFit(dge, design, robust = TRUE)
              qlf <- glmQLFTest(fit, coef = coef_indices)
              top_tags_result <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none")
              tt <- as.data.frame(top_tags_result$table)
              sig_col <- input$sig_metric
              sig_genes <- tt %>%
                filter(!!sym(sig_col) < sig_metric_val) %>%
                pull(Geneid) %>%
                unique()
              list(
                normalized_counts = cpm(dge, log = TRUE), analysis_method = "edgeR", analysis_type = "lrt",
                qlf = qlf, top_tags = top_tags_result, comparison = comparison_label, groups = input$lrt_groups,
                significant_genes = sig_genes, significant_up_genes = character(0), significant_down_genes = character(0),
                background_genes_original = dge$genes$Geneid, metric_at_run = sig_col,
                sig_threshold_at_run = sig_metric_val, logfc_threshold_at_run = logfc_val, status = "analysis_completed"
              )
            } else { # DESeq2 LRT
              colData <- data.frame(row.names = metadata_for_deg$current_name, group = group_factor)
              if (input$use_batch && !is.null(input$batch_col)) {
                colData$batch <- factor(metadata_for_deg[[input$batch_col]])
                design_formula <- ~ batch + group
                reduced_formula <- ~ batch
              } else {
                design_formula <- ~ group
                reduced_formula <- ~ 1
              }
              dds <- DESeqDataSetFromMatrix(countData = round(count_matrix_for_deg), colData = colData, design = design_formula)
              dds <- dds[rowSums(counts(dds)) >= 10, ]
              shiny::validate(shiny::need(nrow(dds) > 0, "フィルタリングの結果、遺伝子が0になりました。"))
              dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
              res <- results(dds)
              res_df <- as.data.frame(res) %>%
                rownames_to_column("Geneid") %>%
                rename(PValue = pvalue, FDR = padj)
              res_df$FDR[is.na(res_df$FDR)] <- 1
              sig_col <- input$sig_metric
              sig_genes <- res_df %>%
                filter(!is.na(!!sym(sig_col)), !!sym(sig_col) < sig_metric_val) %>%
                pull(Geneid) %>%
                unique()
              list(
                normalized_counts = log2(counts(dds, normalized = TRUE) + 1), analysis_method = "DESeq2", analysis_type = "lrt",
                deseq_res = res, top_tags = list(table = res_df), comparison = comparison_label, groups = input$lrt_groups,
                significant_genes = sig_genes, significant_up_genes = character(0), significant_down_genes = character(0),
                background_genes_original = rownames(dds), metric_at_run = sig_col,
                sig_threshold_at_run = sig_metric_val, logfc_threshold_at_run = logfc_val, status = "analysis_completed"
              )
            }
          } else if (input$deg_method == "edgeR") {
            dge <- DGEList(counts = count_matrix_for_deg, group = group_factor, genes = data.frame(Geneid = rownames(count_matrix_for_deg)))
            
            if (input$use_batch && !is.null(input$batch_col)) {
                batch <- factor(metadata_for_deg[[input$batch_col]])
                design <- model.matrix(~ 0 + batch + group_factor)
                # group_factorの最後のレベル（target_group_name）に対応する係数とreferenceに対応する係数の差を見る
                # しかし ~ 0 + batch + group_factor の場合、group_factorのベースラインはbatchの各レベルに吸収される
                # より安全なのは、batchをブロック因子として扱う DESeq2的なやり方か、
                # または ~ batch + group_factor としておいて coef = ncol(design) で指定するか。
                design <- model.matrix(~ batch + group_factor)
                coef_idx <- ncol(design) 
                
                dge <- dge[filterByExpr(dge, design = design), , keep.lib.sizes = FALSE]
                shiny::validate(shiny::need(nrow(dge) > 0, "filterByExprの結果、遺伝子が0になりました。"))
                dge <- calcNormFactors(dge)
                dge <- estimateDisp(dge, design, robust = TRUE)
                fit <- glmQLFit(dge, design, robust = TRUE)
                # ~ batch + group_factor の場合、最後の列が target vs reference になる
                qlf <- glmQLFTest(fit, coef = coef_idx)
            } else {
                design_filter <- model.matrix(~group_factor)
                dge <- dge[filterByExpr(dge, design = design_filter), , keep.lib.sizes = FALSE]
                shiny::validate(shiny::need(nrow(dge) > 0, "filterByExprの結果、遺伝子が0になりました。フィルタリング条件を緩めてください。"))
                dge <- calcNormFactors(dge)
                design <- model.matrix(~ 0 + group, data = dge$samples)
                colnames(design) <- levels(dge$samples$group)
                contrast_str <- paste(target_group_name, "-", reference_group_name)
                my_contrast <- makeContrasts(contrasts = contrast_str, levels = design)
                dge <- estimateDisp(dge, design, robust = TRUE)
                fit <- glmQLFit(dge, design, robust = TRUE)
                qlf <- glmQLFTest(fit, contrast = my_contrast)
            }
            top_tags_result <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none")
            top_tags_table <- as.data.frame(top_tags_result$table)

            sig_col_name <- input$sig_metric
            up_genes <- top_tags_table %>%
              filter(!!sym(sig_col_name) < sig_metric_val, logFC > logfc_val) %>%
              pull(Geneid) %>%
              unique()
            down_genes <- top_tags_table %>%
              filter(!!sym(sig_col_name) < sig_metric_val, logFC < -logfc_val) %>%
              pull(Geneid) %>%
              unique()

            list(
              normalized_counts = cpm(dge, log = TRUE),
              analysis_method = "edgeR", analysis_type = "pairwise",
              qlf = qlf, top_tags = top_tags_result,
              comparison = c(target_group_name, reference_group_name),
              significant_up_genes = up_genes, significant_down_genes = down_genes,
              significant_genes = c(up_genes, down_genes),
              background_genes_original = dge$genes$Geneid,
              metric_at_run = sig_col_name,
              sig_threshold_at_run = sig_metric_val,
              logfc_threshold_at_run = logfc_val,
              status = "analysis_completed"
            )
          } else { # DESeq2
            colData <- data.frame(row.names = metadata_for_deg$current_name, group = group_factor)
            if (input$use_batch && !is.null(input$batch_col)) {
              colData$batch <- factor(metadata_for_deg[[input$batch_col]])
              design_formula <- ~ batch + group
            } else {
              design_formula <- ~ group
            }
            dds <- DESeqDataSetFromMatrix(countData = round(count_matrix_for_deg), colData = colData, design = design_formula)
            keep <- rowSums(counts(dds)) >= 10
            dds <- dds[keep, ]
            shiny::validate(shiny::need(nrow(dds) > 0, "フィルタリングの結果、遺伝子が0になりました。"))

            dds <- DESeq(dds)
            res <- results(dds, contrast = c("group", target_group_name, reference_group_name))

            res_df <- as.data.frame(res) %>%
              rownames_to_column("Geneid") %>%
              rename(logFC = log2FoldChange, PValue = pvalue, FDR = padj)
            res_df$FDR[is.na(res_df$FDR)] <- 1

            top_tags_like_object <- list(table = res_df)

            sig_col_name <- input$sig_metric
            up_genes <- res_df %>%
              filter(!is.na(!!sym(sig_col_name)), !!sym(sig_col_name) < sig_metric_val, !is.na(logFC), logFC > logfc_val) %>%
              pull(Geneid) %>%
              unique()
            down_genes <- res_df %>%
              filter(!is.na(!!sym(sig_col_name)), !!sym(sig_col_name) < sig_metric_val, !is.na(logFC), logFC < -logfc_val) %>%
              pull(Geneid) %>%
              unique()

            list(
              normalized_counts = log2(counts(dds, normalized = TRUE) + 1),
              analysis_method = "DESeq2", analysis_type = "pairwise",
              deseq_res = res,
              top_tags = top_tags_like_object,
              comparison = c(target_group_name, reference_group_name),
              significant_up_genes = up_genes, significant_down_genes = down_genes,
              significant_genes = c(up_genes, down_genes),
              background_genes_original = rownames(dds),
              metric_at_run = sig_col_name,
              sig_threshold_at_run = sig_metric_val,
              logfc_threshold_at_run = logfc_val,
              status = "analysis_completed"
            )
          }
        },
        error = function(e) {
          call_str <- tryCatch(paste(deparse(conditionCall(e)), collapse = " "), error = function(e2) "(unknown call)")
          message("DEG解析エラー詳細 - call: ", call_str, " / message: ", e$message)
          showNotification(paste0("DEG解析エラー: ", e$message, "\n[発生箇所: ", call_str, "]"), type = "error", duration = 20)
          NULL
        }
      )
      rv$deg_results <- deg_results_list
    })

    # === アップロードした DE 結果表を読み込み rv$deg_results に流し込む =========
    # カウント再計算なしで DESeq2_all / edgeR の結果表を直接受け取り、
    # 既存の Volcano/MA・結果表・下流タブ(GSEA/GO)で利用できる形に整形する。
    # 列名は .fe_find_col で大小無視の自動検出 (module_figure_enrichment.R のヘルパを再利用)。
    observeEvent(input$loadUploadedDE, {
      req(input$uploadDE)
      df <- .fe_read_table(input$uploadDE$datapath)
      shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "ファイルを読み込めませんでした (空または不正な形式)。"))

      id_col   <- .fe_find_col(df, c("Geneid", "gene", "gene_id", "ID", "GeneSymbol"))
      lfc_col  <- .fe_find_col(df, c("log2FoldChange", "logFC", "log2FC"))
      p_col    <- .fe_find_col(df, c("pvalue", "PValue", "p.value", "P.Value"))
      fdr_col  <- .fe_find_col(df, c("padj", "FDR", "adj.P.Val", "p.adjust", "qvalue"))
      base_col <- .fe_find_col(df, c("baseMean", "logCPM", "AveExpr"))
      stat_col <- .fe_find_col(df, c("stat", "t", "F", "LR"))
      lfcse_col <- .fe_find_col(df, c("lfcSE"))

      shiny::validate(
        shiny::need(!is.null(id_col),  "遺伝子ID列 (Geneid 等) が見つかりません。"),
        shiny::need(!is.null(lfc_col), "log2FoldChange / logFC 列が見つかりません。"),
        shiny::need(!is.null(p_col) || !is.null(fdr_col), "pvalue または padj 列が見つかりません。")
      )

      res_df <- data.frame(Geneid = as.character(df[[id_col]]), stringsAsFactors = FALSE)
      res_df$logFC <- suppressWarnings(as.numeric(df[[lfc_col]]))
      if (!is.null(base_col))  res_df$baseMean <- suppressWarnings(as.numeric(df[[base_col]]))
      if (!is.null(lfcse_col)) res_df$lfcSE    <- suppressWarnings(as.numeric(df[[lfcse_col]]))
      if (!is.null(stat_col))  res_df$stat     <- suppressWarnings(as.numeric(df[[stat_col]]))
      res_df$PValue <- if (!is.null(p_col)) suppressWarnings(as.numeric(df[[p_col]])) else NA_real_
      res_df$FDR    <- if (!is.null(fdr_col)) suppressWarnings(as.numeric(df[[fdr_col]])) else NA_real_
      # padj列が無い場合は p値から BH 補正で補完
      if (is.null(fdr_col) && !is.null(p_col)) res_df$FDR <- p.adjust(res_df$PValue, method = "BH")
      res_df$FDR[is.na(res_df$FDR)] <- 1
      res_df <- res_df[!is.na(res_df$Geneid) & nzchar(res_df$Geneid), , drop = FALSE]
      shiny::validate(shiny::need(nrow(res_df) > 0, "有効な遺伝子行がありません。"))

      # --- 比較ラベル: ファイル名から自動抽出し、テキスト欄があれば上書き ---
      base_name <- tools::file_path_sans_ext(input$uploadDE$name)
      base_name <- sub("_(deseq2|edger)?_?(all|results|res)$", "", base_name, ignore.case = TRUE)
      tgt <- "Group1"; ref <- "Group2"
      if (grepl("_vs_", base_name, ignore.case = TRUE)) {
        parts <- strsplit(base_name, "(?i)_vs_", perl = TRUE)[[1]]
        if (length(parts) >= 2) {
          tgt <- sub("^[Cc][0-9]+_", "", parts[1])  # 先頭の contrast プレフィックス (C1_ 等) 除去
          ref <- parts[2]
        }
      }
      if (!is.null(input$up_target)    && nzchar(trimws(input$up_target)))    tgt <- trimws(input$up_target)
      if (!is.null(input$up_reference) && nzchar(trimws(input$up_reference))) ref <- trimws(input$up_reference)

      # --- 現在の閾値で有意遺伝子を算出 ---
      sig_metric_val <- if (input$sig_metric == "FDR") input$degFDR else input$degPValue
      logfc_val <- input$degLogFC
      sig_col_name <- input$sig_metric
      up_genes <- res_df %>%
        filter(!is.na(!!sym(sig_col_name)), !!sym(sig_col_name) < sig_metric_val, !is.na(logFC), logFC > logfc_val) %>%
        pull(Geneid) %>% unique()
      down_genes <- res_df %>%
        filter(!is.na(!!sym(sig_col_name)), !!sym(sig_col_name) < sig_metric_val, !is.na(logFC), logFC < -logfc_val) %>%
        pull(Geneid) %>% unique()

      rv$deg_results <- list(
        normalized_counts = NULL,
        analysis_method = "DESeq2", analysis_type = "pairwise",
        deseq_res = NULL,
        top_tags = list(table = res_df),
        comparison = c(tgt, ref),
        significant_up_genes = up_genes, significant_down_genes = down_genes,
        significant_genes = c(up_genes, down_genes),
        background_genes_original = res_df$Geneid,
        metric_at_run = sig_col_name,
        sig_threshold_at_run = sig_metric_val,
        logfc_threshold_at_run = logfc_val,
        status = "analysis_completed", source = "upload"
      )
      rv$background_genes_original <- res_df$Geneid

      showNotification(
        sprintf("DE結果を読み込みました: %d 遺伝子 / 比較 '%s vs %s' / Up %d, Down %d",
                nrow(res_df), tgt, ref, length(up_genes), length(down_genes)),
        type = "message", duration = 8
      )
    })

    output$summary_boxes_ui <- renderUI({
      req(rv$deg_results, rv$deg_results$status == "analysis_completed")
      res <- rv$deg_results
      is_lrt <- !is.null(res$analysis_type) && res$analysis_type == "lrt"

      if (is_lrt) {
        res_table <- as.data.frame(res$top_tags$table)
        n_sig <- sum(res_table[[res$metric_at_run]] < res$sig_threshold_at_run, na.rm = TRUE)
        layout_columns(
          fill = FALSE,
          value_box(title = "Significant Genes", value = n_sig,
                    showcase = icon("dna"), theme = "primary"),
          value_box(title = "Tested Genes", value = nrow(res_table),
                    showcase = icon("list"), theme = "secondary")
        )
      } else {
        layout_columns(
          fill = FALSE,
          value_box(title = "Up-regulated", value = length(res$significant_up_genes),
                    showcase = icon("arrow-trend-up"), theme = "danger"),
          value_box(title = "Down-regulated", value = length(res$significant_down_genes),
                    showcase = icon("arrow-trend-down"), theme = "info"),
          value_box(title = "Tested Genes", value = length(res$background_genes_original),
                    showcase = icon("dna"), theme = "secondary")
        )
      }
    })

    output$degSummary <- renderPrint({
      req(rv$deg_results, rv$deg_results$status == "analysis_completed")
      res <- rv$deg_results
      cat("解析手法:", res$analysis_method, "\n")
      is_lrt <- !is.null(res$analysis_type) && res$analysis_type == "lrt"
      is_int <- !is.null(res$analysis_type) && res$analysis_type == "interaction"
      if (is_lrt) {
        cat("解析タイプ: LRT (多群比較)\n")
        cat("比較グループ:", res$comparison, "\n\n")
      } else if (is_int) {
        cat("解析タイプ: Interaction (相互作用モデル)\n")
        cat("評価項:", res$comparison[1], "\n\n")
      } else {
        cat("比較:", paste0(res$comparison[1], "_vs_", res$comparison[2]), "\n\n")
      }
      current_sig_metric <- input$sig_metric
      current_sig_threshold <- if (current_sig_metric == "FDR") input$degFDR else input$degPValue
      current_logfc_threshold <- input$degLogFC

      if (is_lrt) {
        cat(" 有意変動遺伝子サマリー (", current_sig_metric, " < ", current_sig_threshold, "):\n")
        res_table <- as.data.frame(res$top_tags$table)
        n_sig <- sum(res_table[[current_sig_metric]] < current_sig_threshold, na.rm = TRUE)
        cat("  有意遺伝子:", n_sig, "/", nrow(res_table), "\n")
        cat("\n背景遺伝子数 (解析対象):", length(res$background_genes_original), "\n")
      } else {
        cat(" 有意変動遺伝子サマリー (", current_sig_metric, " < ", current_sig_threshold, ", |LogFC| > ", current_logfc_threshold, "):\n")

        res_table <- as.data.frame(res$top_tags$table)
        up <- res_table %>%
          filter(!is.na(!!sym(current_sig_metric)), !!sym(current_sig_metric) < current_sig_threshold, !is.na(logFC), logFC > current_logfc_threshold) %>%
          nrow()
        down <- res_table %>%
          filter(!is.na(!!sym(current_sig_metric)), !!sym(current_sig_metric) < current_sig_threshold, !is.na(logFC), logFC < -current_logfc_threshold) %>%
          nrow()
        notsig <- nrow(res_table) - up - down

        summary_mat <- matrix(c(down, notsig, up), ncol = 1, dimnames = list(c("Down", "NotSig", "Up"), paste(res$comparison[1], "vs", res$comparison[2])))
        print(summary_mat)

        cat("\n背景遺伝子数 (解析対象):", length(res$background_genes_original), "\n")

        cat("\n[注意] プロットやテーブルは左メニューの閾値に連動して即座に更新されますが、\n下流タブ(GO解析など)には「DEG解析実行」を押した時点のリストが送信されます。\n")
        cat("(下流への送信時点の閾値: ", res$metric_at_run, " < ", res$sig_threshold_at_run, " & |LogFC| > ", res$logfc_threshold_at_run, ")\n", sep = "")
      }
    })

    filtered_and_converted_table <- reactive({
      req(rv$deg_results, input$deg_id_display_type)
      res_table <- as.data.frame(rv$deg_results$top_tags$table)

      sig_metric_filter <- input$sig_metric
      sig_threshold_filter <- if (sig_metric_filter == "FDR") input$degFDR else input$degPValue
      logfc_threshold_filter <- input$degLogFC

      is_lrt <- !is.null(rv$deg_results$analysis_type) && rv$deg_results$analysis_type == "lrt"
      is_int <- !is.null(rv$deg_results$analysis_type) && rv$deg_results$analysis_type == "interaction"
      if (is_lrt) {
        filtered_table <- res_table %>%
          filter(!is.na(!!sym(sig_metric_filter)), !!sym(sig_metric_filter) < sig_threshold_filter) %>%
          arrange(!!sym(sig_metric_filter))
      } else {
        filtered_table <- res_table %>%
          filter(!is.na(!!sym(sig_metric_filter)), !!sym(sig_metric_filter) < sig_threshold_filter, !is.na(logFC), abs(logFC) > logfc_threshold_filter) %>%
          arrange(!!sym(sig_metric_filter))
      }

      if (nrow(filtered_table) == 0) {
        return(data.frame(メッセージ = "指定した閾値で有意な遺伝子はありません。"))
      }

      display_table <- data.table::copy(filtered_table)
      display_table[[1]] <- convert_entrez_ids_for_display(display_table$Geneid, input$deg_id_display_type, rv$selected_species)
      display_col_name <- names(c("Gene Symbol" = "SYMBOL", "Entrez ID" = "ENTREZID", "Gene Name" = "GENENAME"))[c("SYMBOL", "ENTREZID", "GENENAME") == input$deg_id_display_type]
      colnames(display_table)[1] <- if (length(display_col_name) > 0) display_col_name else input$deg_id_display_type

      display_table %>% select(any_of(c(colnames(display_table)[1], "logFC", "logCPM", "baseMean", "PValue", "FDR", "LR", "F", "stat")))
    })

    output$degResultTable <- renderDT({
      final_table <- filtered_and_converted_table()
      req(final_table)
      if ("メッセージ" %in% colnames(final_table)) {
        return(datatable(final_table, options = list(dom = "t")))
      }
      datatable(final_table, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE)) %>%
        formatRound(columns = intersect(c("logFC", "logCPM", "baseMean"), colnames(final_table)), digits = 3) %>%
        formatSignif(columns = intersect(c("PValue", "FDR"), colnames(final_table)), digits = 3)
    })

    output$degMDPlot <- renderPlot({
      req(rv$deg_results)
      res <- rv$deg_results
      is_lrt <- !is.null(res$analysis_type) && res$analysis_type == "lrt"
      is_int <- !is.null(res$analysis_type) && res$analysis_type == "interaction"
      if (is_lrt) {
        plot.new()
        text(0.5, 0.5, "LRTモードではMAプロットは利用できません。\nヒートマップとテーブルをご確認ください。", cex = 1.3)
        return(NULL)
      }
      plot_title <- if (is_int) res$comparison[1] else paste0(res$comparison[1], "_vs_", res$comparison[2])

      if (res$analysis_method == "edgeR") {
        req(res$qlf)
        sig_metric_plot <- input$sig_metric
        sig_threshold_plot <- if (sig_metric_plot == "FDR") input$degFDR else input$degPValue
        logfc_threshold_plot <- input$degLogFC
        status_vec <- decideTests(res$qlf, p.value = sig_threshold_plot, lfc = logfc_threshold_plot, adjust.method = if (sig_metric_plot == "FDR") "BH" else "none")
        plotMD(res$qlf, status = status_vec, values = c(1, -1), col = c("red", "blue"), legend = "topright", main = plot_title)
        abline(h = c(-logfc_threshold_plot, logfc_threshold_plot), col = "dodgerblue", lty = 2)
      } else { # DESeq2
        req(res$top_tags$table)
        res_df <- as.data.frame(res$top_tags$table)
        sig_metric_plot <- input$sig_metric
        sig_threshold_plot <- if (sig_metric_plot == "FDR") input$degFDR else input$degPValue
        logfc_threshold_plot <- input$degLogFC

        is_sig_up <- !is.na(res_df[[sig_metric_plot]]) & res_df[[sig_metric_plot]] < sig_threshold_plot & !is.na(res_df$logFC) & res_df$logFC > logfc_threshold_plot
        is_sig_down <- !is.na(res_df[[sig_metric_plot]]) & res_df[[sig_metric_plot]] < sig_threshold_plot & !is.na(res_df$logFC) & res_df$logFC < -logfc_threshold_plot

        col_vec <- rep(rgb(0, 0, 0, 0.3), nrow(res_df))
        col_vec[is_sig_up] <- "red"
        col_vec[is_sig_down] <- "blue"

        # 外れ値への対応と描画
        y_val <- res_df$logFC
        y_val_capped <- pmax(-5, pmin(5, y_val))
        pch_vec <- ifelse(y_val > 5 | y_val < -5, 17, 16) # 枠外は三角形

        # log=xを指定するとき、x<=0があると描画がバグるため除外または処理
        x_val <- pmax(res_df$baseMean, 0.1)

        plot(x_val, y_val_capped,
          pch = pch_vec, cex = 0.5, col = col_vec, log = "x",
          xlab = "Mean of normalized counts (baseMean)", ylab = "Log2 Fold Change (capped at ±5)",
          main = paste("MA Plot (DESeq2) -", plot_title), ylim = c(-5, 5)
        )

        abline(h = c(-logfc_threshold_plot, 0, logfc_threshold_plot), col = c("dodgerblue", "grey", "dodgerblue"), lty = c(2, 1, 2))
        legend("topright", legend = c("Up", "Down", "Not Sig"), col = c("red", "blue", "black"), pch = c(16, 16, 1))
      }
    })

    volcano_data_reactive <- reactive({
      req(rv$deg_results, input$deg_id_display_type)
      is_lrt <- !is.null(rv$deg_results$analysis_type) && rv$deg_results$analysis_type == "lrt"
      is_int <- !is.null(rv$deg_results$analysis_type) && rv$deg_results$analysis_type == "interaction"
      if (is_lrt) {
        return(NULL)
      }
      res_table <- as.data.frame(rv$deg_results$top_tags$table)

      sig_metric_volcano <- input$sig_metric
      sig_threshold_volcano <- if (sig_metric_volcano == "FDR") input$degFDR else input$degPValue
      logfc_threshold_volcano <- input$degLogFC

      res_table %>%
        filter(!is.na(logFC), !is.na(.data[[sig_metric_volcano]])) %>%
        mutate(
          DisplayGeneid = convert_entrez_ids_for_display(Geneid, input$deg_id_display_type, rv$selected_species),
          SignificanceValue = .data[[sig_metric_volcano]],
          negLog10Sig = -log10(pmax(SignificanceValue, .Machine$double.xmin)),
          SignificanceCategory = case_when(
            .data[[sig_metric_volcano]] < sig_threshold_volcano & logFC > logfc_threshold_volcano ~ "Up-regulated",
            .data[[sig_metric_volcano]] < sig_threshold_volcano & logFC < -logfc_threshold_volcano ~ "Down-regulated",
            TRUE ~ "Not Significant"
          ) %>%
            factor(levels = c("Up-regulated", "Down-regulated", "Not Significant")),
          Highlight = ifelse(DisplayGeneid %in% (input$highlight_genes_select %||% ""), "Highlight", "Normal"),
          tooltip_text = paste(
            "Gene:", DisplayGeneid, "<br>logFC:", round(logFC, 3),
            "<br>PValue:", format(PValue, scientific = T, digits = 3),
            "<br>FDR:", format(FDR, scientific = T, digits = 3)
          )
        )
    })

    output$degVolcanoPlot <- renderPlotly({
      plot_data <- volcano_data_reactive()
      if (is.null(plot_data)) {
        p <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "LRTモードではボルケーノプロットは\n利用できません。", size = 6) +
          theme_void()
        return(ggplotly(p))
      }
      req(plot_data)

      sig_metric_volcano <- input$sig_metric
      sig_threshold_volcano <- if (sig_metric_volcano == "FDR") input$degFDR else input$degPValue
      logfc_lines_volcano <- c(-input$degLogFC, input$degLogFC)

      y_axis_label <- if (sig_metric_volcano == "FDR") "-Log10 FDR" else "-Log10 P-value"
      sig_line_y_intercept <- -log10(sig_threshold_volcano)

      p <- ggplot(plot_data, aes(x = logFC, y = negLog10Sig, key = Geneid, text = tooltip_text))

      if (any(plot_data$Highlight == "Highlight")) {
        p <- p + geom_point(color = "grey", alpha = 0.6) +
          geom_point(data = . %>% filter(Highlight == "Highlight"), color = "black", size = 3, shape = 18)
      } else {
        p <- p + geom_point(aes(color = SignificanceCategory), alpha = 0.6) +
          scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "grey"))
      }

      p <- p + labs(x = "Log2 Fold Change", y = y_axis_label, title = "Volcano Plot") +
        theme_bw() +
        geom_hline(yintercept = sig_line_y_intercept, linetype = "dashed", color = "darkgrey") +
        geom_vline(xintercept = logfc_lines_volcano, linetype = "dashed", color = "darkgrey")

      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center"))
    })

    heatmap_obj <- reactive({
      req(rv$deg_results)
      if (is.null(rv$deg_results$normalized_counts)) {
        return(NULL)
      }

      res <- rv$deg_results
      top_n <- input$heatmapTopN
      res_table <- as.data.frame(res$top_tags$table)
      sort_col <- input$sig_metric

      if (!sort_col %in% colnames(res_table)) {
        sort_col <- "PValue"
      }

      res_table_sorted <- res_table %>% arrange(!!sym(sort_col))

      if (nrow(res_table_sorted) < 2) {
        return(NULL)
      }

      top_genes <- head(res_table_sorted$Geneid, top_n)
      common_genes <- intersect(top_genes, rownames(res$normalized_counts))
      if (length(common_genes) < 2) {
        return(NULL)
      }

      mat <- res$normalized_counts[common_genes, , drop = FALSE]

      display_ids <- convert_entrez_ids_for_display(rownames(mat), input$deg_id_display_type, rv$selected_species)
      if (any(duplicated(display_ids))) {
        display_ids <- make.unique(as.character(display_ids))
      }
      rownames(mat) <- display_ids

      sample_names <- colnames(mat)
      active_metadata <- rv$sample_metadata %>% filter(current_name %in% sample_names, active == TRUE)

      if (nrow(active_metadata) == 0) {
        return(NULL)
      }

      active_metadata <- active_metadata %>% arrange(group, current_name)
      sorted_samples <- active_metadata$current_name

      common_samples <- intersect(sorted_samples, colnames(mat))
      mat <- mat[, common_samples, drop = FALSE]

      if (nrow(mat) > 0) {
        var_vec <- apply(mat, 1, var)
        mat <- mat[var_vec > 0, , drop = FALSE]
      }

      if (nrow(mat) < 2) {
        return(NULL)
      }

      annotation_df <- data.frame(Group = active_metadata$group)
      rownames(annotation_df) <- active_metadata$current_name
      annotation_df <- annotation_df[colnames(mat), , drop = FALSE]

      ph <- pheatmap::pheatmap(mat,
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        annotation_col = annotation_df,
        show_rownames = TRUE,
        show_colnames = TRUE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        border_color = NA,
        fontsize = input$heatmapFontSize,
        fontsize_row = ifelse(nrow(mat) > 50, input$heatmapFontSize * 0.6, input$heatmapFontSize * 0.8),
        fontsize_col = input$heatmapFontSize,
        main = paste("Top", nrow(mat), "Differentially Expressed Genes"),
        angle_col = 45,
        silent = TRUE
      )
      return(ph)
    })

    output$degHeatmap <- renderPlot({
      ph <- heatmap_obj()
      if (is.null(ph)) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "Not enough data or variable genes to display.", cex = 1.5)
        return(NULL)
      }
      grid::grid.draw(ph$gtable)
    })

    output$downloadHeatmapPDF <- downloadHandler(
      filename = function() {
        paste0("DEG_Heatmap_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ph <- heatmap_obj()
        if (!is.null(ph)) {
          pdf(file, width = 10, height = max(8, input$heatmapTopN * 0.15 + 3))
          grid::grid.draw(ph$gtable)
          dev.off()
        }
      }
    )


    kmeans_data_shared <- reactive({
      req(rv$deg_results)
      res <- rv$deg_results
      # 外部DE結果アップロード時はカウント行列(normalized_counts)が無いため
      # K-meansクラスタリングは不可。NULLを返してプレースホルダ表示にする。
      if (is.null(res$normalized_counts)) {
        return(NULL)
      }
      sig_genes <- res$significant_genes
      if (length(sig_genes) < 2) {
        return(NULL)
      }
      k <- input$kmeans_k
      if (length(sig_genes) < k) {
        return(NULL)
      }
      common_genes <- intersect(sig_genes, rownames(res$normalized_counts))
      mat <- res$normalized_counts[common_genes, , drop = FALSE]
      sample_names <- colnames(mat)
      active_metadata <- rv$sample_metadata %>% filter(current_name %in% sample_names, active == TRUE)
      if (nrow(active_metadata) == 0) {
        return(NULL)
      }
      mat <- mat[, active_metadata$current_name, drop = FALSE]
      mat_scaled <- t(scale(t(mat)))
      mat_scaled <- na.omit(mat_scaled)
      if (nrow(mat_scaled) < k) {
        return(NULL)
      }

      list(mat_scaled = mat_scaled, active_metadata = active_metadata, k = k)
    })

    output$kmeansPlot <- renderPlot({
      km_data <- kmeans_data_shared()
      if (is.null(km_data)) {
        plot.new()
        text(0.5, 0.5, "有効な遺伝子数が不足しています。", cex = 1.2)
        return(NULL)
      }
      mat_scaled <- km_data$mat_scaled
      active_metadata <- km_data$active_metadata
      k <- km_data$k

      if (isTRUE(input$show_elbow)) {
        max_k <- min(15, nrow(mat_scaled) - 1)
        if (max_k < 2) {
          plot.new()
          text(0.5, 0.5, "遺伝子数が少なすぎるためエルボー法は実行できません。", cex = 1.2)
          return(NULL)
        }
        set.seed(42)
        wss <- sapply(1:max_k, function(k_test) {
          kmeans(mat_scaled, centers = k_test, nstart = 10)$tot.withinss
        })
        elbow_df <- data.frame(k = 1:max_k, WSS = wss)
        p <- ggplot2::ggplot(elbow_df, ggplot2::aes(x = k, y = WSS)) +
          ggplot2::geom_line(color = "steelblue", linewidth = 1) +
          ggplot2::geom_point(color = "steelblue", size = 3) +
          ggplot2::scale_x_continuous(breaks = 1:max_k) +
          ggplot2::labs(title = "Elbow Method for optimal k", x = "Number of Clusters (k)", y = "Total Within Sum of Squares (WSS)") +
          ggplot2::theme_bw(base_size = 14)
        print(p)
        return(NULL)
      }

      set.seed(42)
      km <- kmeans(mat_scaled, centers = k, nstart = 25)

      df_long <- as.data.frame(mat_scaled) %>%
        tibble::rownames_to_column("Geneid") %>%
        dplyr::mutate(Cluster = paste("Cluster", km$cluster[Geneid])) %>%
        tidyr::pivot_longer(cols = -c(Geneid, Cluster), names_to = "current_name", values_to = "Z_score") %>%
        dplyr::left_join(active_metadata %>% dplyr::select(current_name, group), by = "current_name")

      df_summary <- df_long %>%
        dplyr::group_by(Cluster, group) %>%
        dplyr::summarise(
          mean_Z = mean(Z_score, na.rm = TRUE),
          sd_Z = sd(Z_score, na.rm = TRUE),
          .groups = "drop"
        )

      cluster_counts <- df_long %>%
        dplyr::group_by(Cluster) %>%
        dplyr::summarise(n = length(unique(Geneid)))
      df_summary <- df_summary %>%
        dplyr::left_join(cluster_counts, by = "Cluster") %>%
        dplyr::mutate(Cluster_Label = paste0(Cluster, " (n=", n, ")"))

      df_long <- df_long %>%
        dplyr::left_join(cluster_counts, by = "Cluster") %>%
        dplyr::mutate(Cluster_Label = paste0(Cluster, " (n=", n, ")"))

      unique_groups <- unique(active_metadata$group)
      df_long$group <- factor(df_long$group, levels = unique_groups)
      df_summary$group <- factor(df_summary$group, levels = unique_groups)

      p <- ggplot2::ggplot() +
        ggplot2::geom_line(data = df_long, ggplot2::aes(x = group, y = Z_score, group = Geneid), color = "grey80", alpha = 0.3) +
        ggplot2::geom_line(data = df_summary, ggplot2::aes(x = group, y = mean_Z, group = 1), color = "firebrick3", linewidth = 1.2) +
        ggplot2::geom_point(data = df_summary, ggplot2::aes(x = group, y = mean_Z), color = "firebrick3", size = 2) +
        ggplot2::facet_wrap(~Cluster_Label, scales = "free_y") +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::labs(
          title = "K-means Clustering of Significant Genes (Expression Trends)",
          x = "Group",
          y = "Expression (Z-score)"
        ) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          strip.background = ggplot2::element_rect(fill = "#f0f0f0")
        )
      print(p)
    })

    observeEvent(input$send_to_go, {
      km_data <- kmeans_data_shared()
      req(km_data)
      set.seed(42)
      km <- kmeans(km_data$mat_scaled, centers = km_data$k, nstart = 25)

      cluster_list <- split(rownames(km_data$mat_scaled), paste("Cluster", km$cluster))
      res <- isolate(rv$deg_results)
      if (!is.null(res)) {
        res$kmeans_clusters <- cluster_list
        rv$deg_results <- res
        showNotification("現在のクラスター一覧をGO解析/KEGGモジュールに送信しました！", type = "message", duration = 5)
      }
    })

    output$downloadKmeansCsv <- downloadHandler(
      filename = function() {
        paste0("Kmeans_Clusters_k", input$kmeans_k, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        km_data <- kmeans_data_shared()
        req(km_data)
        set.seed(42)
        km <- kmeans(km_data$mat_scaled, centers = km_data$k, nstart = 25)
        df_export <- data.frame(Geneid = rownames(km_data$mat_scaled), Cluster = paste("Cluster", km$cluster))
        display_ids <- convert_entrez_ids_for_display(df_export$Geneid, input$deg_id_display_type, rv$selected_species)
        df_export$Display_ID <- display_ids

        # メタデータをコメント行として先頭に追加
        res <- rv$deg_results
        comment_lines <- c(
          paste0("# Analysis Date: ", Sys.time()),
          paste0("# K-means k: ", km_data$k),
          paste0("# DEG Method: ", res$analysis_method %||% "N/A"),
          paste0("# Comparison: ", paste(res$comparison, collapse = " vs ")),
          paste0("# Species: ", rv$selected_species %||% "N/A"),
          paste0("# Gene ID Display Type: ", input$deg_id_display_type %||% "N/A"),
          paste0("# Total Genes in Clustering: ", nrow(km_data$mat_scaled)),
          ""
        )
        writeLines(comment_lines, con = file, useBytes = TRUE)
        write.table(df_export, file, row.names = FALSE, quote = TRUE,
                     sep = ",", append = TRUE, fileEncoding = "UTF-8")
      }
    )

    display_gene_list_converted <- function(entrez_gene_ids) {
      req(input$deg_id_display_type, rv$selected_species)
      if (length(entrez_gene_ids) == 0) {
        return("該当する遺伝子はありません。")
      }
      convert_entrez_ids_for_display(entrez_gene_ids, input$deg_id_display_type, rv$selected_species) %>%
        paste(collapse = ",
")
    }
    output$upGenesList <- renderPrint({
      req(rv$deg_results)
      cat(display_gene_list_converted(rv$deg_results$significant_up_genes))
    })
    output$downGenesList <- renderPrint({
      req(rv$deg_results)
      cat(display_gene_list_converted(rv$deg_results$significant_down_genes))
    })

    download_table_reactive <- reactive({
      req(rv$deg_results)
      res_table <- as.data.frame(rv$deg_results$top_tags$table)
      # FDR/PValueでソート
      sort_col <- if ("FDR" %in% colnames(res_table)) "FDR" else "PValue"
      res_table <- res_table %>% arrange(!!sym(sort_col))

      display_table <- data.table::copy(res_table)
      display_table[[1]] <- convert_entrez_ids_for_display(display_table$Geneid, input$deg_id_display_type, rv$selected_species)
      display_col_name <- names(c("Gene Symbol" = "SYMBOL", "Entrez ID" = "ENTREZID", "Gene Name" = "GENENAME"))[c("SYMBOL", "ENTREZID", "GENENAME") == input$deg_id_display_type]
      colnames(display_table)[1] <- if (length(display_col_name) > 0) display_col_name else input$deg_id_display_type

      # 手法に応じて列を選択
      cols_to_select <- if (rv$deg_results$analysis_method == "edgeR") {
        c(colnames(display_table)[1], "logFC", "logCPM", "PValue", "FDR", "LR", "F")
      } else {
        c(colnames(display_table)[1], "logFC", "baseMean", "PValue", "FDR", "stat")
      }
      display_table %>% select(any_of(cols_to_select))
    })

    # DEG解析パラメータのメタデータをdata.frameとして生成
    deg_metadata_df <- reactive({
      req(rv$deg_results)
      res <- rv$deg_results
      is_lrt <- !is.null(res$analysis_type) && res$analysis_type == "lrt"
      is_int <- !is.null(res$analysis_type) && res$analysis_type == "interaction"

      comparison_str <- if (is_lrt) {
        paste("LRT:", paste(res$comparison, collapse = " / "))
      } else if (is_int) {
        paste("Interaction:", res$comparison[1])
      } else {
        paste0(res$comparison[1], " vs ", res$comparison[2])
      }

      # サンプル-グループ対応を構築
      active_meta <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      group_cols <- setdiff(colnames(active_meta), c("original_name", "current_name", "active"))
      sample_group_lines <- paste0(active_meta$current_name, " -> ", apply(active_meta[, group_cols, drop = FALSE], 1, paste, collapse = " / "))

      params <- data.frame(
        Parameter = c(
          "Analysis Date",
          "Analysis Method",
          "Analysis Type",
          "Comparison",
          "Significance Metric",
          "Significance Threshold",
          "LogFC Threshold",
          "Species",
          "Gene ID Display Type",
          "Total Genes Tested",
          "Significant Up Genes",
          "Significant Down Genes",
          paste0("Sample_Group_", seq_along(sample_group_lines))
        ),
        Value = c(
          as.character(Sys.time()),
          res$analysis_method,
          res$analysis_type,
          comparison_str,
          res$metric_at_run %||% "FDR",
          as.character(res$sig_threshold_at_run %||% "N/A"),
          as.character(res$logfc_threshold_at_run %||% "N/A"),
          rv$selected_species %||% "N/A",
          input$deg_id_display_type %||% "N/A",
          as.character(nrow(res$top_tags$table)),
          as.character(length(res$significant_up_genes)),
          as.character(length(res$significant_down_genes)),
          sample_group_lines
        ),
        stringsAsFactors = FALSE
      )
      params
    })

    common_filename_parts_deg <- reactive({
      comp_name <- "DEG_Results"
      if (!is.null(rv$deg_results)) {
        is_lrt <- !is.null(rv$deg_results$analysis_type) && rv$deg_results$analysis_type == "lrt"
        comp_name <- if (is_lrt) paste0("DEG_LRT_", gsub(" / ", "_", rv$deg_results$comparison)) else paste0("DEG_results_", rv$deg_results$comparison[1], "_vs_", rv$deg_results$comparison[2])
      }
      method_part <- rv$deg_results$analysis_method %||% "method"
      id_type_part <- input$deg_id_display_type %||% "ID"
      paste0(comp_name, "_", method_part, "_", id_type_part, "_all-genes_", Sys.Date())
    })
    output$downloadExcelResults <- downloadHandler(
      filename = function() paste0(common_filename_parts_deg(), ".xlsx"),
      content = function(file) {
        sheets <- list(
          "DEG_Results" = download_table_reactive(),
          "Analysis_Parameters" = deg_metadata_df()
        )
        writexl::write_xlsx(sheets, file)
      }
    )
    output$downloadCsvResults <- downloadHandler(
      filename = function() paste0(common_filename_parts_deg(), ".csv"),
      content = function(file) {
        meta <- deg_metadata_df()
        comment_lines <- paste0("# ", meta$Parameter, ": ", meta$Value)
        writeLines(c(comment_lines, ""), con = file, useBytes = TRUE)
        write.table(download_table_reactive(), file, row.names = FALSE, quote = TRUE,
                     sep = ",", append = TRUE, fileEncoding = "UTF-8")
      }
    )

    `%||%` <- function(a, b) if (!is.null(a)) a else b
  })
}
