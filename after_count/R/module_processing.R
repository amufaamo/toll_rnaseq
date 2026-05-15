# R/module_processing.R
# (内部処理をEntrez IDに統一し、表示時のみ変換する修正版)
# (★★★ 複数遺伝子・複数グループのプロットに対応する改修版 ★★★)

library(shiny)
library(edgeR) # DGEList, calcNormFactors, cpm
library(dplyr) # inner_join, group_by, summarise
library(matrixStats) # rowVars
library(DT) # DTOutput, renderDT, datatable
library(shinycssloaders) # withSpinner
library(data.table) # data.table::copy
library(ggplot2) # ggplot, geom_bar, etc.
library(tidyr) # pivot_longer
library(AnnotationDbi) # mapIds, columns, keytypes
library(pheatmap) # pheatmap

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
  "Myxococcus_xanthus_DK1622" = "org.Mxanthus.db"
)

# --- Helper function: Detect Gene ID type (他のモジュールと共通) ---
detect_gene_id_type <- function(ids) {
  ids_clean <- na.omit(ids[ids != "" & !is.na(ids)])
  if (length(ids_clean) == 0) {
    return("UNKNOWN")
  }
  n_total <- length(ids_clean)
  n_sample <- min(n_total, 1000)
  ids_sample <- sample(ids_clean, n_sample)

  ensembl_pattern <- "^ENS[A-Z0-9]*[FPTG]\\d{10,}(\\.\\d+)?$"
  ensembl_match_rate <- mean(grepl(ensembl_pattern, ids_sample, ignore.case = TRUE))
  if (ensembl_match_rate > 0.8) {
    return("ENSEMBL")
  }

  refseq_pattern <- "^[NX][CMPWTZ]_[0-9]+(\\.\\d+)?$"
  refseq_match_rate <- mean(grepl(refseq_pattern, ids_sample, ignore.case = TRUE))
  if (refseq_match_rate > 0.8) {
    return("REFSEQ")
  }

  entrez_pattern <- "^[0-9]+$"
  entrez_match_rate <- mean(grepl(entrez_pattern, ids_sample))
  if (entrez_match_rate > 0.9) {
    return("ENTREZID")
  }

  symbol_pattern <- "^([A-Za-z][A-Za-z0-9-]*[A-Za-z0-9]|[A-Za-z])$"
  symbol_match_rate <- mean(grepl(symbol_pattern, ids_sample))
  if (symbol_match_rate > 0.7 && ensembl_match_rate < 0.2 && refseq_match_rate < 0.2 && entrez_match_rate < 0.2) {
    return("SYMBOL")
  }

  return("UNKNOWN")
}

processingUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4(icon("cogs"), "表示形式を選択"),
      tags$div(class = "alert alert-info", 
               icon("info-circle"), " 「フィルタリング」タブで設定された基準をクリアしたデータに対して正規化・変換を行います。内部処理はすべてEntrezIDで行われます。"),
      radioButtons(ns("normalization_method"),
        label = "手法:",
        choices = c(
          "Raw Counts (EntrezID)" = "raw",
          "CPM (EntrezID)" = "cpm",
          "TPM (EntrezID)" = "tpm",
          "FPKM (EntrezID)" = "fpkm",
          "logCPM (Scaled, Z-score, EntrezID)" = "logcpm_scaled"
        ),
        selected = "raw"
      ),
      hr(),
      selectInput(ns("id_display_type"), "表示する遺伝子IDタイプ:",
        choices = c(
          "Gene Symbol" = "SYMBOL",
          "Entrez ID (内部ID)" = "ENTREZID",
          "Gene Name" = "GENENAME"
        ),
        selected = "SYMBOL"
      ),
      hr(),
      tags$div(class = "well well-sm", style = "background-color: #fcfcfc;",
        HTML("<ul style='padding-left: 20px; font-size: 0.9em;'>
          <li><b>Raw Counts:</b> フィルタリング後の整数カウント</li>
          <li><b>CPM:</b> Counts Per Million (ライブラリサイズ補正)</li>
          <li><b>TPM/FPKM:</b> 遺伝子長による補正 (※遺伝子長データが必要)</li>
          <li><b>logCPM (Scaled):</b> 対数変換後にZ-score化 (ヒートマップ等に最適)</li>
        </ul>")
      ),
      hr(),
      h4("バッチ補正 (任意)"),
      checkboxInput(ns("use_batch"), "バッチ補正を適用して表示 (Batch Correction)", value = FALSE),
      conditionalPanel(
        condition = paste0("input['", ns("use_batch"), "'] == true"),
        selectInput(ns("batch_col"), "バッチ項の列 (Batch column):", choices = NULL)
      ),
      hr(),
      h4("全データダウンロード"),
      p("表示されているテーブル全体をダウンロードします。"),
      downloadButton(ns("downloadXLSX"), "Excel形式でダウンロード (.xlsx)", icon = icon("file-excel")),
      br(), br(),
      downloadButton(ns("downloadTSV"), "TSV形式でダウンロード (.tsv)", icon = icon("file-csv"))
    ),
    mainPanel(
      width = 8,
      h4(textOutput(ns("processed_table_title"))),
      withSpinner(DTOutput(ns("processedDataTable")), type = 6),
      hr(),
      h4(icon("chart-bar"), "個別遺伝子の発現量プロット"),
      helpText(icon("question-circle"), " 下のテーブルに表示されている遺伝子名を選択/入力して、発現量を可視化します。(複数選択可)"),
      fluidRow(
        column(6, selectizeInput(ns("genes_to_plot"), "遺伝子名を選択:", choices = NULL, multiple = TRUE, options = list(placeholder = "遺伝子名を入力または選択..."))),
        column(6, selectInput(ns("plot_type"), "プロットタイプ:",
          choices = c("サンプルごと" = "sample", "グループごと" = "group", "ヒートマップ" = "heatmap")
        ))
      ),
      conditionalPanel(
        condition = paste0("input['", ns("plot_type"), "'] == 'heatmap'"),
        tags$div(
          style = "background-color: #f9f9f9; padding: 15px; border-radius: 5px; margin-bottom: 15px;",
          h5("ヒートマップ設定"),
          fluidRow(
            column(4, checkboxInput(ns("hm1_cluster_rows"), "行(遺伝子)のクラスタリング", TRUE)),
            column(4, checkboxInput(ns("hm1_cluster_cols"), "列(サンプル)のクラスタリング", FALSE)),
            column(4, checkboxInput(ns("hm1_display_numbers"), "セル内に数値を表示", FALSE))
          ),
          fluidRow(
            column(6, numericInput(ns("hm1_fontsize_row"), "行(遺伝子)フォントサイズ:", value = 10, min = 1, max = 50)),
            column(6, numericInput(ns("hm1_fontsize_col"), "列(サンプル)フォントサイズ:", value = 10, min = 1, max = 50))
          ),
          fluidRow(
            column(6, selectInput(ns("hm1_scale"), "値のスケーリング:", choices = c("Z-score (行ごとに標準化)" = "row", "そのままの値" = "none"))),
            column(6, uiOutput(ns("hm1_sample_order_ui")))
          ),
          checkboxInput(ns("hm1_show_group"), "サンプル上のグループ名(カラーバー)を表示", value = TRUE)
        )
      ),
      withSpinner(plotOutput(ns("expression_barplot"), width = "800px", height = "600px"), type = 6),
      downloadButton(ns("download_expression_plot"), "プロットをダウンロード (.png)", icon = icon("download")),
      hr(),
      h4(icon("fire"), "カスタム遺伝子ヒートマップ"),
      helpText(icon("question-circle"), " カンマ区切りで複数の遺伝子を入力してヒートマップを描画します。（重複や存在しない遺伝子は自動的に除かれます）"),
      textInput(ns("custom_heatmap_genes"), "遺伝子 (カンマ区切り):", width = "100%", placeholder = "例: TP53, BRCA1, MYC"),
      tags$div(
        style = "background-color: #f9f9f9; padding: 15px; border-radius: 5px; margin-bottom: 15px;",
        h5("ヒートマップ設定"),
        fluidRow(
          column(4, checkboxInput(ns("hm2_cluster_rows"), "行(遺伝子)のクラスタリング", TRUE)),
          column(4, checkboxInput(ns("hm2_cluster_cols"), "列(サンプル)のクラスタリング", FALSE)),
          column(4, checkboxInput(ns("hm2_display_numbers"), "セル内に数値を表示", FALSE))
        ),
        fluidRow(
          column(6, numericInput(ns("hm2_fontsize_row"), "行(遺伝子)フォントサイズ:", value = 10, min = 1, max = 50)),
          column(6, numericInput(ns("hm2_fontsize_col"), "列(サンプル)フォントサイズ:", value = 10, min = 1, max = 50))
        ),
        fluidRow(
          column(6, selectInput(ns("hm2_scale"), "値のスケーリング:", choices = c("Z-score (行ごとに標準化)" = "row", "そのままの値" = "none"))),
          column(6, uiOutput(ns("hm2_sample_order_ui")))
        ),
        checkboxInput(ns("hm2_show_group"), "サンプル上のグループ名(カラーバー)を表示", value = TRUE)
      ),
      actionButton(ns("draw_custom_heatmap"), "ヒートマップを生成", icon = icon("paint-brush")),
      br(), br(),
      withSpinner(plotOutput(ns("custom_heatmap_plot"), width = "800px", height = "800px"), type = 6),
      downloadButton(ns("download_custom_heatmap"), "ヒートマップをPDFで保存 (.pdf)", icon = icon("file-pdf"))
    )
  )
}

processingServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    detected_input_id_type <- reactiveVal("ENTREZID")

    observeEvent(rv$merged_data,
      {
        req(rv$merged_data)
        if (nrow(rv$merged_data) > 0 && "Geneid" %in% colnames(rv$merged_data)) {
          detected_input_id_type("ENTREZID")
          message("[Processing Module] Assumed internal GeneID type: ENTREZID")

          display_choices_map <- c(
            "Gene Symbol" = "SYMBOL",
            "Entrez ID (内部ID)" = "ENTREZID",
            "Gene Name" = "GENENAME"
          )

          updateSelectInput(session, "id_display_type",
            choices = display_choices_map,
            selected = "SYMBOL"
          )
        } else {
          detected_input_id_type("UNKNOWN")
          updateSelectInput(session, "id_display_type",
            choices = c("Gene Symbol" = "SYMBOL", "Entrez ID" = "ENTREZID"),
            selected = "SYMBOL"
          )
        }
        
        # バッチ選択用の列を更新
        if (!is.null(rv$sample_metadata)) {
          factor_cols <- setdiff(colnames(rv$sample_metadata), c("id", "current_name", "active", "time"))
          updateSelectInput(session, "batch_col", choices = factor_cols, selected = if (length(factor_cols) > 1) factor_cols[2] else factor_cols[1])
        }
      },
      ignoreNULL = FALSE,
      ignoreInit = TRUE
    )

    processed_data <- reactive({
      req(rv$merged_data, rv$sample_metadata)

      active_samples_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      shiny::validate(shiny::need(nrow(active_samples_metadata) > 0, "解析対象のサンプルが選択されていません。データ入力タブで少なくとも1つのサンプルをチェックしてください。"))
      active_sample_names <- active_samples_metadata$current_name

      counts_df <- rv$merged_data
      counts_to_process_full <- counts_df[, c("Geneid", intersect(colnames(counts_df), active_sample_names)), drop = FALSE]

      filter_status_label <- "フィルタリング未適用"
      n_genes_before_filter <- nrow(counts_to_process_full)
      counts_to_process <- counts_to_process_full

      if (!is.null(rv$filtered_keep)) {
        shiny::validate(
          shiny::need(
            length(rv$filtered_keep) == n_genes_before_filter,
            paste(
              "フィルタリング情報の不整合:",
              "rv$filtered_keepの長さ (", length(rv$filtered_keep), ") が",
              "現在のデータ行数 (", n_genes_before_filter, ") と一致しません。",
              "フィルタリングタブを再確認するか、データを再アップロードしてください。"
            )
          )
        )
        counts_to_process <- counts_to_process_full[rv$filtered_keep, , drop = FALSE]
        filter_status_label <- "フィルタリング適用済み"
        shiny::validate(shiny::need(nrow(counts_to_process) > 0, "フィルタリングの結果、遺伝子が残りませんでした。"))
      }
      n_genes_after_processing <- nrow(counts_to_process)
      method <- input$normalization_method

      if (method == "raw") {
        sample_cols_raw <- setdiff(colnames(counts_to_process), "Geneid")
        if (any(!sapply(counts_to_process[, sample_cols_raw, drop = FALSE], is.numeric))) {
          stop("Raw countsが選択されましたが、一部のサンプル列が数値型ではありません。")
        }
        raw_df <- data.frame(
          Geneid = counts_to_process$Geneid,
          counts_to_process[, sample_cols_raw, drop = FALSE],
          check.names = FALSE
        )
        return(list(processed_table = raw_df, filter_status = filter_status_label, n_genes = n_genes_after_processing))
      } else if (method == "logcpm_scaled") {
        count_matrix_logcpm <- as.matrix(counts_to_process[, -which(colnames(counts_to_process) == "Geneid"), drop = FALSE])
        rownames(count_matrix_logcpm) <- counts_to_process$Geneid

        if (!is.numeric(count_matrix_logcpm)) {
          stop("logCPM用のカウント行列がGeneid列除去後に数値型ではありません。")
        }

        samples_metadata_ordered_logcpm <- active_samples_metadata[match(colnames(count_matrix_logcpm), active_samples_metadata$current_name), ]
        shiny::validate(shiny::need(nrow(samples_metadata_ordered_logcpm) == ncol(count_matrix_logcpm), "サンプルメタデータの不整合。"))

        y <- edgeR::DGEList(counts = count_matrix_logcpm, group = factor(samples_metadata_ordered_logcpm$group))
        y <- edgeR::calcNormFactors(y)
        shiny::validate(shiny::need(nrow(y) >= 1, "logCPM計算対象の遺伝子がありません。"))

        logcpm_matrix_calc <- tryCatch(edgeR::cpm(y, log = TRUE, prior.count = 2), error = function(e) shiny::validate(paste("logCPM計算エラー:", e$message)))
        req(logcpm_matrix_calc)

        # バッチ補正の適用 (Scaled logCPMの場合のみ反映)
        if (input$use_batch && !is.null(input$batch_col) && input$batch_col %in% colnames(samples_metadata_ordered_logcpm)) {
            batch <- factor(samples_metadata_ordered_logcpm[[input$batch_col]])
            if (nlevels(batch) >= 2) {
                logcpm_matrix_calc <- limma::removeBatchEffect(logcpm_matrix_calc, batch = batch)
                message("[Processing Module] Batch correction applied to logCPM matrix.")
            }
        }

        gene_vars <- matrixStats::rowVars(logcpm_matrix_calc)
        if (sum(gene_vars < 1e-8, na.rm = TRUE) > 0) {
          showNotification(paste("スケーリング警告:", sum(gene_vars < 1e-8, na.rm = TRUE), "個の遺伝子は分散がほぼゼロです (NaNになります)。"), type = "warning", duration = 10)
        }
        scaled_logcpm_matrix_calc <- t(scale(t(logcpm_matrix_calc), center = TRUE, scale = TRUE))
        scaled_df_logcpm <- data.frame(Geneid = rownames(scaled_logcpm_matrix_calc), round(scaled_logcpm_matrix_calc, 3), stringsAsFactors = FALSE, check.names = FALSE)
        return(list(processed_table = scaled_df_logcpm, filter_status = filter_status_label, n_genes = n_genes_after_processing))
      } else {
        req(rv$gene_lengths)
        gene_lengths_df_entrez <- rv$gene_lengths

        data_for_norm <- dplyr::inner_join(
          mutate(counts_to_process, Geneid = as.character(Geneid)),
          mutate(gene_lengths_df_entrez, Geneid = as.character(Geneid)),
          by = "Geneid"
        )

        shiny::validate(shiny::need(
          nrow(data_for_norm) > 0,
          if (n_genes_after_processing > 0) {
            "フィルタリング後の遺伝子に長さ情報が見つかりません (EntrezIDでの結合失敗)。"
          } else {
            "処理対象の遺伝子がありません。"
          }
        ))

        count_matrix_norm <- as.matrix(data_for_norm[, intersect(colnames(data_for_norm), active_sample_names), drop = FALSE])
        if (!is.numeric(count_matrix_norm)) stop("正規化用のカウント行列が数値型ではありません。")

        library_sizes_norm <- colSums(count_matrix_norm, na.rm = TRUE)
        library_sizes_norm[library_sizes_norm == 0] <- 1

        norm_matrix_result <- switch(method,
          "cpm" = sweep(count_matrix_norm, 2, library_sizes_norm, "/") * 1e6,
          "fpkm" = {
            req(data_for_norm$Length)
            shiny::validate(shiny::need(any(!is.na(data_for_norm$Length)), "FPKM計算のための有効な遺伝子の長さ情報がありません。"))
            length_kb <- data_for_norm$Length / 1000
            length_kb[is.na(length_kb)] <- Inf
            sweep(sweep(count_matrix_norm, 2, library_sizes_norm, "/") * 1e6, 1, length_kb, "/")
          },
          "tpm" = {
            req(data_for_norm$Length)
            shiny::validate(shiny::need(any(!is.na(data_for_norm$Length)), "TPM計算のための有効な遺伝子の長さ情報がありません。"))
            length_kb <- data_for_norm$Length / 1000
            length_kb[is.na(length_kb)] <- 0
            rpk_matrix <- count_matrix_norm
            non_zero_len_indices <- which(length_kb > 0)
            if (length(non_zero_len_indices) > 0) {
              rpk_matrix[non_zero_len_indices, ] <- rpk_matrix[non_zero_len_indices, ] / length_kb[non_zero_len_indices]
            }
            if (length(non_zero_len_indices) < nrow(rpk_matrix)) rpk_matrix[-non_zero_len_indices, ] <- 0
            rpk_sums <- colSums(rpk_matrix, na.rm = TRUE)
            rpk_sums[rpk_sums == 0] <- 1
            sweep(rpk_matrix, 2, rpk_sums, "/") * 1e6
          }
        )

        norm_matrix_result[is.infinite(norm_matrix_result) | is.na(norm_matrix_result)] <- 0
        shiny::validate(shiny::need(!is.null(norm_matrix_result), "正規化計算エラー。"))
        norm_df_result <- data.frame(Geneid = data_for_norm$Geneid, round(norm_matrix_result, 3), stringsAsFactors = FALSE, check.names = FALSE)
        return(list(processed_table = norm_df_result, filter_status = filter_status_label, n_genes = nrow(norm_df_result)))
      }
    })

    output$processed_table_title <- renderText({
      res <- processed_data()
      req(res)
      method_label_full <- switch(input$normalization_method,
        "raw" = "Raw Counts",
        "cpm" = "CPM",
        "tpm" = "TPM",
        "fpkm" = "FPKM",
        "logcpm_scaled" = "logCPM (Scaled)",
        ""
      )
      paste0("データ (", res$filter_status, ", ", method_label_full, ", ", res$n_genes, " 遺伝子, Active Samples Only)")
    })

    display_table_reactive <- reactive({
      res <- processed_data()
      req(res, res$processed_table)
      table_for_display <- data.table::copy(as.data.table(res$processed_table))
      req(nrow(table_for_display) > 0)
      selected_display_type_ui <- input$id_display_type
      ids_to_display <- as.character(table_for_display[[1]])
      if (selected_display_type_ui != "ENTREZID") {
        req(rv$selected_species)
        orgdb_pkg_name_display <- orgdb_species_map[[rv$selected_species]]
        if (!is.null(orgdb_pkg_name_display) && requireNamespace(orgdb_pkg_name_display, quietly = TRUE)) {
          require(orgdb_pkg_name_display, character.only = TRUE, quietly = TRUE)
          org_db_display <- get(orgdb_pkg_name_display)
          if (selected_display_type_ui %in% columns(org_db_display) && "ENTREZID" %in% keytypes(org_db_display)) {
            original_entrez_ids <- as.character(table_for_display[[1]])
            converted_map <- tryCatch(suppressMessages(mapIds(org_db_display, keys = unique(original_entrez_ids), column = selected_display_type_ui, keytype = "ENTREZID", multiVals = "first")), error = function(e) NULL)
            if (!is.null(converted_map)) {
              mapped_values <- converted_map[original_entrez_ids]
              na_indices <- is.na(mapped_values)
              if (any(na_indices)) mapped_values[na_indices] <- paste0(original_entrez_ids[na_indices], " (変換不可)")
              ids_to_display <- mapped_values
            }
          }
        }
      }
      table_for_display[[1]] <- ids_to_display
      choices_map <- c("Gene Symbol" = "SYMBOL", "Entrez ID (内部ID)" = "ENTREZID", "Gene Name" = "GENENAME")
      display_col_name <- names(choices_map)[choices_map == selected_display_type_ui]
      if (length(display_col_name) == 0) display_col_name <- selected_display_type_ui
      colnames(table_for_display)[1] <- display_col_name
      return(table_for_display)
    })

    output$processedDataTable <- renderDT({
      final_table <- display_table_reactive()
      req(final_table)
      display_col_name <- colnames(final_table)[1]
      dt_options <- list(pageLength = 15, scrollX = TRUE, searching = TRUE, columnDefs = list(list(className = "dt-right", targets = "_all")))
      datatable(final_table, rownames = FALSE, options = dt_options, filter = "top") %>%
        formatRound(columns = if (input$normalization_method != "raw" && ncol(final_table) > 1) base::setdiff(colnames(final_table), display_col_name) else NULL, digits = 3)
    })

    # ★★★ ここからプロット機能のサーバーロジック（全面改修） ★★★

    output$hm1_sample_order_ui <- renderUI({
      req(input$hm1_cluster_cols == FALSE)
      active_samples <- rv$sample_metadata[rv$sample_metadata$active == TRUE, , drop = FALSE]
      if(nrow(active_samples) > 0 && "group" %in% colnames(active_samples) && "current_name" %in% colnames(active_samples)) {
        active_samples <- active_samples[order(active_samples$group, active_samples$current_name), ]
      }
      
      selectizeInput(ns("hm1_sample_order"), "サンプルの順序 (削除・選択して並び替え):",
                     choices = active_samples$current_name,
                     selected = active_samples$current_name,
                     multiple = TRUE)
    })

    output$hm2_sample_order_ui <- renderUI({
      req(input$hm2_cluster_cols == FALSE)
      active_samples <- rv$sample_metadata[rv$sample_metadata$active == TRUE, , drop = FALSE]
      if(nrow(active_samples) > 0 && "group" %in% colnames(active_samples) && "current_name" %in% colnames(active_samples)) {
        active_samples <- active_samples[order(active_samples$group, active_samples$current_name), ]
      }
      
      selectizeInput(ns("hm2_sample_order"), "サンプルの順序 (削除・選択して並び替え):",
                     choices = active_samples$current_name,
                     selected = active_samples$current_name,
                     multiple = TRUE)
    })

    # 遺伝子選択UIをサーバーサイドselectizeで動的に更新
    # [Phase 1] wrapped via execute_command (op = "update_gene_choices")
    observeEvent(display_table_reactive(), {
      display_table <- display_table_reactive()
      if (is.data.frame(display_table) && nrow(display_table) > 0) {
        gene_id_col_name <- colnames(display_table)[1]
        available_genes <- sort(unique(na.omit(as.character(display_table[[gene_id_col_name]]))))
        # Route UI update through the command layer so it is logged like any other op.
        execute_command(
          list(
            op = "update_gene_choices",
            input_id = "genes_to_plot",
            gene_choices = as.character(available_genes),
            source = "processing/display_table_reactive"
          ),
          rv, session
        )
      }
    })

    # ★★★ グループ選択UIは不要になったため、サーバー側のUI生成ロジックを削除 ★★★

    # プロット用データを作成するReactive
    plot_data_reactive <- reactive({
      shiny::validate(shiny::need(input$genes_to_plot, "プロットする遺伝子を1つ以上選択してください。"))

      full_display_table <- display_table_reactive()
      shiny::validate(shiny::need(is.data.frame(full_display_table) && nrow(full_display_table) > 0, "表示データがありません。"))

      target_gene_names <- input$genes_to_plot
      gene_id_col_name <- colnames(full_display_table)[1]

      gene_rows <- full_display_table[full_display_table[[gene_id_col_name]] %in% target_gene_names, ]
      shiny::validate(shiny::need(nrow(gene_rows) > 0, "選択された遺伝子がデータ内に見つかりません。"))

      plot_df <- gene_rows %>%
        tidyr::pivot_longer(cols = -!!sym(gene_id_col_name), names_to = "current_name", values_to = "Expression") %>%
        rename(Gene = !!sym(gene_id_col_name))

      active_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      plot_df <- plot_df %>%
        left_join(active_metadata, by = "current_name")

      return(plot_df)
    })

    # プロットオブジェクトを作成するReactive
    expression_plot_object <- reactive({
      plot_df <- plot_data_reactive()
      req(plot_df)

      y_axis_label <- paste("Expression Level (", input$normalization_method, ")", sep = "")

      if (input$plot_type == "sample") {
        p <- ggplot(plot_df, aes(x = current_name, y = Expression, fill = group)) +
          geom_bar(stat = "identity", color = "black") +
          facet_wrap(~Gene, scales = "free_y") +
          labs(title = "Gene Expression by Sample", x = "Sample", y = y_axis_label) +
          theme_bw(base_size = 12) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top")

        return(p)
      } else { # グループごとのプロット
        active_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
        all_active_groups <- unique(active_metadata$group)

        plot_df_group <- plot_df %>%
          filter(group %in% all_active_groups)

        summary_df <- plot_df_group %>%
          group_by(Gene, group) %>%
          summarise(
            mean_expr = mean(Expression, na.rm = TRUE),
            sd_expr = sd(Expression, na.rm = TRUE),
            .groups = "drop"
          )

        # ★★★ 統計検定をt検定からANOVA（3群以上）またはt検定（2群）に自動切り替え ★★★
        p_values_df <- plot_df_group %>%
          group_by(Gene) %>%
          summarise(
            p_value_text = {
              n_groups <- length(unique(group))
              if (n_groups >= 3) {
                # ANOVA
                aov_res <- aov(Expression ~ group, data = cur_data())
                p_val <- summary(aov_res)[[1]][["Pr(>F)"]][1]
                paste("ANOVA p =", format(p_val, digits = 2, scientific = TRUE))
              } else if (n_groups == 2) {
                # t-test
                group_names <- unique(group)
                group1_vals <- Expression[group == group_names[1]]
                group2_vals <- Expression[group == group_names[2]]
                if (length(group1_vals) >= 2 && length(group2_vals) >= 2) {
                  paste("t-test p =", format(t.test(group1_vals, group2_vals)$p.value, digits = 2, scientific = TRUE))
                } else {
                  "N/A"
                }
              } else {
                "N/A"
              }
            },
            .groups = "drop"
          ) %>%
          filter(p_value_text != "N/A")

        p <- ggplot(summary_df, aes(x = group, y = mean_expr, fill = group)) +
          geom_bar(stat = "identity", color = "black") +
          geom_errorbar(aes(ymin = pmax(0, mean_expr - sd_expr), ymax = mean_expr + sd_expr), width = 0.2) +
          facet_wrap(~Gene, scales = "free_y") +
          labs(title = "Gene Expression by Group", x = "Group", y = y_axis_label) +
          theme_bw(base_size = 14) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

        if (nrow(p_values_df) > 0) {
          p <- p + geom_text(
            data = p_values_df,
            aes(x = Inf, y = Inf, label = p_value_text),
            hjust = 1.05, vjust = 1.5, inherit.aes = FALSE, size = 4
          )
        }

        return(p)
      }
    })

    output$expression_barplot <- renderPlot({
      if (input$plot_type == "heatmap") {
        target_genes <- input$genes_to_plot
        shiny::validate(shiny::need(length(target_genes) >= 2, "Heatmap requires at least 2 genes."))

        df <- display_table_reactive()
        shiny::validate(shiny::need(nrow(df) > 0, "No data available."))

        gene_col <- colnames(df)[1]
        df_sub <- df[df[[gene_col]] %in% target_genes, ]

        shiny::validate(shiny::need(nrow(df_sub) > 0, "Selected genes not found in data."))

        mat <- as.matrix(df_sub[, -1, drop = FALSE])
        mat <- apply(mat, 2, as.numeric)
        rownames(mat) <- make.unique(as.character(df_sub[[1]]))

        # Active samples only & Sort by group
        meta <- rv$sample_metadata[rv$sample_metadata$active, ]
        meta <- meta[order(meta$group, meta$current_name), ]
        valid_samples <- intersect(meta$current_name, colnames(mat))

        if (!is.null(input$hm1_cluster_cols) && !input$hm1_cluster_cols) {
          if (!is.null(input$hm1_sample_order) && length(input$hm1_sample_order) > 0) {
            valid_samples <- intersect(input$hm1_sample_order, valid_samples)
          }
        }

        mat <- mat[, valid_samples, drop = FALSE]
        shiny::validate(shiny::need(ncol(mat) > 0, "No active samples."))

        # Annotation
        anno <- data.frame(Group = meta$group)
        rownames(anno) <- meta$current_name
        anno <- anno[valid_samples, , drop = FALSE]

        # Config
        scale_val <- if (!is.null(input$hm1_scale)) input$hm1_scale else "row"
        clus_rows <- if (!is.null(input$hm1_cluster_rows)) input$hm1_cluster_rows else TRUE
        clus_cols <- if (!is.null(input$hm1_cluster_cols)) input$hm1_cluster_cols else FALSE
        disp_num  <- if (!is.null(input$hm1_display_numbers)) input$hm1_display_numbers else FALSE
        f_row <- if (!is.null(input$hm1_fontsize_row) && !is.na(input$hm1_fontsize_row)) as.numeric(input$hm1_fontsize_row) else 10
        f_col <- if (!is.null(input$hm1_fontsize_col) && !is.na(input$hm1_fontsize_col)) as.numeric(input$hm1_fontsize_col) else 10
        show_grp <- if (!is.null(input$hm1_show_group)) input$hm1_show_group else TRUE
        main_title <- if (scale_val == "row") "Expression Heatmap (Z-score)" else "Expression Heatmap"

        # Plot
        ph <- pheatmap::pheatmap(mat,
          scale = scale_val, cluster_rows = clus_rows, cluster_cols = clus_cols,
          display_numbers = disp_num,
          annotation_col = if(show_grp) anno else NA, show_colnames = TRUE,
          color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          border_color = NA, main = main_title,
          fontsize_row = f_row, fontsize_col = f_col,
          silent = TRUE
        )
        return(ggplotify::as.ggplot(ph$gtable))
      } else {
        expression_plot_object()
      }
    })

    output$download_expression_plot <- downloadHandler(
      filename = function() {
        gene_names <- "selected_genes"
        if (!is.null(input$genes_to_plot) && length(input$genes_to_plot) > 0) {
          gene_names <- paste(input$genes_to_plot, collapse = "_")
        }
        paste0("expression_plot_", gene_names, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        p <- expression_plot_object()
        req(p)
        # Shiny UIで表示されているサイズと同じサイズでダウンロード
        # plotOutputで width="800px", height="600px" と指定しているため、
        # 800px / 75ppi ≈ 10.67 inches, 600px / 75ppi = 8 inches
        plot_width <- 10.67
        plot_height <- 8
        ggsave(file, plot = p, device = "png", width = plot_width, height = plot_height, dpi = 75, limitsize = FALSE)
      }
    )

    custom_heatmap_data <- eventReactive(input$draw_custom_heatmap, {
      req(input$custom_heatmap_genes)
      gene_text <- input$custom_heatmap_genes
      genes <- strsplit(gene_text, ",")[[1]]
      genes <- trimws(genes)
      genes <- genes[nzchar(genes)]
      target_genes <- unique(genes)

      shiny::validate(shiny::need(length(target_genes) >= 2, "ヒートマップを描画するには、少なくとも2つ以上の遺伝子を入力してください。"))

      df <- display_table_reactive()
      shiny::validate(shiny::need(nrow(df) > 0, "データがありません。"))

      gene_col <- colnames(df)[1]

      target_genes_lower <- tolower(target_genes)
      df_genes_lower <- tolower(as.character(df[[gene_col]]))

      df_sub <- df[df_genes_lower %in% target_genes_lower, ]

      found_genes <- as.character(df_sub[[gene_col]])
      missing_idx <- which(!(target_genes_lower %in% tolower(found_genes)))
      if (length(missing_idx) > 0) {
        missing_genes <- target_genes[missing_idx]
        showNotification(paste("以下の遺伝子はフィルタリング後データ内に見つかりませんでした:", paste(missing_genes, collapse = ", ")), type = "warning", duration = 12)
      }

      shiny::validate(shiny::need(nrow(df_sub) >= 2, "データ内に見つかった指定遺伝子が2つ未満のため、ヒートマップを描画できません。"))

      mat <- as.matrix(df_sub[, -1, drop = FALSE])
      mat <- apply(mat, 2, as.numeric)
      rownames(mat) <- make.unique(as.character(df_sub[[1]]))

      meta <- rv$sample_metadata[rv$sample_metadata$active, ]
      meta <- meta[order(meta$group, meta$current_name), ]
      valid_samples <- intersect(meta$current_name, colnames(mat))

      if (!is.null(input$hm2_cluster_cols) && !input$hm2_cluster_cols) {
        if (!is.null(input$hm2_sample_order) && length(input$hm2_sample_order) > 0) {
            valid_samples <- intersect(input$hm2_sample_order, valid_samples)
        }
      }

      mat <- mat[, valid_samples, drop = FALSE]
      shiny::validate(shiny::need(ncol(mat) > 0, "有効なサンプルがありません。"))

      anno <- data.frame(Group = meta$group)
      rownames(anno) <- meta$current_name
      anno <- anno[valid_samples, , drop = FALSE]

      list(mat = mat, anno = anno)
    })

    output$custom_heatmap_plot <- renderPlot({
      res <- custom_heatmap_data()
      scale_val <- if (!is.null(input$hm2_scale)) input$hm2_scale else "row"
      clus_rows <- if (!is.null(input$hm2_cluster_rows)) input$hm2_cluster_rows else TRUE
      clus_cols <- if (!is.null(input$hm2_cluster_cols)) input$hm2_cluster_cols else FALSE
      disp_num  <- if (!is.null(input$hm2_display_numbers)) input$hm2_display_numbers else FALSE
      f_row <- if (!is.null(input$hm2_fontsize_row) && !is.na(input$hm2_fontsize_row)) as.numeric(input$hm2_fontsize_row) else 10
      f_col <- if (!is.null(input$hm2_fontsize_col) && !is.na(input$hm2_fontsize_col)) as.numeric(input$hm2_fontsize_col) else 10
      show_grp <- if (!is.null(input$hm2_show_group)) input$hm2_show_group else TRUE
      main_title <- if (scale_val == "row") "Custom Gene Expression Heatmap (Z-score)" else "Custom Gene Expression Heatmap"

      ph <- pheatmap::pheatmap(res$mat,
        scale = scale_val, cluster_rows = clus_rows, cluster_cols = clus_cols,
        display_numbers = disp_num,
        annotation_col = if(show_grp) res$anno else NA, show_colnames = TRUE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        border_color = NA, main = main_title,
        fontsize_row = f_row, fontsize_col = f_col,
        silent = TRUE
      )
      return(ggplotify::as.ggplot(ph$gtable))
    })

    output$download_custom_heatmap <- downloadHandler(
      filename = function() {
        paste0("Custom_Heatmap_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        res <- custom_heatmap_data()
        pdf(file, width = 8, height = 8)
        scale_val <- if (!is.null(input$hm2_scale)) input$hm2_scale else "row"
        clus_rows <- if (!is.null(input$hm2_cluster_rows)) input$hm2_cluster_rows else TRUE
        clus_cols <- if (!is.null(input$hm2_cluster_cols)) input$hm2_cluster_cols else FALSE
        disp_num  <- if (!is.null(input$hm2_display_numbers)) input$hm2_display_numbers else FALSE
        main_title <- if (scale_val == "row") "Custom Gene Expression Heatmap (Z-score)" else "Custom Gene Expression Heatmap"

        pheatmap::pheatmap(res$mat,
          scale = scale_val, cluster_rows = clus_rows, cluster_cols = clus_cols,
          display_numbers = disp_num,
          annotation_col = res$anno, show_colnames = TRUE,
          color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          border_color = NA, main = main_title
        )
        dev.off()
      }
    )

    common_filename_parts <- reactive({
      res_dl <- processed_data()
      req(res_dl)
      filter_part_dl <- gsub(" ", "_", res_dl$filter_status)
      method_part_dl <- input$normalization_method
      id_type_part_dl <- input$id_display_type %||% "data"

      id_type_part_dl_safe <- gsub("[^A-Za-z0-9_.-]", "_", id_type_part_dl)

      paste0("data_", id_type_part_dl_safe, "_", filter_part_dl, "_", method_part_dl, "_", res_dl$n_genes, "genes_", Sys.Date())
    })

    # 前処理パラメータのメタデータを生成
    processing_metadata_df <- reactive({
      res_dl <- processed_data()
      req(res_dl)
      active_meta <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      group_cols <- setdiff(colnames(active_meta), c("original_name", "current_name", "active"))
      sample_group_lines <- paste0(active_meta$current_name, " -> ", apply(active_meta[, group_cols, drop = FALSE], 1, paste, collapse = " / "))

      params <- data.frame(
        Parameter = c(
          "Export Date",
          "Normalization Method",
          "Filtering Status",
          "Gene ID Display Type",
          "Total Genes",
          "Species",
          paste0("Sample_Group_", seq_along(sample_group_lines))
        ),
        Value = c(
          as.character(Sys.time()),
          input$normalization_method,
          res_dl$filter_status,
          input$id_display_type %||% "N/A",
          as.character(res_dl$n_genes),
          rv$selected_species %||% "N/A",
          sample_group_lines
        ),
        stringsAsFactors = FALSE
      )
      params
    })

    output$downloadXLSX <- downloadHandler(
      filename = function() {
        paste0(common_filename_parts(), ".xlsx")
      },
      content = function(file) {
        table_to_download <- display_table_reactive()
        if (is.null(table_to_download) || nrow(table_to_download) == 0) {
          df_error <- data.frame(Message = "ダウンロードするデータがありません。")
          writexl::write_xlsx(df_error, file)
          return()
        }
        sheets <- list(
          "Expression_Data" = as.data.frame(table_to_download),
          "Analysis_Parameters" = processing_metadata_df()
        )
        writexl::write_xlsx(sheets, file)
      }
    )

    output$downloadTSV <- downloadHandler(
      filename = function() {
        paste0(common_filename_parts(), ".tsv")
      },
      content = function(file) {
        table_to_download <- display_table_reactive()
        if (is.null(table_to_download) || nrow(table_to_download) == 0) {
          write("ダウンロードするデータがありません。", file)
          return()
        }
        # メタデータをコメント行として先頭に追加
        meta <- processing_metadata_df()
        comment_lines <- paste0("# ", meta$Parameter, ": ", meta$Value)
        writeLines(c(comment_lines, ""), con = file, useBytes = TRUE)
        data.table::fwrite(table_to_download, file, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
      }
    )

    `%||%` <- function(a, b) if (!is.null(a)) a else b
  })
}
