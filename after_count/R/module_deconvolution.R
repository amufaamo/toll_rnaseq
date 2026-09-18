# R/module_deconvolution.R

library(shiny)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(shinycssloaders)
library(pheatmap)
library(AnnotationDbi)

detect_deconv_gene_id_type <- function(ids) {
  ids_clean <- na.omit(as.character(ids[ids != "" & !is.na(ids)]))
  if (length(ids_clean) == 0) return("UNKNOWN")
  ids_sample <- sample(ids_clean, min(length(ids_clean), 1000))
  if (mean(grepl("^ENS[A-Z0-9]*[FPTG]\\d{10,}(\\.\\d+)?$", ids_sample, ignore.case = TRUE)) > 0.8) return("ENSEMBL")
  if (mean(grepl("^[NX][CMPWTZ]_[0-9]+(\\.\\d+)?$", ids_sample, ignore.case = TRUE)) > 0.8) return("REFSEQ")
  if (mean(grepl("^[0-9]+$", ids_sample)) > 0.9) return("ENTREZID")
  if (mean(grepl("^([A-Za-z][A-Za-z0-9-]*[A-Za-z0-9]|[A-Za-z])$", ids_sample)) > 0.7) return("SYMBOL")
  "UNKNOWN"
}

get_current_gene_id_type <- function(rv, gene_ids) {
  id_type <- NULL
  if (!is.null(rv$current_gene_id_type)) {
    id_type <- if (is.function(rv$current_gene_id_type)) rv$current_gene_id_type() else rv$current_gene_id_type
  }
  if (is.null(id_type) || !nzchar(id_type) || grepl("ORIGINAL|MIXED|UNKNOWN", id_type, ignore.case = TRUE)) {
    id_type <- detect_deconv_gene_id_type(gene_ids)
  }
  id_type
}

prepare_deconv_expression <- function(expr_matrix, gene_ids, selected_species, id_type) {
  rownames(expr_matrix) <- as.character(gene_ids)
  rownames(expr_matrix) <- sub("\\.\\d+$", "", rownames(expr_matrix))

  if (id_type == "SYMBOL" || selected_species == "Others_Original") {
    symbols <- rownames(expr_matrix)
  } else {
    orgdb_species_map <- c(
      "Homo_sapiens" = "org.Hs.eg.db",
      "Mus_musculus" = "org.Mm.eg.db"
    )
    orgdb_pkg_name <- orgdb_species_map[[selected_species]]
    if (is.null(orgdb_pkg_name) || !requireNamespace(orgdb_pkg_name, quietly = TRUE)) {
      stop("Gene Symbolへの変換に必要なOrgDbパッケージが見つかりません。生物種設定を確認してください。")
    }
    require(orgdb_pkg_name, character.only = TRUE, quietly = TRUE)
    orgDb <- get(orgdb_pkg_name)
    if (!id_type %in% keytypes(orgDb)) {
      detected_type <- detect_deconv_gene_id_type(rownames(expr_matrix))
      if (detected_type == "SYMBOL") {
        symbols <- rownames(expr_matrix)
      } else {
        stop(paste0("Gene ID type '", id_type, "' は ", orgdb_pkg_name, " でSYMBOLへ変換できません。"))
      }
    } else {
      symbols <- suppressMessages(mapIds(orgDb, keys = rownames(expr_matrix), column = "SYMBOL", keytype = id_type, multiVals = "first"))
    }
  }

  valid_idx <- !is.na(symbols) & nzchar(symbols)
  if (sum(valid_idx) == 0) {
    stop("Gene Symbolへ変換できた遺伝子がありません。アップロード時の生物種とGene ID形式を確認してください。")
  }

  expr_matrix <- expr_matrix[valid_idx, , drop = FALSE]
  symbols <- toupper(as.character(symbols[valid_idx]))
  rowsum(expr_matrix, group = symbols, na.rm = TRUE)
}

prepare_deconv_method <- function(method) {
  if (method == "epic") {
    if (!requireNamespace("EPIC", quietly = TRUE)) {
      stop("EPICパッケージが見つかりません。install_pkg.R を再実行してください。")
    }
    if (!"package:EPIC" %in% search()) {
      suppressPackageStartupMessages(require("EPIC", character.only = TRUE))
    }
  }

  if (method == "xcell") {
    if (!requireNamespace("xCell", quietly = TRUE)) {
      stop("xCellパッケージが見つかりません。install_pkg.R を再実行してください。")
    }
    if (!"package:xCell" %in% search()) {
      suppressPackageStartupMessages(require("xCell", character.only = TRUE))
    }
    if (!exists("xCell.data", envir = .GlobalEnv, inherits = FALSE)) {
      data("xCell.data", package = "xCell", envir = .GlobalEnv)
    }
  }

  if (method == "mcp_counter" && !requireNamespace("MCPcounter", quietly = TRUE)) {
    stop("MCPcounterパッケージが見つかりません。install_pkg.R を再実行してください。")
  }
}

format_deconv_result <- function(res) {
  res_df <- as.data.frame(res, check.names = FALSE)
  if (!"cell_type" %in% colnames(res_df)) {
    res_df <- tibble::rownames_to_column(res_df, "cell_type")
  }
  score_cols <- setdiff(colnames(res_df), "cell_type")
  if (length(score_cols) == 0) {
    stop("Deconvolution結果にサンプル列がありません。サンプルメタデータのactive設定とサンプル名を確認してください。")
  }
  res_df[, score_cols] <- lapply(res_df[, score_cols, drop = FALSE], as.numeric)
  res_df
}

deconvolutionUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4(icon("layer-group"), "Deconvolution"),
      helpText("Bulk RNA-seqデータから各細胞タイプの割合や存在量を推定します。(主としてヒトの免疫・間質細胞向け)"),
      selectInput(ns("deconv_method"), "推定手法:",
                  choices = c("quanTIseq" = "quantiseq", 
                              "EPIC" = "epic", 
                              "xCell" = "xcell", 
                              "MCP-counter" = "mcp_counter")),
      tags$div(class = "well well-sm", style = "background-color: #fcfcfc;",
        HTML("<ul style='padding-left: 20px; font-size: 0.9em;'>
          <li><b>quanTIseq:</b> 免疫細胞の絶対割合を推定 (TPM推奨)</li>
          <li><b>EPIC:</b> 免疫・間質細胞の絶対割合を推定 (TPM推奨)</li>
          <li><b>xCell:</b> 64種類の細胞のエンリッチメントスコア</li>
          <li><b>MCP-counter:</b> 免疫・間質細胞の存在量スコア</li>
        </ul>")
      ),
      actionButton(ns("run_deconv"), "Deconvolutionを実行", class = "btn-primary", icon = icon("play")),
      hr(),
      downloadButton(ns("download_results"), "結果をダウンロード (.csv)", icon = icon("download"))
    ),
    mainPanel(
      width = 8,
      h4("細胞割合 / スコア プロット"),
      withSpinner(plotOutput(ns("barplot"), width = "100%", height = "500px"), type = 6),
      hr(),
      h4("ヒートマップ (Row Z-score)"),
      withSpinner(plotOutput(ns("heatmap"), width = "100%", height = "600px"), type = 6),
      hr(),
      h4("推定マトリクス"),
      withSpinner(DTOutput(ns("result_table")), type = 6)
    )
  )
}

deconvolutionServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    deconv_results <- reactiveVal(NULL)

    observeEvent(input$run_deconv, {
      req(rv$merged_data, rv$sample_metadata)

      if (!requireNamespace("immunedeconv", quietly = TRUE)) {
        showNotification("immunedeconvパッケージがインストールされていません。install_pkg.R を実行してください。", type = "error", duration = 10)
        return()
      }

      showNotification("Deconvolutionを実行中...", type = "message", duration = 5)

      tryCatch({
        # 1. データの準備 (Filtered raw counts)
        counts_df <- rv$merged_data
        if (!is.null(rv$filtered_keep)) {
            counts_df <- counts_df[rv$filtered_keep, , drop = FALSE]
        }
        
        # 2. TPMの計算 (immunedeconvはTPMを推奨)
        count_matrix <- as.matrix(counts_df[, setdiff(colnames(counts_df), "Geneid"), drop = FALSE])
        
        if (!is.null(rv$gene_lengths) && "Length" %in% colnames(rv$gene_lengths)) {
           len_df <- rv$gene_lengths
           idx <- match(counts_df$Geneid, len_df$Geneid)
           lengths <- len_df$Length[idx]
           lengths[is.na(lengths) | lengths == 0] <- 1
           
           rpk <- count_matrix / (lengths / 1000)
           expr_matrix <- sweep(rpk, 2, colSums(rpk, na.rm=TRUE), "/") * 1e6
        } else {
           # 長さがない場合はCPM
           expr_matrix <- sweep(count_matrix, 2, colSums(count_matrix, na.rm=TRUE), "/") * 1e6
        }

        # 3. Gene ID -> Symbol変換 (immunedeconvは主にHGNC Symbolを要求)
        req(rv$selected_species)
        id_type <- get_current_gene_id_type(rv, counts_df$Geneid)
        expr_matrix <- prepare_deconv_expression(expr_matrix, counts_df$Geneid, rv$selected_species, id_type)

        # 4. Active samples only & run deconvolution
        active_meta <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
        expr_matrix <- expr_matrix[, intersect(colnames(expr_matrix), active_meta$current_name), drop = FALSE]
        if (ncol(expr_matrix) == 0) {
          stop("activeなサンプル名と発現データの列名が一致していません。サンプルメタデータを確認してください。")
        }

        if (!is.null(rv$selected_species) && rv$selected_species != "Homo_sapiens") {
          showNotification(
            "注意: EPIC/xCell/MCPcounter/quanTIseqはヒト向けツールです。非ヒトデータの場合、結果の生物学的解釈に注意してください。",
            type = "warning", duration = 10
          )
        }

        prepare_deconv_method(input$deconv_method)
        run_deconv <- function(mat) {
          prepare_deconv_method(input$deconv_method)
          immunedeconv::deconvolute(mat, input$deconv_method)
        }
        res <- tryCatch(
          run_deconv(expr_matrix),
          error = function(e) {
            gene_match_error <- grepl(
              "No match found between signature genes|0 signature genes|non-numeric argument to binary operator",
              e$message
            )
            has_lower <- any(grepl("[a-z]", rownames(expr_matrix)))
            if (!gene_match_error || !has_lower) stop(e)

            upper_symbols <- toupper(rownames(expr_matrix))
            expr_matrix_upper <- rowsum(expr_matrix, group = upper_symbols, na.rm = TRUE)
            run_deconv(expr_matrix_upper)
          }
        )
        # MCPcounterが遺伝子マッチ0件で空結果を返した場合、大文字で再試行
        if (input$deconv_method == "mcp_counter") {
          res_check <- as.data.frame(res, check.names = FALSE)
          if (!"cell_type" %in% colnames(res_check)) res_check <- tibble::rownames_to_column(res_check, "cell_type")
          if (length(setdiff(colnames(res_check), "cell_type")) == 0 && any(grepl("[a-z]", rownames(expr_matrix)))) {
            upper_symbols <- toupper(rownames(expr_matrix))
            expr_matrix_upper <- rowsum(expr_matrix, group = upper_symbols, na.rm = TRUE)
            res <- run_deconv(expr_matrix_upper)
          }
        }
        res_df <- format_deconv_result(res)
        
        deconv_results(list(data = res_df, method = input$deconv_method))
        showNotification("Deconvolutionが完了しました！", type = "message")

      }, error = function(e) {
        showNotification(paste("エラー:", e$message), type = "error", duration = 10)
      })
    })

    # --- Plot Barplot ---
    output$barplot <- renderPlot({
      res <- deconv_results()
      req(res)
      
      df <- res$data
      shiny::validate(shiny::need(ncol(df) > 1, "Deconvolution結果にプロット可能なサンプル列がありません。"))
      df_long <- df %>%
        tidyr::pivot_longer(cols = -cell_type, names_to = "Sample", values_to = "Score")
      
      meta <- rv$sample_metadata
      df_long <- df_long %>% left_join(meta, by = c("Sample" = "current_name"))
      
      y_label <- if(res$method %in% c("quantiseq", "epic")) "Estimated Cell Fraction" else "Enrichment Score"
      
      ggplot(df_long, aes(x = Sample, y = Score, fill = cell_type)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_grid(~group, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("Deconvolution by", res$method), y = y_label, x = "Sample", fill = "Cell Type")
    })

    # --- Plot Heatmap ---
    output$heatmap <- renderPlot({
      res <- deconv_results()
      req(res)
      
      mat <- as.matrix(res$data[, -1, drop = FALSE])
      shiny::validate(shiny::need(ncol(mat) > 0, "Deconvolution結果にヒートマップ用のサンプル列がありません。"))
      rownames(mat) <- res$data$cell_type
      
      # Hide rows with zero variance
      row_vars <- apply(mat, 1, var, na.rm = TRUE)
      mat <- mat[row_vars > 1e-8, , drop = FALSE]
      req(nrow(mat) > 1)
      
      meta <- rv$sample_metadata[rv$sample_metadata$active, ]
      meta <- meta[match(colnames(mat), meta$current_name), ]
      anno_col <- data.frame(Group = meta$group)
      rownames(anno_col) <- meta$current_name
      
      main_title <- paste("Cell Type Scores (", res$method, ")")
      
      pheatmap::pheatmap(mat,
                         cluster_cols = FALSE,
                         annotation_col = anno_col,
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         scale = "row",
                         main = main_title)
    })

    # --- Table ---
    output$result_table <- renderDT({
      res <- deconv_results()
      req(res)
      datatable(res$data, rownames = FALSE, style = "bootstrap5", class = "table-hover table-sm", options = list(scrollX = TRUE)) %>%
        formatRound(columns = setdiff(colnames(res$data), "cell_type"), digits = 4)
    })

    # --- Download ---
    output$download_results <- downloadHandler(
      filename = function() {
        res <- deconv_results()
        paste0("Deconvolution_", res$method, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        res <- deconv_results()
        write.csv(res$data, file, row.names = FALSE)
      }
    )
  })
}
