# R/module_gene_id_converter.R (修正版)

# --- 1. Load Required Libraries ---
library(shiny)
library(AnnotationDbi)
library(dplyr)
library(shinycssloaders) # For withSpinner loading indicators

# --- ★★★ NULL Coalescing Operator (ファイルのトップレベルで定義) ★★★ ---
# Returns 'a' if 'a' is not NULL, otherwise returns 'b'
# Provides compatibility for R versions < 4.0.0 where `??` is not available.
# If using R >= 4.0.0, you could potentially replace `%||%` with `??`
# but defining it ensures compatibility.
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

# --- 2. UI Definition ---
#' Gene ID Converter UI Function
# ... (UIコードは変更なし) ...
geneIdConverterUI <- function(id) {
  ns <- NS(id) # Namespace function
  tagList(
    h4("遺伝子ID変換 (Ensembl IDへ)"),
    helpText("DEG解析結果の有意変動遺伝子リストをEnsembl IDに変換します。"),
    fluidRow(
      # Species Selection
      column(6,
             selectInput(ns("species"), "生物種:",
                         choices = c("選択してください" = "",
                                     "Human" = "Human",
                                     "Mouse" = "Mouse"
                                     # Add other species here if needed
                                     ),
                         selected = "")
      ),
      # Original Keytype Selection (with auto-detection hint)
      column(6,
             selectInput(ns("original_keytype"), "元の遺伝子IDタイプ (自動推定):",
                         choices = c("選択してください" = "",
                                     "Gene Symbol" = "SYMBOL",
                                     "Entrez ID" = "ENTREZID",
                                     "RefSeq ID" = "REFSEQ",
                                     "Ensembl ID (変換元)" = "ENSEMBL"
                                     # Add other relevant keytypes if needed
                                     ),
                         selected = "")
      )
    ),
    # Action Button to Trigger Conversion
    actionButton(ns("run_conversion"), "Ensembl IDに変換", icon = icon("sync-alt")),
    hr(), # Horizontal rule
    # Output Area for Converted Lists
    h5("変換結果: Up-regulated Genes (Ensembl ID)"),
    withSpinner(verbatimTextOutput(ns("upGenesEnsemblList")), type = 6),
    hr(),
    h5("変換結果: Down-regulated Genes (Ensembl ID)"),
    withSpinner(verbatimTextOutput(ns("downGenesEnsemblList")), type = 6)
  )
}


# --- 3. Server Definition ---
#' Gene ID Converter Server Module with Auto-detection
# ... (関数の説明コメントは変更なし) ...
geneIdConverterServer <- function(id, rv_deg_results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns # Get namespace function

    # --- Helper Function: Detect Original Keytype ---
    # ... (detect_keytype 関数は変更なし) ...
    detect_keytype <- function(gene_ids) {
      if (length(gene_ids) == 0 || all(is.na(gene_ids))) return(NULL)
      sample_ids <- head(unique(na.omit(as.character(gene_ids))), 100)
      if (length(sample_ids) == 0) return(NULL)
      n_samples <- length(sample_ids)
      message(paste("[detect_keytype] Analyzing", n_samples, "unique non-NA IDs."))
      is_ensembl <- grepl("^ENS[A-Z]+G[0-9]+(\\.?[0-9]+)?$", sample_ids, ignore.case = TRUE)
      is_transcript <- grepl("^ENS[A-Z]+T[0-9]+(\\.?[0-9]+)?$", sample_ids, ignore.case = TRUE)
      is_ensembl_gene <- is_ensembl & !is_transcript
      ensembl_ratio <- sum(is_ensembl_gene) / n_samples
      message(paste("[detect_keytype] Ensembl Gene pattern match ratio:", round(ensembl_ratio, 2)))
      if (ensembl_ratio > 0.8) return("ENSEMBL")
      is_refseq <- grepl("^[NX][MP]_[0-9]+(\\.?[0-9]+)?$", sample_ids, ignore.case = TRUE)
      refseq_ratio <- sum(is_refseq) / n_samples
      message(paste("[detect_keytype] RefSeq pattern match ratio:", round(refseq_ratio, 2)))
      if (refseq_ratio > 0.8) return("REFSEQ")
      is_entrez <- grepl("^[0-9]+$", sample_ids)
      entrez_ratio <- sum(is_entrez) / n_samples
      message(paste("[detect_keytype] Numeric pattern match ratio:", round(entrez_ratio, 2)))
      if (entrez_ratio > 0.95) return("ENTREZID")
      message("[detect_keytype] Defaulting to SYMBOL.")
      return("SYMBOL")
    }


    # --- Observer: Auto-detect Keytype when DEG results change ---
    # ... (observe コードは変更なし) ...
     observe({
      deg_results <- rv_deg_results()
      req(deg_results, !is.null(deg_results$significant_up_genes) || !is.null(deg_results$significant_down_genes))
      message("[GeneIdConverter AutoDetect] DEG results updated. Detecting keytype...")
      all_ids <- c( na.omit(deg_results$significant_up_genes %||% character(0)), na.omit(deg_results$significant_down_genes %||% character(0)) )
      if (length(all_ids) > 0) {
        predicted_keytype <- detect_keytype(all_ids)
        message(paste("[GeneIdConverter AutoDetect] Predicted keytype:", predicted_keytype))
        if (!is.null(predicted_keytype)) {
          updateSelectInput(session, "original_keytype", selected = predicted_keytype)
          keytype_choices <- c("Gene Symbol" = "SYMBOL", "Entrez ID" = "ENTREZID", "RefSeq ID" = "REFSEQ", "Ensembl ID (変換元)" = "ENSEMBL")
          keytype_label <- names(keytype_choices)[which(keytype_choices == predicted_keytype)] %||% predicted_keytype
          showNotification( paste("元の遺伝子IDタイプを「", keytype_label, "」と推定しました。異なる場合は選択し直してください。"), type = "message", duration = 8 )
        } else { message("[GeneIdConverter AutoDetect] Could not predict keytype.") }
      } else { message("[GeneIdConverter AutoDetect] No significant genes found to detect keytype from."); updateSelectInput(session, "original_keytype", selected = "") }
    })


    # --- Reactive Value: Store Conversion Results ---
    # ... (変更なし) ...
    conversion_results <- reactiveVal(list( up_ensembl = NULL, down_ensembl = NULL, up_na_count = 0, down_na_count = 0, message = "変換を実行してください。" ))


    # --- Event Observer: Run Conversion on Button Click ---
    # ... (observeEvent コードは変更なし, 内部で %||% を使用) ...
    observeEvent(input$run_conversion, {
      message("[GeneIdConverter] 'Run Conversion' button clicked.")
      conversion_results(list(message="処理中...", up_ensembl=NULL, down_ensembl=NULL))
      deg_results <- rv_deg_results()
      req( deg_results, input$species != "", input$original_keytype != "" )
      message("[GeneIdConverter] Input checks passed for conversion.")
      selected_species <- input$species
      orgdb_pkg_name <- switch(selected_species, "Human" = "org.Hs.eg.db", "Mouse" = "org.Mm.eg.db", NULL)
      if (is.null(orgdb_pkg_name)) { msg <- "..."; conversion_results(list(message = msg)); showNotification(msg, type = "error"); return() }
      if (!requireNamespace(orgdb_pkg_name, quietly = TRUE)) { msg <- "..."; conversion_results(list(message = msg)); showNotification(msg, type = "error", duration = NULL); return() }
      library(orgdb_pkg_name, character.only = TRUE)
      orgDb <- get(orgdb_pkg_name)
      message(paste("[GeneIdConverter] Using OrgDb:", orgdb_pkg_name))
      original_keytype <- input$original_keytype
      up_genes_original <- deg_results$significant_up_genes %||% character(0) # ← ここで %||% が使われる
      down_genes_original <- deg_results$significant_down_genes %||% character(0) # ← ここで %||% が使われる
      message(paste("[GeneIdConverter] Keytype selected for conversion:", original_keytype))
      message(paste("[GeneIdConverter] Original Up genes count:", length(up_genes_original)))
      message(paste("[GeneIdConverter] Original Down genes count:", length(down_genes_original)))
      target_keytype <- "ENSEMBL"
      ensembl_up_genes <- NULL; ensembl_down_genes <- NULL; up_na_count <- 0; down_na_count <- 0
      conversion_message <- "変換完了。"
      if (length(up_genes_original) > 0) {
        tryCatch({
          ensembl_up_genes <- mapIds(orgDb, keys = unique(na.omit(up_genes_original)), column = target_keytype, keytype = original_keytype, multiVals = "first")
          up_na_count <- sum(is.na(ensembl_up_genes))
          message(paste("[GeneIdConverter] Up-genes mapped. NA count:", up_na_count))
          if (up_na_count == length(unique(na.omit(up_genes_original)))) { conversion_message <- paste(conversion_message, "警告: Up-regulated遺伝子のEnsembl IDが一つも見つかりませんでした。"); ensembl_up_genes <- NULL; }
        }, error = function(e) { errmsg <- paste("Up変換エラー:", e$message); message(paste("[GeneIdConverter] Error mapping Up genes:", errmsg)); conversion_message <<- paste(conversion_message, errmsg); ensembl_up_genes <<- NULL; })
      } else { ensembl_up_genes <- character(0) }
       if (length(down_genes_original) > 0) {
         tryCatch({
           ensembl_down_genes <- mapIds(orgDb, keys = unique(na.omit(down_genes_original)), column = target_keytype, keytype = original_keytype, multiVals = "first")
           down_na_count <- sum(is.na(ensembl_down_genes))
           message(paste("[GeneIdConverter] Down-genes mapped. NA count:", down_na_count))
           if (down_na_count == length(unique(na.omit(down_genes_original)))) { conversion_message <- paste(conversion_message, "警告: Down-regulated遺伝子のEnsembl IDが一つも見つかりませんでした。"); ensembl_down_genes <- NULL; }
         }, error = function(e) { errmsg <- paste("Down変換エラー:", e$message); message(paste("[GeneIdConverter] Error mapping Down genes:", errmsg)); conversion_message <<- paste(conversion_message, errmsg); ensembl_down_genes <<- NULL; })
      } else { ensembl_down_genes <- character(0) }
      final_results <- list( up_ensembl = ensembl_up_genes, down_ensembl = ensembl_down_genes, up_na_count = up_na_count, down_na_count = down_na_count, message = conversion_message )
      conversion_results(final_results)
      message("[GeneIdConverter] Conversion finished. Results stored.")
    })


    # --- Render Output: Display Converted Lists ---
    # ... (renderPrint コードは変更なし) ...
     output$upGenesEnsemblList <- renderPrint({
      results <- conversion_results()
      if (is.null(results$up_ensembl) && results$message != "変換完了。") { cat(results$message) }
      else if (is.null(results$up_ensembl) && length(rv_deg_results()$significant_up_genes %||% character(0)) == 0) { cat("元のリストに Up-regulated 遺伝子はありませんでした。") } # Use %||% here too
      else if (is.null(results$up_ensembl) || length(na.omit(results$up_ensembl)) == 0) {
         cat("Ensembl IDに変換できる Up-regulated 遺伝子はありませんでした。\n")
         cat(paste("(元のID数:", length(unique(na.omit(rv_deg_results()$significant_up_genes %||% character(0)))), ", 変換失敗(NA):", results$up_na_count, ")\n")) # Use %||% here too
         if (results$message != "変換完了。" && !is.null(results$up_ensembl)) cat(results$message)
      } else {
        up_genes_ensembl_no_na <- na.omit(results$up_ensembl)
        cat(paste(unique(up_genes_ensembl_no_na), collapse = ",\n"))
        if (results$up_na_count > 0) { cat(paste0("\n(注: ", results$up_na_count, " 個の一意な元のIDはEnsembl IDに変換できませんでした)")) }
      }
    })
    output$downGenesEnsemblList <- renderPrint({
      results <- conversion_results()
      if (is.null(results$down_ensembl) && results$message != "変換完了。") { cat(results$message) }
       else if (is.null(results$down_ensembl) && length(rv_deg_results()$significant_down_genes %||% character(0)) == 0) { cat("元のリストに Down-regulated 遺伝子はありませんでした。") } # Use %||% here too
       else if (is.null(results$down_ensembl) || length(na.omit(results$down_ensembl)) == 0) {
         cat("Ensembl IDに変換できる Down-regulated 遺伝子はありませんでした。\n")
         cat(paste("(元のID数:", length(unique(na.omit(rv_deg_results()$significant_down_genes %||% character(0)))), ", 変換失敗(NA):", results$down_na_count, ")\n")) # Use %||% here too
          if (results$message != "変換完了。" && !is.null(results$down_ensembl)) cat(results$message)
      } else {
        down_genes_ensembl_no_na <- na.omit(results$down_ensembl)
        cat(paste(unique(down_genes_ensembl_no_na), collapse = ",\n"))
        if (results$down_na_count > 0) { cat(paste0("\n(注: ", results$down_na_count, " 個の一意な元のIDはEnsembl IDに変換できませんでした)")) }
      }
    })


    # --- ★★★ Server関数内の %||% 定義は削除 ★★★ ---
    # `%||%` <- function(a, b) if (!is.null(a)) a else b # ← この行を削除する

    # --- Return Value ---
    # Return the reactiveVal itself, allowing parent modules to observe it.
    return(conversion_results)

  }) # End moduleServer
}

# --- End of R/module_gene_id_converter.R ---