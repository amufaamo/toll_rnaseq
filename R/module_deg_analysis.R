# R/module_deg_analysis.R
# (フィルタリング指標としてFDR/P-valueを選択可能にし、プロットも連動させる修正版)

library(shiny)
library(edgeR)
library(DESeq2) # DESeq2パッケージをロード
library(dplyr)
library(DT)
library(plotly)
library(ggplot2)
library(shinycssloaders)
library(AnnotationDbi) # ID変換に必要
library(writexl)       # Excelダウンロードに必要
library(tibble)        # rownames_to_columnのため

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
  if (length(ids_clean) == 0) return("UNKNOWN")
  n_total <- length(ids_clean); n_sample <- min(n_total, 1000)
  ids_sample <- sample(ids_clean, n_sample)
  if (mean(grepl("^ENS[A-Z0-9]*[FPTG]\\d{10,}(\\.\\d+)?$", ids_sample, ignore.case = TRUE)) > 0.8) return("ENSEMBL")
  if (mean(grepl("^[NX][CMPWTZ]_[0-9]+(\\.\\d+)?$", ids_sample, ignore.case = TRUE)) > 0.8) return("REFSEQ")
  if (mean(grepl("^[0-9]+$", ids_sample)) > 0.9) return("ENTREZID")
  if (mean(grepl("^([A-Za-z][A-Za-z0-9-]*[A-Za-z0-9]|[A-Za-z])$", ids_sample)) > 0.7 &&
      mean(grepl("^ENS", ids_sample, ignore.case=TRUE)) < 0.2 &&
      mean(grepl("^[NX][M_]", ids_sample, ignore.case=TRUE)) < 0.2 &&
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
      # ★★★ UIの変更箇所 ★★★
      h4("解析手法の選択"),
      radioButtons(ns("deg_method"), "解析パッケージ:",
                   choices = c("edgeR" = "edgeR", "DESeq2" = "DESeq2"),
                   selected = "edgeR", inline = TRUE),
      helpText("edgeRは迅速、DESeq2はより多くのサンプルで頑健です。"),
      hr(),
      h4("比較グループ選択"),
      uiOutput(ns("degGroupSelectionUI")),
      hr(),
      h4("フィルタリング閾値"),
      helpText("以下の閾値は、表示されるテーブルとプロットに適用されます。"),
      radioButtons(ns("sig_metric"), "有意差の指標:",
                   choices = c("FDR (adjusted P-value)" = "FDR",
                               "P-value" = "PValue"),
                   selected = "FDR", inline = TRUE),
      conditionalPanel(
        condition = paste0("input['", ns("sig_metric"), "'] == 'FDR'"),
        numericInput(ns("degFDR"), "FDR 閾値:", value = 0.05, min = 0, max = 1, step = 0.01)
      ),
      conditionalPanel(
        condition = paste0("input['", ns("sig_metric"), "'] == 'PValue'"),
        numericInput(ns("degPValue"), "P-value 閾値:", value = 0.05, min = 0, max = 1, step = 0.01)
      ),
      numericInput(ns("degLogFC"), "Log2 Fold Change 閾値 (|LogFC| >):", value = 1, min = 0, step = 0.1),
      hr(),
      actionButton(ns("runDEG"), "DEG解析実行", icon = icon("play")),
      hr(),
      h4("ボルケーノプロット設定 (任意)"),
      uiOutput(ns("highlightGenesUI")),
      hr(),
      h4("表示IDタイプ選択"),
      selectInput(ns("deg_id_display_type"), "結果テーブルの遺伝子IDタイプ:",
                  choices = c("Gene Symbol" = "SYMBOL", "Entrez ID (内部ID)" = "ENTREZID"),
                  selected = "SYMBOL"),
      hr(),
      h4("結果ダウンロード"),
      p("DEG解析結果の全遺伝子リスト（フィルタリングなし）をダウンロードします。"),
      downloadButton(ns("downloadExcelResults"), "Excel形式でダウンロード (.xlsx)", icon = icon("file-excel")),
      br(),br(),
      downloadButton(ns("downloadCsvResults"), "CSV形式でダウンロード (.csv)", icon = icon("file-csv"))
    ),
    mainPanel(
      width = 8,
      h4("解析サマリー"),
      withSpinner(verbatimTextOutput(ns("degSummary")), type = 6),
      hr(),
      h4("MAプロット (旧MDプロット)"),
      helpText("赤点・青点が有意に発現変動している遺伝子を示します。"),
      withSpinner(plotOutput(ns("degMDPlot")), type = 6),
      hr(),
      h4("ボルケーノプロット"),
      helpText("有意な発現変動遺伝子をLogFCと選択した指標でプロットします。"),
      withSpinner(plotlyOutput(ns("degVolcanoPlot")), type = 6),
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
    
    detected_deg_input_id_type <- reactiveVal("ENTREZID")
    observeEvent(rv$merged_data, {
      req(rv$merged_data)
      detected_deg_input_id_type("ENTREZID")
      updateSelectInput(session, "deg_id_display_type",
                        choices = c("Gene Symbol" = "SYMBOL", "Entrez ID (内部ID)" = "ENTREZID", "Gene Name" = "GENENAME"),
                        selected = "SYMBOL")
    })
    convert_entrez_ids_for_display <- function(entrez_ids, target_display_type, selected_species_code) {
      if (target_display_type == "ENTREZID" || is.null(target_display_type) || !nzchar(target_display_type)) return(as.character(entrez_ids))
      req(selected_species_code); orgdb_pkg_name <- orgdb_species_map[[selected_species_code]]
      if (is.null(orgdb_pkg_name) || !requireNamespace(orgdb_pkg_name, quietly = TRUE)) return(as.character(entrez_ids))
      org_db <- get(orgdb_pkg_name)
      if (!target_display_type %in% columns(org_db) || !"ENTREZID" %in% keytypes(org_db)) return(as.character(entrez_ids))
      unique_entrez_keys <- unique(as.character(entrez_ids))
      converted_map <- tryCatch(suppressMessages(mapIds(org_db, keys = unique_entrez_keys, column = target_display_type, keytype = "ENTREZID", multiVals = "first")), error = function(e) NULL)
      if (is.null(converted_map)) return(as.character(entrez_ids))
      final_converted_ids <- converted_map[as.character(entrez_ids)]
      na_indices <- is.na(final_converted_ids)
      if (any(na_indices)) final_converted_ids[na_indices] <- paste0(as.character(entrez_ids[na_indices]), " (変換不可)")
      return(as.character(final_converted_ids))
    }
    
    output$degGroupSelectionUI <- renderUI({
      req(rv$sample_metadata)
      # ★★★ 修正箇所: アクティブなサンプルのグループのみを表示 ★★★
      active_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      available_groups <- unique(active_metadata$group)
      
      if (length(na.omit(available_groups)) < 2) {
        return(tags$p("比較可能なグループが2つ以上ありません（アクティブなサンプル内）。"))
      }
      selected_target <- available_groups[1]
      selected_ref <- if(length(available_groups) > 1) available_groups[2] else available_groups[1]
      tagList(
        selectInput(ns("target_group"), "Target:", choices = available_groups, selected = selected_target),
        selectInput(ns("reference_group"), "Reference:", choices = available_groups, selected = selected_ref)
      )
    })
    
    output$highlightGenesUI <- renderUI({
      all_entrez_ids_in_deg <- NULL
      if (!is.null(rv$deg_results) && !is.null(rv$deg_results$top_tags) && "Geneid" %in% colnames(rv$deg_results$top_tags$table)) {
        all_entrez_ids_in_deg <- unique(as.character(rv$deg_results$top_tags$table$Geneid))
      }
      display_type <- input$deg_id_display_type
      id_type_label <- switch(display_type, "SYMBOL" = "Gene Symbol", "ENTREZID" = "Entrez ID", "GENENAME" = "Gene Name", display_type)
      highlight_label <- paste0("ハイライトする遺伝子を選択 (", id_type_label, "で検索):")
      choices_for_ui <- NULL
      if (!is.null(all_entrez_ids_in_deg)) {
        choices_for_ui <- convert_entrez_ids_for_display(all_entrez_ids_in_deg, display_type, rv$selected_species)
        choices_for_ui <- sort(unique(na.omit(choices_for_ui[!grepl("\\(変換不可\\)$", choices_for_ui)])))
      }
      if (is.null(choices_for_ui)) {
        return(tags$p(em("DEG解析を実行すると、ここで遺伝子を選択できます。")))
      }
      selectizeInput(ns("highlight_genes_select"), label = highlight_label, choices = choices_for_ui,
                     selected = NULL, multiple = TRUE, options = list(placeholder = '入力または選択...', plugins = list('remove_button')))
    })
    
    observeEvent(input$runDEG, {
      req(rv$merged_data, rv$sample_metadata, input$target_group, input$reference_group, input$deg_method)
      
      sig_metric_val <- if (input$sig_metric == "FDR") input$degFDR else input$degPValue
      logfc_val <- input$degLogFC
      validate(need(is.numeric(sig_metric_val) && sig_metric_val >= 0 && sig_metric_val <= 1, " 有意差の閾値は0から1の数値で入力してください。"),
               need(is.numeric(logfc_val) && logfc_val >= 0, "LogFC閾値は0以上の数値で入力してください。"))
      
      active_samples_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      validate(need(nrow(active_samples_metadata) > 0, "解析対象のサンプルが選択されていません。"))
      active_sample_names <- active_samples_metadata$current_name
      
      counts_df_iso_entrez <- isolate(rv$merged_data)
      counts_df_iso_entrez_active <- counts_df_iso_entrez[, c("Geneid", intersect(colnames(counts_df_iso_entrez), active_sample_names)), drop = FALSE]
      samples_metadata_iso <- active_samples_metadata
      
      keep_vector_initial_filter <- isolate(rv$filtered_keep)
      validate(need(is.data.frame(counts_df_iso_entrez_active) && "Geneid" %in% colnames(counts_df_iso_entrez_active), "カウントデータが不正です。"))
      
      count_matrix_full_entrez <- as.matrix(counts_df_iso_entrez_active[,-1])
      rownames(count_matrix_full_entrez) <- counts_df_iso_entrez_active$Geneid
      validate(need(length(keep_vector_initial_filter) == nrow(isolate(rv$merged_data)), "フィルタリング情報が古いため、フィルタリングタブを再確認してください。"))
      
      count_matrix_ui_filtered_entrez <- count_matrix_full_entrez[keep_vector_initial_filter, , drop = FALSE]
      validate(need(nrow(count_matrix_ui_filtered_entrez) > 0, "フィルタリングの結果、遺伝子が残りませんでした。"))
      
      target_group_name <- input$target_group; reference_group_name <- input$reference_group
      samples_to_keep <- samples_metadata_iso$current_name[samples_metadata_iso$group %in% c(target_group_name, reference_group_name)]
      valid_samples <- intersect(colnames(count_matrix_ui_filtered_entrez), samples_to_keep)
      
      validate(need(length(valid_samples) >= 2 && length(unique(samples_metadata_iso$group[samples_metadata_iso$current_name %in% valid_samples])) == 2, "各比較グループに最低1サンプル必要です。"))
      
      count_matrix_for_deg <- count_matrix_ui_filtered_entrez[, valid_samples, drop = FALSE]
      metadata_for_deg <- samples_metadata_iso %>% filter(current_name %in% valid_samples)
      group_factor <- factor(metadata_for_deg$group) %>% relevel(ref = reference_group_name)
      
      deg_results_list <- tryCatch({
        # ★★★ 解析手法による分岐 ★★★
        if (input$deg_method == "edgeR") {
          dge <- DGEList(counts = count_matrix_for_deg, group = group_factor, genes = data.frame(Geneid = rownames(count_matrix_for_deg)))
          design_filter <- model.matrix(~group, data=dge$samples)
          dge <- dge[filterByExpr(dge, design=design_filter), , keep.lib.sizes=FALSE]
          validate(need(nrow(dge) > 0, "filterByExprの結果、遺伝子が0になりました。フィルタリング条件を緩めてください。"))
          
          dge <- calcNormFactors(dge)
          design <- model.matrix(~0 + group, data = dge$samples)
          colnames(design) <- levels(dge$samples$group)
          contrast_str <- paste(target_group_name, "-", reference_group_name)
          my_contrast <- makeContrasts(contrasts = contrast_str, levels = design)
          dge <- estimateDisp(dge, design, robust = TRUE)
          fit <- glmQLFit(dge, design, robust = TRUE)
          qlf <- glmQLFTest(fit, contrast = my_contrast)
          top_tags_result <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none")
          top_tags_table <- as.data.frame(top_tags_result$table)
          
          sig_col_name <- input$sig_metric
          up_genes <- top_tags_table %>% filter(!!sym(sig_col_name) < sig_metric_val, logFC > logfc_val) %>% pull(Geneid) %>% unique()
          down_genes <- top_tags_table %>% filter(!!sym(sig_col_name) < sig_metric_val, logFC < -logfc_val) %>% pull(Geneid) %>% unique()
          
          list(
            analysis_method = "edgeR",
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
          dds <- DESeqDataSetFromMatrix(countData = round(count_matrix_for_deg), colData = colData, design = ~ group)
          keep <- rowSums(counts(dds)) >= 10
          dds <- dds[keep,]
          validate(need(nrow(dds) > 0, "フィルタリングの結果、遺伝子が0になりました。"))
          
          dds <- DESeq(dds)
          res <- results(dds, contrast = c("group", target_group_name, reference_group_name))
          
          res_df <- as.data.frame(res) %>%
            rownames_to_column("Geneid") %>%
            rename(logFC = log2FoldChange, PValue = pvalue, FDR = padj)
          res_df$FDR[is.na(res_df$FDR)] <- 1
          
          top_tags_like_object <- list(table = res_df)
          
          sig_col_name <- input$sig_metric
          up_genes <- res_df %>% filter(!is.na(!!sym(sig_col_name)), !!sym(sig_col_name) < sig_metric_val, !is.na(logFC), logFC > logfc_val) %>% pull(Geneid) %>% unique()
          down_genes <- res_df %>% filter(!is.na(!!sym(sig_col_name)), !!sym(sig_col_name) < sig_metric_val, !is.na(logFC), logFC < -logfc_val) %>% pull(Geneid) %>% unique()
          
          list(
            analysis_method = "DESeq2",
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
      }, error = function(e) {
        showNotification(paste("DEG解析エラー:", e$message), type = "error", duration = 15); NULL
      })
      rv$deg_results <- deg_results_list
    })
    
    # (以降のコードは変更なし)
    # ... (output$degSummary から下のコードは元のままでOK) ...
    output$degSummary <- renderPrint({
      req(rv$deg_results, rv$deg_results$status == "analysis_completed")
      res <- rv$deg_results
      cat("解析手法:", res$analysis_method, "\n")
      cat("比較:", paste0(res$comparison[1], "_vs_", res$comparison[2]), "\n\n")
      
      current_sig_metric <- input$sig_metric
      current_sig_threshold <- if (current_sig_metric == "FDR") input$degFDR else input$degPValue
      current_logfc_threshold <- input$degLogFC
      
      cat(" 有意変動遺伝子サマリー (現在のUI設定 ", current_sig_metric, " < ", current_sig_threshold, ", |LogFC| > ", current_logfc_threshold, "):\n")
      
      if (res$analysis_method == "edgeR") {
        dt <- tryCatch(
          decideTests(res$qlf, p.value = current_sig_threshold, lfc = current_logfc_threshold,
                      adjust.method = if(current_sig_metric == "FDR") "BH" else "none"),
          error = function(e) { cat("サマリー計算エラー:", e$message, "\n"); NULL }
        )
        if (!is.null(dt)) print(summary(dt))
      } else { # DESeq2
        res_table <- res$top_tags$table
        up <- res_table %>% filter(!is.na(!!sym(current_sig_metric)), !!sym(current_sig_metric) < current_sig_threshold, !is.na(logFC), logFC > current_logfc_threshold) %>% nrow() 
        down <- res_table %>% filter(!is.na(!!sym(current_sig_metric)), !!sym(current_sig_metric) < current_sig_threshold, !is.na(logFC), logFC < -current_logfc_threshold) %>% nrow() 
        notsig <- nrow(res_table) - up - down
        summary_mat <- matrix(c(down, notsig, up), ncol = 1, dimnames = list(c("Down", "NotSig", "Up"), paste(res$comparison[1], "vs", res$comparison[2])))
        print(summary_mat)
      }
      
      cat("\n 有意変動遺伝子数 (実行時の閾値に基づくリスト):\n")
      cat("  (実行時 ", res$metric_at_run, " < ", res$sig_threshold_at_run, ", |LogFC| > ", res$logfc_threshold_at_run, "):\n")
      cat("  Up-regulated:", length(res$significant_up_genes), "\n")
      cat("  Down-regulated:", length(res$significant_down_genes), "\n")
      cat("  Total (Unique):", length(res$significant_genes), "\n")
      cat("背景遺伝子数 (解析対象):", length(res$background_genes_original), "\n")
    })
    
    filtered_and_converted_table <- reactive({
      req(rv$deg_results, input$deg_id_display_type, rv$selected_species)
      res_table <- as.data.frame(rv$deg_results$top_tags$table)
      
      sig_metric_filter <- input$sig_metric
      sig_threshold_filter <- if (sig_metric_filter == "FDR") input$degFDR else input$degPValue
      logfc_threshold_filter <- input$degLogFC
      
      filtered_table <- res_table %>%
        filter(!is.na(!!sym(sig_metric_filter)), !!sym(sig_metric_filter) < sig_threshold_filter, !is.na(logFC), abs(logFC) > logfc_threshold_filter) %>%
        arrange(!!sym(sig_metric_filter))
      
      if (nrow(filtered_table) == 0) return(data.frame(メッセージ = "指定した閾値で有意な遺伝子はありません。"))
      
      display_table <- data.table::copy(filtered_table)
      display_table[[1]] <- convert_entrez_ids_for_display(display_table$Geneid, input$deg_id_display_type, rv$selected_species)
      display_col_name <- names(c("Gene Symbol"="SYMBOL", "Entrez ID"="ENTREZID", "Gene Name"="GENENAME"))[c("SYMBOL", "ENTREZID", "GENENAME") == input$deg_id_display_type]
      colnames(display_table)[1] <- if(length(display_col_name) > 0) display_col_name else input$deg_id_display_type
      
      display_table %>% select(any_of(c(colnames(display_table)[1], "logFC", "logCPM", "baseMean", "PValue", "FDR", "LR", "F", "stat")))
    })
    
    output$degResultTable <- renderDT({
      final_table <- filtered_and_converted_table(); req(final_table)
      if ("メッセージ" %in% colnames(final_table)) return(datatable(final_table, options = list(dom = 't')))
      datatable(final_table, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE)) %>%
        formatRound(columns = intersect(c("logFC", "logCPM", "baseMean"), colnames(final_table)), digits = 3) %>%
        formatSignif(columns = intersect(c("PValue", "FDR"), colnames(final_table)), digits = 3)
    })
    
    output$degMDPlot <- renderPlot({
      req(rv$deg_results)
      res <- rv$deg_results
      plot_title <- paste0(res$comparison[1], "_vs_", res$comparison[2])
      
      if (res$analysis_method == "edgeR") {
        req(res$qlf)
        sig_metric_plot <- input$sig_metric
        sig_threshold_plot <- if (sig_metric_plot == "FDR") input$degFDR else input$degPValue
        logfc_threshold_plot <- input$degLogFC
        status_vec <- decideTests(res$qlf, p.value = sig_threshold_plot, lfc = logfc_threshold_plot, adjust.method = if(sig_metric_plot == "FDR") "BH" else "none")
        plotMD(res$qlf, status = status_vec, values = c(1, -1), col = c("red", "blue"), legend = "topright", main = plot_title)
        abline(h = 0, col = "grey", lty = 2)
      } else { # DESeq2
        req(res$deseq_res)
        # DESeq2のplotMAはFDR(padj)で色付けするため、alphaにFDR閾値を指定
        plotMA(res$deseq_res, alpha = input$degFDR, main = paste("MA Plot (DESeq2) -", plot_title), ylim = c(-5, 5))
        # UIで設定したLogFC閾値の線を追加
        abline(h = c(-input$degLogFC, input$degLogFC), col = "dodgerblue", lty = 2)
      }
    })
    
    volcano_data_reactive <- reactive({
      req(rv$deg_results, input$deg_id_display_type, rv$selected_species)
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
          tooltip_text = paste("Gene:", DisplayGeneid, "<br>logFC:", round(logFC, 3),
                               "<br>PValue:", format(PValue, scientific=T, digits=3),
                               "<br>FDR:", format(FDR, scientific=T, digits=3))
        )
    })
    
    output$degVolcanoPlot <- renderPlotly({
      plot_data <- volcano_data_reactive(); req(plot_data)
      
      sig_metric_volcano <- input$sig_metric
      sig_threshold_volcano <- if (sig_metric_volcano == "FDR") input$degFDR else input$degPValue
      logfc_lines_volcano <- c(-input$degLogFC, input$degLogFC)
      
      y_axis_label <- if(sig_metric_volcano == "FDR") "-Log10 FDR" else "-Log10 P-value"
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
    
    display_gene_list_converted <- function(entrez_gene_ids) {
      req(input$deg_id_display_type, rv$selected_species)
      if (length(entrez_gene_ids) == 0) return("該当する遺伝子はありません。")
      convert_entrez_ids_for_display(entrez_gene_ids, input$deg_id_display_type, rv$selected_species) %>%
        paste(collapse = ",
")
    }
    output$upGenesList <- renderPrint({ req(rv$deg_results); cat(display_gene_list_converted(rv$deg_results$significant_up_genes)) })
    output$downGenesList <- renderPrint({ req(rv$deg_results); cat(display_gene_list_converted(rv$deg_results$significant_down_genes)) })
    
    download_table_reactive <- reactive({
      req(rv$deg_results)
      res_table <- as.data.frame(rv$deg_results$top_tags$table)
      # FDR/PValueでソート
      sort_col <- if("FDR" %in% colnames(res_table)) "FDR" else "PValue"
      res_table <- res_table %>% arrange(!!sym(sort_col))
      
      display_table <- data.table::copy(res_table)
      display_table[[1]] <- convert_entrez_ids_for_display(display_table$Geneid, input$deg_id_display_type, rv$selected_species)
      display_col_name <- names(c("Gene Symbol"="SYMBOL", "Entrez ID"="ENTREZID", "Gene Name"="GENENAME"))[c("SYMBOL", "ENTREZID", "GENENAME") == input$deg_id_display_type]
      colnames(display_table)[1] <- if(length(display_col_name) > 0) display_col_name else input$deg_id_display_type
      
      # 手法に応じて列を選択
      cols_to_select <- if(rv$deg_results$analysis_method == "edgeR") {
        c(colnames(display_table)[1], "logFC", "logCPM", "PValue", "FDR", "LR", "F")
      } else {
        c(colnames(display_table)[1], "logFC", "baseMean", "PValue", "FDR", "stat")
      }
      display_table %>% select(any_of(cols_to_select))
    })
    common_filename_parts_deg <- reactive({
      comp_name <- "DEG_Results"; if (!is.null(rv$deg_results)) comp_name <- paste0("DEG_results_", rv$deg_results$comparison[1], "_vs_", rv$deg_results$comparison[2])
      method_part <- rv$deg_results$analysis_method %||% "method"
      id_type_part <- input$deg_id_display_type %||% "ID"
      paste0(comp_name, "_", method_part, "_", id_type_part, "_all-genes_", Sys.Date())
    })
    output$downloadExcelResults <- downloadHandler(
      filename = function() paste0(common_filename_parts_deg(), ".xlsx"),
      content = function(file) writexl::write_xlsx(download_table_reactive(), file)
    )
    output$downloadCsvResults <- downloadHandler(
      filename = function() paste0(common_filename_parts_deg(), ".csv"),
      content = function(file) write.csv(download_table_reactive(), file, row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
    )
    
    `%||%` <- function(a, b) if (!is.null(a)) a else b
  })
}