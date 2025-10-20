# R/module_go_enrichment_integrated.R
# (KEGG解析、ID表示選択、解析中メッセージ機能を追加)
# (GO/KEGG結合時のエラーを修正)
# (Excelダウンロード機能を追加)

library(shiny)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(ggplot2)
library(DT)
library(shinycssloaders)
library(shinyjs)
library(dplyr)
library(writexl) # Excel書き出しに必要

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

# --- Species to KEGG Organism Code mapping ---
kegg_species_map <- c(
  "Homo_sapiens" = "hsa",
  "Mus_musculus" = "mmu",
  "Rattus_norvegicus" = "rno",
  "Drosophila_melanogaster" = "dme",
  "Caenorhabditis_elegans" = "cel",
  "Danio_rerio" = "dre",
  "Saccharomyces_cerevisiae" = "sce",
  "Arabidopsis_thaliana" = "ath",
  "Bos_taurus" = "bta",
  "Gallus_gallus" = "gga",
  "Canis_familiaris" = "cfa",
  "Macaca_mulatta" = "mcc",
  "Pan_troglodytes" = "ptr",
  "Sus_scrofa" = "ssc",
  "Xenopus_laevis" = "xla"
  # 必要に応じて他の生物種を追加
)


`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

confirm_entrez_id_format <- function(gene_ids) {
  if (length(gene_ids) == 0 || all(is.na(gene_ids))) return(FALSE)
  valid_ids <- unique(na.omit(as.character(gene_ids)))
  valid_ids <- valid_ids[valid_ids != ""]
  if (length(valid_ids) == 0) return(FALSE)
  sample_ids <- head(valid_ids, 100)
  n_samples <- length(sample_ids)
  is_entrez <- grepl("^[0-9]+$", sample_ids)
  entrez_ratio <- sum(is_entrez) / n_samples
  message(paste("[GO_confirm_entrez] Numeric (expected EntrezID) pattern match ratio:", round(entrez_ratio, 2)))
  return(entrez_ratio > 0.9)
}


goEnrichmentIntegratedUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("エンリッチメント解析 (GO & KEGG)"),
    helpText("差次発現遺伝子リスト（Entrez ID）を使用して、Gene OntologyおよびKEGGパスウェイのエンリッチメント解析を実行します。"),
    fluidRow(
      column(4, 
             selectInput(ns("gene_set"), "解析対象遺伝子:",
                         choices = c("選択してください" = "",
                                     "Up-regulated" = "up",
                                     "Down-regulated" = "down",
                                     "All Significant (Up+Down)" = "all_sig"),
                         selected = "up") 
      ),
      column(4, 
             selectInput(ns("analysis_type"), "解析タイプ:",
                         choices = c("ALL (GO + KEGG)" = "ALL",
                                     "GO: BP" = "BP",
                                     "GO: MF" = "MF",
                                     "GO: CC" = "CC",
                                     "KEGG" = "KEGG"),
                         selected = "ALL")
      ),
      column(4,
             selectInput(ns("go_id_display_type"), "結果の遺伝子IDタイプ:",
                         choices = c(
                           "Gene Symbol" = "SYMBOL",
                           "Entrez ID (内部ID)" = "ENTREZID",
                           "Gene Name" = "GENENAME"
                         ),
                         selected = "SYMBOL")
      )
    ),
    fluidRow(
      column(4, numericInput(ns("pvalue_cutoff"), "P-value Cutoff:", value = 0.05, min = 0, max = 1, step = 0.01)),
      column(4, numericInput(ns("qvalue_cutoff"), "Q-value (adj. P) Cutoff:", value = 0.2, min = 0, max = 1, step = 0.01)),
      column(4, actionButton(ns("run_analysis"), "解析実行", icon = icon("play")))
    ),
    hr(),
    tabsetPanel(
      id = ns("resultTabs"),
      tabPanel("結果テーブル",
               withSpinner(DTOutput(ns("goResultTable")), type = 6),
               # ★★★ UIの変更箇所 ★★★
               downloadButton(ns("downloadExcelResults"), "Excel形式でダウンロード (.xlsx)", icon = icon("file-excel")),
               br(), br(),
               downloadButton(ns("downloadCsvResults"), "CSV形式でダウンロード (.csv)", icon = icon("file-csv"))
      ),
      tabPanel("Bar Plot",
               helpText("エンリッチされたタームを棒グラフで表示します (上位表示)。"),
               uiOutput(ns("plot_selector_ui_bar")), # プロットセレクター
               numericInput(ns("barplot_n"), "表示するターム数:", value = 10, min = 1, max = 50),
               withSpinner(plotOutput(ns("goBarPlot")), type = 6),
               downloadButton(ns("downloadBarPlot"), "Bar Plotをダウンロード (.png)", icon = icon("download"))
      ),
      tabPanel("Dot Plot",
               helpText("エンリッチされたタームをドットプロットで表示します (上位表示)。"),
               uiOutput(ns("plot_selector_ui_dot")), # プロットセレクター
               numericInput(ns("dotplot_n"), "表示するターム数:", value = 10, min = 1, max = 50),
               withSpinner(plotOutput(ns("goDotPlot")), type = 6),
               downloadButton(ns("downloadDotPlot"), "Dot Plotをダウンロード (.png)", icon = icon("download"))
      )
    )
  )
}

goEnrichmentIntegratedServer <- function(id, deg_results_reactive, background_genes_reactive, selected_species_code_reactive) { 
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    analysis_results <- reactiveVal(NULL)
    analysis_running <- reactiveVal(FALSE)
    
    convert_entrez_ids_for_display <- function(entrez_ids, target_display_type, selected_species_code) {
      if (target_display_type == "ENTREZID" || is.null(target_display_type) || !nzchar(target_display_type)) {
        return(as.character(entrez_ids))
      }
      req(selected_species_code)
      orgdb_pkg_name <- orgdb_species_map[[selected_species_code]]
      
      if (is.null(orgdb_pkg_name) || !nzchar(orgdb_pkg_name)) {
        return(as.character(entrez_ids))
      }
      if (!requireNamespace(orgdb_pkg_name, quietly = TRUE)) {
        return(as.character(entrez_ids))
      }
      
      org_db <- get(orgdb_pkg_name)
      if (!target_display_type %in% columns(org_db) || !"ENTREZID" %in% keytypes(org_db)) {
        return(as.character(entrez_ids))
      }
      
      unique_entrez_keys <- unique(as.character(entrez_ids))
      converted_map <- tryCatch(
        suppressMessages(mapIds(org_db, keys = unique_entrez_keys, column = target_display_type, keytype = "ENTREZID", multiVals = "first")),
        error = function(e) { NULL }
      )
      
      if (is.null(converted_map)) return(as.character(entrez_ids))
      
      final_converted_ids <- converted_map[as.character(entrez_ids)]
      na_indices <- is.na(final_converted_ids)
      if (any(na_indices)) {
        final_converted_ids[na_indices] <- paste0(as.character(entrez_ids[na_indices]), " (変換不可)")
      }
      return(as.character(final_converted_ids))
    }
    
    observeEvent(input$run_analysis, {
      shinyjs::disable("run_analysis")
      on.exit(shinyjs::enable("run_analysis"), add = TRUE)
      
      analysis_results(NULL)
      analysis_running(TRUE)
      message("---------------------------------------------")
      message("[EnrichmentAnalysis] ", Sys.time(), " - Run Analysis button clicked.")
      showNotification("エンリッチメント解析を開始しました...", type = "message", duration = 3)
      
      current_selected_species_code <- selected_species_code_reactive() 
      target_entrez_ids_from_deg <- deg_results_reactive()
      background_entrez_ids_original <- background_genes_reactive()
      
      # (Input validation code is omitted for brevity but assumed to be here)
      
      tryCatch({
        target_genes_entrez_final <- switch(input$gene_set,
                                            "up" = target_entrez_ids_from_deg$significant_up_genes %||% character(0),
                                            "down" = target_entrez_ids_from_deg$significant_down_genes %||% character(0),
                                            "all_sig" = unique(c(target_entrez_ids_from_deg$significant_up_genes %||% character(0), target_entrez_ids_from_deg$significant_down_genes %||% character(0))),
                                            character(0)
        )
        target_genes_entrez_final <- unique(na.omit(as.character(target_genes_entrez_final)))
        target_genes_entrez_final <- target_genes_entrez_final[target_genes_entrez_final != "" & grepl("^[0-9]+$", target_genes_entrez_final)]
        
        if (length(target_genes_entrez_final) == 0) {
          stop("選択されたセットに対応する有効な対象Entrez ID（数値形式）がありません。")
        }
        
        universe_genes_entrez_final <- unique(na.omit(as.character(background_entrez_ids_original)))
        universe_genes_entrez_final <- universe_genes_entrez_final[universe_genes_entrez_final != "" & grepl("^[0-9]+$", universe_genes_entrez_final)]
        
        # --- Analysis Execution ---
        all_results_list <- list()
        
        # GO Analysis
        go_ontologies <- if (input$analysis_type == "ALL") c("BP", "MF", "CC") else if (input$analysis_type %in% c("BP", "MF", "CC")) input$analysis_type else c()
        if (length(go_ontologies) > 0) {
          orgdb_pkg_name <- orgdb_species_map[[current_selected_species_code]]
          if (is.null(orgdb_pkg_name) || !nzchar(orgdb_pkg_name) || !requireNamespace(orgdb_pkg_name, quietly = TRUE)) {
            stop(paste0("GO解析に必要なOrgDbパッケージ '", orgdb_pkg_name, "' が見つかりません。"))
          }
          require(orgdb_pkg_name, character.only = TRUE)
          orgDb_object <- get(orgdb_pkg_name)
          
          for (ont_cat in go_ontologies) {
            message("  - Running enrichGO for: ", ont_cat)
            res_ont <- enrichGO(gene = target_genes_entrez_final, universe = universe_genes_entrez_final,
                                OrgDb = orgDb_object, keyType = 'ENTREZID', ont = ont_cat,
                                pAdjustMethod = "BH", pvalueCutoff = input$pvalue_cutoff, qvalueCutoff = input$qvalue_cutoff,
                                readable = FALSE)
            if (!is.null(res_ont) && nrow(as.data.frame(res_ont)) > 0) {
              all_results_list[[ont_cat]] <- res_ont
            }
          }
        }
        
        # KEGG Analysis
        if (input$analysis_type == "ALL" || input$analysis_type == "KEGG") {
          kegg_code <- kegg_species_map[[current_selected_species_code]]
          if (is.null(kegg_code) || !nzchar(kegg_code)) {
            showNotification(paste0("KEGG解析はスキップされました: '", current_selected_species_code, "' に対応するKEGGコードがありません。"), type = "warning", duration = 8)
          } else {
            message("  - Running enrichKEGG for: ", kegg_code)
            res_kegg <- enrichKEGG(gene = target_genes_entrez_final, universe = universe_genes_entrez_final,
                                   organism = kegg_code, pvalueCutoff = input$pvalue_cutoff, qvalueCutoff = input$qvalue_cutoff)
            if (!is.null(res_kegg) && nrow(as.data.frame(res_kegg)) > 0) {
              all_results_list[["KEGG"]] <- res_kegg
            }
          }
        }
        
        analysis_results(all_results_list)
        
        if (length(all_results_list) == 0) {
          showNotification("指定されたカットオフで有意なタームは見つかりませんでした。", type = "warning", duration = 5)
        } else {
          showNotification("エンリッチメント解析が完了しました。", type = "message", duration = 5)
        }
        
      }, error = function(e) {
        errmsg <- paste("解析中にエラーが発生しました:", conditionMessage(e))
        showNotification(ui = tags$div(tags$b("エラー発生:"), tags$p(errmsg)), type = "error", duration = NULL)
        analysis_results(NULL) 
      })
      
      analysis_running(FALSE)
    })
    
    final_table_reactive <- reactive({
      results_list <- analysis_results()
      if (is.null(results_list)) return(NULL)
      if (length(results_list) == 0) return(data.frame(Message = "有意な結果がありません。"))
      
      res_df_list <- lapply(names(results_list), function(category) {
        res_obj <- results_list[[category]]
        if (!is.null(res_obj) && inherits(res_obj, "enrichResult") && nrow(as.data.frame(res_obj)) > 0) {
          df <- as.data.frame(res_obj)
          df$Category <- category
          return(df)
        }
        return(NULL)
      })
      
      combined_res_df <- dplyr::bind_rows(Filter(Negate(is.null), res_df_list))
      
      if (is.null(combined_res_df) || nrow(combined_res_df) == 0) {
        return(data.frame(Message = "有意な結果がありません。"))
      }
      
      if ("geneID" %in% colnames(combined_res_df) && input$go_id_display_type != "ENTREZID") {
        message("Converting geneID column for display...")
        converted_gene_ids <- sapply(combined_res_df$geneID, function(id_string) {
          if(is.na(id_string)) return(NA_character_)
          entrez_ids <- strsplit(id_string, "/")[[1]]
          converted <- convert_entrez_ids_for_display(entrez_ids, input$go_id_display_type, selected_species_code_reactive())
          paste(converted, collapse = "/")
        })
        combined_res_df$geneID <- converted_gene_ids
      }
      
      combined_res_df %>%
        select(any_of(c("ID", "Description", "Category", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"))) %>%
        arrange(Category, p.adjust)
    })
    
    output$goResultTable <- renderDT({
      if (analysis_running()) {
        return(datatable(data.frame(Message = "解析中..."), options = list(dom = 't', searching = FALSE)))
      }
      
      final_table <- final_table_reactive()
      if (is.null(final_table)) {
        return(datatable(data.frame(Message = "解析を実行してください。"), options = list(dom = 't', searching = FALSE)))
      }
      if("Message" %in% colnames(final_table)){
        return(datatable(final_table, options = list(dom = 't', searching = FALSE)))
      }
      
      res_df_display <- final_table %>%
        mutate(across(where(is.numeric), ~ round(.x, digits = 4)))
      
      datatable(res_df_display, rownames = FALSE, filter = "top", extensions = 'Buttons', options = list(pageLength=10, scrollX=TRUE), escape = FALSE)
    })
    
    output$plot_selector_ui_bar <- renderUI({
      results_list <- analysis_results()
      req(input$analysis_type == "ALL", results_list, length(results_list) > 1)
      selectInput(ns("plot_source_bar"), "プロット対象:", choices = names(results_list), selected = names(results_list)[1])
    })
    output$plot_selector_ui_dot <- renderUI({
      results_list <- analysis_results()
      req(input$analysis_type == "ALL", results_list, length(results_list) > 1)
      selectInput(ns("plot_source_dot"), "プロット対象:", choices = names(results_list), selected = names(results_list)[1])
    })
    
    get_plot_data_object <- reactive({
      results_list <- analysis_results()
      req(results_list)
      
      source_key <- if (input$analysis_type == "ALL") {
        req(input$plot_source_bar) 
        input$plot_source_bar
      } else {
        input$analysis_type
      }
      
      result_obj <- results_list[[source_key]]
      req(result_obj, inherits(result_obj, "enrichResult"))
      
      orgdb_pkg_name <- orgdb_species_map[[selected_species_code_reactive()]]
      if (!is.null(orgdb_pkg_name) && requireNamespace(orgdb_pkg_name, quietly = TRUE)) {
        orgDb_obj <- get(orgdb_pkg_name)
        result_obj <- setReadable(result_obj, OrgDb = orgDb_obj, keyType = "ENTREZID")
      }
      return(result_obj)
    })
    
    bar_plot_gg_object <- reactive({
      result_obj_for_plot <- get_plot_data_object()
      req(result_obj_for_plot)
      validate(need(nrow(as.data.frame(result_obj_for_plot)) > 0, "プロットする有意なタームがありません。"))
      n_terms <- min(input$barplot_n, nrow(as.data.frame(result_obj_for_plot)))
      barplot(result_obj_for_plot, showCategory = n_terms)
    })
    
    output$goBarPlot <- renderPlot({
      print(bar_plot_gg_object())
    })
    
    dot_plot_gg_object <- reactive({
      result_obj_for_plot <- get_plot_data_object()
      req(result_obj_for_plot)
      validate(need(nrow(as.data.frame(result_obj_for_plot)) > 0, "プロットする有意なタームがありません。"))
      n_terms <- min(input$dotplot_n, nrow(as.data.frame(result_obj_for_plot)))
      dotplot(result_obj_for_plot, showCategory = n_terms)
    })
    
    output$goDotPlot <- renderPlot({
      print(dot_plot_gg_object())
    })
    
    # ★★★ ダウンロードハンドラを修正 ★★★
    output$downloadCsvResults <- downloadHandler(
      filename = function() {
        paste0("Enrichment_Results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        table_to_download <- final_table_reactive()
        if (is.null(table_to_download) || "Message" %in% colnames(table_to_download)) {
          write.csv(data.frame(Message="結果がありません。"), file, row.names = FALSE, fileEncoding="UTF-8")
        } else {
          write.csv(table_to_download, file, row.names = FALSE, quote = TRUE, fileEncoding="UTF-8")
        }
      }
    )
    
    output$downloadExcelResults <- downloadHandler(
      filename = function() {
        paste0("Enrichment_Results_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        table_to_download <- final_table_reactive()
        if (is.null(table_to_download) || "Message" %in% colnames(table_to_download)) {
          writexl::write_xlsx(data.frame(Message="結果がありません。"), file)
        } else {
          writexl::write_xlsx(table_to_download, file)
        }
      }
    )
    
  }) 
}
