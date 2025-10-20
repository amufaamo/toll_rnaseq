# R/module_gsea.R
# (ランキング指標のデフォルト変更と選択肢追加)

# UI関数
gseaUI <- function(id) {
  ns <- NS(id)
  
  # msigdbrで利用可能な生物種のリストを取得
  # UI表示用の選択肢を作成 (例: "Human (Homo sapiens)" = "Homo sapiens")
  msigdbr_supported_species_choices <- c(
    "Human (Homo sapiens)" = "Homo sapiens",
    "Mouse (Mus musculus)" = "Mus musculus",
    "Rat (Rattus norvegicus)" = "Rattus norvegicus",
    "Zebrafish (Danio rerio)" = "Danio rerio",
    "Fly (Drosophila melanogaster)" = "Drosophila melanogaster",
    "Worm (Caenorhabditis elegans)" = "Caenorhabditis elegans",
    "Yeast (Saccharomyces cerevisiae)" = "Saccharomyces cerevisiae",
    "Bovine (Bos taurus)" = "Bos taurus",
    "Chicken (Gallus gallus)" = "Gallus gallus",
    "Dog (Canis lupus familiaris)" = "Canis lupus familiaris",
    "Rhesus monkey (Macaca mulatta)" = "Macaca mulatta",
    "Chimpanzee (Pan troglodytes)" = "Pan troglodytes",
    "Pig (Sus scrofa)" = "Sus scrofa"
  )
  
  fluidPage(
    h4("GSEA (Gene Set Enrichment Analysis)"),
    helpText("DEG解析結果の全遺伝子ランキングに基づき、指定した遺伝子セットコレクション（パスウェイなど）が、比較条件間で全体として発現上昇または下降しているかを評価します。DEG解析結果の遺伝子ID(EntrezID)はGene Symbolに変換して実行されます。"),
    hr(),
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("GSEA 設定"),
        selectInput(ns("gsea_species"), "生物種 (遺伝子セット用):",
                    choices = msigdbr_supported_species_choices,
                    selected = "Homo sapiens"),
        radioButtons(ns("gsea_selection_method"), "遺伝子セット指定方法:",
                     choices = c("カテゴリから選択" = "category",
                                 "キーワードで検索・選択" = "search"),
                     selected = "category", inline = TRUE),
        conditionalPanel(
          condition = paste0("input['", ns("gsea_selection_method"), "'] == 'category'"),
          selectInput(ns("gsea_collection_category"), "遺伝子セット コレクション (MSigDB):",
                      choices = c(
                        "Hallmark gene sets" = "H",
                        "Curated: KEGG" = "C2_KEGG",
                        "Curated: Reactome" = "C2_REACTOME",
                        "Curated: BioCarta" = "C2_BIOCARTA",
                        "GO Biological Process" = "C5_BP",
                        "GO Cellular Component" = "C5_CC",
                        "GO Molecular Function" = "C5_MF",
                        "Oncogenic Signatures" = "C6",
                        "Immunologic Signatures" = "C7"
                      ), selected = "H")
        ),
        conditionalPanel(
          condition = paste0("input['", ns("gsea_selection_method"), "'] == 'search'"),
          selectizeInput(
            ns("gsea_search_pathway_name"),
            label = "遺伝子セット名を検索・選択:",
            choices = NULL,
            selected = NULL,
            multiple = FALSE,
            options = list(
              placeholder = '例: SAUL_SEN_MAYO または mayo などキーワードを入力',
              maxOptions = 30000
            )
          ),
          helpText("生物種を変更した場合、再度検索してください。利用可能な全遺伝子セットから検索します。")
        ),
        # ★★★ UIの変更箇所 ★★★
        radioButtons(ns("gsea_rank_metric"), "ランキング指標 (DEG結果より):",
                     choices = c("符号付き -log10(PValue)" = "signed_p",
                                 "logFC x (-log10(PValue))" = "logFC_x_signed_p",
                                 "log2 Fold Change" = "logFC",
                                 "符号付き 検定統計量 (F)" = "stat"),
                     selected = "signed_p"), # デフォルトを signed_p に変更
        numericInput(ns("gsea_min_size"), "最小遺伝子セットサイズ:", value = 15, min = 1),
        numericInput(ns("gsea_max_size"), "最大遺伝子セットサイズ:", value = 500, min = 1),
        actionButton(ns("runGSEA"), "GSEA実行", icon = icon("play")),
        hr(),
        h4("結果ダウンロード"),
        downloadButton(ns("downloadGSEAResults"), "GSEA結果テーブルをダウンロード (.csv)", icon = icon("download"))
      ),
      mainPanel(
        width = 8,
        h4("GSEA 結果テーブル"),
        helpText("NES (正): Target群で上昇傾向, NES (負): Target群で下降傾向。padjが有意性の指標。テーブルの行を選択すると下にプロット表示。"),
        withSpinner(DTOutput(ns("gseaResultTable")), type = 6),
        hr(),
        h4("エンリッチメントプロット"),
        helpText("選択されたパスウェイについて、遺伝子ランキングに対するエンリッチメントスコアの推移を示します。"),
        withSpinner(plotOutput(ns("enrichmentPlot"), height = "500px"), type = 6)
      )
    )
  )
}

# Server関数
gseaServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    
    orgdb_pkg_map_gsea <- c(
      "Homo sapiens" = "org.Hs.eg.db",
      "Mus musculus" = "org.Mm.eg.db",
      "Rattus norvegicus" = "org.Rn.eg.db",
      "Danio rerio" = "org.Dr.eg.db",
      "Drosophila melanogaster" = "org.Dm.eg.db",
      "Caenorhabditis elegans" = "org.Ce.eg.db",
      "Saccharomyces cerevisiae" = "org.Sc.sgd.db",
      "Bos taurus" = "org.Bt.eg.db",
      "Gallus gallus" = "org.Gg.eg.db",
      "Canis lupus familiaris" = "org.Cf.eg.db",
      "Macaca mulatta" = "org.Mmu.eg.db",
      "Pan troglodytes" = "org.Pt.eg.db",
      "Sus scrofa" = "org.Ss.eg.db"
    )
    
    all_gene_sets_for_species <- reactive({
      req(input$gsea_species)
      message(paste("GSEA: Loading all gene sets for species:", input$gsea_species))
      tryCatch({
        msigdbr::msigdbr(species = input$gsea_species)
      }, error = function(e) {
        showNotification(paste("選択された生物種 (", input$gsea_species, ") の全遺伝子セット情報の取得に失敗しました:", e$message), type = "error", duration = 10)
        NULL
      })
    })
    
    observeEvent(all_gene_sets_for_species(), {
      all_gs <- all_gene_sets_for_species()
      if (!is.null(all_gs) && nrow(all_gs) > 0) {
        unique_gs_names <- sort(unique(all_gs$gs_name))
        updateSelectizeInput(session, "gsea_search_pathway_name",
                             choices = unique_gs_names,
                             selected = if (length(unique_gs_names) > 0) unique_gs_names[1] else NULL,
                             server = FALSE)
        message(paste("GSEA: Updated pathway search choices for", input$gsea_species, "-", length(unique_gs_names), "unique set names loaded."))
      } else {
        updateSelectizeInput(session, "gsea_search_pathway_name", choices = character(0), selected = NULL)
        message(paste("GSEA: No gene sets found or error loading for species:", input$gsea_species, "- pathway search choices cleared."))
      }
    })
    
    gsea_results_reactive <- eventReactive(input$runGSEA, {
      req(rv$deg_results, rv$deg_results$top_tags, rv$deg_results$top_tags$table)
      deg_table_original_entrez <- as.data.frame(rv$deg_results$top_tags$table)
      validate(
        need("Geneid" %in% colnames(deg_table_original_entrez), "DEG結果テーブルに 'Geneid' 列（EntrezIDを想定）が見つかりません。"),
        need(nrow(deg_table_original_entrez) > 0, "DEG結果テーブルが空です。")
      )
      req(input$gsea_species)
      message("GSEA: ID変換開始 (EntrezID -> Gene Symbol)...")
      entrez_ids <- deg_table_original_entrez$Geneid
      selected_gsea_species_name <- input$gsea_species
      orgdb_pkg_name <- orgdb_pkg_map_gsea[[selected_gsea_species_name]]
      validate(need(!is.null(orgdb_pkg_name) && nzchar(orgdb_pkg_name),
                    paste0("選択された生物種 '", selected_gsea_species_name, "' に対応するOrgDbパッケージが見つかりません。")))
      if (!requireNamespace(orgdb_pkg_name, quietly = TRUE)) {
        stop(paste0("GSEAのためのID変換には '", orgdb_pkg_name, "' パッケージが必要です。\n",
                    "インストールしてください: BiocManager::install('", orgdb_pkg_name, "')"))
      }
      orgdb_pkg <- get(orgdb_pkg_name)
      gene_symbols_map <- tryCatch({
        suppressMessages(AnnotationDbi::mapIds(orgdb_pkg,
                                               keys = unique(as.character(entrez_ids)),
                                               keytype = "ENTREZID",
                                               column = "SYMBOL",
                                               multiVals = "first"))
      }, error = function(e) {
        showNotification(paste("ID変換エラー(Entrez->Symbol):", e$message), type="error", duration=10)
        NULL
      })
      validate(need(!is.null(gene_symbols_map),
                    paste("EntrezIDからGene Symbolへの変換に失敗しました。")))
      deg_table_symbol <- deg_table_original_entrez
      deg_table_symbol$GeneSymbol <- gene_symbols_map[as.character(deg_table_symbol$Geneid)]
      n_before_filter_symbol <- nrow(deg_table_symbol)
      deg_table_symbol <- deg_table_symbol[!is.na(deg_table_symbol$GeneSymbol) & deg_table_symbol$GeneSymbol != "", ]
      n_na_removed_symbol <- n_before_filter_symbol - nrow(deg_table_symbol)
      if(n_na_removed_symbol > 0) {
        message(paste("GSEA情報:", n_na_removed_symbol, "個のEntrezIDはGene Symbolに変換できず除外。"))
      }
      validate(need(nrow(deg_table_symbol) > 0, "Gene SymbolへのID変換後、有効な遺伝子が残りませんでした。"))
      has_duplicates_symbol <- any(duplicated(deg_table_symbol$GeneSymbol))
      if(has_duplicates_symbol) {
        message("GSEA情報: 重複Gene Symbolあり。ランキング指標を平均化。")
      }
      metric_type <- input$gsea_rank_metric
      ranks_df <- NULL
      # ★★★ Serverの変更箇所 ★★★
      if (metric_type == "logFC") {
        validate(need("logFC" %in% colnames(deg_table_symbol), "logFC列なし"))
        ranks_df <- deg_table_symbol[, c("GeneSymbol", "logFC"), drop=FALSE]; colnames(ranks_df)[2] <- "rank_val"
      } else if (metric_type == "signed_p") {
        validate(need(all(c("PValue", "logFC") %in% colnames(deg_table_symbol)), "PValueまたはlogFC列なし"))
        pvals <- deg_table_symbol$PValue; pvals[pvals == 0] <- .Machine$double.eps
        ranks_df <- data.frame(GeneSymbol = deg_table_symbol$GeneSymbol, rank_val = -log10(pvals) * sign(deg_table_symbol$logFC))
      } else if (metric_type == "logFC_x_signed_p") { # 新しい指標の計算
        validate(need(all(c("PValue", "logFC") %in% colnames(deg_table_symbol)), "PValueまたはlogFC列なし"))
        pvals <- deg_table_symbol$PValue; pvals[pvals == 0] <- .Machine$double.eps
        ranks_df <- data.frame(GeneSymbol = deg_table_symbol$GeneSymbol, rank_val = deg_table_symbol$logFC * -log10(pvals))
      } else if (metric_type == "stat") {
        f_stat_col <- if ("F" %in% colnames(deg_table_symbol)) "F" else if ("LR" %in% colnames(deg_table_symbol)) "LR" else NULL
        validate(need(!is.null(f_stat_col), "検定統計量列 (FまたはLR) がDEG結果にありません。"))
        message(paste("GSEA: 検定統計量として", f_stat_col, "列を使用します。"))
        ranks_df <- data.frame(GeneSymbol = deg_table_symbol$GeneSymbol, rank_val = deg_table_symbol[[f_stat_col]] * sign(deg_table_symbol$logFC))
      }
      ranks_df <- na.omit(ranks_df)
      validate(need(nrow(ranks_df) > 0, "有効なランキング指標を持つ遺伝子なし(NA除去後)。"))
      if(has_duplicates_symbol) { ranks_agg <- aggregate(rank_val ~ GeneSymbol, data = ranks_df, FUN = mean) } else { ranks_agg <- ranks_df }
      ranks <- setNames(ranks_agg$rank_val, ranks_agg$GeneSymbol); ranks <- sort(ranks, decreasing = TRUE)
      validate(need(length(ranks) > 0, "最終ランキングメトリクス作成失敗。"))
      message(paste("GSEA: ランキングに使用する遺伝子数 (GeneSymbolベース, 重複集約/NA除去後):", length(ranks)))
      
      message("GSEA: 遺伝子セット準備中...")
      pathways <- list()
      
      if (input$gsea_selection_method == "category") {
        req(input$gsea_collection_category, input$gsea_species)
        gsea_collection_parts <- strsplit(input$gsea_collection_category, "_")[[1]]
        msig_cat <- gsea_collection_parts[1]
        msig_subcat <- if(length(gsea_collection_parts) > 1) paste(gsea_collection_parts[-1], collapse="_") else NA
        
        m_df_category <- tryCatch({
          msigdbr::msigdbr(species = input$gsea_species,
                           category = msig_cat,
                           subcategory = if (!is.na(msig_subcat) && nzchar(msig_subcat)) msig_subcat else NULL)
        }, error = function(e) { NULL })
        
        validate(need(!is.null(m_df_category) && nrow(m_df_category) > 0,
                      paste("カテゴリ指定で遺伝子セットを取得できませんでした:", input$gsea_collection_category)))
        validate(need("gene_symbol" %in% colnames(m_df_category), "msigdbr(category)結果に 'gene_symbol' 列なし。"))
        pathways <- split(x = m_df_category$gene_symbol, f = m_df_category$gs_name)
        
      } else if (input$gsea_selection_method == "search") {
        req(input$gsea_search_pathway_name, nzchar(input$gsea_search_pathway_name))
        selected_gs_name <- input$gsea_search_pathway_name
        
        all_gs_data <- all_gene_sets_for_species()
        validate(need(!is.null(all_gs_data) && nrow(all_gs_data) > 0,
                      "検索対象の全遺伝子セット情報がロードされていません。生物種を確認してください。"))
        
        target_pathway_genes <- all_gs_data[all_gs_data$gs_name == selected_gs_name, "gene_symbol", drop = TRUE]
        
        validate(need(length(target_pathway_genes) > 0,
                      paste("選択された遺伝子セット名 '", selected_gs_name, "' に対応する遺伝子が見つかりません。")))
        
        pathways[[selected_gs_name]] <- unique(target_pathway_genes)
      } else {
        validate("遺伝子セットの指定方法が選択されていません。")
      }
      
      validate(need(length(pathways) > 0, "解析対象のパスウェイがありませんでした。遺伝子セットの指定を確認してください。"))
      message(paste("GSEA: 使用パスウェイ数:", length(pathways)))
      if(length(pathways) == 1) message(paste("GSEA: 解析対象パスウェイ:", names(pathways)[1]))
      
      message("GSEA: fgsea実行中...")
      set.seed(42)
      fgsea_res <- tryCatch({
        fgsea::fgsea(pathways = pathways,
                     stats = ranks,
                     minSize = input$gsea_min_size,
                     maxSize = input$gsea_max_size,
                     eps = 0.0,
                     nproc = 0)
      }, error = function(e) {
        showNotification(paste("fgsea解析エラー:", e$message), type="error", duration=10)
        NULL
      })
      validate(need(!is.null(fgsea_res) && nrow(fgsea_res) > 0,
                    "GSEA解析結果が得られませんでした。パラメータや入力データ、遺伝子セットを確認してください。"))
      message("GSEA: fgsea完了")
      fgsea_res_ordered <- fgsea_res[order(padj, -abs(NES)), ]
      fgsea_res_ordered$leadingEdge <- sapply(fgsea_res_ordered$leadingEdge, paste, collapse=",")
      return(list( results_table = fgsea_res_ordered, ranks = ranks, pathways = pathways ))
    })
    
    output$gseaResultTable <- renderDT({
      gsea_output <- gsea_results_reactive()
      req(gsea_output, gsea_output$results_table)
      dt <- datatable(gsea_output$results_table,
                      rownames = FALSE,
                      selection = 'single',
                      filter = 'top',
                      extensions = 'Buttons',
                      options = list(pageLength = 10,
                                     scrollX = TRUE,
                                     searching = TRUE,
                                     search = list(regex = TRUE, caseInsensitive = TRUE),
                                     dom = 'Blfrtip',
                                     buttons = list(
                                       list(extend = 'copy', text = 'コピー'),
                                       list(extend = 'csv', text = 'CSV', filename = "gsea_results"),
                                       list(extend = 'excel', text = 'Excel', filename = "gsea_results"),
                                       list(extend = 'print', text = '印刷')
                                     ),
                                     order = list(list(which(colnames(gsea_output$results_table) == "padj") -1, 'asc'))
                      )
      ) %>%
        formatSignif(columns = c('pval', 'padj'), digits = 3) %>%
        formatRound(columns = c('NES', 'ES', 'log2err'), digits = 3)
      return(dt)
    })
    
    output$enrichmentPlot <- renderPlot({
      gsea_output <- gsea_results_reactive()
      selected_row_indices <- input$gseaResultTable_rows_selected
      req(gsea_output, gsea_output$results_table, gsea_output$ranks, gsea_output$pathways)
      validate(need(length(selected_row_indices) == 1, "結果テーブルから1行選択してください。"))
      selected_pathway_name <- gsea_output$results_table$pathway[selected_row_indices]
      validate(need(selected_pathway_name %in% names(gsea_output$pathways),
                    paste0("選択されたパスウェイ '", selected_pathway_name, "' が遺伝子セットリストに見つかりません。")))
      plot_obj <- fgsea::plotEnrichment(gsea_output$pathways[[selected_pathway_name]],
                                        gsea_output$ranks) +
        labs(title = selected_pathway_name) +
        theme_bw(base_size = 12)
      print(plot_obj)
    })
    
    output$downloadGSEAResults <- downloadHandler(
      filename = function() {
        comparison_name_dl <- "GSEA_comparison"
        if(!is.null(rv$deg_results) && !is.null(rv$deg_results$comparison) && length(rv$deg_results$comparison) == 2){
          comparison_name_dl <- paste0(rv$deg_results$comparison[1], "_vs_", rv$deg_results$comparison[2])
        }
        
        collection_info_dl <- ""
        if (input$gsea_selection_method == "category") {
          collection_info_dl <- gsub("[: /]", "_", input$gsea_collection_category)
        } else if (input$gsea_selection_method == "search") {
          collection_info_dl <- gsub("[^A-Za-z0-9_.-]", "_", input$gsea_search_pathway_name) 
          if (nchar(collection_info_dl) > 50) collection_info_dl <- substr(collection_info_dl, 1, 50)
        }
        
        species_name_dl <- gsub(" ", "_", input$gsea_species)
        paste0("GSEA_results_", comparison_name_dl, "_", species_name_dl, "_", collection_info_dl, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        gsea_output_dl <- gsea_results_reactive()
        if (is.null(gsea_output_dl) || is.null(gsea_output_dl$results_table)) {
          write.csv(data.frame(Message="GSEA結果がありません。"), file, row.names = FALSE, quote = TRUE)
          return()
        }
        write.csv(gsea_output_dl$results_table, file, row.names = FALSE, quote = TRUE)
      }
    )
  })
}