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
  "Xenopus_laevis" = "xla",
  "Lotus_japonicus" = "lja"
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
                         choices = c("ALL (GO + KEGG + Reactome)" = "ALL",
                                     "GO: BP" = "BP",
                                     "GO: MF" = "MF",
                                     "GO: CC" = "CC",
                                     "KEGG" = "KEGG",
                                     "Reactome" = "REACTOME"),
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
      ),
      tabPanel("Network Plot",
               helpText("遺伝子とエンリッチされたタームの関係をネットワーク図で表示します。"),
               uiOutput(ns("plot_selector_ui_net")), # プロットセレクター
               numericInput(ns("netplot_n"), "表示するターム数:", value = 5, min = 1, max = 20),
               withSpinner(plotOutput(ns("goNetPlot")), type = 6),
               downloadButton(ns("downloadNetPlot"), "Network Plotをダウンロード (.png)", icon = icon("download"))
      )
    )
  )
}

goEnrichmentIntegratedServer <- function(id, rv, deg_results_reactive, background_genes_reactive, selected_species_code_reactive, gene_annotation_reactive = reactive(NULL), gtf_id_choices_reactive = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    analysis_results <- reactiveVal(NULL)
    analysis_running <- reactiveVal(FALSE)

    # GTFアップロード時は、その中身に応じて表示する遺伝子IDタイプの選択肢・ラベルを追従させる
    observeEvent(gtf_id_choices_reactive(), {
      ch <- gtf_id_choices_reactive()
      if (is.null(ch)) {
        ch <- list(choices = c("Gene Symbol" = "SYMBOL", "Entrez ID (内部ID)" = "ENTREZID"), selected = "SYMBOL")
      }
      updateSelectInput(session, "go_id_display_type", choices = ch$choices, selected = ch$selected)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)
    
    # 共通ヘルパー annotate_display_ids に委譲 (GTFアノテーション優先 -> OrgDb -> 生ID)
    convert_entrez_ids_for_display <- function(entrez_ids, target_display_type, selected_species_code) {
      if (target_display_type == "ENTREZID" || is.null(target_display_type) || !nzchar(target_display_type)) {
        return(as.character(entrez_ids))
      }
      req(selected_species_code)
      annotate_display_ids(entrez_ids, target_display_type, selected_species_code,
                           gene_annotation = gene_annotation_reactive(), orgdb_map = orgdb_species_map)
    }
    
    observeEvent(deg_results_reactive(), {
      deg <- deg_results_reactive()
      req(deg)
      clusters <- deg$kmeans_clusters
      base_choices <- c("選択してください" = "", "Up-regulated" = "up", "Down-regulated" = "down", "All Significant (Up+Down)" = "all_sig")
      
      if (!is.null(clusters) && length(clusters) > 0) {
        cluster_choices <- setNames(names(clusters), names(clusters))
        new_choices <- c(base_choices, cluster_choices)
      } else {
        new_choices <- base_choices
      }
      
      curr_val <- isolate(input$gene_set)
      if (!(curr_val %in% new_choices)) {
         curr_val <- "up"
      }
      updateSelectInput(session, "gene_set", choices = new_choices, selected = curr_val)
    })

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
        if (input$gene_set %in% names(target_entrez_ids_from_deg$kmeans_clusters)) {
            target_genes_entrez_final <- target_entrez_ids_from_deg$kmeans_clusters[[input$gene_set]]
        }
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
          if (!is.null(rv$custom_annotations) && !is.null(rv$custom_annotations$go_t2g)) {
            message("  - Running enricher (Custom GO)")
            t2g <- rv$custom_annotations$go_t2g
            res_ont <- tryCatch(
              clusterProfiler::enricher(
                gene = target_genes_entrez_final,
                universe = universe_genes_entrez_final,
                TERM2GENE = t2g,
                pAdjustMethod = "BH",
                pvalueCutoff = input$pvalue_cutoff,
                qvalueCutoff = input$qvalue_cutoff
              ),
              error = function(e) { message("enricher (custom GO) error: ", e$message); NULL }
            )
            if (!is.null(res_ont) && nrow(as.data.frame(res_ont)) > 0) {
              all_results_list[["Custom_GO"]] <- res_ont
            }
          } else if (!is.null(rv$gene2go_annotation) || !is.null(resolve_bundled_gene2go(current_selected_species_code))) {
            # OrgDb不在の生物種: gene2go (NCBI gene2go形式: GeneID, GO, [Term], [Category]) ベースの
            # TERM2GENE で GO エンリッチメントを代替する。
            #   優先1: ユーザーがアップロードした rv$gene2go_annotation (任意の非モデル生物)
            #   優先2: アプリに同梱した種別 gene2go (bundled_gene2go_registry; 例 Lotus=lja)
            #         → ユーザーは種を選ぶだけ。アップロード不要。
            gene2go_df <- NULL
            if (!is.null(rv$gene2go_annotation) && is.data.frame(rv$gene2go_annotation) && nrow(rv$gene2go_annotation) > 0) {
              gene2go_df <- rv$gene2go_annotation
              message("  - GO: using uploaded gene2go annotation (", nrow(gene2go_df), " rows)")
            } else {
              bundled_path <- resolve_bundled_gene2go(current_selected_species_code)
              if (is.null(bundled_path)) {
                showNotification(paste0("同梱 gene2go が見つかりません (", current_selected_species_code,
                                        ")。GO解析をスキップします。"), type = "warning", duration = 8)
              } else {
                message("  - GO: using bundled gene2go for ", current_selected_species_code, ": ", bundled_path)
                gene2go_df <- read.table(gzfile(bundled_path), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
              }
            }
            if (!is.null(gene2go_df) && all(c("GeneID", "GO") %in% colnames(gene2go_df))) {
              gene2go_df$GeneID <- as.character(gene2go_df$GeneID)
              has_cat  <- "Category" %in% colnames(gene2go_df)
              has_term <- "Term" %in% colnames(gene2go_df)
              # Category列があれば NCBI形式 (Process/Function/Component) で ontology 別に集計、
              # 無ければ ontology で分けず一括 (1回) で解析する。
              ont_loop <- if (has_cat) go_ontologies else "ALL_GO"
              for (ont_cat in ont_loop) {
                if (has_cat) {
                  ont_label <- switch(ont_cat, "BP" = "Process", "MF" = "Function", "CC" = "Component", ont_cat)
                  sub_df <- gene2go_df[gene2go_df$Category == ont_label, , drop = FALSE]
                } else {
                  sub_df <- gene2go_df
                }
                t2g <- sub_df[, c("GO", "GeneID")]
                colnames(t2g) <- c("term", "gene")
                message("  - Running enricher (non-model GO) for: ", ont_cat, " (", nrow(t2g), " annotations)")
                res_ont <- tryCatch(
                  clusterProfiler::enricher(
                    gene = target_genes_entrez_final,
                    universe = universe_genes_entrez_final,
                    TERM2GENE = t2g,
                    pAdjustMethod = "BH",
                    pvalueCutoff = input$pvalue_cutoff,
                    qvalueCutoff = input$qvalue_cutoff
                  ),
                  error = function(e) { message("enricher error: ", e$message); NULL }
                )
                if (!is.null(res_ont) && nrow(as.data.frame(res_ont)) > 0) {
                  # Descriptionをterm名(Term列)で補完
                  if (has_term) {
                    t2d <- unique(gene2go_df[, c("GO", "Term")])
                    idx <- match(res_ont@result$ID, t2d$GO)
                    res_ont@result$Description <- ifelse(is.na(idx), res_ont@result$ID, t2d$Term[idx])
                  }
                  res_key <- if (has_cat) ont_cat else "GO"
                  all_results_list[[res_key]] <- res_ont
                }
              }
            }
          } else {
            orgdb_pkg_name <- orgdb_species_map[[current_selected_species_code]]
            if (is.null(orgdb_pkg_name) || !nzchar(orgdb_pkg_name) || !requireNamespace(orgdb_pkg_name, quietly = TRUE)) {
              stop(paste0("GO解析に必要なOrgDbパッケージ '", orgdb_pkg_name,
                          "' が見つかりません。非モデル生物の場合は、データアップロードタブで ",
                          "gene2go アノテーション (GeneID, GO, Category 列) または eggNOG-mapper アノテーションを",
                          "アップロードしてください。"))
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
        }
        
        # KEGG Analysis
        if (input$analysis_type == "ALL" || input$analysis_type == "KEGG") {
          if (!is.null(rv$custom_annotations) && (!is.null(rv$custom_annotations$ko_t2g) || !is.null(rv$custom_annotations$pathway_t2g))) {
            if (!is.null(rv$custom_annotations$pathway_t2g)) {
              t2g <- rv$custom_annotations$pathway_t2g
              message("  - Running enricher (Custom KEGG Pathway)")
            } else {
              t2g <- rv$custom_annotations$ko_t2g
              message("  - Running enricher (Custom KEGG KO)")
            }
            res_kegg <- tryCatch(
              clusterProfiler::enricher(
                gene = target_genes_entrez_final,
                universe = universe_genes_entrez_final,
                TERM2GENE = t2g,
                pAdjustMethod = "BH",
                pvalueCutoff = input$pvalue_cutoff,
                qvalueCutoff = input$qvalue_cutoff
              ),
              error = function(e) { message("enricher (custom KEGG) error: ", e$message); NULL }
            )
            if (!is.null(res_kegg) && nrow(as.data.frame(res_kegg)) > 0) {
              all_results_list[["Custom_KEGG"]] <- res_kegg
            }
          } else {
            # KEGG生物種コード: ユーザー指定(rv$kegg_organism_code, 例 lja)を優先、無ければ内蔵マップ
            kegg_code <- rv$kegg_organism_code
            if (is.null(kegg_code) || !nzchar(kegg_code)) kegg_code <- kegg_species_map[[current_selected_species_code]]
            if (is.null(kegg_code) || !nzchar(kegg_code)) {
              showNotification(paste0("KEGG解析はスキップされました: '", current_selected_species_code, "' に対応するKEGGコードがありません。データアップロードタブでKEGG生物種コード(例: lja)を指定できます。"), type = "warning", duration = 8)
            } else {
              message("  - Running enrichKEGG for: ", kegg_code)
              # enrichKEGG はオンラインKEGG取得に失敗すると "長さ0の変数名" 等で落ちるため、
              # tryCatch で隔離し失敗してもGO等の結果を巻き込まないようにする。
              res_kegg <- tryCatch(
                enrichKEGG(gene = target_genes_entrez_final, universe = universe_genes_entrez_final,
                           organism = kegg_code, pvalueCutoff = input$pvalue_cutoff, qvalueCutoff = input$qvalue_cutoff),
                error = function(e) {
                  message("enrichKEGG error: ", e$message)
                  showNotification(paste0("KEGG解析をスキップしました (", conditionMessage(e), ")"), type = "warning", duration = 8)
                  NULL
                }
              )
              if (!is.null(res_kegg) && nrow(as.data.frame(res_kegg)) > 0) {
                all_results_list[["KEGG"]] <- res_kegg
              }
            }
          }
        }
        
        # Reactome Analysis
        if (input$analysis_type == "ALL" || input$analysis_type == "REACTOME") {
          msig_species_name <- switch(current_selected_species_code,
            "Homo_sapiens" = "Homo sapiens",
            "Mus_musculus" = "Mus musculus",
            "Rattus_norvegicus" = "Rattus norvegicus",
            "Drosophila_melanogaster" = "Drosophila melanogaster",
            "Caenorhabditis_elegans" = "Caenorhabditis elegans",
            "Danio_rerio" = "Danio rerio",
            "Saccharomyces_cerevisiae" = "Saccharomyces cerevisiae",
            "Bos_taurus" = "Bos taurus",
            "Gallus_gallus" = "Gallus gallus",
            "Canis_familiaris" = "Canis lupus familiaris",
            "Macaca_mulatta" = "Macaca mulatta",
            "Pan_troglodytes" = "Pan troglodytes",
            "Sus_scrofa" = "Sus scrofa",
            NULL
          )
          if (!is.null(msig_species_name)) {
            message("  - Running enrichReactome (via msigdbr) for: ", msig_species_name)
            reactome_df <- tryCatch({
              msigdbr::msigdbr(species = msig_species_name, category = "C2", subcategory = "CP:REACTOME")
            }, error = function(e) {
              message("Reactome load error: ", e$message)
              NULL
            })
            if (!is.null(reactome_df) && nrow(reactome_df) > 0) {
              m_t2g <- reactome_df[, c("gs_name", "entrez_gene")]
              res_reactome <- tryCatch({
                clusterProfiler::enricher(
                  gene = target_genes_entrez_final,
                  universe = universe_genes_entrez_final,
                  TERM2GENE = m_t2g,
                  pvalueCutoff = input$pvalue_cutoff,
                  qvalueCutoff = input$qvalue_cutoff
                )
              }, error = function(e) {
                message("Reactome enricher error: ", e$message)
                NULL
              })
              if (!is.null(res_reactome) && nrow(as.data.frame(res_reactome)) > 0) {
                all_results_list[["REACTOME"]] <- res_reactome
              }
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
        return(datatable(data.frame(Message = "解析中..."), style = "bootstrap5", options = list(dom = 't', searching = FALSE)))
      }

      final_table <- final_table_reactive()
      if (is.null(final_table)) {
        return(datatable(data.frame(Message = "解析を実行してください。"), style = "bootstrap5", options = list(dom = 't', searching = FALSE)))
      }
      if("Message" %in% colnames(final_table)){
        return(datatable(final_table, style = "bootstrap5", options = list(dom = 't', searching = FALSE)))
      }

      res_df_display <- final_table %>%
        mutate(across(where(is.numeric), ~ round(.x, digits = 4)))

      datatable(res_df_display, rownames = FALSE, style = "bootstrap5", class = "table-hover table-sm", filter = "top", extensions = 'Buttons', options = list(pageLength=10, scrollX=TRUE), escape = FALSE)
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
    output$plot_selector_ui_net <- renderUI({
      results_list <- analysis_results()
      req(input$analysis_type == "ALL", results_list, length(results_list) > 1)
      selectInput(ns("plot_source_net"), "プロット対象:", choices = names(results_list), selected = names(results_list)[1])
    })
    
    .apply_set_readable <- function(result_obj, species_code) {
      if (species_code == "Lotus_japonicus") return(result_obj)
      orgdb_pkg_name <- orgdb_species_map[[species_code]]
      if (!is.null(orgdb_pkg_name) && requireNamespace(orgdb_pkg_name, quietly = TRUE)) {
        orgDb_obj <- get(orgdb_pkg_name)
        result_obj <- setReadable(result_obj, OrgDb = orgDb_obj, keyType = "ENTREZID")
      }
      return(result_obj)
    }

    get_plot_data_object_bar <- reactive({
      results_list <- analysis_results()
      req(results_list)
      source_key <- if (input$analysis_type == "ALL") { req(input$plot_source_bar); input$plot_source_bar } else { input$analysis_type }
      result_obj <- results_list[[source_key]]
      req(result_obj, inherits(result_obj, "enrichResult"))
      .apply_set_readable(result_obj, selected_species_code_reactive())
    })

    get_plot_data_object_dot <- reactive({
      results_list <- analysis_results()
      req(results_list)
      source_key <- if (input$analysis_type == "ALL") { req(input$plot_source_dot); input$plot_source_dot } else { input$analysis_type }
      result_obj <- results_list[[source_key]]
      req(result_obj, inherits(result_obj, "enrichResult"))
      .apply_set_readable(result_obj, selected_species_code_reactive())
    })

    get_plot_data_object_net <- reactive({
      results_list <- analysis_results()
      req(results_list)
      source_key <- if (input$analysis_type == "ALL") { req(input$plot_source_net); input$plot_source_net } else { input$analysis_type }
      result_obj <- results_list[[source_key]]
      req(result_obj, inherits(result_obj, "enrichResult"))
      .apply_set_readable(result_obj, selected_species_code_reactive())
    })
    
    bar_plot_gg_object <- reactive({
      result_obj_for_plot <- get_plot_data_object_bar()
      req(result_obj_for_plot)
      shiny::validate(shiny::need(nrow(as.data.frame(result_obj_for_plot)) > 0, "プロットする有意なタームがありません。"))
      n_terms <- min(input$barplot_n, nrow(as.data.frame(result_obj_for_plot)))
      barplot(result_obj_for_plot, showCategory = n_terms)
    })
    
    output$goBarPlot <- renderPlot({
      print(bar_plot_gg_object())
    })
    
    dot_plot_gg_object <- reactive({
      result_obj_for_plot <- get_plot_data_object_dot()
      req(result_obj_for_plot)
      shiny::validate(shiny::need(nrow(as.data.frame(result_obj_for_plot)) > 0, "プロットする有意なタームがありません。"))
      n_terms <- min(input$dotplot_n, nrow(as.data.frame(result_obj_for_plot)))
      dotplot(result_obj_for_plot, showCategory = n_terms)
    })
    
    output$goDotPlot <- renderPlot({
      print(dot_plot_gg_object())
    })
    
    net_plot_gg_object <- reactive({
      result_obj_for_plot <- get_plot_data_object_net()
      req(result_obj_for_plot)
      shiny::validate(shiny::need(nrow(as.data.frame(result_obj_for_plot)) > 0, "プロットする有意なタームがありません。"))
      n_terms <- min(input$netplot_n, nrow(as.data.frame(result_obj_for_plot)))
      enrichplot::cnetplot(result_obj_for_plot, showCategory = n_terms)
    })
    
    output$goNetPlot <- renderPlot({
      print(net_plot_gg_object())
    })
    
    output$downloadNetPlot <- downloadHandler(
      filename = function() {
        paste0("NetworkPlot_", Sys.Date(), ".png")
      },
      content = function(file) {
        ggsave(file, plot = net_plot_gg_object(), device = "png", width = 10, height = 8, dpi = 300)
      }
    )
    
    # エンリッチメント解析パラメータのメタデータを生成
    enrichment_metadata_df <- reactive({
      deg_res <- deg_results_reactive()
      comparison_str <- if (!is.null(deg_res) && !is.null(deg_res$comparison)) {
        paste(deg_res$comparison, collapse = " vs ")
      } else { "N/A" }

      params <- data.frame(
        Parameter = c(
          "Analysis Date",
          "Analysis Type",
          "Gene Set",
          "P-value Cutoff",
          "Q-value Cutoff",
          "Species",
          "DEG Comparison",
          "DEG Method",
          "DEG Significance Metric",
          "DEG Significance Threshold",
          "DEG LogFC Threshold",
          "Gene ID Display Type"
        ),
        Value = c(
          as.character(Sys.time()),
          input$analysis_type %||% "N/A",
          input$gene_set %||% "N/A",
          as.character(input$pvalue_cutoff %||% "N/A"),
          as.character(input$qvalue_cutoff %||% "N/A"),
          selected_species_code_reactive() %||% "N/A",
          comparison_str,
          deg_res$analysis_method %||% "N/A",
          deg_res$metric_at_run %||% "N/A",
          as.character(deg_res$sig_threshold_at_run %||% "N/A"),
          as.character(deg_res$logfc_threshold_at_run %||% "N/A"),
          input$go_id_display_type %||% "N/A"
        ),
        stringsAsFactors = FALSE
      )
      params
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
          meta <- enrichment_metadata_df()
          comment_lines <- paste0("# ", meta$Parameter, ": ", meta$Value)
          writeLines(c(comment_lines, ""), con = file, useBytes = TRUE)
          write.table(table_to_download, file, row.names = FALSE, quote = TRUE,
                       sep = ",", append = TRUE, fileEncoding = "UTF-8")
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
          sheets <- list(
            "Enrichment_Results" = table_to_download,
            "Analysis_Parameters" = enrichment_metadata_df()
          )
          writexl::write_xlsx(sheets, file)
        }
      }
    )
    return(analysis_results)
  }) 
}
