#' Enrichment Module — GSEA (fgsea) + GO/KEGG (clusterProfiler)
#' @import shiny
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box input_task_button
#' @importFrom bsicons bs_icon
#' @noRd
mod_enrichment_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 350,
      # ── Contrast selector ────────────────────────────────────────────────
      card(
        card_header(bsicons::bs_icon("diagram-3"), " Settings"),
        card_body(
          selectInput(ns("contrast"), "DEG contrast", choices = NULL),
          selectInput(ns("organism"), "Organism",
                      choices = c("Human (hsa)" = "human",
                                  "Mouse (mmu)" = "mouse")),
          radioButtons(ns("analysis_type"), "Analysis type",
                       choices = c("GSEA (ranked gene list)" = "gsea",
                                   "ORA (DEG gene set)"      = "ora"),
                       selected = "gsea")
        )
      ),
      # ── GSEA settings ────────────────────────────────────────────────────
      conditionalPanel(
        condition = paste0("input['", ns("analysis_type"), "'] == 'gsea'"),
        card(
          card_header("GSEA Options"),
          card_body(
            selectInput(ns("gene_set_db"), "Gene set database",
                        choices = c("Hallmark (MSigDB)"      = "H",
                                    "GO Biological Process"  = "GO_BP",
                                    "GO Molecular Function"  = "GO_MF",
                                    "KEGG Pathway"           = "KEGG",
                                    "Reactome"               = "C2_CP_REACTOME")),
            selectInput(ns("rank_by"), "Rank metric",
                        choices = c("stat (Wald statistic)" = "stat",
                                    "-log10(padj) × sign(LFC)" = "signed_log10p")),
            numericInput(ns("n_perm"), "Permutations", 1000, min = 100, step = 500),
            numericInput(ns("gsea_padj"), "FDR cutoff", 0.25, min = 0, max = 1, step = 0.05)
          )
        )
      ),
      # ── ORA settings ─────────────────────────────────────────────────────
      conditionalPanel(
        condition = paste0("input['", ns("analysis_type"), "'] == 'ora'"),
        card(
          card_header("ORA Options"),
          card_body(
            selectInput(ns("ora_db"), "Database",
                        choices = c("GO Biological Process" = "BP",
                                    "GO Molecular Function" = "MF",
                                    "KEGG Pathway"          = "KEGG")),
            numericInput(ns("ora_padj_deg"),  "DEG FDR cutoff",  0.05, min=0, max=1, step=0.01),
            numericInput(ns("ora_lfc_deg"),   "DEG LFC cutoff",  1,    min=0, step=0.5),
            numericInput(ns("ora_padj_enr"),  "Enrichment FDR",  0.05, min=0, max=1, step=0.05)
          )
        )
      ),
      input_task_button(ns("run_enr"), "Run Analysis", icon = icon("play"),
                        class = "btn-primary w-100 mt-2")
    ),
    uiOutput(ns("main_content"))
  )
}

#' @noRd
mod_enrichment_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Populate contrast selector ────────────────────────────────────────────
    observe({
      results <- state$deg_results()
      if (length(results) > 0)
        updateSelectInput(session, "contrast", choices = names(results))
    })

    # ── OrgDb helper ─────────────────────────────────────────────────────────
    .orgdb <- function(org) {
      pkg <- if (org == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
      if (!requireNamespace(pkg, quietly = TRUE))
        stop(sprintf("%s not installed. Run: BiocManager::install('%s')", pkg, pkg))
      getExportedValue(pkg, pkg)
    }

    # ── ExtendedTask ──────────────────────────────────────────────────────────
    enr_task <- ExtendedTask$new(function(deg_df, contrast, analysis_type, organism,
                                          gene_set_db, rank_by, n_perm, gsea_padj,
                                          ora_db, ora_padj_deg, ora_lfc_deg, ora_padj_enr) {
      promises::future_promise({
        if (!requireNamespace("fgsea", quietly = TRUE))
          stop("fgsea not installed.")
        if (!requireNamespace("clusterProfiler", quietly = TRUE))
          stop("clusterProfiler not installed.")

        orgdb <- if (organism == "human") {
          if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
            stop("org.Hs.eg.db not installed.")
          org.Hs.eg.db::org.Hs.eg.db
        } else {
          if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
            stop("org.Mm.eg.db not installed.")
          org.Mm.eg.db::org.Mm.eg.db
        }

        gene_col <- if ("gene" %in% colnames(deg_df)) "gene" else "rownames"
        genes    <- if (gene_col == "gene") deg_df$gene else rownames(deg_df)

        if (analysis_type == "gsea") {
          # ── GSEA via fgsea ────────────────────────────────────────────────
          # Build ranked vector
          rank_vec <- if (rank_by == "stat" && "stat" %in% colnames(deg_df)) {
            setNames(deg_df$stat, genes)
          } else {
            sign_lfc <- sign(deg_df$log2FoldChange)
            sign_lfc[is.na(sign_lfc)] <- 0
            log10p <- -log10(pmax(deg_df$pvalue, 1e-300))
            log10p[is.na(log10p)] <- 0
            setNames(sign_lfc * log10p, genes)
          }
          rank_vec <- sort(rank_vec[!is.na(rank_vec)], decreasing = TRUE)

          # Get gene sets
          if (gene_set_db %in% c("H","C2_CP_REACTOME")) {
            if (!requireNamespace("msigdbr", quietly = TRUE))
              stop("msigdbr not installed. Run: install.packages('msigdbr')")
            category <- if (gene_set_db == "H") "H" else "C2"
            subcategory <- if (gene_set_db == "C2_CP_REACTOME") "CP:REACTOME" else NULL
            species_msig <- if (organism == "human") "Homo sapiens" else "Mus musculus"
            m_df <- msigdbr::msigdbr(species = species_msig,
                                      category = category,
                                      subcategory = subcategory)
            pathways <- split(m_df$gene_symbol, m_df$gs_name)
          } else if (gene_set_db == "KEGG") {
            kegg_org <- if (organism == "human") "hsa" else "mmu"
            # Use clusterProfiler KEGG data via enrichKEGG path list
            kegg_data <- clusterProfiler::download_KEGG(kegg_org)
            pathways  <- split(kegg_data$KEGGPATHID2EXTID$to,
                                kegg_data$KEGGPATHID2EXTID$from)
            # Convert ENTREZ to SYMBOL for ranking
            sym2eg <- clusterProfiler::bitr(names(rank_vec), fromType="SYMBOL",
                                            toType="ENTREZID", OrgDb=orgdb)
            rank_vec <- rank_vec[sym2eg$SYMBOL]
            names(rank_vec) <- sym2eg$ENTREZID
            rank_vec <- sort(rank_vec, decreasing = TRUE)
          } else {
            # GO_BP / GO_MF via clusterProfiler
            go_ont <- gsub("GO_", "", gene_set_db)
            go_df  <- clusterProfiler::bitr_kegg(NULL)  # just to load pkg
            # Build pathways from GO annotation
            sym2eg <- clusterProfiler::bitr(names(rank_vec), fromType="SYMBOL",
                                            toType="ENTREZID", OrgDb=orgdb)
            rank_eg <- rank_vec[sym2eg$SYMBOL]
            names(rank_eg) <- sym2eg$ENTREZID
            rank_eg <- sort(rank_eg, decreasing=TRUE)

            gse_res <- clusterProfiler::gseGO(
              geneList  = rank_eg,
              OrgDb     = orgdb,
              ont       = go_ont,
              nPerm     = n_perm,
              pvalueCutoff = gsea_padj,
              verbose   = FALSE
            )
            return(list(type="gsea_go", result=gse_res, db=gene_set_db))
          }

          res <- fgsea::fgsea(pathways = pathways,
                              stats    = rank_vec,
                              nperm    = n_perm)
          res <- res[res$padj < gsea_padj & !is.na(res$padj), ]
          res <- res[order(res$padj), ]
          return(list(type="gsea_fgsea", result=as.data.frame(res), db=gene_set_db))

        } else {
          # ── ORA via clusterProfiler ───────────────────────────────────────
          sig_genes <- genes[!is.na(deg_df$padj) &
                               deg_df$padj  < ora_padj_deg &
                               abs(deg_df$log2FoldChange) > ora_lfc_deg]
          if (length(sig_genes) < 5)
            stop("Too few DEGs (< 5) to run ORA. Adjust thresholds.")

          eg <- clusterProfiler::bitr(sig_genes, fromType="SYMBOL",
                                      toType="ENTREZID", OrgDb=orgdb,
                                      drop=TRUE)$ENTREZID

          if (ora_db == "KEGG") {
            kegg_org <- if (organism == "human") "hsa" else "mmu"
            res <- clusterProfiler::enrichKEGG(gene = eg, organism = kegg_org,
                                               pvalueCutoff = ora_padj_enr)
          } else {
            res <- clusterProfiler::enrichGO(gene = eg, OrgDb = orgdb,
                                             ont  = ora_db,
                                             pAdjustMethod = "BH",
                                             pvalueCutoff  = ora_padj_enr,
                                             readable = TRUE)
          }
          return(list(type="ora", result=res, db=ora_db))
        }
      })
    }) |> bslib::bind_task_button("run_enr")

    observeEvent(input$run_enr, {
      req(length(state$deg_results()) > 0, input$contrast)
      state$set_status("enrichment", "running")
      deg_df <- state$deg_results()[[input$contrast]]
      enr_task$invoke(
        deg_df       = deg_df,
        contrast     = input$contrast,
        analysis_type = input$analysis_type,
        organism     = input$organism,
        gene_set_db  = input$gene_set_db,
        rank_by      = input$rank_by,
        n_perm       = input$n_perm,
        gsea_padj    = input$gsea_padj,
        ora_db       = input$ora_db,
        ora_padj_deg = input$ora_padj_deg,
        ora_lfc_deg  = input$ora_lfc_deg,
        ora_padj_enr = input$ora_padj_enr
      )
    })

    enr_results <- reactiveVal(NULL)

    observe({
      req(enr_task$status() == "success")
      res <- enr_task$result()
      enr_results(res)
      state$log_param("enrichment", "contrast",      isolate(input$contrast))
      state$log_param("enrichment", "type",          res$type)
      state$log_param("enrichment", "db",            res$db)
      n_sig <- if (res$type == "gsea_fgsea") nrow(res$result)
               else if (res$type == "gsea_go") nrow(as.data.frame(res$result))
               else nrow(as.data.frame(res$result))
      state$log_param("enrichment", "n_sig_pathways", n_sig)
      state$set_status("enrichment", "done")
    })

    observe({
      req(enr_task$status() == "error")
      showNotification(paste("Enrichment error:", enr_task$result()$message),
                       type = "error", duration = 10)
      state$set_status("enrichment", "pending")
    })

    # ── Dot plot (clusterProfiler / fgsea) ────────────────────────────────────
    dot_gg <- reactive({
      req(enr_results())
      res <- enr_results()

      if (res$type %in% c("gsea_go", "ora")) {
        # enrichplot dotplot
        if (!requireNamespace("enrichplot", quietly = TRUE))
          return(.placeholder_plot("enrichplot not installed."))
        enrichplot::dotplot(res$result, showCategory = 20) +
          ggplot2::scale_color_distiller(palette = "RdBu", direction = -1) +
          ggplot2::theme_minimal(base_size = 11)

      } else {
        # fgsea barplot
        df <- res$result
        req(nrow(df) > 0)
        df <- head(df[order(df$padj), ], 20)
        df$pathway <- factor(df$pathway, levels = rev(df$pathway))
        ggplot2::ggplot(df, ggplot2::aes(NES, pathway, fill = padj)) +
          ggplot2::geom_col() +
          ggplot2::scale_fill_distiller(palette = "Blues", direction = -1,
                                        name = "FDR") +
          ggplot2::labs(x = "NES", y = NULL,
                        title = paste("GSEA:", res$db)) +
          ggplot2::theme_minimal(base_size = 11)
      }
    }) |> bindCache(enr_results())

    # ── UI ────────────────────────────────────────────────────────────────────
    output$main_content <- renderUI({
      if (!state$step_ready("deg")) {
        return(ui_empty_state(
          "lock", "Enrichment Locked",
          "Run DEG analysis first.",
          actionButton("jump_deg", "Go to DEG",
                       onclick = "Shiny.setInputValue('main_nav','deg')")
        ))
      }
      res <- enr_results()
      tagList(
        if (!is.null(res)) {
          n_sig <- if (res$type == "gsea_fgsea") nrow(res$result)
                   else nrow(as.data.frame(res$result))
          fluidRow(
            column(4, value_box("Sig. Pathways", n_sig,
                                bsicons::bs_icon("diagram-3"), theme = "primary")),
            column(4, value_box("Database", res$db,
                                bsicons::bs_icon("database"), theme = "secondary")),
            column(4, value_box("Analysis", toupper(res$type),
                                bsicons::bs_icon("bar-chart"), theme = "light"))
          )
        },
        card(
          card_header(
            "Enrichment Plot",
            actionButton(ns("open_export_enr"), bsicons::bs_icon("gear"),
                         class = "btn btn-sm btn-outline-secondary float-end",
                         title = "Export options")
          ),
          card_body(plotOutput(ns("dot_plot"), height = "500px"))
        ),
        card(
          card_header(
            "Result Table",
            downloadButton(ns("dl_enr_csv"), "CSV",
                           class = "btn btn-sm btn-outline-secondary float-end")
          ),
          card_body(DT::DTOutput(ns("enr_table")))
        )
      )
    })

    output$dot_plot <- renderPlot({ req(dot_gg()); dot_gg() })

    output$enr_table <- DT::renderDT({
      req(enr_results())
      res <- enr_results()
      df <- if (res$type == "gsea_fgsea") res$result
            else as.data.frame(res$result)
      # Keep key columns
      keep_cols <- intersect(
        c("pathway","Description","NES","padj","pvalue","size","ES",
          "GeneRatio","BgRatio","Count","leadingEdge"),
        colnames(df)
      )
      df_show <- df[, keep_cols, drop = FALSE]
      # Truncate leadingEdge list column if present
      if ("leadingEdge" %in% colnames(df_show))
        df_show$leadingEdge <- sapply(df_show$leadingEdge,
                                      function(x) paste(head(x, 5), collapse=", "))
      DT::datatable(df_show, filter = "top", rownames = FALSE,
                    options = list(scrollX = TRUE, pageLength = 15)) |>
        DT::formatSignif(intersect(c("padj","pvalue","NES"), colnames(df_show)), 3)
    }, server = TRUE)

    output$dl_enr_csv <- downloadHandler(
      filename = function() paste0("enrichment_", isolate(input$contrast), ".csv"),
      content  = function(file) {
        req(enr_results())
        df <- if (enr_results()$type == "gsea_fgsea") enr_results()$result
              else as.data.frame(enr_results()$result)
        if ("leadingEdge" %in% colnames(df))
          df$leadingEdge <- sapply(df$leadingEdge, paste, collapse=",")
        write.csv(df, file, row.names = FALSE)
      }
    )

    observeEvent(input$open_export_enr, {
      req(dot_gg())
      open_export_modal(session, plot_fn = function() dot_gg(),
                        base_name = paste0("enrichment_", isolate(input$contrast)))
    })
  })
}

# ── Helper ────────────────────────────────────────────────────────────────────
.placeholder_plot <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x=0.5, y=0.5, label=msg, size=5, color="#666") +
    ggplot2::theme_void()
}
