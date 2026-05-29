#' GSVA Module — Single-sample pathway scoring
#' @import shiny
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box input_task_button
#' @importFrom bsicons bs_icon
#' @noRd
mod_gsva_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 350,
      card(
        card_header(bsicons::bs_icon("diagram-2"), " GSVA Settings"),
        card_body(
          selectInput(ns("organism"), "Organism",
                      choices = c("Human" = "human", "Mouse" = "mouse")),
          selectInput(ns("db"), "Gene set database",
                      choices = c("Hallmark (MSigDB)"     = "H",
                                  "KEGG Pathway"          = "C2_CP_KEGG",
                                  "GO Biological Process" = "C5_GO_BP")),
          numericInput(ns("top_n"), "Top pathways in heatmap", 20, min = 5, max = 50),
          input_task_button(ns("run_gsva"), "Run GSVA", icon = icon("play"),
                            class = "btn-primary w-100 mt-2")
        )
      )
    ),
    uiOutput(ns("main_content"))
  )
}

#' @noRd
mod_gsva_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns       <- session$ns
    gsva_res <- reactiveVal(NULL)

    gsva_task <- ExtendedTask$new(function(mat, organism, db) {
      promises::future_promise({
        if (!requireNamespace("GSVA",    quietly = TRUE)) stop("GSVA not installed. Run: BiocManager::install('GSVA')")
        if (!requireNamespace("msigdbr", quietly = TRUE)) stop("msigdbr not installed.")

        species <- if (organism == "human") "Homo sapiens" else "Mus musculus"
        cat_sub <- switch(db,
          "H"          = list(cat = "H",  sub = NULL),
          "C2_CP_KEGG" = list(cat = "C2", sub = "CP:KEGG"),
          "C5_GO_BP"   = list(cat = "C5", sub = "GO:BP")
        )
        m_df     <- msigdbr::msigdbr(species = species, category = cat_sub$cat, subcategory = cat_sub$sub)
        pathways <- split(m_df$gene_symbol, m_df$gs_name)
        pathways <- lapply(pathways, function(g) intersect(g, rownames(mat)))
        pathways <- Filter(function(g) length(g) >= 10, pathways)
        if (length(pathways) == 0) stop("No gene sets with ≥10 overlapping genes.")

        param <- GSVA::gsvaParam(exprData = mat, geneSets = pathways)
        GSVA::gsva(param, verbose = FALSE)
      })
    }) |> bslib::bind_task_button("run_gsva")

    observeEvent(input$run_gsva, {
      mat <- isolate(state$counts_filtered() %||% state$counts())
      req(mat)
      state$set_status("gsva", "running")
      gsva_task$invoke(mat, isolate(input$organism), isolate(input$db))
    })

    observe({
      req(gsva_task$status() == "success")
      scores <- gsva_task$result()
      gsva_res(scores)
      state$log_param("gsva", "db",       isolate(input$db))
      state$log_param("gsva", "pathways", nrow(scores))
      state$set_status("gsva", "done")
    })

    observe({
      req(gsva_task$status() == "error")
      showNotification(paste("GSVA error:", gsva_task$result()$message), type = "error", duration = 10)
      state$set_status("gsva", "pending")
    })

    heatmap_gg <- reactive({
      req(gsva_res())
      scores <- gsva_res()
      n_top  <- min(input$top_n, nrow(scores))
      vars   <- apply(scores, 1, var)
      top_pw <- names(sort(vars, decreasing = TRUE)[seq_len(n_top)])
      mat_t  <- scores[top_pw, , drop = FALSE]

      df <- data.frame(
        Pathway = rep(rownames(mat_t), times = ncol(mat_t)),
        Sample  = rep(colnames(mat_t), each  = nrow(mat_t)),
        Score   = as.vector(mat_t)
      )
      df$Pathway <- factor(df$Pathway, levels = rev(top_pw))

      ggplot2::ggplot(df, ggplot2::aes(Sample, Pathway, fill = Score)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.3) +
        ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                                      midpoint = 0, name = "GSVA\nScore") +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }) |> bindCache(gsva_res(), input$top_n)

    output$main_content <- renderUI({
      if (!state$step_ready("upload")) {
        return(ui_empty_state("lock", "GSVA Locked", "Upload data first."))
      }
      res <- gsva_res()
      tagList(
        if (!is.null(res)) fluidRow(
          column(4, value_box("Pathways", nrow(res), bsicons::bs_icon("diagram-2"), theme = "primary")),
          column(4, value_box("Samples",  ncol(res), bsicons::bs_icon("people"),    theme = "secondary")),
          column(4, value_box("Shown",    input$top_n, bsicons::bs_icon("filter"),  theme = "light"))
        ),
        card(
          card_header(
            "GSVA Score Heatmap",
            actionButton(ns("open_export"), bsicons::bs_icon("gear"),
                         class = "btn btn-sm btn-outline-secondary float-end",
                         title = "Export"),
            downloadButton(ns("dl_scores"), "Scores CSV",
                           class = "btn btn-sm btn-outline-secondary float-end me-1")
          ),
          card_body(plotOutput(ns("heatmap"), height = "500px"))
        )
      )
    })

    output$heatmap <- renderPlot({ req(heatmap_gg()); heatmap_gg() })

    output$dl_scores <- downloadHandler(
      filename = "GSVA_scores.csv",
      content  = function(file) { req(gsva_res()); write.csv(gsva_res(), file) }
    )

    observeEvent(input$open_export, {
      req(heatmap_gg())
      open_export_modal(session, plot_fn = function() heatmap_gg(),
                        base_name = paste0("GSVA_", isolate(input$db)))
    })
  })
}
