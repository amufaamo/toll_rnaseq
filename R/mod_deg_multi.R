#' Multi-contrast DEG Module — UpSet plot
#' @import shiny
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box
#' @importFrom bsicons bs_icon
#' @noRd
mod_deg_multi_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 350,
      card(
        card_header(bsicons::bs_icon("list-check"), " UpSet Settings"),
        card_body(
          uiOutput(ns("contrast_selector")),
          numericInput(ns("fdr_thr"), "FDR cutoff",    0.05, min = 0, max = 1, step = 0.01),
          numericInput(ns("lfc_thr"), "Log2FC cutoff", 1,    min = 0, step = 0.5),
          actionButton(ns("run_upset"), "Draw UpSet",
                       class = "btn-primary w-100 mt-2",
                       icon  = icon("chart-bar"))
        )
      )
    ),
    uiOutput(ns("main_content"))
  )
}

#' @noRd
mod_deg_multi_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns        <- session$ns
    upset_data <- reactiveVal(NULL)

    # ── Populate contrast selector ────────────────────────────────────────────
    observe({
      cnts <- names(state$deg_results())
      if (length(cnts) > 0)
        updateCheckboxGroupInput(session, "contrasts",
                                 choices = cnts, selected = cnts)
    })

    output$contrast_selector <- renderUI({
      cnts <- names(state$deg_results())
      if (length(cnts) == 0)
        return(tags$p("No DEG results yet.", class = "text-muted small"))
      checkboxGroupInput(ns("contrasts"), "Contrasts",
                         choices = cnts, selected = cnts)
    })

    # ── Build membership matrix ───────────────────────────────────────────────
    observeEvent(input$run_upset, {
      req(input$contrasts, length(input$contrasts) >= 2)
      res_list <- state$deg_results()[input$contrasts]

      sig_list <- lapply(res_list, function(df) {
        gene_col <- if ("gene" %in% colnames(df)) df$gene else rownames(df)
        gene_col[!is.na(df$padj) & df$padj < input$fdr_thr &
                   abs(df$log2FoldChange) > input$lfc_thr]
      })
      sig_list <- Filter(function(g) length(g) > 0, sig_list)

      if (length(sig_list) < 2) {
        showNotification("Need ≥2 contrasts with DEGs.", type = "warning")
        return()
      }

      all_genes <- unique(unlist(sig_list))
      mat_df <- as.data.frame(
        lapply(sig_list, function(g) as.integer(all_genes %in% g)),
        check.names = FALSE
      )
      rownames(mat_df) <- all_genes
      upset_data(list(mat = mat_df, contrasts = names(sig_list)))
    })

    # ── UpSet plot (ComplexUpset or ggplot2 fallback) ─────────────────────────
    upset_gg <- reactive({
      req(upset_data())
      d    <- upset_data()
      mat  <- d$mat
      cnts <- d$contrasts

      if (requireNamespace("ComplexUpset", quietly = TRUE)) {
        ComplexUpset::upset(mat, intersect = cnts,
                            name = "Gene intersection",
                            min_size = 1,
                            width_ratio = 0.2)
      } else {
        # Fallback: intersection count bar chart
        mat$comb <- apply(mat[cnts], 1, function(x) {
          hits <- cnts[x == 1]
          if (length(hits) == 0) NA else paste(sort(hits), collapse = " ∩ ")
        })
        mat <- mat[!is.na(mat$comb), ]
        counts <- sort(table(mat$comb), decreasing = TRUE)
        df_c <- data.frame(Intersection = names(counts), Count = as.integer(counts))
        df_c$Intersection <- factor(df_c$Intersection, levels = rev(df_c$Intersection))

        ggplot2::ggplot(df_c, ggplot2::aes(Count, Intersection)) +
          ggplot2::geom_col(fill = "#2C3E50") +
          ggplot2::labs(x = "Intersection size", y = NULL,
                        caption = "Install ComplexUpset for full UpSet plot") +
          ggplot2::theme_minimal(base_size = 12)
      }
    }) |> bindCache(upset_data(), input$fdr_thr, input$lfc_thr)

    # ── UI ────────────────────────────────────────────────────────────────────
    output$main_content <- renderUI({
      if (!state$step_ready("deg")) {
        return(ui_empty_state("lock", "UpSet Locked", "Run DEG first.",
          actionButton("jump_deg", "Go to DEG",
                       onclick = "Shiny.setInputValue('main_nav','deg')")))
      }
      d <- upset_data()
      tagList(
        if (!is.null(d)) fluidRow(
          column(4, value_box("Contrasts", length(d$contrasts),
                              bsicons::bs_icon("layers"), theme = "primary")),
          column(4, value_box("Total DEGs", nrow(d$mat),
                              bsicons::bs_icon("asterisk"), theme = "secondary")),
          column(4, value_box("Shared DEGs",
                              sum(rowSums(d$mat[d$contrasts]) == length(d$contrasts)),
                              bsicons::bs_icon("intersect"), theme = "success"))
        ),
        card(
          card_header(
            "UpSet Plot",
            actionButton(ns("open_export"), bsicons::bs_icon("gear"),
                         class = "btn btn-sm btn-outline-secondary float-end",
                         title = "Export"),
            downloadButton(ns("dl_pdf"), "PDF",
                           class = "btn btn-sm btn-outline-secondary float-end me-1")
          ),
          card_body(plotOutput(ns("upset_plot"), height = "500px"))
        )
      )
    })

    output$upset_plot <- renderPlot({ req(upset_gg()); upset_gg() })

    output$dl_pdf <- downloadHandler(
      filename = function() paste0("UpSet_", Sys.Date(), ".pdf"),
      content  = function(file) {
        req(upset_gg())
        ggplot2::ggsave(file, plot = upset_gg(),
                        width = 183, height = 120, units = "mm",
                        device = grDevices::cairo_pdf)
      }
    )

    observeEvent(input$open_export, {
      req(upset_gg())
      open_export_modal(session, plot_fn = function() upset_gg(),
                        base_name = "UpSet_DEG")
    })
  })
}
