#' EDA Module — PCA + Sample Distance Heatmap
#' @import shiny
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box
#' @importFrom bsicons bs_icon
#' @noRd
mod_eda_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 350,
      card(
        card_header(bsicons::bs_icon("sliders"), " PCA Options"),
        card_body(
          selectInput(ns("color_by"),  "Color by",           choices = NULL),
          selectInput(ns("shape_by"),  "Shape by (optional)", choices = NULL),
          numericInput(ns("top_genes"), "Top variable genes", 500,
                       min = 100, max = 5000, step = 100),
          checkboxInput(ns("show_labels"), "Show sample labels", FALSE),
          actionButton(ns("run_pca"), "Run PCA",
                       class = "btn-primary w-100 mt-2",
                       icon  = icon("play"))
        )
      )
    ),
    uiOutput(ns("main_content"))
  )
}

#' @noRd
mod_eda_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns      <- session$ns
    pca_data <- reactiveVal(NULL)

    # ── Populate selectors ────────────────────────────────────────────────────
    observe({
      meta <- state$metadata()
      if (!is.null(meta)) {
        choices <- c("none", colnames(meta)[-1])
        updateSelectInput(session, "color_by", choices = choices)
        updateSelectInput(session, "shape_by", choices = choices)
      }
    })

    # ── PCA computation ───────────────────────────────────────────────────────
    observeEvent(input$run_pca, {
      req(state$counts_filtered() %||% state$counts())
      mat   <- state$counts_filtered() %||% state$counts()
      n_top <- min(input$top_genes, nrow(mat))
      lvars <- apply(log1p(mat), 1, var)
      top   <- mat[order(lvars, decreasing = TRUE)[seq_len(n_top)], ]
      pca   <- prcomp(t(log1p(top)), scale. = TRUE)
      pct   <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

      df    <- as.data.frame(pca$x[, seq_len(min(4, ncol(pca$x)))])
      df$sample <- rownames(df)

      meta <- state$metadata()
      if (!is.null(meta)) {
        df <- merge(df, meta, by.x = "sample",
                    by.y = colnames(meta)[1], all.x = TRUE)
      }

      pca_data(list(df = df, pct = pct, loadings = pca$rotation[, 1:2],
                    n_top = n_top))
      state$set_status("eda", "done")
    })

    # ── PCA ggplot (cached) ───────────────────────────────────────────────────
    pca_gg <- reactive({
      req(pca_data())
      d   <- pca_data()
      df  <- d$df
      pct <- d$pct
      p   <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2))

      color_col <- input$color_by
      shape_col <- input$shape_by

      if (!is.null(color_col) && color_col != "none" && color_col %in% colnames(df))
        p <- p + ggplot2::aes(color = .data[[color_col]])
      if (!is.null(shape_col) && shape_col != "none" && shape_col %in% colnames(df))
        p <- p + ggplot2::aes(shape = .data[[shape_col]])

      p <- p +
        ggplot2::geom_point(size = 3.5, alpha = 0.85) +
        ggplot2::scale_color_manual(values = .okabe_ito(), na.value = "#AAAAAA") +
        ggplot2::labs(
          x = sprintf("PC1 (%.1f%%)", pct[1]),
          y = sprintf("PC2 (%.1f%%)", pct[2])
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(legend.position = "right")

      if (isTRUE(input$show_labels) && "sample" %in% colnames(df)) {
        if (!requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggplot2::geom_text(ggplot2::aes(label = sample),
                                      size = 3, vjust = -0.8)
        } else {
          p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 3)
        }
      }
      p
    }) |> bindCache(pca_data(), input$color_by, input$shape_by, input$show_labels)

    # ── Sample distance heatmap ───────────────────────────────────────────────
    dist_gg <- reactive({
      req(state$counts_filtered() %||% state$counts())
      mat   <- state$counts_filtered() %||% state$counts()
      vst   <- log1p(mat)
      dmat  <- as.matrix(dist(t(vst)))

      # Melt to data.frame for ggplot
      samples <- colnames(dmat)
      df_dist <- data.frame(
        s1   = rep(samples, each  = length(samples)),
        s2   = rep(samples, times = length(samples)),
        dist = as.vector(dmat)
      )
      df_dist$s1 <- factor(df_dist$s1, levels = samples)
      df_dist$s2 <- factor(df_dist$s2, levels = rev(samples))

      ggplot2::ggplot(df_dist, ggplot2::aes(s1, s2, fill = dist)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_distiller(palette = "Blues", direction = -1,
                                      name = "Euclidean\ndistance") +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }) |> bindCache(state$counts_filtered(), state$counts())

    # ── UI ────────────────────────────────────────────────────────────────────
    output$main_content <- renderUI({
      if (!state$step_ready("upload")) {
        return(ui_empty_state("lock", "EDA Locked", "Upload data first.",
          actionButton("jump_upload", "Go to Upload",
                       onclick = "Shiny.setInputValue('main_nav','upload')")))
      }
      tagList(
        fluidRow(
          column(7,
            card(
              card_header(
                "PCA",
                actionButton(ns("open_export_pca"), bsicons::bs_icon("gear"),
                             class = "btn btn-sm btn-outline-secondary float-end",
                             title = "Export options")
              ),
              card_body(plotOutput(ns("pca_plot"), height = "400px"))
            )
          ),
          column(5,
            card(
              card_header(
                "Sample Distance",
                actionButton(ns("open_export_dist"), bsicons::bs_icon("gear"),
                             class = "btn btn-sm btn-outline-secondary float-end",
                             title = "Export options")
              ),
              card_body(plotOutput(ns("dist_plot"), height = "400px"))
            )
          )
        ),
        if (!is.null(pca_data())) {
          d <- pca_data()
          card(
            card_header("Variance Explained"),
            card_body(plotOutput(ns("scree_plot"), height = "200px"))
          )
        }
      )
    })

    output$pca_plot  <- renderPlot({ req(pca_gg());  pca_gg()  })
    output$dist_plot <- renderPlot({ req(dist_gg()); dist_gg() })

    output$scree_plot <- renderPlot({
      req(pca_data())
      pct <- pca_data()$pct
      df  <- data.frame(PC = paste0("PC", seq_along(pct)), var = pct,
                        cumvar = cumsum(pct))
      df$PC <- factor(df$PC, levels = df$PC)
      ggplot2::ggplot(df[seq_len(min(10, nrow(df))), ],
                      ggplot2::aes(PC, var)) +
        ggplot2::geom_col(fill = "#18BC9C") +
        ggplot2::geom_line(ggplot2::aes(y = cumvar, group = 1),
                           color = "#2C3E50", linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes(y = cumvar), color = "#2C3E50", size = 2) +
        ggplot2::labs(x = NULL, y = "% Variance") +
        ggplot2::theme_minimal(base_size = 11)
    })

    # ── Export modal hooks ────────────────────────────────────────────────────
    observeEvent(input$open_export_pca, {
      req(pca_gg())
      open_export_modal(session, plot_fn = function() pca_gg(),
                        base_name = "PCA")
    })
    observeEvent(input$open_export_dist, {
      req(dist_gg())
      open_export_modal(session, plot_fn = function() dist_gg(),
                        base_name = "sample_distance")
    })
  })
}

.okabe_ito <- function() {
  c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
