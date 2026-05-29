#' DEG Module — DESeq2 + apeglm LFC shrinkage + ComBat-seq batch correction
#' @import shiny
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box input_task_button
#' @importFrom bsicons bs_icon
#' @noRd
mod_deg_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 350,
      # ── Analysis settings ────────────────────────────────────────────────
      card(
        card_header(bsicons::bs_icon("graph-up"), " DESeq2"),
        card_body(
          selectInput(ns("condition_col"), "Condition column", choices = NULL),
          selectInput(ns("ref_level"),     "Reference level",  choices = NULL),
          selectInput(ns("test_level"),    "Test level",       choices = NULL),
          checkboxInput(ns("use_apeglm"), "LFC shrinkage (apeglm)",         value = TRUE),
          checkboxInput(ns("use_ihw"),    "IHW multiple testing correction", value = FALSE)
        )
      ),
      # ── Batch correction ─────────────────────────────────────────────────
      card(
        card_header(bsicons::bs_icon("layers"), " Batch Correction"),
        card_body(
          checkboxInput(ns("use_combat"), "ComBat-seq (sva)", value = FALSE),
          conditionalPanel(
            condition = paste0("input['", ns("use_combat"), "'] == true"),
            selectInput(ns("batch_col"), "Batch column", choices = NULL)
          )
        )
      ),
      # ── Thresholds ───────────────────────────────────────────────────────
      card(
        card_header(bsicons::bs_icon("sliders"), " Thresholds"),
        card_body(
          numericInput(ns("padj_cutoff"), "FDR cutoff",    0.05, min = 0, max = 1, step = 0.01),
          numericInput(ns("lfc_cutoff"),  "Log2FC cutoff", 1,    min = 0, step = 0.5),
          input_task_button(ns("run_deg"), "Run DEG", icon = icon("play"),
                            class = "btn-primary w-100 mt-2")
        )
      )
    ),
    uiOutput(ns("main_content"))
  )
}

#' @noRd
mod_deg_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Populate selectors from metadata ─────────────────────────────────────
    observe({
      meta <- state$metadata()
      if (!is.null(meta) && ncol(meta) >= 2) {
        cond_choices <- colnames(meta)[-1]
        updateSelectInput(session, "condition_col", choices = cond_choices)
        updateSelectInput(session, "batch_col",     choices = cond_choices)
      }
    })

    observeEvent(input$condition_col, {
      req(state$metadata(), input$condition_col)
      lvls <- unique(state$metadata()[[input$condition_col]])
      updateSelectInput(session, "ref_level",  choices = lvls)
      updateSelectInput(session, "test_level", choices = lvls,
                        selected = lvls[min(2, length(lvls))])
    })

    # ── ExtendedTask: DESeq2 async ────────────────────────────────────────────
    deg_task <- ExtendedTask$new(function(mat, meta, cond_col, ref, test,
                                          use_apeglm, use_ihw,
                                          use_combat, batch_col) {
      promises::future_promise({
        if (!requireNamespace("DESeq2", quietly = TRUE))
          stop("DESeq2 not installed. Run: BiocManager::install('DESeq2')")

        # ── Align samples ────────────────────────────────────────────────────
        sample_col <- colnames(meta)[1]
        samples    <- intersect(meta[[sample_col]], colnames(mat))
        if (length(samples) == 0)
          stop("No sample names match between count matrix and metadata.")
        mat  <- mat[, samples, drop = FALSE]
        meta <- meta[meta[[sample_col]] %in% samples, , drop = FALSE]
        rownames(meta) <- meta[[sample_col]]
        meta <- meta[samples, , drop = FALSE]

        # ── ComBat-seq batch correction ──────────────────────────────────────
        if (use_combat && !is.null(batch_col) && batch_col %in% colnames(meta)) {
          if (!requireNamespace("sva", quietly = TRUE))
            stop("sva not installed. Run: BiocManager::install('sva')")
          batch <- meta[[batch_col]]
          group <- meta[[cond_col]]
          mat   <- sva::ComBat_seq(mat, batch = batch, group = group)
        }

        # ── DESeq2 ───────────────────────────────────────────────────────────
        meta[[cond_col]] <- relevel(factor(meta[[cond_col]]), ref = ref)
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = mat,
          colData   = meta,
          design    = as.formula(paste("~", cond_col))
        )
        dds <- DESeq2::DESeq(dds)

        # ── IHW (optional) ───────────────────────────────────────────────────
        if (use_ihw && requireNamespace("IHW", quietly = TRUE)) {
          res <- DESeq2::results(dds, contrast = c(cond_col, test, ref),
                                 filterFun = IHW::ihw)
        } else {
          res <- DESeq2::results(dds, contrast = c(cond_col, test, ref))
        }

        # ── apeglm LFC shrinkage ─────────────────────────────────────────────
        if (use_apeglm && requireNamespace("apeglm", quietly = TRUE)) {
          coef_names <- DESeq2::resultsNames(dds)
          # Find the coefficient matching test vs ref
          coef_pat <- paste0(cond_col, "_", gsub("[^A-Za-z0-9]", ".", test),
                             "_vs_", gsub("[^A-Za-z0-9]", ".", ref))
          coef_hit <- grep(coef_pat, coef_names, value = TRUE)
          if (length(coef_hit) == 1) {
            res <- DESeq2::lfcShrink(dds, coef = coef_hit, type = "apeglm")
          }
        }

        df <- as.data.frame(res)
        df$gene <- rownames(df)
        df[order(df$padj, na.last = TRUE), ]
      })
    }) |> bslib::bind_task_button("run_deg")

    observeEvent(input$run_deg, {
      req(state$counts_filtered() %||% state$counts(), state$metadata())
      state$set_status("deg", "running")
      deg_task$invoke(
        mat        = state$counts_filtered() %||% state$counts(),
        meta       = state$metadata(),
        cond_col   = input$condition_col,
        ref        = input$ref_level,
        test       = input$test_level,
        use_apeglm = input$use_apeglm,
        use_ihw    = input$use_ihw,
        use_combat = input$use_combat,
        batch_col  = input$batch_col
      )
    })

    observe({
      req(deg_task$status() == "success")
      result_df <- deg_task$result()
      existing  <- isolate(state$deg_results())
      key       <- paste0(isolate(input$test_level), "_vs_", isolate(input$ref_level))
      existing[[key]] <- result_df
      state$deg_results(existing)
      state$log_param("deg", "contrast", key)
      state$log_param("deg", "n_sig",
        sum(!is.na(result_df$padj) & result_df$padj < isolate(input$padj_cutoff)))
      state$set_status("deg", "done")
    })

    observe({
      req(deg_task$status() == "error")
      showNotification(paste("DEG error:", deg_task$result()$message),
                       type = "error", duration = 10)
      state$set_status("deg", "pending")
    })

    # ── Cached volcano plot ───────────────────────────────────────────────────
    current_res <- reactive({
      req(length(state$deg_results()) > 0)
      state$deg_results()[[length(state$deg_results())]]
    })

    volcano_plot <- reactive({
      res <- current_res()
      req(nrow(res) > 0)
      res <- res[!is.na(res$padj), ]
      padj_cut <- input$padj_cutoff
      lfc_cut  <- input$lfc_cutoff
      res$significance <- "NS"
      res$significance[res$padj < padj_cut & res$log2FoldChange >  lfc_cut] <- "Up"
      res$significance[res$padj < padj_cut & res$log2FoldChange < -lfc_cut] <- "Down"
      ggplot2::ggplot(res, ggplot2::aes(log2FoldChange, -log10(padj),
                                        color = significance)) +
        ggplot2::geom_point(alpha = 0.5, size = 1.2) +
        ggplot2::scale_color_manual(
          values = c(Up = "#D55E00", Down = "#0072B2", NS = "#AAAAAA")) +
        ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut),
                            linetype = "dashed", color = "#666666") +
        ggplot2::geom_hline(yintercept = -log10(padj_cut),
                            linetype = "dashed", color = "#666666") +
        ggplot2::labs(x = "Log2 Fold Change", y = expression(-log[10](FDR)),
                      color = NULL) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(legend.position = "top")
    }) |> bindCache(state$deg_results(), input$padj_cutoff, input$lfc_cutoff)

    # ── MA plot ───────────────────────────────────────────────────────────────
    ma_plot <- reactive({
      res <- current_res()
      req("baseMean" %in% colnames(res))
      res <- res[!is.na(res$padj), ]
      padj_cut <- input$padj_cutoff
      res$sig <- res$padj < padj_cut
      ggplot2::ggplot(res, ggplot2::aes(log2(baseMean + 1), log2FoldChange,
                                        color = sig)) +
        ggplot2::geom_point(alpha = 0.4, size = 0.9) +
        ggplot2::scale_color_manual(
          values = c("TRUE" = "#D55E00", "FALSE" = "#AAAAAA"),
          labels = c("TRUE" = "Significant", "FALSE" = "NS")) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::labs(x = "Log2(mean expression + 1)", y = "Log2 Fold Change",
                      color = NULL) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(legend.position = "top")
    }) |> bindCache(state$deg_results(), input$padj_cutoff)

    # ── UI ────────────────────────────────────────────────────────────────────
    output$main_content <- renderUI({
      if (!state$step_ready("upload")) {
        return(ui_empty_state("lock", "DEG Locked", "Upload data first.",
          actionButton("jump_upload", "Go to Upload",
                       onclick = "Shiny.setInputValue('main_nav','upload')")))
      }
      results <- state$deg_results()
      tagList(
        if (length(results) > 0) {
          r       <- results[[length(results)]]
          padj    <- input$padj_cutoff
          lfc     <- input$lfc_cutoff
          n_up    <- sum(!is.na(r$padj) & r$padj < padj & r$log2FoldChange >  lfc)
          n_down  <- sum(!is.na(r$padj) & r$padj < padj & r$log2FoldChange < -lfc)
          fluidRow(
            column(3, value_box("Contrasts",     length(results),
                                bsicons::bs_icon("layers"),    theme = "primary")),
            column(3, value_box("Up",            n_up,
                                bsicons::bs_icon("arrow-up"),  theme = "danger")),
            column(3, value_box("Down",          n_down,
                                bsicons::bs_icon("arrow-down"), theme = "info")),
            column(3, value_box("Total DEG",     n_up + n_down,
                                bsicons::bs_icon("asterisk"), theme = "secondary"))
          )
        },
        # ── Plot cards ──────────────────────────────────────────────────────
        fluidRow(
          column(6,
            card(
              card_header(
                "Volcano Plot",
                actionButton(ns("open_export_volcano"), bsicons::bs_icon("gear"),
                             class = "btn btn-sm btn-outline-secondary float-end",
                             title = "Export options")
              ),
              card_body(
                plotOutput(ns("volcano_plot"), height = "400px")
              )
            )
          ),
          column(6,
            card(
              card_header("MA Plot"),
              card_body(
                plotOutput(ns("ma_plot"), height = "400px")
              )
            )
          )
        ),
        # ── DEG table ───────────────────────────────────────────────────────
        card(
          card_header(
            "DEG Table",
            downloadButton(ns("dl_deg_csv"), "CSV",
                           class = "btn btn-sm btn-outline-secondary float-end ms-1"),
            downloadButton(ns("dl_volcano_pdf"), "Volcano PDF",
                           class = "btn btn-sm btn-outline-secondary float-end")
          ),
          card_body(DT::DTOutput(ns("deg_table")))
        )
      )
    })

    output$volcano_plot <- renderPlot({ req(volcano_plot()); volcano_plot() })
    output$ma_plot      <- renderPlot({ req(ma_plot());      ma_plot()      })

    output$deg_table <- DT::renderDT({
      req(length(state$deg_results()) > 0)
      res <- current_res()
      cols_show <- intersect(c("gene","baseMean","log2FoldChange","lfcSE",
                                "stat","pvalue","padj"), colnames(res))
      tbl <- res[, cols_show, drop = FALSE]
      tbl[, sapply(tbl, is.numeric)] <- round(
        tbl[, sapply(tbl, is.numeric)], 4)
      DT::datatable(
        tbl, filter = "top", rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15,
                       order = list(list(which(cols_show == "padj") - 1, "asc")))
      ) |>
        DT::formatStyle(
          "log2FoldChange",
          background = DT::styleInterval(
            c(-input$lfc_cutoff, input$lfc_cutoff),
            c("#CCE5FF", "white", "#FFD7CC")
          )
        )
    }, server = TRUE)

    # ── Downloads ─────────────────────────────────────────────────────────────
    output$dl_volcano_pdf <- downloadHandler(
      filename = function() paste0("volcano_", names(state$deg_results())[length(state$deg_results())], ".pdf"),
      content  = function(file) {
        req(volcano_plot())
        ggplot2::ggsave(file, plot = volcano_plot(),
                        width = 89, height = 89, units = "mm",
                        device = cairo_pdf)
      }
    )

    output$dl_deg_csv <- downloadHandler(
      filename = function() paste0("DEG_", names(state$deg_results())[length(state$deg_results())], ".csv"),
      content  = function(file) {
        req(length(state$deg_results()) > 0)
        write.csv(current_res(), file)
      }
    )

    # ── Journal-Ready Export modal ────────────────────────────────────────────
    observeEvent(input$open_export_volcano, {
      req(volcano_plot())
      key <- names(state$deg_results())[length(state$deg_results())]
      open_export_modal(session,
        plot_fn   = function() volcano_plot(),
        base_name = paste0("volcano_", key))
    })
  })
}
