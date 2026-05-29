#' Time-series Module — maSigPro temporal expression analysis
#' @import shiny
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box input_task_button
#' @importFrom bsicons bs_icon
#' @noRd
mod_timeseries_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 350,
      card(
        card_header(bsicons::bs_icon("activity"), " Experimental Design"),
        card_body(
          helpText("Metadata must have numeric time, integer replicate, and 0/1 group columns."),
          selectInput(ns("time_col"),  "Time column",      choices = NULL),
          selectInput(ns("rep_col"),   "Replicate column", choices = NULL),
          checkboxGroupInput(ns("group_cols"), "Group columns", choices = NULL)
        )
      ),
      card(
        card_header(bsicons::bs_icon("sliders"), " maSigPro Options"),
        card_body(
          numericInput(ns("degree"),      "Polynomial degree",  2,    min=1, max=3),
          numericInput(ns("q_thr"),       "Q-value cutoff",     0.05, min=0, max=1, step=0.01),
          numericInput(ns("rsq"),         "R² cutoff",          0.7,  min=0, max=1, step=0.05),
          numericInput(ns("top_n_plot"),  "Top genes to plot",  9,    min=1, max=25),
          input_task_button(ns("run_masigpro"), "Run maSigPro",
                            icon  = icon("play"),
                            class = "btn-primary w-100 mt-2")
        )
      )
    ),
    uiOutput(ns("main_content"))
  )
}

#' @noRd
mod_timeseries_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns      <- session$ns
    ts_res  <- reactiveVal(NULL)

    # ── Populate selectors from metadata ─────────────────────────────────────
    observe({
      meta <- state$metadata()
      if (is.null(meta) || ncol(meta) < 2) return()
      cols <- colnames(meta)[-1]
      time_guess  <- cols[grepl("time|week|day|hour", tolower(cols))][1]
      rep_guess   <- cols[grepl("rep|replicate|batch", tolower(cols))][1]
      grp_guess   <- cols[grepl("group|condition|trt|treat", tolower(cols))]
      updateSelectInput(session, "time_col",  choices=cols, selected=time_guess %||% cols[1])
      updateSelectInput(session, "rep_col",   choices=cols, selected=rep_guess  %||% cols[min(2,length(cols))])
      updateCheckboxGroupInput(session, "group_cols", choices=cols, selected=grp_guess)
    })

    # ── ExtendedTask ──────────────────────────────────────────────────────────
    ts_task <- ExtendedTask$new(function(mat, meta, sample_col,
                                          time_col, rep_col, group_cols,
                                          degree, q_thr, rsq) {
      promises::future_promise({
        if (!requireNamespace("maSigPro", quietly=TRUE))
          stop("maSigPro not installed. Run: BiocManager::install('maSigPro')")

        samples  <- meta[[sample_col]]
        common   <- intersect(colnames(mat), samples)
        if (length(common) < 4)
          stop(sprintf("Only %d samples match metadata.", length(common)))
        mat     <- mat[, common, drop=FALSE]
        meta    <- meta[match(common, meta[[sample_col]]), , drop=FALSE]

        edesign <- data.frame(
          Time       = as.numeric(meta[[time_col]]),
          Replicates = as.integer(meta[[rep_col]]),
          row.names  = common
        )
        for (g in group_cols) {
          val <- meta[[g]]
          edesign[[g]] <- if (is.numeric(val)) val else as.integer(factor(val)) - 1L
        }
        if (length(group_cols) == 0) edesign$Group1 <- 1L

        # log2 normalise for maSigPro
        expr   <- log2(mat + 1)
        design <- maSigPro::make.design.matrix(edesign, degree=degree,
                                                time.col=1, repl.col=2,
                                                group.cols=seq(3, ncol(edesign)))
        pvec  <- maSigPro::p.vector(expr, design, Q=q_thr,
                                     MT.adjust="BH", min.obs=max(4, ncol(mat)/4))
        if (pvec$i == 0)
          return(list(n_sig=0, pvec=pvec, tstep=NULL, sigs=NULL, edesign=edesign, expr=expr))

        tstep <- maSigPro::T.fit(pvec, step.method="backward", alfa=q_thr)
        sigs  <- maSigPro::get.siggenes(tstep, rsq=rsq, var="groups")
        list(n_sig=pvec$i, pvec=pvec, tstep=tstep, sigs=sigs, edesign=edesign, expr=expr)
      })
    }) |> bslib::bind_task_button("run_masigpro")

    observeEvent(input$run_masigpro, {
      req(state$counts_filtered() %||% state$counts(), state$metadata())
      req(input$time_col, input$rep_col)
      state$set_status("timeseries", "running")
      ts_task$invoke(
        mat        = state$counts_filtered() %||% state$counts(),
        meta       = state$metadata(),
        sample_col = colnames(state$metadata())[1],
        time_col   = input$time_col,
        rep_col    = input$rep_col,
        group_cols = input$group_cols %||% character(0),
        degree     = input$degree,
        q_thr      = input$q_thr,
        rsq        = input$rsq
      )
    })

    observe({
      req(ts_task$status() == "success")
      res <- ts_task$result()
      ts_res(res)
      state$log_param("timeseries", "n_sig",  res$n_sig)
      state$log_param("timeseries", "degree", isolate(input$degree))
      state$set_status("timeseries", "done")
    })

    observe({
      req(ts_task$status() == "error")
      showNotification(paste("maSigPro:", ts_task$result()$message),
                       type="error", duration=10)
      state$set_status("timeseries", "pending")
    })

    # ── Time-series line plot ─────────────────────────────────────────────────
    ts_gg <- reactive({
      req(ts_res())
      res   <- ts_res()
      if (res$n_sig == 0 || is.null(res$sigs))
        return(.ts_msg("No significant genes. Try relaxing Q-value or R² cutoff."))

      sigs    <- res$sigs
      expr    <- res$expr
      edesign <- res$edesign
      n_top   <- min(isolate(input$top_n_plot), res$n_sig, nrow(expr))

      sig_names <- rownames(sigs$summary)[sigs$summary[,1] != ""]
      if (length(sig_names) == 0) return(.ts_msg("No genes passed the R² filter."))
      top_genes <- head(sig_names, n_top)

      df <- do.call(rbind, lapply(top_genes, function(g) {
        data.frame(gene=g, time=edesign$Time,
                   expr=as.numeric(expr[g, rownames(edesign)]),
                   stringsAsFactors=FALSE)
      }))

      ggplot2::ggplot(df, ggplot2::aes(time, expr)) +
        ggplot2::stat_summary(fun=mean, geom="line", color="#0072B2", linewidth=0.9) +
        ggplot2::geom_point(color="#0072B2", alpha=0.6, size=1.5) +
        ggplot2::facet_wrap(~gene, scales="free_y",
                            ncol=min(3, ceiling(sqrt(n_top)))) +
        ggplot2::labs(x="Time", y="Log2 expression") +
        ggplot2::theme_minimal(base_size=10) +
        ggplot2::theme(strip.text=ggplot2::element_text(size=7, face="bold"))
    }) |> bindCache(ts_res(), input$top_n_plot)

    # ── UI ────────────────────────────────────────────────────────────────────
    output$main_content <- renderUI({
      if (!state$step_ready("upload")) {
        return(ui_empty_state("lock", "Time-series Locked",
          "Upload data and set metadata first.",
          actionButton("jump_upload", "Go to Upload",
                       onclick="Shiny.setInputValue('main_nav','upload')")))
      }
      if (!requireNamespace("maSigPro", quietly=TRUE)) {
        return(ui_empty_state(
          "box-arrow-in-down", "maSigPro Not Installed",
          "Run in R console: BiocManager::install('maSigPro')",
          tags$code("BiocManager::install('maSigPro')")
        ))
      }
      res <- ts_res()
      tagList(
        if (!is.null(res)) fluidRow(
          column(4, value_box("Sig. genes", res$n_sig,
                              bsicons::bs_icon("activity"), theme="primary")),
          column(4, value_box("Degree", isolate(input$degree),
                              bsicons::bs_icon("graph-up"), theme="secondary")),
          column(4, value_box("Q cutoff", isolate(input$q_thr),
                              bsicons::bs_icon("filter"), theme="light"))
        ),
        card(
          card_header(
            "Top Significant Gene Profiles",
            actionButton(ns("open_export"), bsicons::bs_icon("gear"),
                         class="btn btn-sm btn-outline-secondary float-end",
                         title="Export"),
            downloadButton(ns("dl_csv"), "Sig. CSV",
                           class="btn btn-sm btn-outline-secondary float-end me-1")
          ),
          card_body(plotOutput(ns("ts_plot"), height="500px"))
        )
      )
    })

    output$ts_plot <- renderPlot({ req(ts_gg()); ts_gg() })

    output$dl_csv <- downloadHandler(
      filename = "maSigPro_significant_genes.csv",
      content  = function(file) {
        req(ts_res()); write.csv(ts_res()$sigs$summary, file)
      }
    )

    observeEvent(input$open_export, {
      req(ts_gg())
      open_export_modal(session, plot_fn=function() ts_gg(),
                        base_name="maSigPro_timeseries")
    })
  })
}

.ts_msg <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x=0.5, y=0.5, label=msg, size=5, color="#666") +
    ggplot2::theme_void()
}
