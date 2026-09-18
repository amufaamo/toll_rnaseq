#' Multiverse DEG Module
#' @import shiny
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box input_task_button
#' @importFrom bsicons bs_icon
#' @noRd
mod_deg_multiverse_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 350,
      card(card_header(bsicons::bs_icon("diagram-3"), " Contrast"), card_body(
        selectInput(ns("condition_col"), "Condition column", choices = NULL),
        selectInput(ns("ref_level"), "Reference level", choices = NULL),
        selectInput(ns("test_level"), "Test level", choices = NULL)
      )),
      card(card_header(bsicons::bs_icon("layers"), " Multiverse choices"), card_body(
        checkboxInput(ns("use_covariate"), "Include covariate/no-covariate paths", TRUE),
        conditionalPanel(paste0("input['", ns("use_covariate"), "'] == true"),
                         selectInput(ns("covariate"), "Covariate", choices = NULL)),
        checkboxInput(ns("use_apeglm"), "Include apeglm path when available", TRUE),
        selectInput(ns("bootstrap_B"), "Null bootstraps", choices = c("20 (exploratory)" = 20,
          "100 (default)" = 100), selected = 100),
        numericInput(ns("seed"), "Random seed", 20260918, min = 1, step = 1)
      )),
      card(card_header(bsicons::bs_icon("sliders"), " Calling rule"), card_body(
        selectInput(ns("padj_cutoff"), "Per-specification FDR", choices = c("0.01" = .01, "0.05" = .05, "0.10" = .10), selected = .05),
        selectInput(ns("lfc_cutoff"), "Absolute log2FC", choices = c("0" = 0, "0.5" = .5, "1" = 1), selected = 1),
        selectInput(ns("target_efdr"), "Target effective FDR", choices = c("0.10" = .10, "0.20" = .20), selected = .10),
        input_task_button(ns("run_multiverse"), "Run Multiverse DEG", icon = icon("play"), class = "btn-primary w-100 mt-2")
      )),
      tags$p("Each column is one reasonable analysis. Stability requires agreement in direction, not merely many analyses.", class = "small text-muted")
    ), uiOutput(ns("main_content"))
  )
}

#' Payload for one null-bootstrap replicate
#'
#' Only the fitted observed run and the call rule are shipped to the worker; the
#' accumulated bootstrap results stay in the main session.
#' @noRd
.mv_boot_payload <- function(job) {
  list(run = job$run, null_fit = job$null_fit, b = job$b,
       seed = job$payload$seed, padj_cutoff = job$payload$padj_cutoff,
       lfc_cutoff = job$payload$lfc_cutoff)
}

#' @noRd
mod_deg_multiverse_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    observe({
      meta <- state$metadata()
      if (!is.null(meta) && ncol(meta) >= 2) {
        choices <- colnames(meta)[-1]
        updateSelectInput(session, "condition_col", choices = choices)
        updateSelectInput(session, "covariate", choices = choices)
      }
    })
    observeEvent(input$condition_col, {
      req(state$metadata(), input$condition_col)
      levels <- unique(state$metadata()[[input$condition_col]])
      updateSelectInput(session, "ref_level", choices = levels)
      updateSelectInput(session, "test_level", choices = levels, selected = levels[min(2, length(levels))])
      updateSelectInput(session, "covariate", choices = setdiff(colnames(state$metadata())[-1], input$condition_col))
    })

    # ── ExtendedTask: each completed future queues the next null replicate ──
    job <- reactiveVal(NULL)
    mv_task <- ExtendedTask$new(function(stage, payload) {
      promises::future_promise({
        if (identical(stage, "observed")) {
          run <- mv_run_observed(payload$mat, payload$meta, payload$condition_col, payload$ref,
                                 payload$test, payload$covariate, payload$use_apeglm)
          run$stability <- mv_compute_stability(run$stats, run$tested, payload$padj_cutoff,
                                                payload$lfc_cutoff)
          list(run = run, null_fit = mv_fit_full_null(run$counts, run$metadata,
            payload$condition_col, payload$ref, payload$test, payload$covariate))
        } else {
          mv_bootstrap_one(payload$run, payload$null_fit, payload$b, payload$seed,
                           payload$padj_cutoff, payload$lfc_cutoff)
        }
      }, seed = TRUE)
    }) |> bslib::bind_task_button("run_multiverse")

    observeEvent(input$run_multiverse, {
      req(state$counts(), state$metadata(), input$condition_col, input$ref_level, input$test_level)
      cov <- if (isTRUE(input$use_covariate)) input$covariate else NULL
      if (identical(cov, input$condition_col)) {
        showNotification("Choose a covariate different from the condition, or turn it off.", type = "error")
        return()
      }
      state$set_status("deg_multiverse", "running")
      payload <- list(mat = state$counts(), meta = state$metadata(), condition_col = input$condition_col,
        ref = input$ref_level, test = input$test_level, covariate = cov, use_apeglm = input$use_apeglm,
        B = as.integer(input$bootstrap_B), seed = as.integer(input$seed),
        padj_cutoff = as.numeric(input$padj_cutoff), lfc_cutoff = as.numeric(input$lfc_cutoff),
        target = as.numeric(input$target_efdr))
      job(list(stage = "observed", payload = payload, boot = list()))
      state$deg_multiverse_progress(list(stage = "observed", completed = 0, total = payload$B,
                                         message = "Running observed specification grid\u2026"))
      mv_task$invoke("observed", payload)
    })

    observe({
      req(mv_task$status() == "success")
      active <- job(); ans <- mv_task$result()
      if (identical(active$stage, "observed")) {
        active$run <- ans$run; active$null_fit <- ans$null_fit; active$stage <- "bootstrap"; active$b <- 1L
        job(active)
        state$deg_multiverse_progress(list(stage = "bootstrap", completed = 0, total = active$payload$B,
                                           message = "Running null bootstrap 1\u2026"))
        mv_task$invoke("bootstrap", .mv_boot_payload(active))
      } else if (identical(active$stage, "bootstrap")) {
        active$boot[[active$b]] <- ans
        completed <- active$b
        if (completed < active$payload$B) {
          active$b <- completed + 1L; job(active)
          state$deg_multiverse_progress(list(stage = "bootstrap", completed = completed, total = active$payload$B,
            message = paste0("Running null bootstrap ", active$b, "\u2026")))
          mv_task$invoke("bootstrap", .mv_boot_payload(active))
        } else {
          calibration <- .mv_efdr_curve(active$run, active$boot, active$payload$target)
          stability <- mv_attach_efdr(active$run$stability, calibration$curve, active$payload$target)
          mv_apply_run(state, active$run, stability, calibration$curve, active$payload$target,
                       active$payload$seed, active$payload$B)
          job(list(stage = "done"))
          state$deg_multiverse_progress(list(stage = "done", completed = completed, total = completed, message = "Complete"))
        }
      }
    })
    observe({
      req(mv_task$status() == "error")
      msg <- tryCatch({ mv_task$result(); "unknown error" }, error = conditionMessage)
      showNotification(paste("Multiverse DEG error:", msg), type = "error", duration = 12)
      state$set_status("deg_multiverse", "pending")
      state$deg_multiverse_progress(list(stage = "error", completed = 0, total = 1, message = "Failed"))
    })

    current <- reactive({ state$deg_multiverse() })
    output$main_content <- renderUI({
      if (!state$step_ready("upload")) return(ui_empty_state("lock", "Multiverse DEG Locked", "Upload data first."))
      x <- current()
      if (is.null(x)) return(card(card_header("Multiverse DEG"), card_body("Choose one contrast and run the grid of defensible analyses.")))
      fluidRow(
        column(3, value_box("Specifications", nrow(x$run$specifications), bsicons::bs_icon("layers"), theme = "primary")),
        column(3, value_box("Stable calls", sum(x$stability$called), bsicons::bs_icon("check-circle"), theme = "success")),
        column(3, value_box("Bootstraps", x$metadata$B, bsicons::bs_icon("repeat"), theme = "secondary")),
        column(3, value_box("Target eFDR", format(x$metadata$target, digits = 2), bsicons::bs_icon("sliders"), theme = "light")),
        card(card_header("Specification Curve"), card_body(selectInput(ns("gene"), "Gene", choices = x$stability$gene), plotOutput(ns("curve"), height = "650px"))),
        card(card_header("Stability-ranked DEG table", downloadButton(ns("download"), "CSV", class = "btn btn-sm btn-outline-secondary float-end")), card_body(DT::DTOutput(ns("table"))))
      )
    })
    output$curve <- renderPlot({ req(current(), input$gene); mv_plot_specification_curve(current()$run, input$gene) })
    output$table <- DT::renderDT({
      req(current())
      DT::datatable(current()$stability[order(current()$stability$efdr, -current()$stability$stability, na.last = TRUE), ],
                    rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
    })
    output$download <- downloadHandler(filename = function() "multiverse_deg_stability.csv",
      content = function(file) write.csv(current()$stability, file, row.names = FALSE))
  })
}
