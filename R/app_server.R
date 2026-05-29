#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}. DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # ── Centralised state ──────────────────────────────────────────────────────
  state <- AppState$new()

  # ── Species propagation ────────────────────────────────────────────────────
  observeEvent(input$species, {
    state$species(input$species)
    state$log_param("global", "species", input$species)
  })

  # ── Reset handler ──────────────────────────────────────────────────────────
  observeEvent(input$btn_reset, {
    showModal(modalDialog(
      title  = "Reset Analysis",
      "All data and results will be cleared. Continue?",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("btn_reset_confirm", "Reset", class = "btn-danger")
      )
    ))
  })
  observeEvent(input$btn_reset_confirm, {
    state$reset()
    removeModal()
    nav_select("main_nav", "upload")
  })

  # ── Pipeline nav badges ────────────────────────────────────────────────────
  observe({
    status <- state$pipeline_status()
    steps  <- c("upload", "qc", "eda", "deg", "enrichment", "gsva", "report")
    for (s in steps) {
      badge_id  <- paste0("badge_", s)
      badge_val <- switch(status[[s]],
        done    = tags$span("Done",     class = "badge bg-success ms-1"),
        running = tags$span("Running…", class = "badge bg-warning ms-1"),
        NULL
      )
      output[[badge_id]] <- renderUI(badge_val)
    }
  })

  # ── Module servers ─────────────────────────────────────────────────────────
  mod_upload_server("upload",     state = state)
  mod_qc_server("qc",             state = state)
  mod_eda_server("eda",           state = state)
  mod_deg_server("deg",           state = state)
  mod_enrichment_server("enrichment", state = state)
  mod_gsva_server("gsva",             state = state)
  mod_deg_multi_server("upset",       state = state)
  mod_report_server("report",         state = state)
}
