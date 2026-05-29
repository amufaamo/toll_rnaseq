#' QC & Filter Module
#' @import shiny
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box
#' @importFrom bsicons bs_icon
#' @noRd
mod_qc_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 350,
      card(
        card_header(bsicons::bs_icon("funnel"), " Filtering"),
        card_body(
          numericInput(ns("min_count"), "Min total count per gene", 10, min = 0),
          numericInput(ns("min_samples"), "Min samples with count > 0", 2, min = 1),
          actionButton(ns("run_filter"), "Apply Filters",
                       class = "btn-primary w-100 mt-2",
                       icon  = icon("play"))
        )
      )
    ),
    uiOutput(ns("main_content"))
  )
}

#' @noRd
mod_qc_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$run_filter, {
      req(state$counts())
      mat <- state$counts()
      keep <- rowSums(mat) >= input$min_count &
              rowSums(mat > 0) >= input$min_samples
      filtered <- mat[keep, , drop = FALSE]
      state$counts_filtered(filtered)
      state$log_param("qc", "min_count",   input$min_count)
      state$log_param("qc", "min_samples", input$min_samples)
      state$log_param("qc", "genes_kept",  sum(keep))
      state$set_status("qc", "done")
      showNotification(
        sprintf("Filtered: %s / %s genes retained", sum(keep), nrow(mat)),
        type = "message", duration = 4
      )
    })

    output$main_content <- renderUI({
      if (!state$step_ready("upload")) {
        return(ui_empty_state(
          "lock", "QC Locked",
          "Upload a count matrix first.",
          actionButton("jump_upload", "Go to Upload",
                       onclick = "Shiny.setInputValue('main_nav','upload')")
        ))
      }
      counts   <- state$counts()
      filtered <- state$counts_filtered()
      tagList(
        fluidRow(
          column(4, value_box("Before filter", format(nrow(counts), big.mark=","),
                              bsicons::bs_icon("grid"), theme = "light")),
          column(4, value_box("After filter",
                              if (!is.null(filtered)) format(nrow(filtered), big.mark=",") else "—",
                              bsicons::bs_icon("funnel-fill"),
                              theme = if (!is.null(filtered)) "success" else "light")),
          column(4, value_box("Removed",
                              if (!is.null(filtered)) format(nrow(counts)-nrow(filtered), big.mark=",") else "—",
                              bsicons::bs_icon("trash"), theme = "light"))
        ),
        card(
          card_header("Library Size Distribution"),
          card_body(plotOutput(ns("libsize_plot"), height = "300px"))
        )
      )
    })

    output$libsize_plot <- renderPlot({
      req(state$counts())
      mat <- state$counts()
      df  <- data.frame(
        sample = factor(colnames(mat), levels = colnames(mat)),
        libsize = colSums(mat) / 1e6
      )
      ggplot2::ggplot(df, ggplot2::aes(sample, libsize)) +
        ggplot2::geom_col(fill = "#18BC9C") +
        ggplot2::labs(x = NULL, y = "Library size (M reads)") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })
  })
}
