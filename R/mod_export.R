#' Journal-Ready Export Modal Module
#'
#' Shared module providing a modal dialog with journal presets
#' (Nature / Cell / Science / Custom), format selector, font controls,
#' and live preview. Any plot module can open this modal via
#' `open_export_modal(session, plot_reactive, filename_base)`.
#'
#' @import shiny
#' @importFrom bslib card card_header card_body layout_columns sidebar layout_sidebar
#' @importFrom bsicons bs_icon
#' @noRd

# ── Journal presets ────────────────────────────────────────────────────────────
.journal_presets <- list(
  "Nature 1-col (89mm)"  = list(w = 89,  h = 89,  unit = "mm"),
  "Nature 2-col (183mm)" = list(w = 183, h = 120, unit = "mm"),
  "Cell 1-col (85mm)"    = list(w = 85,  h = 85,  unit = "mm"),
  "Cell 2-col (174mm)"   = list(w = 174, h = 120, unit = "mm"),
  "Science 1-col (57mm)" = list(w = 57,  h = 57,  unit = "mm"),
  "Science 2-col (120mm)"= list(w = 120, h = 90,  unit = "mm"),
  "Custom"               = list(w = 120, h = 90,  unit = "mm")
)

#' Open the export modal from any module server
#'
#' @param session Shiny session
#' @param plot_fn Zero-arg function that returns the ggplot2 object
#' @param base_name Default filename base (no extension)
#' @noRd
open_export_modal <- function(session, plot_fn, base_name = "figure") {
  ns <- session$ns
  showModal(modalDialog(
    title = tagList(bsicons::bs_icon("download"), " Export Figure"),
    size  = "xl",
    easyClose = TRUE,
    footer = tagList(
      modalButton("Cancel"),
      downloadButton(ns("dl_export_final"), "Download",
                     class = "btn-primary")
    ),
    fluidRow(
      # ── Settings panel ─────────────────────────────────────────────────────
      column(4,
        h6("Journal Preset"),
        selectInput(ns("export_preset"), NULL,
                    choices  = names(.journal_presets),
                    selected = "Nature 1-col (89mm)"),
        h6("Format"),
        radioButtons(ns("export_fmt"), NULL,
                     choices = c("PDF (cairo)" = "pdf",
                                 "SVG"          = "svg",
                                 "TIFF 300 DPI" = "tiff",
                                 "PNG 300 DPI"  = "png"),
                     inline  = FALSE),
        h6("Dimensions"),
        fluidRow(
          column(6, numericInput(ns("export_w"), "Width",  89, min = 20, max = 500)),
          column(6, numericInput(ns("export_h"), "Height", 89, min = 20, max = 500))
        ),
        selectInput(ns("export_unit"), "Unit",
                    choices = c("mm", "cm", "in"),
                    selected = "mm"),
        h6("Typography"),
        selectInput(ns("export_font"), "Font family",
                    choices  = c("Helvetica", "Arial", "Times New Roman",
                                 "Roboto", "Inter"),
                    selected = "Helvetica"),
        numericInput(ns("export_base_size"), "Base font size (pt)", 7,
                     min = 5, max = 20, step = 0.5),
        checkboxInput(ns("export_transparent"), "Transparent background", FALSE),
        tags$small("Preview is approximate; exact rendering may differ.", class="text-muted")
      ),
      # ── Live preview ────────────────────────────────────────────────────────
      column(8,
        plotOutput(ns("export_preview"), height = "480px")
      )
    )
  ))

  # Sync preset → dimensions
  observeEvent(session$input$export_preset, {
    p <- .journal_presets[[session$input$export_preset]]
    if (!is.null(p)) {
      updateNumericInput(session, "export_w",    value = p$w)
      updateNumericInput(session, "export_h",    value = p$h)
      updateSelectInput(session,  "export_unit", selected = p$unit)
    }
  }, ignoreInit = TRUE)

  # Live preview
  output <- session$output
  output[[ns("export_preview")]] <- renderPlot({
    p <- plot_fn()
    req(inherits(p, "ggplot"))
    size <- session$input$export_base_size %||% 7
    p + ggplot2::theme_minimal(base_size = size + 4) +
      ggplot2::theme(
        text = ggplot2::element_text(family = session$input$export_font %||% "Helvetica"),
        plot.background = if (isTRUE(session$input$export_transparent))
          ggplot2::element_blank() else ggplot2::element_rect(fill = "white")
      )
  }, bg = "transparent")

  # Final download
  output[[ns("dl_export_final")]] <- downloadHandler(
    filename = function() {
      ext <- switch(session$input$export_fmt %||% "pdf",
        pdf = "pdf", svg = "svg", tiff = "tiff", png = "png", "pdf")
      paste0(base_name, ".", ext)
    },
    content = function(file) {
      p    <- plot_fn()
      req(inherits(p, "ggplot"))
      w    <- session$input$export_w    %||% 89
      h    <- session$input$export_h    %||% 89
      unit <- session$input$export_unit %||% "mm"
      size <- session$input$export_base_size %||% 7
      font <- session$input$export_font %||% "Helvetica"
      bg   <- if (isTRUE(session$input$export_transparent)) "transparent" else "white"

      p_out <- p + ggplot2::theme_minimal(base_size = size + 4) +
        ggplot2::theme(
          text = ggplot2::element_text(family = font),
          plot.background = if (bg == "transparent")
            ggplot2::element_blank() else ggplot2::element_rect(fill = "white")
        )

      fmt <- session$input$export_fmt %||% "pdf"
      switch(fmt,
        pdf  = ggplot2::ggsave(file, plot = p_out, width = w, height = h,
                               units = unit, device = grDevices::cairo_pdf,
                               dpi = 300),
        svg  = {
          if (!requireNamespace("svglite", quietly = TRUE))
            stop("SVG export requires the svglite package: install.packages('svglite')")
          ggplot2::ggsave(file, plot = p_out, width = w, height = h,
                          units = unit, device = svglite::svglite)
        },
        tiff = ggplot2::ggsave(file, plot = p_out, width = w, height = h,
                               units = unit, dpi = 300,
                               device = grDevices::tiff, compression = "lzw"),
        png  = ggplot2::ggsave(file, plot = p_out, width = w, height = h,
                               units = unit, dpi = 300)
      )
    }
  )
}
