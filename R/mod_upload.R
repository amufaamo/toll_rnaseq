#' Upload Module
#'
#' Handles count matrix upload (CSV/TSV/RDS), format auto-detection,
#' and sample metadata input. Sets state$counts and state$metadata.
#' @import shiny
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box
#' @importFrom bsicons bs_icon
#'
#' @param id Module namespace id
#' @param state AppState R6 object
#' @noRd
mod_upload_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 350,
      open  = TRUE,
      # ── File input ──────────────────────────────────────────────────────
      card(
        card_header(bsicons::bs_icon("cloud-upload"), " Count Matrix"),
        card_body(
          fileInput(ns("count_file"), NULL,
                    accept  = c(".csv", ".tsv", ".txt", ".rds", ".rdata"),
                    placeholder = "CSV / TSV / RDS"),
          helpText("Rows = genes, Columns = samples. Header required.")
        )
      ),
      # ── Metadata ────────────────────────────────────────────────────────
      card(
        card_header(bsicons::bs_icon("table"), " Sample Metadata"),
        card_body(
          fileInput(ns("meta_file"), NULL,
                    accept = c(".csv", ".tsv", ".txt"),
                    placeholder = "CSV / TSV (optional)"),
          helpText("Column 1 = sample names matching count matrix columns."),
          actionButton(ns("use_colnames"), "Auto-generate from column names",
                       class = "btn-outline-secondary btn-sm w-100 mt-1")
        )
      )
    ),
    # ── Main panel ──────────────────────────────────────────────────────────
    uiOutput(ns("main_content"))
  )
}

#' @noRd
mod_upload_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Auto-detect format and load counts ──────────────────────────────────
    observeEvent(input$count_file, {
      req(input$count_file)
      path <- input$count_file$datapath
      ext  <- tolower(tools::file_ext(input$count_file$name))

      mat <- tryCatch({
        if (ext %in% c("rds", "rdata")) {
          obj <- readRDS(path)
          if (is.data.frame(obj) || is.matrix(obj)) obj else stop("RDS must contain a matrix or data.frame")
        } else {
          sep <- if (ext == "csv") "," else "\t"
          read.table(path, header = TRUE, sep = sep, row.names = 1, check.names = FALSE)
        }
      }, error = function(e) {
        showNotification(paste("Upload error:", e$message), type = "error")
        NULL
      })

      if (!is.null(mat)) {
        state$counts(as.matrix(mat))
        state$log_param("upload", "file_name", input$count_file$name)
        state$log_param("upload", "n_genes",   nrow(mat))
        state$log_param("upload", "n_samples", ncol(mat))
        state$set_status("upload", "done")
        showNotification(
          sprintf("Loaded %s genes × %s samples", nrow(mat), ncol(mat)),
          type = "message", duration = 4
        )
      }
    })

    # ── Auto-generate metadata from column names ─────────────────────────────
    observeEvent(input$use_colnames, {
      req(state$counts())
      meta <- data.frame(
        sample    = colnames(state$counts()),
        condition = "untreated",
        row.names = NULL,
        stringsAsFactors = FALSE
      )
      state$metadata(meta)
      showNotification("Metadata auto-generated. Edit the condition column.", type = "message")
    })

    # ── Load metadata file ────────────────────────────────────────────────────
    observeEvent(input$meta_file, {
      req(input$meta_file)
      path <- input$meta_file$datapath
      ext  <- tolower(tools::file_ext(input$meta_file$name))
      sep  <- if (ext == "csv") "," else "\t"
      meta <- tryCatch(
        read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE),
        error = function(e) { showNotification(paste("Metadata error:", e$message), type = "error"); NULL }
      )
      if (!is.null(meta)) {
        state$metadata(meta)
        state$log_param("upload", "meta_cols", paste(colnames(meta), collapse = ","))
      }
    })

    # ── Main panel rendering ──────────────────────────────────────────────────
    output$main_content <- renderUI({
      counts <- state$counts()
      meta   <- state$metadata()

      if (is.null(counts)) {
        ui_empty_state(
          "cloud-arrow-up", "No Data Loaded",
          "Upload a count matrix to begin."
        )
      } else {
        tagList(
          fluidRow(
            column(4, value_box(
              title = "Genes",    value = format(nrow(counts), big.mark = ","),
              showcase = bsicons::bs_icon("dna"), theme = "primary"
            )),
            column(4, value_box(
              title = "Samples",  value = ncol(counts),
              showcase = bsicons::bs_icon("people"), theme = "secondary"
            )),
            column(4, value_box(
              title = "Metadata", value = if (!is.null(meta)) "Loaded" else "None",
              showcase = bsicons::bs_icon("table"),
              theme = if (!is.null(meta)) "success" else "light"
            ))
          ),
          card(
            card_header("Count Matrix Preview"),
            card_body(
              DT::DTOutput(ns("count_preview"))
            )
          ),
          if (!is.null(meta)) card(
            card_header("Metadata Preview"),
            card_body(DT::DTOutput(ns("meta_preview")))
          )
        )
      }
    })

    output$count_preview <- DT::renderDT({
      req(state$counts())
      DT::datatable(
        head(state$counts(), 10),
        options = list(scrollX = TRUE, pageLength = 10, dom = "t"),
        class   = "compact"
      )
    }, server = TRUE)

    output$meta_preview <- DT::renderDT({
      req(state$metadata())
      DT::datatable(state$metadata(),
                    options = list(pageLength = 10, dom = "t"),
                    class   = "compact")
    }, server = TRUE)
  })
}
