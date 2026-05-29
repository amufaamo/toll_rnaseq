#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`. DO NOT REMOVE.
#' @import shiny
#' @importFrom bslib bs_theme font_google bs_add_rules page_navbar nav_panel nav_spacer nav_item
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box input_dark_mode
#' @importFrom bsicons bs_icon
#' @noRd
app_ui <- function(request) {
  tagList(
    golem_add_external_resources(),
    page_navbar(
      title = tags$span(
        tags$img(src = "www/favicon.ico", height = "24px", class = "me-2"),
        "EasyRNA-Seq"
      ),
      id = "main_nav",
      theme = .app_theme(),
      # ── Sticky top-right controls ──────────────────────────────────────────
      nav_spacer(),
      nav_item(
        selectInput(
          "species", NULL,
          choices  = c("Human" = "human", "Mouse" = "mouse", "Rat" = "rat"),
          selected = "human",
          width    = "130px"
        )
      ),
      nav_item(input_dark_mode(id = "dark_mode")),

      # ── Pipeline tabs ──────────────────────────────────────────────────────
      nav_panel(
        title = .nav_label("1", "Upload", "upload"),
        value = "upload",
        mod_upload_ui("upload")
      ),
      nav_panel(
        title = .nav_label("2", "QC & Filter", "funnel"),
        value = "qc",
        mod_qc_ui("qc")
      ),
      nav_panel(
        title = .nav_label("3", "EDA", "bar-chart"),
        value = "eda",
        mod_eda_ui("eda")
      ),
      nav_panel(
        title = .nav_label("4", "DEG", "graph-up"),
        value = "deg",
        mod_deg_ui("deg")
      ),
      nav_panel(
        title = .nav_label("5", "Enrichment", "diagram-3"),
        value = "enrichment",
        mod_enrichment_ui("enrichment")
      ),
      nav_panel(
        title = .nav_label("6", "GSVA", "diagram-2"),
        value = "gsva",
        mod_gsva_ui("gsva")
      ),
      nav_panel(
        title = .nav_label("7", "UpSet", "list-check"),
        value = "upset",
        mod_deg_multi_ui("upset")
      ),
      nav_panel(
        title = .nav_label("8", "Report", "file-earmark-text"),
        value = "report",
        mod_report_ui("report")
      ),

      # ── Footer items ───────────────────────────────────────────────────────
      nav_item(
        actionButton(
          "btn_export_report", "Export Report",
          icon  = icon("download"),
          class = "btn btn-outline-secondary btn-sm"
        )
      ),
      nav_item(
        actionButton(
          "btn_reset", "Reset",
          icon  = icon("rotate-left"),
          class = "btn btn-outline-danger btn-sm"
        )
      )
    )
  )
}

# ── Helpers ──────────────────────────────────────────────────────────────────

.app_theme <- function() {
  bslib::bs_theme(
    version     = 5,
    primary     = "#2C3E50",
    secondary   = "#18BC9C",
    base_font   = bslib::font_google("Inter"),
    heading_font = bslib::font_google("Roboto Slab")
  ) |>
    bslib::bs_add_rules("
      .card {
        border: none !important;
        box-shadow: 0 4px 6px -1px rgba(0,0,0,.10);
        border-radius: .5rem;
      }
      .nav-link { font-size: .88rem; }
      .value-box .value-box-showcase { font-size: 1.4rem; }
    ")
}

#' Build a nav_panel title with step badge + Bootstrap icon
.nav_label <- function(step, label, icon_name) {
  tags$span(
    tags$span(step, class = "badge bg-secondary me-1"),
    bsicons::bs_icon(icon_name, class = "me-1"),
    label
  )
}

#' Empty-state card shown when a prerequisite step is not done yet
#' @noRd
ui_empty_state <- function(icon_name, title, message, action = NULL) {
  card(
    class = "text-center py-5",
    bsicons::bs_icon(icon_name, size = "3em", class = "text-muted mb-3"),
    tags$h5(title, class = "text-muted"),
    tags$p(message, class = "text-muted"),
    if (!is.null(action)) action
  )
}

#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path("www", app_sys("app/www"))
  tags$head(
    favicon(),
    bundle_resources(path = app_sys("app/www"), app_title = "EasyRNA-Seq")
  )
}
