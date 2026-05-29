#' R6 AppState: Centralized Reactive State Manager
#'
#' Replaces scattered reactiveValues with a single R6 object.
#' All Shiny reactivity is handled via reactiveVal wrappers so this
#' class can also be instantiated and tested outside a Shiny session.
#'
#' @import R6
#' @import shiny
#' @noRd
AppState <- R6::R6Class(
  "AppState",
  public = list(

    # ── Reactive fields ──────────────────────────────────────────────────────

    #' @field counts Raw count matrix (genes x samples), NULL until uploaded
    counts = NULL,

    #' @field metadata Sample metadata data.frame
    metadata = NULL,

    #' @field counts_filtered Filtered count matrix after QC
    counts_filtered = NULL,

    #' @field deg_results List of DESeq2 results (one per contrast)
    deg_results = NULL,

    #' @field gsea_results List of fgsea/enrichplot results
    gsea_results = NULL,

    #' @field params_log Named list recording every parameter change for
    #'   reproducibility report generation
    params_log = NULL,

    #' @field species Organism identifier ("human", "mouse", "rat")
    species = NULL,

    #' @field pipeline_status Named list: step name → "pending"/"running"/"done"
    pipeline_status = NULL,

    # ── Constructor ──────────────────────────────────────────────────────────

    initialize = function() {
      self$counts          <- reactiveVal(NULL)
      self$metadata        <- reactiveVal(NULL)
      self$counts_filtered <- reactiveVal(NULL)
      self$deg_results     <- reactiveVal(list())
      self$gsea_results    <- reactiveVal(list())
      self$params_log      <- reactiveVal(list())
      self$species         <- reactiveVal("human")
      self$pipeline_status <- reactiveVal(list(
        upload      = "pending",
        qc          = "pending",
        eda         = "pending",
        deg         = "pending",
        enrichment  = "pending",
        gsva        = "pending",
        timeseries  = "pending",
        report      = "pending"
      ))
    },

    # ── Helpers ───────────────────────────────────────────────────────────────

    #' Record a parameter into params_log for reproducibility tracking
    log_param = function(step, key, value) {
      current <- self$params_log()
      entry <- list(
        step      = step,
        key       = key,
        value     = value,
        timestamp = Sys.time()
      )
      current[[paste0(step, ".", key)]] <- entry
      self$params_log(current)
      invisible(self)
    },

    #' Update pipeline_status for a given step
    set_status = function(step, status) {
      current <- self$pipeline_status()
      current[[step]] <- status
      self$pipeline_status(current)
      invisible(self)
    },

    #' TRUE if a given step's data is available
    step_ready = function(step) {
      status <- isolate(self$pipeline_status())
      isTRUE(status[[step]] == "done")
    },

    #' Reset all state (used by "Reset Analysis" button)
    reset = function() {
      self$counts(NULL)
      self$metadata(NULL)
      self$counts_filtered(NULL)
      self$deg_results(list())
      self$gsea_results(list())
      self$params_log(list())
      self$species("human")
      self$pipeline_status(list(
        upload     = "pending",
        qc         = "pending",
        eda        = "pending",
        deg        = "pending",
        enrichment = "pending",
        gsva       = "pending",
        report     = "pending"
      ))
      invisible(self)
    }
  )
)
