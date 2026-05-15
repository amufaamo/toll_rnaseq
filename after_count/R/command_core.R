# command_core.R -- Phase 1 unified command layer + logging
# Provides:
#   execute_command(cmd, rv, session = NULL)
#   log_command(cmd, status, message, result_summary, rv, session_id = NULL)
#   replay_commands(file, rv, session = NULL)
#   COMMAND_HANDLERS  (named list, extensible)
#   `%||%`            (nullish fallback)
#
# History persistence:
#   Default dir: /tmp/toll_history (or Sys.getenv("TOLL_HISTORY_DIR"))
#   File format: {session_id}.jsonl (one JSON object per line)

# -------- helpers --------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

.toll_history_dir <- function() {
  d <- Sys.getenv("TOLL_HISTORY_DIR", unset = "/tmp/toll_history")
  if (!nzchar(d)) d <- "/tmp/toll_history"
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  d
}

# -------- handler registry --------
# Extend by: COMMAND_HANDLERS$my_op <- function(cmd, rv, session) { ... }
COMMAND_HANDLERS <- list(
  noop = function(cmd, rv, session) {
    list(
      status = "success",
      message = as.character(cmd$message %||% "noop executed"),
      result_summary = list()
    )
  }
)

# -------- log_command --------
# Appends one row to rv$command_log (tibble) and one line to history JSONL.
log_command <- function(cmd, status, message, result_summary, rv, session_id = NULL) {
  if (!requireNamespace("tibble", quietly = TRUE)) return(invisible(NULL))
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(invisible(NULL))

  # rv may be a reactiveValues (needs isolate) or plain list
  safe_get <- function(x, name) {
    tryCatch(
      if (inherits(x, "reactivevalues")) shiny::isolate(x[[name]]) else x[[name]],
      error = function(e) NULL
    )
  }
  sid <- session_id %||% safe_get(rv, "session_id") %||% "unknown"
  ts  <- Sys.time()

  params_json <- tryCatch(
    as.character(jsonlite::toJSON(cmd, auto_unbox = TRUE, null = "null", force = TRUE)),
    error = function(e) "{}"
  )
  result_json <- tryCatch(
    as.character(jsonlite::toJSON(result_summary, auto_unbox = TRUE, null = "null", force = TRUE)),
    error = function(e) "{}"
  )

  row <- tibble::tibble(
    timestamp   = ts,
    session_id  = as.character(sid),
    op          = as.character(cmd$op %||% ""),
    status      = as.character(status %||% ""),
    message     = as.character(message %||% ""),
    params_json = params_json,
    result_json = result_json
  )

  # append to in-memory log (reactiveValues or plain list both OK)
  cur <- safe_get(rv, "command_log")
  if (is.null(cur) || (is.data.frame(cur) && nrow(cur) == 0 && ncol(cur) == 0)) {
    rv$command_log <- row
  } else {
    rv$command_log <- if (requireNamespace("dplyr", quietly = TRUE)) {
      dplyr::bind_rows(cur, row)
    } else {
      rbind(cur, row)
    }
  }

  # append to JSONL file (fire-and-forget; warn but never stop)
  tryCatch({
    d <- .toll_history_dir()
    fp <- file.path(d, paste0(sid, ".jsonl"))
    line <- jsonlite::toJSON(
      list(
        timestamp   = format(ts, "%Y-%m-%dT%H:%M:%OS3%z"),
        session_id  = as.character(sid),
        op          = as.character(cmd$op %||% ""),
        status      = as.character(status %||% ""),
        message     = as.character(message %||% ""),
        params_json = params_json,
        result_json = result_json
      ),
      auto_unbox = TRUE, force = TRUE
    )
    cat(line, "\n", sep = "", file = fp, append = TRUE)
  }, error = function(e) {
    warning("log_command: failed to write history file: ", conditionMessage(e))
  })

  invisible(row)
}

# -------- execute_command --------
# Dispatches cmd$op to COMMAND_HANDLERS[[op]]; always logs outcome.
execute_command <- function(cmd, rv, session = NULL) {
  if (!is.list(cmd)) {
    res <- list(status = "error",
                message = "cmd must be a list",
                result_summary = list())
    log_command(list(op = "<invalid>"), res$status, res$message, res$result_summary, rv)
    return(res)
  }

  op <- as.character(cmd$op %||% "")
  if (!nzchar(op)) {
    res <- list(status = "error",
                message = "missing op",
                result_summary = list())
    log_command(cmd, res$status, res$message, res$result_summary, rv)
    return(res)
  }

  handler <- COMMAND_HANDLERS[[op]]
  if (is.null(handler)) {
    res <- list(status = "error",
                message = paste0("Unknown op: ", op),
                result_summary = list())
    log_command(cmd, res$status, res$message, res$result_summary, rv)
    return(res)
  }

  res <- tryCatch(
    handler(cmd, rv, session),
    error = function(e) list(
      status = "error",
      message = conditionMessage(e),
      result_summary = list()
    )
  )

  if (!is.list(res) || is.null(res$status)) {
    res <- list(status = "error",
                message = "handler returned invalid result",
                result_summary = list())
  }
  res$message        <- res$message        %||% ""
  res$result_summary <- res$result_summary %||% list()

  log_command(cmd, res$status, res$message, res$result_summary, rv)
  res
}

# -------- replay_commands --------
# Accepts .yaml/.yml/.jsonl/.json files. Returns summary list.
replay_commands <- function(file, rv, session = NULL) {
  if (!file.exists(file)) stop("replay_commands: file not found: ", file)
  ext <- tolower(tools::file_ext(file))

  cmds <- if (ext %in% c("yaml", "yml")) {
    if (!requireNamespace("yaml", quietly = TRUE))
      stop("yaml package required for YAML replay")
    y <- yaml::read_yaml(file)
    # single command or list-of-commands
    if (is.list(y) && !is.null(y$op)) list(y) else y
  } else if (ext %in% c("jsonl", "json")) {
    if (!requireNamespace("jsonlite", quietly = TRUE))
      stop("jsonlite package required for JSONL replay")
    lines <- readLines(file, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    lapply(lines, function(l) jsonlite::fromJSON(l, simplifyVector = FALSE))
  } else {
    stop("replay_commands: unsupported extension: ", ext)
  }

  if (!is.list(cmds)) cmds <- list(cmds)

  results <- vector("list", length(cmds))
  n_success <- 0L
  n_error   <- 0L
  for (i in seq_along(cmds)) {
    r <- tryCatch(
      execute_command(cmds[[i]], rv, session),
      error = function(e) list(
        status = "error",
        message = conditionMessage(e),
        result_summary = list()
      )
    )
    results[[i]] <- r
    if (isTRUE(r$status == "success")) {
      n_success <- n_success + 1L
    } else {
      n_error <- n_error + 1L
    }
  }

  list(
    n_total   = length(cmds),
    n_success = n_success,
    n_error   = n_error,
    results   = results
  )
}

# -------- Phase 1 demo handler: update_gene_choices --------
# Used from module_processing.R observeEvent(display_table_reactive(), ...) to
# demonstrate that UI-side effects can also flow through the command layer.
# Parameters via cmd:
#   cmd$gene_choices  (character vector)
#   cmd$input_id      (character; selectize input id, e.g. "genes_to_plot")
COMMAND_HANDLERS$update_gene_choices <- function(cmd, rv, session) {
  choices <- cmd$gene_choices %||% character()
  input_id <- as.character(cmd$input_id %||% "genes_to_plot")
  if (is.null(session)) {
    return(list(
      status = "skipped",
      message = "no session; choices not applied",
      result_summary = list(n_choices = length(choices), input_id = input_id)
    ))
  }
  tryCatch({
    shiny::updateSelectizeInput(session, input_id,
      choices = choices,
      server = TRUE
    )
    list(
      status = "success",
      message = sprintf("updated %s with %d choices", input_id, length(choices)),
      result_summary = list(n_choices = length(choices), input_id = input_id)
    )
  }, error = function(e) list(
    status = "error",
    message = conditionMessage(e),
    result_summary = list(input_id = input_id)
  ))
}

# ---- end command_core.R ----
