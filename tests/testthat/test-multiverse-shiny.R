# Integration tests for the Multiverse DEG Shiny module.
#
# shinytest2/chromote are unavailable in this environment, so the module is driven
# through shiny::testServer() instead: this exercises the real observeEvent wiring,
# the ExtendedTask scheduler that chains the null-bootstrap replicates, and the
# writes back into AppState -- none of which the engine unit tests touch.

# testServer never runs the event loop, so promises resolve only when we pump
# later::run_now() and flush reactives ourselves.
pump_until <- function(session, done, timeout = 900) {
  deadline <- Sys.time() + timeout
  repeat {
    if (isTRUE(done())) return(TRUE)
    if (Sys.time() > deadline) return(FALSE)
    later::run_now(0.05)
    session$flushReact()
  }
}

test_that("module server runs a multiverse end-to-end and populates AppState", {
  skip_mv_deps()
  d <- make_mv_counts(seed = 21, n_genes = 120, n_de = 20)
  state <- AppState$new()

  shiny::testServer(mod_deg_multiverse_server, args = list(state = state), {
    state$counts(d$counts)
    state$metadata(d$metadata)
    state$set_status("upload", "done")
    session$flushReact()

    session$setInputs(
      condition_col = "condition", ref_level = "C", test_level = "T",
      use_covariate = FALSE, use_apeglm = FALSE,
      bootstrap_B = 2, seed = 5,
      padj_cutoff = 0.05, lfc_cutoff = 1, target_efdr = 0.10
    )
    session$setInputs(run_multiverse = 1)

    # The run is asynchronous: observed grid first, then B chained bootstraps.
    expect_equal(state$pipeline_status()$deg_multiverse, "running")
    finished <- pump_until(session, function() !is.null(state$deg_multiverse()))
    expect_true(finished)

    res <- state$deg_multiverse()
    expect_equal(state$pipeline_status()$deg_multiverse, "done")
    expect_equal(res$metadata$B, 2L)
    expect_equal(res$metadata$target, 0.10)
    expect_equal(nrow(res$stability), nrow(d$counts))
    expect_true(all(c("gene", "stability", "efdr", "called") %in% names(res$stability)))
    expect_true(nrow(res$run$specifications) > 0)

    # Progress record reaches the terminal stage with every replicate accounted for.
    prog <- state$deg_multiverse_progress()
    expect_equal(prog$stage, "done")
    expect_equal(prog$completed, 2L)

    # Downstream hand-off contract for Enrichment / UpSet.
    key <- "multiverse:T_vs_C"
    expect_true(key %in% names(state$deg_results()))
    expect_equal(state$pipeline_status()$deg, "done")

    # Outputs render without error.
    ui <- output$main_content
    expect_true(is.list(ui) && "html" %in% names(ui))
    expect_match(as.character(ui$html), "Specification Curve", fixed = TRUE)
    expect_match(as.character(ui$html), "Stability-ranked DEG table", fixed = TRUE)

    session$setInputs(gene = res$stability$gene[1])
    plt <- output$curve
    expect_true(is.list(plt) && nzchar(plt$src))
    expect_s3_class(output$table, "json")
  })
})

test_that("bootstrap payload carries the call rule and seed of the run", {
  # Regression guard: the scheduler used to hand the whole job object to the
  # worker, so seed/padj_cutoff/lfc_cutoff (which live under job$payload)
  # arrived as NULL and set.seed() aborted every bootstrap replicate.
  job <- list(
    run = "RUN", null_fit = "FIT", b = 3L, boot = list("already done"),
    payload = list(seed = 77L, padj_cutoff = 0.10, lfc_cutoff = 0, B = 5L)
  )
  out <- .mv_boot_payload(job)
  expect_identical(out$run, "RUN")
  expect_identical(out$null_fit, "FIT")
  expect_identical(out$b, 3L)
  expect_identical(out$seed, 77L)
  expect_identical(out$padj_cutoff, 0.10)
  expect_identical(out$lfc_cutoff, 0)
  # Accumulated replicates must not be shipped to the worker.
  expect_false("boot" %in% names(out))
})

test_that("module runs are reproducible and respond to the call rule", {
  skip_mv_deps()
  d <- make_mv_counts(seed = 22, n_genes = 120, n_de = 20)

  run_with <- function(padj_cutoff, lfc_cutoff, seed) {
    state <- AppState$new()
    shiny::testServer(mod_deg_multiverse_server, args = list(state = state), {
      state$counts(d$counts)
      state$metadata(d$metadata)
      state$set_status("upload", "done")
      session$flushReact()
      session$setInputs(
        condition_col = "condition", ref_level = "C", test_level = "T",
        use_covariate = FALSE, use_apeglm = FALSE,
        bootstrap_B = 2, seed = seed,
        padj_cutoff = padj_cutoff, lfc_cutoff = lfc_cutoff, target_efdr = 0.10
      )
      session$setInputs(run_multiverse = 1)
      expect_true(pump_until(session, function() !is.null(state$deg_multiverse())))
    })
    shiny::isolate(state$deg_multiverse())
  }

  a <- run_with(0.05, 1, 5)
  b <- run_with(0.05, 1, 5)
  expect_equal(a$efdr_curve, b$efdr_curve)
  expect_equal(a$stability, b$stability)

  # Loosening the rule must change the stability actually stored for the user.
  loose <- run_with(0.10, 0, 5)
  expect_false(identical(loose$stability$stability, a$stability$stability))
  expect_gte(sum(loose$stability$stability), sum(a$stability$stability))
})

test_that("a covariate equal to the condition is rejected before any compute", {
  skip_mv_deps()
  d <- make_mv_counts(seed = 23, n_genes = 60, n_de = 10)
  state <- AppState$new()

  shiny::testServer(mod_deg_multiverse_server, args = list(state = state), {
    state$counts(d$counts)
    state$metadata(d$metadata)
    state$set_status("upload", "done")
    session$flushReact()
    session$setInputs(
      condition_col = "condition", ref_level = "C", test_level = "T",
      use_covariate = TRUE, covariate = "condition", use_apeglm = FALSE,
      bootstrap_B = 2, seed = 5,
      padj_cutoff = 0.05, lfc_cutoff = 1, target_efdr = 0.10
    )
    session$setInputs(run_multiverse = 1)
    expect_equal(state$pipeline_status()$deg_multiverse, "pending")
    expect_null(state$deg_multiverse())
  })
})

test_that("a failing run surfaces an error and leaves the step not done", {
  skip_mv_deps()
  d <- make_mv_counts(seed = 24, n_genes = 60, n_de = 10)
  state <- AppState$new()

  shiny::testServer(mod_deg_multiverse_server, args = list(state = state), {
    state$counts(d$counts)
    state$metadata(d$metadata)
    state$set_status("upload", "done")
    session$flushReact()
    session$setInputs(
      condition_col = "condition", ref_level = "C", test_level = "NOT_A_LEVEL",
      use_covariate = FALSE, use_apeglm = FALSE,
      bootstrap_B = 2, seed = 5,
      padj_cutoff = 0.05, lfc_cutoff = 1, target_efdr = 0.10
    )
    suppressWarnings({
      session$setInputs(run_multiverse = 1)
      finished <- pump_until(
        session,
        function() identical(state$deg_multiverse_progress()$stage, "error"),
        timeout = 300
      )
    })
    expect_true(finished)
    expect_equal(state$pipeline_status()$deg_multiverse, "pending")
    expect_null(state$deg_multiverse())
  })
})

test_that("reset leaves a status for every step the nav badges render", {
  # app_server switch()es on pipeline_status()[[step]] for each of these, so a
  # missing entry crashes the badge observer as soon as Reset is confirmed.
  steps <- c("upload", "qc", "eda", "deg", "deg_multiverse", "enrichment",
             "gsva", "timeseries", "report")
  shiny::isolate({
    state <- AppState$new()
    expect_true(all(steps %in% names(state$pipeline_status())))
    state$set_status("timeseries", "done")
    state$reset()
    expect_true(all(steps %in% names(state$pipeline_status())))
    expect_true(all(vapply(state$pipeline_status()[steps],
                           function(x) identical(x, "pending"), logical(1))))
  })
})
