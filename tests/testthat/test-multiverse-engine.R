make_efdr_fixture <- function() {
  # Synthetic multiverse outcomes: observed calls 10, 5, 2 and null calls 8, 3, 0.
  # This gives strictly decreasing raw eFDR 0.45, 0.40, 0.25 at tau 0, .5, 1.
  mk <- function(stability) data.frame(stability = stability, passes_direction = TRUE)
  observed <- list(tested = matrix(TRUE, 10, 2), stability = mk(c(rep(1, 2), rep(.5, 3), rep(0, 5))))
  bootstrap <- list(mk(c(rep(.5, 3), rep(0, 5))))
  list(observed = observed, bootstrap = bootstrap)
}

test_that("eFDR q-value sweep preserves an already monotone raw curve", {
  f <- make_efdr_fixture()
  curve <- .mv_efdr_curve(f$observed, f$bootstrap, target = .30)$curve
  expect_equal(curve$efdr, curve$efdr_raw)
})

test_that("eFDR declines with increasing stability and selects a nonzero tau", {
  f <- make_efdr_fixture()
  ans <- .mv_efdr_curve(f$observed, f$bootstrap, target = .30)
  curve <- ans$curve
  good <- !is.na(curve$efdr)
  expect_true(all(diff(curve$efdr[good]) <= 0))
  expect_gt(curve$efdr[1], tail(curve$efdr[good], 1))
  expect_gt(ans$selected_tau, 0)
  expect_gt(curve$efdr[1], .10)
})

test_that("per-gene eFDR varies and float-safe indexing covers all stable genes", {
  f <- make_efdr_fixture()
  curve <- .mv_efdr_curve(f$observed, f$bootstrap, target = .30)$curve
  stability <- data.frame(gene = c("g1", "g2", "g3"), stability_up = c(0, .5, 1),
    stability_down = 0, stability = c(0 + 2e-16, .5, 1), direction_consistency = 1,
    direction = "Up", passes_direction = TRUE)
  attached <- mv_attach_efdr(stability, curve, target = .30)
  expect_true(all(!is.na(attached$efdr[attached$passes_direction])))
  expect_gt(length(unique(attached$efdr)), 1)
  expect_true(all(diff(attached$efdr[order(attached$stability)]) <= 0))
  stats <- array(NA_real_, dim = c(3, 2, 2), dimnames = list(stability$gene, c("a", "b"), c("baseMean", "log2FoldChange")))
  stats[, , "baseMean"] <- 100; stats[, , "log2FoldChange"] <- 2
  handoff <- mv_make_handoff(list(stats = stats), attached)
  expect_gt(length(unique(handoff$padj)), 1)
})

test_that("stability algebra is signed and monotone in tau", {
  stats <- array(NA_real_, dim = c(3, 2, 2), dimnames = list(c("all", "none", "mixed"), c("a", "b"), c("log2FoldChange", "padj")))
  stats["all", , "log2FoldChange"] <- 2; stats["all", , "padj"] <- .01
  stats["none", , "log2FoldChange"] <- 0; stats["none", , "padj"] <- 1
  stats["mixed", , "log2FoldChange"] <- c(2, -2); stats["mixed", , "padj"] <- .01
  x <- mv_compute_stability(stats, matrix(TRUE, 3, 2, dimnames = dimnames(stats)[1:2]))
  expect_equal(x$stability[x$gene == "all"], 1)
  expect_equal(x$stability[x$gene == "none"], 0)
  expect_false(x$passes_direction[x$gene == "mixed"])
  expect_lte(sum(x$stability >= 1), sum(x$stability >= .5))
})

test_that("full-model null reconstruction is exact on the DESeq2 count scale", {
  skip_mv_deps()
  d <- make_mv_counts(n_genes = 160, n_de = 20)
  fit <- mv_fit_full_null(d$counts, d$metadata, "condition", "C", "T")
  expect_true(fit$reconstruction_ok)
  expect_equal(dim(fit$mu0), dim(d$counts))
  expect_true(all(is.finite(fit$dispersions)))
})

test_that("an all-zero-count gene does not poison the null bootstrap", {
  skip_mv_deps()
  d <- make_mv_counts(seed = 5, n_genes = 60, n_de = 10)
  d$counts["g1", ] <- 0L
  fit <- mv_fit_full_null(d$counts, d$metadata, "condition", "C", "T")
  expect_false(anyNA(fit$mu0))
  expect_false(anyNA(fit$dispersions))
  expect_equal(unname(fit$mu0["g1", ]), rep(0, ncol(fit$mu0)))
  sim <- mv_simulate_null(fit, seed = 1)
  expect_false(anyNA(sim))
  expect_true(all(sim["g1", ] == 0))
  # The full engine, not just the two building blocks above, must also complete.
  run <- mv_run_multiverse(d$counts, d$metadata, "condition", "C", "T", B = 2, seed = 3,
                           include_apeglm = FALSE)
  expect_false(anyNA(run$efdr_curve$efdr_raw[run$efdr_curve$observed_calls > 0]))
})

test_that("mv_bootstrap_efdr()/mv_run_multiverse() select tau against the requested target, not a hardcoded 0.10", {
  skip_mv_deps()
  d <- make_mv_counts(seed = 7, n_genes = 200, n_de = 20)
  obs <- mv_run_observed(d$counts, d$metadata, "condition", "C", "T", include_apeglm = FALSE)
  obs$stability <- mv_compute_stability(obs$stats, obs$tested, 0.05, 1, 0.90)
  cal_loose <- mv_bootstrap_efdr(obs, B = 5, seed = 1, target_efdr = 0.50)
  cal_strict <- mv_bootstrap_efdr(obs, B = 5, seed = 1, target_efdr = 0.01)
  # A stricter target can only select an equal-or-higher stability threshold.
  expect_gte(cal_strict$selected_tau, cal_loose$selected_tau)
  expect_gt(cal_strict$selected_tau, cal_loose$selected_tau)
})

test_that("multiverse engine is deterministic and produces a downstream contract", {
  skip_mv_deps()
  d <- make_mv_counts(seed = 11, n_genes = 180, n_de = 30)
  a <- mv_run_multiverse(d$counts, d$metadata, "condition", "C", "T", B = 2, seed = 99,
                         include_apeglm = FALSE)
  b <- mv_run_multiverse(d$counts, d$metadata, "condition", "C", "T", B = 2, seed = 99,
                         include_apeglm = FALSE)
  expect_equal(a$efdr_curve, b$efdr_curve)
  expect_true(all(c("gene", "efdr", "called") %in% names(a$stability)))
  handoff <- mv_make_handoff(a$run, a$stability)
  expect_true(all(c("gene", "baseMean", "log2FoldChange", "stat", "pvalue", "padj") %in% names(handoff)))
  expect_true(all(is.na(handoff$pvalue)))
  expect_true(all(is.finite(handoff$stat)))
})

test_that("pure null is conservative across fixed seeds", {
  skip_mv_deps()
  # Three seeds allow rare discrete null excursions; <=2 calls is a practical bound for 180 genes/B=5.
  calls <- numeric(3); selected <- numeric(3); means <- numeric(3)
  for (i in seq_along(c(31, 32, 33))) {
    d <- make_mv_counts(seed = c(31, 32, 33)[i], n_genes = 180, n_de = 0)
    z <- mv_run_multiverse(d$counts, d$metadata, "condition", "C", "T", B = 5,
                           seed = 100 + i, include_apeglm = FALSE)
    calls[i] <- sum(z$stability$called)
    selected[i] <- if (all(is.na(z$efdr_curve$efdr))) 1 else min(z$efdr_curve$efdr, na.rm = TRUE)
    means[i] <- mean(z$stability$stability)
  }
  expect_gte(sum(calls == 0), 2)
  expect_lte(max(calls), 2)
  expect_gte(min(selected), .10)
  expect_lt(max(means), .05)
})

test_that("AppState can receive a multiverse result outside a Shiny server", {
  skip_mv_deps()
  d <- make_mv_counts(seed = 44, n_genes = 160, n_de = 20)
  z <- mv_run_multiverse(d$counts, d$metadata, "condition", "C", "T", B = 1, seed = 7,
                         include_apeglm = FALSE)
  state <- AppState$new()
  shiny::isolate({
    mv_apply_run(state, z$run, z$stability, z$efdr_curve, z$target_efdr, z$seed, z$B)
    expect_equal(state$pipeline_status()$deg, "done")
    expect_equal(state$pipeline_status()$deg_multiverse, "done")
    expect_true(grepl("^multiverse:", names(state$deg_results())[1]))
  })
})

test_that("known-truth eFDR calibration regression", {
  skip_mv_deps()
  # The FDP <= 0.20 tolerance did not expose the flattened-curve bug: it assesses
  # realised calls, not whether tau/eFDR ranking is internally calibrated.
  skip_if(Sys.getenv("EASYRNASEQ_RUN_SLOW") == "", "set EASYRNASEQ_RUN_SLOW=1 for B=50 validation")
  d <- make_mv_counts(n_genes = 1000, n_de = 120)
  z <- mv_run_multiverse(d$counts, d$metadata, "condition", "C", "T", B = 50,
                         seed = 91, include_apeglm = FALSE)
  called <- z$stability$gene[z$stability$called]
  fdp <- if (length(called)) mean(called %in% d$null_genes) else 0
  expect_lte(fdp, .20)
})
