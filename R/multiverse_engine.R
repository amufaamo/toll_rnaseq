# Multiverse DEG engine -- pure compute functions

.mv_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("%s is required for multiverse DEG analysis.", pkg), call. = FALSE)
}

.mv_align_inputs <- function(counts, metadata, condition_col, ref_level, test_level) {
  if (!is.matrix(counts) || any(counts < 0) || any(counts != round(counts)))
    stop("Counts must be a non-negative integer matrix.", call. = FALSE)
  if (!condition_col %in% colnames(metadata))
    stop("Condition column is not present in metadata.", call. = FALSE)
  sample_col <- colnames(metadata)[1]
  samples <- intersect(as.character(metadata[[sample_col]]), colnames(counts))
  if (length(samples) == 0) stop("No sample names match between counts and metadata.", call. = FALSE)
  meta <- metadata[match(samples, metadata[[sample_col]]), , drop = FALSE]
  rownames(meta) <- samples
  mat <- counts[, samples, drop = FALSE]
  meta[[condition_col]] <- stats::relevel(factor(meta[[condition_col]]), ref = ref_level)
  if (!test_level %in% levels(meta[[condition_col]]))
    stop("Test level is not present in the condition column.", call. = FALSE)
  if (sum(meta[[condition_col]] == ref_level) < 2 || sum(meta[[condition_col]] == test_level) < 2)
    stop("Multiverse DEG requires at least two samples in each contrast group.", call. = FALSE)
  list(counts = mat, metadata = meta)
}

.mv_design <- function(metadata, condition_col, covariate = NULL) {
  if (!is.null(covariate) && !is.na(covariate) && nzchar(covariate)) {
    if (!covariate %in% colnames(metadata)) stop("Covariate is not in metadata.", call. = FALSE)
    if (anyNA(metadata[[covariate]]) || length(unique(metadata[[covariate]])) < 2)
      stop("Covariate must have at least two non-missing values.", call. = FALSE)
    metadata[[covariate]] <- factor(metadata[[covariate]])
    form <- stats::as.formula(paste("~", covariate, "+", condition_col))
  } else {
    form <- stats::as.formula(paste("~", condition_col))
  }
  mm <- stats::model.matrix(form, metadata)
  if (qr(mm)$rank < ncol(mm)) stop("Condition and covariate are confounded or rank deficient.", call. = FALSE)
  list(formula = form, matrix = mm, metadata = metadata)
}

#' Build the reduced MVP specification grid
#' @noRd
mv_build_specifications <- function(counts, metadata, condition_col, ref_level, test_level,
                                    covariate = NULL, include_apeglm = TRUE) {
  x <- .mv_align_inputs(counts, metadata, condition_col, ref_level, test_level)
  covariates <- c(NA_character_, if (!is.null(covariate) && nzchar(covariate)) covariate)
  paths <- c("deseq_wald", "edger_ql", "edger_lrt")
  if (isTRUE(include_apeglm) && requireNamespace("apeglm", quietly = TRUE))
    paths <- c(paths, "deseq_wald_apeglm")
  specs <- expand.grid(
    filter = c("min_count", "filterByExpr"),
    covariate = covariates,
    method = paths,
    stringsAsFactors = FALSE
  )
  specs$id <- sprintf("spec_%03d", seq_len(nrow(specs)))
  specs$normalization <- ifelse(grepl("deseq", specs$method), "median_ratio", "TMM")
  specs$shrinkage <- specs$method == "deseq_wald_apeglm"
  # Validate the covariate arms now so invalid combinations are pruned before compute.
  keep <- vapply(seq_len(nrow(specs)), function(i) {
    tryCatch({ .mv_design(x$metadata, condition_col, specs$covariate[i]); TRUE }, error = function(e) FALSE)
  }, logical(1))
  skipped <- specs[!keep, , drop = FALSE]
  if (nrow(skipped)) skipped$skip_reason <- "rank-deficient covariate design"
  list(specifications = specs[keep, , drop = FALSE], skipped = skipped)
}

.mv_keep_genes <- function(counts, metadata, condition_col, spec) {
  group <- metadata[[condition_col]]
  if (spec$filter == "min_count") {
    n_min <- min(table(group))
    rowSums(counts >= 10) >= n_min
  } else {
    .mv_require("edgeR")
    des <- .mv_design(metadata, condition_col, spec$covariate)$matrix
    edgeR::filterByExpr(edgeR::DGEList(counts = counts), design = des)
  }
}

.mv_condition_coef <- function(coef_names, condition_col, test_level, ref_level) {
  wanted <- paste0(condition_col, "_", gsub("[^A-Za-z0-9]", ".", test_level),
                   "_vs_", gsub("[^A-Za-z0-9]", ".", ref_level))
  hit <- grep(wanted, coef_names, value = TRUE, fixed = TRUE)
  if (length(hit) != 1) stop("Could not identify the condition coefficient.", call. = FALSE)
  hit
}

.mv_empty_result <- function(genes) {
  data.frame(gene = genes, baseMean = NA_real_, log2FoldChange = NA_real_,
             lfcSE = NA_real_, stat = NA_real_, pvalue = NA_real_, padj = NA_real_,
             stringsAsFactors = FALSE)
}

#' Run one legal MVP DEG specification
#' @noRd
mv_run_specification <- function(counts, metadata, condition_col, ref_level, test_level, spec) {
  keep <- .mv_keep_genes(counts, metadata, condition_col, spec)
  out <- .mv_empty_result(rownames(counts))
  if (sum(keep) < 10) return(list(result = out, tested = keep, error = "fewer than 10 genes retained"))
  mat <- counts[keep, , drop = FALSE]
  des <- .mv_design(metadata, condition_col, spec$covariate)
  fit <- tryCatch({
    if (grepl("^deseq", spec$method)) {
      .mv_require("DESeq2")
      dds <- DESeq2::DESeqDataSetFromMatrix(mat, des$metadata, des$formula)
      dds <- DESeq2::DESeq(dds, test = "Wald", quiet = TRUE)
      contrast <- c(condition_col, test_level, ref_level)
      res <- DESeq2::results(dds, contrast = contrast)
      if (identical(spec$method, "deseq_wald_apeglm")) {
        coef_name <- .mv_condition_coef(DESeq2::resultsNames(dds), condition_col, test_level, ref_level)
        res <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
      }
      data.frame(gene = rownames(res), as.data.frame(res), check.names = FALSE)
    } else {
      .mv_require("edgeR")
      y <- edgeR::DGEList(counts = mat)
      y <- edgeR::calcNormFactors(y, method = "TMM")
      y <- edgeR::estimateDisp(y, des$matrix, robust = TRUE)
      coef_index <- ncol(des$matrix)
      if (identical(spec$method, "edger_ql")) {
        qfit <- edgeR::glmQLFit(y, des$matrix, robust = TRUE)
        tt <- edgeR::glmQLFTest(qfit, coef = coef_index)$table
        data.frame(gene = rownames(tt), baseMean = rowMeans(edgeR::cpm(y, log = FALSE)),
                   log2FoldChange = tt$logFC, lfcSE = NA_real_, stat = tt$F,
                   pvalue = tt$PValue, padj = stats::p.adjust(tt$PValue, "BH"))
      } else {
        gfit <- edgeR::glmFit(y, des$matrix)
        tt <- edgeR::glmLRT(gfit, coef = coef_index)$table
        data.frame(gene = rownames(tt), baseMean = rowMeans(edgeR::cpm(y, log = FALSE)),
                   log2FoldChange = tt$logFC, lfcSE = NA_real_, stat = tt$LR,
                   pvalue = tt$PValue, padj = stats::p.adjust(tt$PValue, "BH"))
      }
    }
  }, error = function(e) e)
  if (inherits(fit, "error")) return(list(result = out, tested = keep, error = conditionMessage(fit)))
  cols <- intersect(colnames(out), colnames(fit))
  out[match(fit$gene, out$gene), cols] <- fit[, cols, drop = FALSE]
  list(result = out, tested = keep, error = NULL)
}

#' Run observed multiverse and retain per-specification statistics
#' @noRd
mv_run_observed <- function(counts, metadata, condition_col, ref_level, test_level,
                            covariate = NULL, include_apeglm = TRUE) {
  aligned <- .mv_align_inputs(counts, metadata, condition_col, ref_level, test_level)
  grid <- mv_build_specifications(aligned$counts, aligned$metadata, condition_col, ref_level,
                                  test_level, covariate, include_apeglm)
  specs <- grid$specifications
  if (!nrow(specs)) stop("No legal multiverse specifications remain.", call. = FALSE)
  metrics <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  stats_array <- array(NA_real_, dim = c(nrow(aligned$counts), nrow(specs), length(metrics)),
                       dimnames = list(rownames(aligned$counts), specs$id, metrics))
  tested <- matrix(FALSE, nrow(aligned$counts), nrow(specs),
                   dimnames = list(rownames(aligned$counts), specs$id))
  errors <- character(nrow(specs))
  for (i in seq_len(nrow(specs))) {
    ans <- mv_run_specification(aligned$counts, aligned$metadata, condition_col, ref_level,
                                test_level, specs[i, , drop = FALSE])
    stats_array[, i, metrics] <- as.matrix(ans$result[, metrics])
    tested[, i] <- ans$tested & is.finite(ans$result$pvalue)
    errors[i] <- if (is.null(ans$error)) "" else ans$error
  }
  specs$error <- errors
  list(counts = aligned$counts, metadata = aligned$metadata, specifications = specs,
       skipped = grid$skipped, stats = stats_array, tested = tested,
       contrast = list(condition_col = condition_col, ref = ref_level, test = test_level,
                       covariate = covariate))
}

#' Calculate signed stability from stored statistics
#' @noRd
mv_compute_stability <- function(stats, tested, padj_cutoff = 0.05, lfc_cutoff = 1,
                                 direction_cutoff = 0.90) {
  lfc <- stats[, , "log2FoldChange"]
  padj <- stats[, , "padj"]
  up <- tested & !is.na(padj) & padj <= padj_cutoff & lfc >= lfc_cutoff
  down <- tested & !is.na(padj) & padj <= padj_cutoff & lfc <= -lfc_cutoff
  up_s <- rowMeans(up); down_s <- rowMeans(down)
  stability <- pmax(up_s, down_s)
  consistency <- ifelse(up_s + down_s > 0, stability / (up_s + down_s), NA_real_)
  direction <- ifelse(up_s >= down_s, "Up", "Down")
  data.frame(gene = rownames(stats), stability_up = up_s, stability_down = down_s,
             stability = stability, direction_consistency = consistency,
             direction = direction, passes_direction = !is.na(consistency) & consistency >= direction_cutoff,
             stringsAsFactors = FALSE)
}

#' Fit the full DESeq2 model and reconstruct a condition-null mean matrix
#' @noRd
mv_fit_full_null <- function(counts, metadata, condition_col, ref_level, test_level,
                             covariate = NULL) {
  .mv_require("DESeq2")
  aligned <- .mv_align_inputs(counts, metadata, condition_col, ref_level, test_level)
  des <- .mv_design(aligned$metadata, condition_col, covariate)
  dds <- DESeq2::DESeqDataSetFromMatrix(aligned$counts, des$metadata, des$formula)
  dds <- DESeq2::DESeq(dds, test = "Wald", quiet = TRUE)
  beta <- stats::coef(dds, SE = FALSE)
  x <- stats::model.matrix(des$formula, des$metadata)
  if (ncol(beta) != ncol(x)) stop("DESeq2 coefficient/model matrix mismatch.", call. = FALSE)
  # DESeq2 renames factor columns (for example conditionT -> condition_T_vs_C),
  # but preserves their model-matrix order.
  colnames(x) <- colnames(beta)
  factors <- DESeq2::normalizationFactors(dds)
  if (is.null(factors)) factors <- matrix(DESeq2::sizeFactors(dds), nrow(beta), nrow(x), byrow = TRUE)
  mu_full <- 2^(beta %*% t(x)) * factors
  observed_mu <- SummarizedExperiment::assays(dds)[["mu"]]
  if (!isTRUE(all.equal(unname(mu_full), unname(observed_mu), tolerance = 1e-6, check.attributes = FALSE)))
    stop("DESeq2 mu reconstruction check failed; refusing to simulate a mis-scaled null.", call. = FALSE)
  condition_coef <- .mv_condition_coef(colnames(beta), condition_col, test_level, ref_level)
  beta0 <- beta
  beta0[, condition_coef] <- 0
  mu0 <- 2^(beta0 %*% t(x)) * factors
  # Genes with all-zero observed counts get NA beta/dispersion from DESeq2 (no
  # information to fit a coefficient from). Such a gene has a legitimate null:
  # it stays at zero under any specification, so its simulated null mean is 0,
  # not "unknown". Left as NA, it silently poisons the whole bootstrap replicate
  # (rnbinom(mu = NA) -> NA counts -> the alignment check in the next bootstrap
  # fit trips on any(counts < 0) being NA rather than TRUE/FALSE).
  mu0[is.na(mu0)] <- 0
  dispersions <- DESeq2::dispersions(dds)
  dispersions[is.na(dispersions)] <- 1e-8
  list(mu0 = mu0, dispersions = dispersions, full_mu = mu_full,
       reconstruction_ok = TRUE, metadata = des$metadata)
}

#' Simulate one count matrix from a fitted full-model condition-null
#' @noRd
mv_simulate_null <- function(null_fit, seed) {
  set.seed(seed)
  alpha <- pmax(null_fit$dispersions, 1e-8)
  out <- matrix(stats::rnbinom(length(null_fit$mu0), mu = as.vector(null_fit$mu0),
                               size = rep(1 / alpha, ncol(null_fit$mu0))),
                nrow = nrow(null_fit$mu0), dimnames = dimnames(null_fit$mu0))
  storage.mode(out) <- "integer"
  out
}

.mv_efdr_curve <- function(observed, bootstrap_stability, target = 0.10) {
  m <- ncol(observed$tested)
  tau <- seq.int(0, m) / m
  obs_r <- vapply(tau, function(t) sum(observed$stability$stability >= t & observed$stability$passes_direction), numeric(1))
  null_r <- vapply(tau, function(t) mean(vapply(bootstrap_stability, function(x)
    sum(x$stability >= t & x$passes_direction), numeric(1))), numeric(1))
  raw <- ifelse(obs_r > 0, pmin(1, (1 + length(bootstrap_stability) * null_r) /
                                  ((length(bootstrap_stability) + 1) * obs_r)), NA_real_)
  # A gene of stability T is called at tau <= T, so q-value logic takes the
  # cumulative minimum in ascending tau (not the reverse, stricter direction).
  monotone <- cummin(replace(raw, is.na(raw), Inf))
  monotone[obs_r == 0] <- NA_real_
  curve <- data.frame(tau = tau, observed_calls = obs_r, expected_null_calls = null_r,
                      efdr_raw = raw, efdr = monotone)
  # tau = 0 imposes no cross-specification support and is never a stability call.
  candidates <- which(curve$tau > 0 & !is.na(curve$efdr) & curve$efdr <= target)
  selected_tau <- if (length(candidates)) curve$tau[min(candidates)] else NA_real_
  list(curve = curve, selected_tau = selected_tau)
}

#' Run the MVP bootstrap calibration
#' @noRd
mv_bootstrap_efdr <- function(observed, B = 100, seed = 1, padj_cutoff = 0.05,
                              lfc_cutoff = 1, direction_cutoff = 0.90, target_efdr = 0.10) {
  if (B < 1) stop("B must be positive.", call. = FALSE)
  cc <- observed$contrast
  null_fit <- mv_fit_full_null(observed$counts, observed$metadata, cc$condition_col,
                               cc$ref, cc$test, cc$covariate)
  boot <- vector("list", B)
  for (b in seq_len(B)) {
    sim <- mv_simulate_null(null_fit, seed + b)
    run <- mv_run_observed(sim, observed$metadata, cc$condition_col, cc$ref, cc$test,
                           cc$covariate, include_apeglm = any(observed$specifications$method == "deseq_wald_apeglm"))
    boot[[b]] <- mv_compute_stability(run$stats, run$tested, padj_cutoff, lfc_cutoff, direction_cutoff)
  }
  .mv_efdr_curve(observed, boot, target_efdr)
}

#' Run one independently schedulable null-bootstrap replicate
#' @noRd
mv_bootstrap_one <- function(observed, null_fit, b, seed, padj_cutoff = 0.05,
                             lfc_cutoff = 1, direction_cutoff = 0.90) {
  cc <- observed$contrast
  sim <- mv_simulate_null(null_fit, seed + b)
  run <- mv_run_observed(sim, observed$metadata, cc$condition_col, cc$ref, cc$test,
                         cc$covariate, include_apeglm = any(observed$specifications$method == "deseq_wald_apeglm"))
  mv_compute_stability(run$stats, run$tested, padj_cutoff, lfc_cutoff, direction_cutoff)
}

#' Run a complete MVP multiverse analysis
#' @noRd
mv_run_multiverse <- function(counts, metadata, condition_col, ref_level, test_level,
                              covariate = NULL, B = 100, seed = 1,
                              padj_cutoff = 0.05, lfc_cutoff = 1,
                              direction_cutoff = 0.90, include_apeglm = TRUE,
                              target_efdr = 0.10) {
  if (B < 100 && target_efdr < 0.10)
    stop("B < 100 is exploratory; eFDR targets below 0.10 are unavailable.", call. = FALSE)
  observed <- mv_run_observed(counts, metadata, condition_col, ref_level, test_level,
                              covariate, include_apeglm)
  observed$stability <- mv_compute_stability(observed$stats, observed$tested,
                                             padj_cutoff, lfc_cutoff, direction_cutoff)
  calibration <- mv_bootstrap_efdr(observed, B, seed, padj_cutoff, lfc_cutoff,
                                   direction_cutoff, target_efdr)
  stability <- mv_attach_efdr(observed$stability, calibration$curve, target_efdr)
  list(run = observed, stability = stability, efdr_curve = calibration$curve,
       selected_tau = calibration$selected_tau, B = B, seed = seed,
       target_efdr = target_efdr)
}

#' Add per-gene eFDR values at each gene's attainable stability
#' @noRd
mv_attach_efdr <- function(stability, curve, target = 0.10) {
  m <- length(curve$tau) - 1L
  idx <- as.integer(round(stability$stability * m)) + 1L
  if (any(idx < 1L | idx > nrow(curve)))
    stop("Stability is outside the attainable eFDR grid.", call. = FALSE)
  stability$efdr <- curve$efdr[idx]
  # A direction-inconsistent/untestable gene is never a stability discovery.
  stability$efdr[!stability$passes_direction] <- 1
  stability$called <- stability$passes_direction & !is.na(stability$efdr) & stability$efdr <= target
  stability
}

#' Make a compatibility DEG table for downstream modules
#' @noRd
mv_make_handoff <- function(run, stability) {
  lfc <- run$stats[, , "log2FoldChange"]
  base_mean <- run$stats[, , "baseMean"]
  signed <- ifelse(stability$direction == "Up", stability$stability, -stability$stability)
  keep_dir <- lfc >= 0
  down_rows <- stability$direction == "Down"
  keep_dir[down_rows, ] <- lfc[down_rows, , drop = FALSE] <= 0
  med_lfc <- vapply(seq_len(nrow(lfc)), function(i) stats::median(lfc[i, keep_dir[i, ]], na.rm = TRUE), numeric(1))
  med_lfc[!is.finite(med_lfc)] <- NA_real_
  out <- data.frame(gene = stability$gene, baseMean = rowMeans(base_mean, na.rm = TRUE),
                    log2FoldChange = med_lfc, lfcSE = NA_real_, stat = signed,
                    pvalue = NA_real_, padj = stability$efdr, stability = stability$stability,
                    direction_consistency = stability$direction_consistency,
                    stringsAsFactors = FALSE)
  attr(out, "result_type") <- "multiverse_stability"
  out[order(out$padj, -out$stability, na.last = TRUE), , drop = FALSE]
}

#' Apply a completed multiverse result to AppState-compatible state
#' @noRd
mv_apply_run <- function(state, run, stability, curve, target = 0.10, seed = NA_integer_, B = NA_integer_) {
  handoff <- mv_make_handoff(run, stability)
  key <- paste0("multiverse:", run$contrast$test, "_vs_", run$contrast$ref)
  results <- state$deg_results(); results[[key]] <- handoff; state$deg_results(results)
  state$deg_multiverse(list(run = run, stability = stability, efdr_curve = curve,
                             metadata = list(seed = seed, B = B, target = target)))
  state$log_param("multiverse_deg", "contrast", key)
  state$log_param("multiverse_deg", "specifications", nrow(run$specifications))
  state$log_param("multiverse_deg", "bootstrap_B", B)
  state$log_param("multiverse_deg", "target_efdr", target)
  state$log_param("multiverse_deg", "compatibility_deg_status", "done")
  state$set_status("deg_multiverse", "done")
  # Compatibility flag: existing Enrichment and UpSet gates require the deg step.
  state$set_status("deg", "done")
  invisible(handoff)
}
