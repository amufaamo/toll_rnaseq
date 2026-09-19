#!/usr/bin/env Rscript
# Experiment D: why is the calibrated stability threshold tau not binding?
#
#   conda run -n toll Rscript scripts/analyze_spec_concordance.R [run ...]
#
# The automatic rule picks the smallest attainable positive threshold (tau = 1/m)
# whenever the eFDR curve is already below the target there. Two mechanisms could
# produce that: (a) the m specifications agree so strongly that false positives
# co-occur across specifications and raising tau removes signal and noise in the
# same proportion, or (b) false positives are simply scarce in absolute terms,
# so the ratio is already small at tau = 1/m.
#
# This script measures both directly, without refitting any bootstrap
# calibration and without modifying the engine:
#   * observed run: pairwise Jaccard and overlap of the per-specification call
#     sets, broken down by whether the pair shares a method family, a filter or
#     a covariate arm;
#   * null runs: the same quantities on `n_null` matrices simulated from the
#     fitted condition-null, plus the distribution of null stability. Under
#     independent specifications almost every null false positive would sit at
#     stability 1/m; concentration above that is the signature of mechanism (a).
#     The independence reference is computed by permuting, within each
#     specification, which genes it calls, which destroys cross-specification
#     coupling while preserving each specification's own call count.
#
# Outputs go to validation/results/spec_concordance_*.csv.

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
})
MV_REPO <- "/mnt/g/マイドライブ/toll_rnaseq"
source(file.path(MV_REPO, "scripts", "validate_common.R"))

N_NULL <- 5L
SEED   <- 1L

.call_matrix <- function(run, padj_cutoff = MV_PADJ, lfc_cutoff = MV_LFC) {
  lfc  <- run$stats[, , "log2FoldChange"]
  padj <- run$stats[, , "padj"]
  run$tested & !is.na(padj) & padj <= padj_cutoff & abs(lfc) >= lfc_cutoff
}

.pairwise <- function(calls, specs, label) {
  m <- ncol(calls)
  if (m < 2) return(NULL)
  cmb <- utils::combn(m, 2)
  do.call(rbind, lapply(seq_len(ncol(cmb)), function(k) {
    i <- cmb[1, k]; j <- cmb[2, k]
    a <- calls[, i]; b <- calls[, j]
    inter <- sum(a & b); uni <- sum(a | b)
    data.frame(source = label, spec_i = specs$id[i], spec_j = specs$id[j],
               method_i = specs$method[i], method_j = specs$method[j],
               same_method = specs$method[i] == specs$method[j],
               same_filter = specs$filter[i] == specs$filter[j],
               same_covariate = identical(specs$covariate[i], specs$covariate[j]),
               n_i = sum(a), n_j = sum(b), n_intersect = inter,
               jaccard = if (uni > 0) inter / uni else NA_real_,
               overlap_min = if (min(sum(a), sum(b)) > 0) inter / min(sum(a), sum(b)) else NA_real_,
               stringsAsFactors = FALSE)
  }))
}

# Stability distribution implied by specifications that call the same number of
# genes each but choose them independently.
.independence_reference <- function(calls, seed) {
  set.seed(seed)
  perm <- calls
  for (j in seq_len(ncol(calls))) perm[, j] <- sample(calls[, j])
  rowMeans(perm)
}

.stability_table <- function(stab_vec, m, source, replicate = NA_integer_) {
  k <- round(stab_vec * m)
  tab <- table(factor(k, levels = seq_len(m)))
  data.frame(source = source, replicate = replicate, n_specs_calling = as.integer(names(tab)),
             stability = as.integer(names(tab)) / m, n_genes = as.integer(tab),
             stringsAsFactors = FALSE)
}

analyse_run <- function(key) {
  d <- MV_RUNS[[key]]()
  z <- mv_drop_all_zero(d$counts)
  counts <- z$counts
  cat(sprintf("\n=== %s (%s) : %d genes, %d samples ===\n", key, d$dataset, nrow(counts), ncol(counts)))

  t0 <- Sys.time()
  obs <- mv_run_observed(counts, d$metadata, d$condition_col, d$ref_level, d$test_level,
                         covariate = d$covariate, include_apeglm = TRUE)
  specs <- obs$specifications
  m <- nrow(specs)
  obs_calls <- .call_matrix(obs)
  cat(sprintf("specifications: %d; per-spec call counts: %s\n", m,
              paste(colSums(obs_calls), collapse = ", ")))

  pair_rows <- list(.pairwise(obs_calls, specs, "observed"))
  stab_rows <- list(.stability_table(rowMeans(obs_calls), m, "observed"),
                    .stability_table(.independence_reference(obs_calls, SEED), m,
                                     "observed_independence_reference"))

  null_fit <- mv_fit_full_null(counts, d$metadata, d$condition_col, d$ref_level,
                               d$test_level, d$covariate)
  for (b in seq_len(N_NULL)) {
    sim <- mv_simulate_null(null_fit, SEED + b)
    nrun <- mv_run_observed(sim, d$metadata, d$condition_col, d$ref_level, d$test_level,
                            covariate = d$covariate, include_apeglm = TRUE)
    ncalls <- .call_matrix(nrun)
    cat(sprintf("null replicate %d: per-spec call counts: %s\n", b,
                paste(colSums(ncalls), collapse = ", ")))
    pair_rows[[length(pair_rows) + 1]] <- transform(.pairwise(ncalls, nrun$specifications, "null"),
                                                    replicate = b)
    stab_rows[[length(stab_rows) + 1]] <- .stability_table(rowMeans(ncalls), m, "null", b)
    stab_rows[[length(stab_rows) + 1]] <-
      .stability_table(.independence_reference(ncalls, SEED + 100L + b), m,
                       "null_independence_reference", b)
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  pairs <- do.call(rbind, lapply(pair_rows, function(x) {
    if (is.null(x$replicate)) x$replicate <- NA_integer_
    x
  }))
  pairs$run <- key
  stabs <- do.call(rbind, stab_rows)
  stabs$run <- key

  # Headline numbers for the manuscript.
  obs_pairs <- pairs[pairs$source == "observed", ]
  null_pairs <- pairs[pairs$source == "null", ]
  null_obs <- stabs[stabs$source == "null", ]
  null_ind <- stabs[stabs$source == "null_independence_reference", ]
  frac_ge2 <- function(tab) {
    tot <- sum(tab$n_genes); ge2 <- sum(tab$n_genes[tab$n_specs_calling >= 2])
    if (tot > 0) ge2 / tot else NA_real_
  }
  summ <- data.frame(
    run = key, dataset = d$dataset, n_specs = m, n_genes = nrow(counts),
    obs_calls_union = sum(rowSums(obs_calls) > 0),
    obs_calls_unanimous = sum(rowSums(obs_calls) == m),
    obs_jaccard_median = stats::median(obs_pairs$jaccard, na.rm = TRUE),
    obs_jaccard_min = min(obs_pairs$jaccard, na.rm = TRUE),
    obs_jaccard_same_method = stats::median(obs_pairs$jaccard[obs_pairs$same_method], na.rm = TRUE),
    obs_jaccard_diff_method = stats::median(obs_pairs$jaccard[!obs_pairs$same_method], na.rm = TRUE),
    obs_jaccard_same_cov = stats::median(obs_pairs$jaccard[obs_pairs$same_covariate], na.rm = TRUE),
    obs_jaccard_diff_cov = stats::median(obs_pairs$jaccard[!obs_pairs$same_covariate], na.rm = TRUE),
    null_calls_per_spec_mean = mean(c(null_pairs$n_i, null_pairs$n_j), na.rm = TRUE),
    null_jaccard_median = stats::median(null_pairs$jaccard, na.rm = TRUE),
    null_union_mean = mean(tapply(null_obs$n_genes, null_obs$replicate, sum)),
    null_frac_stability_ge_2m = frac_ge2(null_obs),
    null_frac_stability_ge_2m_if_independent = frac_ge2(null_ind),
    n_null_replicates = N_NULL, elapsed_sec = elapsed, stringsAsFactors = FALSE)
  print(summ, row.names = FALSE, digits = 3)

  list(pairs = pairs, stability = stabs, summary = summ)
}

args <- commandArgs(trailingOnly = TRUE)
keys <- if (length(args)) args else names(MV_RUNS)
res <- lapply(keys, analyse_run)

out <- file.path(MV_REPO, "validation", "results")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
# One file per run so that the four runs can be executed as parallel processes;
# scripts/summarize_tau_experiments.R combines the summaries.
for (i in seq_along(keys)) {
  for (part in c("pairs", "stability", "summary"))
    utils::write.csv(res[[i]][[part]],
                     file.path(out, sprintf("spec_concordance_%s_%s.csv", keys[i], part)),
                     row.names = FALSE)
}
cat(sprintf("\nwrote %s/spec_concordance_<run>_{pairs,stability,summary}.csv for: %s\n",
            out, paste(keys, collapse = ", ")))
