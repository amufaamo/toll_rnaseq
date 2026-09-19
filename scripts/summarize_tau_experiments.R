#!/usr/bin/env Rscript
# Collect the per-run outputs of the threshold-binding experiments into two
# tables for the manuscript.
#
#   conda run -n toll Rscript scripts/summarize_tau_experiments.R
#
#   validation/results/tau_conditions_summary.csv
#       one row per (call rule / eFDR target / dataset) run, from the
#       *_summary.csv files written by mv_validate(), with the default-condition
#       runs included for comparison.
#   validation/results/spec_concordance_summary.csv
#       one row per dataset, from scripts/analyze_spec_concordance.R.

MV_REPO <- "/mnt/g/マイドライブ/toll_rnaseq"
out <- file.path(MV_REPO, "validation", "results")

bind_rows <- function(files) {
  tabs <- lapply(files, utils::read.csv, stringsAsFactors = FALSE)
  cols <- Reduce(union, lapply(tabs, names))
  do.call(rbind, lapply(tabs, function(x) {
    for (cn in setdiff(cols, names(x))) x[[cn]] <- NA
    x[, cols, drop = FALSE]
  }))
}

summ <- list.files(out, pattern = "_summary\\.csv$", full.names = TRUE)
summ <- setdiff(summ, list.files(out, pattern = "^(all_datasets|spec_concordance|tau_conditions)",
                                 full.names = TRUE))
tau <- bind_rows(summ)
if (is.null(tau$experiment)) tau$experiment <- "default"
tau$experiment[is.na(tau$experiment)] <- "default"
tau$padj_cutoff[is.na(tau$padj_cutoff)] <- 0.05
tau$lfc_cutoff[is.na(tau$lfc_cutoff)] <- 1
tau$direction_cutoff[is.na(tau$direction_cutoff)] <- 0.90
tau$target_efdr[is.na(tau$target_efdr)] <- 0.10
tau$min_positive_tau[is.na(tau$min_positive_tau)] <- 1 / tau$n_specifications[is.na(tau$min_positive_tau)]
# mv_run_multiverse() reports a selected_tau that is always computed at the
# default 0.10 eFDR target: mv_bootstrap_efdr() has no target argument and calls
# .mv_efdr_curve() with its default, so the target_efdr argument reaches
# mv_attach_efdr() (which decides the per-gene calls, and is correct) but not the
# threshold-selection step. The Shiny module is not affected -- it calls
# .mv_efdr_curve() with the target directly. The saved eFDR curves are correct,
# so the threshold the selection rule would have chosen is recomputed here from
# the curve, by the same rule .mv_efdr_curve() applies: the smallest strictly
# positive tau whose eFDR is at or below the target. This is an analysis-side
# correction; the engine is not modified.
curve_of <- function(row) {
  f <- file.path(out, sprintf("%s_%s%s%s_efdr_curve.csv", row$dataset, row$contrast,
                              if (is.na(row$covariate)) "" else paste0("_cov_", row$covariate),
                              if (row$experiment == "default") "" else paste0("_", row$experiment)))
  if (!file.exists(f)) stop(sprintf("missing eFDR curve: %s", f), call. = FALSE)
  utils::read.csv(f, stringsAsFactors = FALSE)
}
selected_from_curve <- function(row) {
  cv <- curve_of(row)
  ok <- which(cv$tau > 0 & !is.na(cv$efdr) & cv$efdr <= row$target_efdr)
  if (length(ok)) cv$tau[min(ok)] else NA_real_
}
tau$selected_tau_at_target <- vapply(seq_len(nrow(tau)),
                                     function(i) selected_from_curve(tau[i, ]), numeric(1))
# The oldest runs predate the extra summary columns; fill them from the curves.
for (i in seq_len(nrow(tau))) {
  cv <- curve_of(tau[i, ])
  if (is.na(tau$efdr_at_min_tau[i])) tau$efdr_at_min_tau[i] <- cv$efdr[2]
  if (is.na(tau$observed_calls_at_min_tau[i])) tau$observed_calls_at_min_tau[i] <- cv$observed_calls[2]
  if (is.na(tau$null_calls_at_min_tau[i])) tau$null_calls_at_min_tau[i] <- cv$expected_null_calls[2]
  tau$efdr_at_tau[i] <- cv$efdr[which.min(abs(cv$tau - tau$selected_tau_at_target[i]))]
  tau$calls_at_min_tau_vs_selected[i] <- cv$observed_calls[2] -
    cv$observed_calls[which.min(abs(cv$tau - tau$selected_tau_at_target[i]))]
}
tau$engine_tau_matches <- isTRUE(all.equal(tau$selected_tau, tau$selected_tau_at_target)) |
  (!is.na(tau$selected_tau) & !is.na(tau$selected_tau_at_target) &
     abs(tau$selected_tau - tau$selected_tau_at_target) < 1e-9)
tau$tau_binding <- !is.na(tau$selected_tau_at_target) &
  tau$selected_tau_at_target > tau$min_positive_tau + 1e-9
tau <- tau[order(tau$experiment, tau$dataset), , drop = FALSE]
utils::write.csv(tau, file.path(out, "tau_conditions_summary.csv"), row.names = FALSE)

key <- c("experiment", "dataset", "n_specifications", "B", "padj_cutoff",
         "lfc_cutoff", "target_efdr", "min_positive_tau", "selected_tau_at_target",
         "tau_binding", "selected_tau", "engine_tau_matches", "n_calls",
         "efdr_at_min_tau", "observed_calls_at_min_tau", "null_calls_at_min_tau",
         "n_single_deg", "elapsed_sec")
print(tau[, intersect(key, names(tau)), drop = FALSE], row.names = FALSE, digits = 4)

conc <- list.files(out, pattern = "^spec_concordance_.*_summary\\.csv$", full.names = TRUE)
if (length(conc)) {
  cs <- bind_rows(conc)
  utils::write.csv(cs, file.path(out, "spec_concordance_summary.csv"), row.names = FALSE)
  cat("\n")
  print(cs, row.names = FALSE, digits = 3)
}
# Concordance broken down by which grid axis a specification pair differs on.
# "filter twin" = same method and same covariate arm, differing only in the
# gene filter; these are the pairs that make the grid's effective width smaller
# than its nominal m.
pairfiles <- list.files(out, pattern = "^spec_concordance_.*_pairs\\.csv$", full.names = TRUE)
if (length(pairfiles)) {
  pr <- bind_rows(pairfiles)
  axis <- ifelse(!pr$same_covariate, "covariate differs",
                 ifelse(!pr$same_method, "method differs", "filter twin"))
  agg <- stats::aggregate(list(median_jaccard = pr$jaccard, n_pairs = pr$jaccard),
                          by = list(run = pr$run, source = pr$source, axis = axis),
                          FUN = function(x) c(stats::median(x, na.rm = TRUE), sum(!is.na(x))))
  axis_tab <- data.frame(run = agg$run, source = agg$source, axis = agg$axis,
                         median_jaccard = agg$median_jaccard[, 1],
                         n_pairs = agg$n_pairs[, 2], stringsAsFactors = FALSE)
  axis_tab <- axis_tab[axis_tab$source %in% c("observed", "null"), , drop = FALSE]
  axis_tab <- axis_tab[order(axis_tab$run, axis_tab$source, axis_tab$axis), , drop = FALSE]
  utils::write.csv(axis_tab, file.path(out, "spec_concordance_axis.csv"), row.names = FALSE)
  cat("\n"); print(axis_tab, row.names = FALSE, digits = 3)
}

# Parity of the number of calling specifications: with a redundant filter axis,
# genes are called by an even number of specifications almost always.
stabfiles <- list.files(out, pattern = "^spec_concordance_.*_stability\\.csv$", full.names = TRUE)
if (length(stabfiles)) {
  st <- bind_rows(stabfiles)
  st <- st[st$source == "observed", , drop = FALSE]
  parity <- do.call(rbind, lapply(split(st, st$run), function(x) {
    ev <- sum(x$n_genes[x$n_specs_calling %% 2 == 0])
    od <- sum(x$n_genes[x$n_specs_calling %% 2 == 1])
    data.frame(run = x$run[1], genes_even_k = ev, genes_odd_k = od,
               fraction_even = ev / (ev + od), stringsAsFactors = FALSE)
  }))
  utils::write.csv(parity, file.path(out, "spec_concordance_parity.csv"), row.names = FALSE)
  cat("\n"); print(parity, row.names = FALSE, digits = 3)
}

cat(sprintf("\nwrote %s/{tau_conditions_summary.csv,spec_concordance_summary.csv,spec_concordance_axis.csv,spec_concordance_parity.csv}\n", out))
