# Shared harness for the public-dataset validation of the multiverse DEG engine.
#
# Sourced by scripts/validate_airway.R, scripts/validate_pasilla.R,
# scripts/validate_gse60450.R and scripts/validate_tau_conditions.R. It calls
# the same engine API as the Shiny module (R/mod_deg_multiverse.R):
# mv_run_multiverse() with the default MVP grid. The defaults below are the
# module defaults (call rule padj <= 0.05, |log2FC| >= 1, direction consistency
# >= 0.90; 0.10 eFDR target; B = 20); mv_validate() accepts overrides so that
# alternative call rules and eFDR targets can be swept without touching the
# engine.

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
})

MV_REPO <- "/mnt/g/マイドライブ/toll_rnaseq"
source(file.path(MV_REPO, "R", "multiverse_engine.R"))

MV_B            <- 20L
MV_SEED         <- 1L
MV_PADJ         <- 0.05
MV_LFC          <- 1
MV_DIRECTION    <- 0.90
MV_TARGET_EFDR  <- 0.10

# Genes with zero counts in every sample have an undefined fitted mean, so
# mv_fit_full_null() returns NA there and mv_simulate_null() propagates NA into
# the simulated matrix, which aborts the bootstrap. No filter in the grid would
# ever retain such a gene, so they are removed before the engine is called and
# the number removed is reported.
mv_drop_all_zero <- function(counts) {
  keep <- rowSums(counts) > 0
  list(counts = counts[keep, , drop = FALSE], n_dropped = sum(!keep))
}

mv_validate <- function(dataset, counts, metadata, condition_col, ref_level, test_level,
                        covariate = NULL, symbols = NULL, markers = character(),
                        B = MV_B, seed = MV_SEED,
                        padj_cutoff = MV_PADJ, lfc_cutoff = MV_LFC,
                        direction_cutoff = MV_DIRECTION, target_efdr = MV_TARGET_EFDR,
                        condition_label = "default", tag_suffix = "") {
  stopifnot(is.matrix(counts), storage.mode(counts) == "integer",
            colnames(metadata)[1] == "sample")

  z <- mv_drop_all_zero(counts)
  counts <- z$counts

  cat(sprintf("\n=== %s: %s %s vs %s%s ===\n", dataset, condition_col, test_level, ref_level,
              if (is.null(covariate)) " (no covariate)" else sprintf(" | covariate: %s", covariate)))
  cat(sprintf("condition: %s | call rule: padj <= %g, |log2FC| >= %g, C >= %g | B = %d, target eFDR = %g\n",
              condition_label, padj_cutoff, lfc_cutoff, direction_cutoff, B, target_efdr))
  cat(sprintf("genes: %d (dropped %d all-zero), samples: %d\n", nrow(counts), z$n_dropped, ncol(counts)))
  print(table(metadata[[condition_col]]))
  if (!is.null(covariate)) print(table(metadata[[condition_col]], metadata[[covariate]]))

  t0 <- Sys.time()
  mv <- mv_run_multiverse(counts, metadata, condition_col, ref_level, test_level,
                          covariate = covariate, B = B, seed = seed,
                          padj_cutoff = padj_cutoff, lfc_cutoff = lfc_cutoff,
                          direction_cutoff = direction_cutoff, include_apeglm = TRUE,
                          target_efdr = target_efdr)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  specs <- mv$run$specifications
  stab  <- mv$stability
  curve <- mv$efdr_curve
  tau   <- mv$selected_tau

  # Single default pipeline = the grid member a conventional analysis would use:
  # count >= 10 filter, DESeq2 Wald, no shrinkage, same design as the multiverse.
  single_id <- specs$id[specs$filter == "min_count" & specs$method == "deseq_wald" &
                          (if (is.null(covariate)) is.na(specs$covariate) else
                             (!is.na(specs$covariate) & specs$covariate == covariate))]
  stopifnot(length(single_id) == 1)
  s_padj <- mv$run$stats[, single_id, "padj"]
  s_lfc  <- mv$run$stats[, single_id, "log2FoldChange"]
  single_sig <- !is.na(s_padj) & s_padj <= padj_cutoff & abs(s_lfc) >= lfc_cutoff

  called <- stab$called
  names(single_sig) <- rownames(mv$run$stats)
  single_genes <- names(single_sig)[single_sig]
  called_genes <- stab$gene[called]

  stab$single_padj <- s_padj[match(stab$gene, names(s_padj))]
  stab$single_lfc  <- s_lfc[match(stab$gene, names(s_lfc))]
  stab$single_sig  <- single_sig[match(stab$gene, names(single_sig))]
  if (!is.null(symbols)) stab$symbol <- unname(symbols[stab$gene]) else stab$symbol <- NA_character_

  cat(sprintf("\nspecifications        : %d (pruned %d)\n", nrow(specs), nrow(mv$run$skipped)))
  cat(sprintf("bootstrap replicates  : %d\n", B))
  cat(sprintf("wall-clock            : %.1f s (%.1f min)\n", elapsed, elapsed / 60))
  # mv_run_multiverse() computes selected_tau against the default 0.10 eFDR
  # target regardless of target_efdr (mv_bootstrap_efdr() has no target
  # argument), so this line is only meaningful at target_efdr = 0.10. The eFDR
  # curve written below is correct; scripts/summarize_tau_experiments.R
  # recomputes the threshold from it at the requested target.
  cat(sprintf("selected tau          : %s (minimum attainable positive tau = %.4f; binding = %s)\n",
              format(tau), 1 / nrow(specs), !is.na(tau) && tau > 1 / nrow(specs) + 1e-9))
  cat(sprintf("multiverse calls      : %d\n", length(called_genes)))
  cat(sprintf("distinct eFDR values  : %d\n", length(unique(stats::na.omit(curve$efdr)))))
  d <- diff(stats::na.omit(curve$efdr))
  cat(sprintf("eFDR monotone non-incr: %s\n", all(d <= 1e-12)))
  cat("\neFDR curve:\n"); print(curve, row.names = FALSE, digits = 4)

  cat(sprintf("\nsingle pipeline (%s, DESeq2 Wald, count>=10 filter): %d DEGs\n",
              single_id, length(single_genes)))
  cat(sprintf("  overlap with multiverse calls      : %d\n", length(intersect(single_genes, called_genes))))
  cat(sprintf("  single-only (not called)           : %d\n", length(setdiff(single_genes, called_genes))))
  cat(sprintf("  multiverse-only (not single-sig)   : %d\n", length(setdiff(called_genes, single_genes))))
  cat("\nstability of single-pipeline DEGs:\n")
  print(table(round(stab$stability[stab$single_sig %in% TRUE], 4)))
  cat(sprintf("  fragile (stability < 1, i.e. lost in >=1 specification): %d of %d\n",
              sum(stab$single_sig %in% TRUE & stab$stability < 1), length(single_genes)))
  cat(sprintf("  unanimous (stability = 1)                              : %d\n",
              sum(stab$single_sig %in% TRUE & stab$stability >= 1)))
  cat(sprintf("  direction-inconsistent (C < cutoff)                    : %d\n",
              sum(stab$single_sig %in% TRUE & !stab$passes_direction)))

  top <- stab[stab$called, , drop = FALSE]
  top <- top[order(-top$stability, top$efdr, -abs(top$single_lfc)), , drop = FALSE]
  cat("\ntop 25 genes by stability:\n")
  print(utils::head(top[, c("gene", "symbol", "stability", "direction",
                            "direction_consistency", "efdr", "single_lfc", "single_padj")], 25),
        row.names = FALSE, digits = 3)

  if (length(markers)) {
    cat("\nface validity - pre-specified marker genes from the original publication:\n")
    key <- if (!is.null(symbols)) stab$symbol else stab$gene
    for (mk in markers) {
      i <- which(key == mk)
      if (!length(i)) { cat(sprintf("  %-12s not present in the count matrix\n", mk)); next }
      i <- i[1]
      rk <- if (stab$called[i]) which(top$gene == stab$gene[i]) else NA_integer_
      cat(sprintf("  %-12s stability %.3f  %-4s  eFDR %.4f  called %-5s  rank %s\n",
                  mk, stab$stability[i], stab$direction[i], stab$efdr[i],
                  stab$called[i], if (is.na(rk)) "-" else as.character(rk)))
    }
  }

  summary_row <- data.frame(
    dataset = dataset, condition = condition_col, contrast = sprintf("%s_vs_%s", test_level, ref_level),
    covariate = if (is.null(covariate)) NA_character_ else covariate,
    experiment = condition_label, padj_cutoff = padj_cutoff, lfc_cutoff = lfc_cutoff,
    direction_cutoff = direction_cutoff, target_efdr = target_efdr,
    n_samples = ncol(counts), n_genes = nrow(counts), n_genes_dropped_all_zero = z$n_dropped,
    n_specifications = nrow(specs), B = B, seed = seed,
    min_positive_tau = 1 / nrow(specs),
    tau_binding = !is.na(tau) && tau > 1 / nrow(specs) + 1e-9,
    efdr_at_tau = if (is.na(tau)) NA_real_ else curve$efdr[which.min(abs(curve$tau - tau))],
    efdr_at_min_tau = curve$efdr[2],
    null_calls_at_min_tau = curve$expected_null_calls[2],
    observed_calls_at_min_tau = curve$observed_calls[2],
    selected_tau = tau, n_calls = length(called_genes),
    n_distinct_efdr = length(unique(stats::na.omit(curve$efdr))),
    efdr_monotone = all(d <= 1e-12),
    single_pipeline_spec = single_id, n_single_deg = length(single_genes),
    n_overlap = length(intersect(single_genes, called_genes)),
    n_single_only = length(setdiff(single_genes, called_genes)),
    n_multiverse_only = length(setdiff(called_genes, single_genes)),
    n_single_fragile = sum(stab$single_sig %in% TRUE & stab$stability < 1),
    n_single_unanimous = sum(stab$single_sig %in% TRUE & stab$stability >= 1),
    n_single_dir_inconsistent = sum(stab$single_sig %in% TRUE & !stab$passes_direction),
    elapsed_sec = elapsed, stringsAsFactors = FALSE)

  out <- file.path(MV_REPO, "validation", "results")
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  tag <- sprintf("%s_%s_vs_%s%s%s", dataset, test_level, ref_level,
                 if (is.null(covariate)) "" else paste0("_cov_", covariate),
                 if (nzchar(tag_suffix)) paste0("_", tag_suffix) else "")
  saveRDS(list(summary = summary_row, stability = stab, efdr_curve = curve,
               specifications = specs, skipped = mv$run$skipped,
               session = utils::sessionInfo()),
          file.path(out, sprintf("%s_summary.rds", tag)))
  utils::write.csv(summary_row, file.path(out, sprintf("%s_summary.csv", tag)), row.names = FALSE)
  utils::write.csv(curve, file.path(out, sprintf("%s_efdr_curve.csv", tag)), row.names = FALSE)
  utils::write.csv(top[, c("gene", "symbol", "stability", "direction", "direction_consistency",
                           "efdr", "single_lfc", "single_padj")],
                   file.path(out, sprintf("%s_calls.csv", tag)), row.names = FALSE)
  cat(sprintf("\nwrote %s/%s_{summary.rds,summary.csv,efdr_curve.csv,calls.csv}\n", out, tag))
  invisible(list(mv = mv, stability = stab, summary = summary_row))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

source(file.path(MV_REPO, "scripts", "validate_datasets.R"))
