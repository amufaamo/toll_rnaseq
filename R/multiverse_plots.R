#' Multiverse DEG plots
#' @noRd

#' Build an aligned per-gene specification curve
#' @noRd
mv_plot_specification_curve <- function(run, gene, padj_cutoff = 0.05, lfc_cutoff = 1) {
  i <- match(gene, dimnames(run$stats)[[1]])
  if (is.na(i)) stop("Selected gene is not in the multiverse result.", call. = FALSE)
  specs <- run$specifications
  lfc <- run$stats[i, , "log2FoldChange"]
  padj <- run$stats[i, , "padj"]
  ord <- order(lfc, na.last = TRUE)
  specs <- specs[ord, , drop = FALSE]
  lfc <- lfc[ord]; padj <- padj[ord]
  top <- data.frame(spec_rank = seq_along(lfc), value = sign(lfc) * -log10(pmax(padj, 1e-300)),
                    significance = ifelse(!is.na(padj) & padj <= padj_cutoff & lfc >= lfc_cutoff, "Up",
                    ifelse(!is.na(padj) & padj <= padj_cutoff & lfc <= -lfc_cutoff, "Down", "NS")),
                    panel = "Signed -log10(FDR)")
  choices <- c("filter", "normalization", "method", "shrinkage", "covariate")
  lower <- do.call(rbind, lapply(choices, function(choice) data.frame(
    spec_rank = seq_len(nrow(specs)), value = 1, choice = as.character(specs[[choice]]),
    panel = choice, stringsAsFactors = FALSE
  )))
  top$choice <- top$significance
  all <- rbind(top[, c("spec_rank", "value", "choice", "panel")], lower)
  ggplot2::ggplot(all, ggplot2::aes(spec_rank, value)) +
    ggplot2::geom_point(data = top, ggplot2::aes(color = significance), size = 1.8, na.rm = TRUE) +
    ggplot2::geom_tile(data = lower, ggplot2::aes(fill = choice), height = 0.8) +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_color_manual(values = c(Up = "#D55E00", Down = "#0072B2", NS = "#999999")) +
    ggplot2::scale_fill_manual(values = c("min_count" = "#E69F00", "filterByExpr" = "#009E73",
      "median_ratio" = "#0072B2", "TMM" = "#CC79A7", "deseq_wald" = "#0072B2",
      "deseq_wald_apeglm" = "#CC79A7", "edger_ql" = "#009E73", "edger_lrt" = "#E69F00",
      "TRUE" = "#D55E00", "FALSE" = "#999999", "NA" = "#999999"), na.value = "#999999") +
    ggplot2::labs(x = "Specification (ordered by log2 fold change)", y = NULL, color = NULL, fill = NULL,
                  title = paste("Specification curve:", gene)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
}

#' Build a global multiverse stability summary
#' @noRd
mv_plot_global_summary <- function(stability) {
  ggplot2::ggplot(stability, ggplot2::aes(stability, fill = direction)) +
    ggplot2::geom_histogram(binwidth = 0.05, boundary = 0, color = "white") +
    ggplot2::scale_fill_manual(values = c(Up = "#D55E00", Down = "#0072B2")) +
    ggplot2::labs(x = "Maximum signed stability", y = "Genes", fill = "Direction") +
    ggplot2::theme_minimal(base_size = 12)
}
