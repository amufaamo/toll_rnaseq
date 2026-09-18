# Base-package functions used across modules, plus the column names that
# ggplot2's non-standard evaluation makes invisible to R CMD check.

#' @importFrom grDevices cairo_pdf tiff
#' @importFrom stats as.formula dist prcomp relevel setNames time var
#' @importFrom utils head read.table write.csv
#' @noRd
NULL

utils::globalVariables(c(
  ".data",
  "Count", "Intersection", "NES", "PC", "PC1", "PC2", "Pathway", "Sample",
  "Score", "baseMean", "choice", "cumvar", "direction", "libsize",
  "log2FoldChange", "padj", "pathway", "s1", "s2", "sig", "significance",
  "spec_rank", "value"
))
