# Shared fixtures for the multiverse DEG engine and Shiny module tests.

skip_mv_deps <- function() {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("edgeR")
}

make_mv_counts <- function(seed = 20260918, n_genes = 1000, n_de = 120) {
  set.seed(seed)
  n_per_group <- 6L
  means <- rlnorm(n_genes, log(100), 0.45)
  lib <- rep(seq(0.8, 1.2, length.out = n_per_group), 2)
  lfc <- numeric(n_genes)
  if (n_de > 0) {
    lfc[seq_len(n_de / 2)] <- 1.5
    lfc[seq.int(n_de / 2 + 1, n_de)] <- -1.5
  }
  condition <- factor(rep(c("C", "T"), each = n_per_group))
  mu <- outer(means, lib)
  mu[, condition == "T"] <- mu[, condition == "T", drop = FALSE] * 2^lfc
  counts <- matrix(rnbinom(length(mu), mu = as.vector(mu), size = 1 / 0.15), nrow = n_genes,
                   dimnames = list(paste0("g", seq_len(n_genes)), paste0("s", seq_len(12))))
  list(counts = counts,
       metadata = data.frame(sample = colnames(counts), condition = condition),
       null_genes = paste0("g", seq.int(n_de + 1, n_genes)))
}
