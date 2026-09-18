#!/usr/bin/env Rscript
# Generate manuscript Figure 3: a per-gene specification curve from a real
# multiverse DEG run on the repository's synthetic fixture.
#
#   conda run -n toll Rscript scripts/generate_fig3_specification_curve.R
#
# Fixture and grid match the synthetic calibration reported in the manuscript:
# make_mv_counts() from tests/testthat/helper-multiverse.R with 400 genes,
# 40 true DE genes and n = 6 per group; the 8-path grid (2 filters x
# {DESeq2 Wald, DESeq2 Wald + apeglm, edgeR QL/TMM, edgeR LRT/TMM}); B = 10
# bootstrap replicates; call rule padj <= 0.05 and |log2FC| >= 1; target
# eFDR 0.10. Seeds are fixed so the figure is reproducible.

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
})

repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
source(file.path(repo, "R", "multiverse_engine.R"))
source(file.path(repo, "R", "multiverse_plots.R"))
source(file.path(repo, "tests", "testthat", "helper-multiverse.R"))

DATA_SEED <- 20260918L   # helper-multiverse.R default
RUN_SEED  <- 1L          # mv_run_multiverse() default
N_GENES   <- 400L
N_DE      <- 40L
B         <- 10L

d <- make_mv_counts(seed = DATA_SEED, n_genes = N_GENES, n_de = N_DE)
run <- mv_run_multiverse(d$counts, d$metadata, "condition", ref_level = "C", test_level = "T",
                         B = B, seed = RUN_SEED, include_apeglm = TRUE,
                         padj_cutoff = 0.05, lfc_cutoff = 1, target_efdr = 0.10)

stability <- run$stability
called <- stability$gene[stability$called]
true_de <- setdiff(stability$gene, d$null_genes)

cat("\n--- Multiverse run summary (Figure 3 source) ---\n")
cat(sprintf("specifications        : %d\n", nrow(run$run$specifications)))
cat(sprintf("bootstrap replicates  : %d\n", B))
cat(sprintf("selected tau          : %s\n", format(run$selected_tau)))
cat(sprintf("calls at selected tau : %d\n", length(called)))
cat(sprintf("realised FDP          : %.4f\n", if (length(called)) mean(called %in% d$null_genes) else NA_real_))
cat(sprintf("true DE recovered     : %d of %d\n", sum(called %in% true_de), length(true_de)))
cat(sprintf("distinct eFDR values  : %d\n", length(unique(stats::na.omit(run$efdr_curve$efdr)))))
cat("\neFDR curve:\n"); print(run$efdr_curve, row.names = FALSE)

# Gene selection rule (fixed in advance, not hand-picked): among the genes
# called at the selected threshold, take the highest-ranked gene on which the
# multiverse is not unanimous (stability < 1), ordered by eFDR, then stability,
# then gene id. A gene called in every specification produces a degenerate curve
# with no visible disagreement, which is not what the display exists to show.
ranked <- stability[stability$called & stability$stability < 1, , drop = FALSE]
ranked <- ranked[order(ranked$efdr, -ranked$stability, ranked$gene), , drop = FALSE]
gene <- ranked$gene[1]
cat("\nstability of called genes:\n"); print(table(round(stability$stability[stability$called], 3)))
cat(sprintf("\nplotted gene          : %s (stability %.3f, consistency %.3f, eFDR %.4f, true DE: %s)\n",
            gene, stability$stability[stability$gene == gene],
            stability$direction_consistency[stability$gene == gene],
            stability$efdr[stability$gene == gene], gene %in% true_de))

p <- mv_plot_specification_curve(run$run, gene, padj_cutoff = 0.05, lfc_cutoff = 1)

# Nature double-column width (183 mm) via the journal-ready cairo_pdf device.
width_in <- 183 / 25.4
height_in <- 130 / 25.4
pdf_path <- file.path(repo, "fig3_specification_curve.pdf")
png_path <- file.path(repo, "fig3_specification_curve.png")
ggplot2::ggsave(pdf_path, p, device = grDevices::cairo_pdf, width = width_in, height = height_in)
ggplot2::ggsave(png_path, p, device = grDevices::png, type = "cairo",
                width = width_in, height = height_in, dpi = 600)
cat(sprintf("\nwrote %s\nwrote %s\n", pdf_path, png_path))
