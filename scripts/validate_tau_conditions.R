#!/usr/bin/env Rscript
# Sweep of call rules and eFDR targets, to find conditions under which the
# automatically selected stability threshold tau is binding (i.e. strictly
# larger than the smallest attainable positive value 1/m).
#
#   conda run -n toll Rscript scripts/validate_tau_conditions.R <condition> <run>
#
# <condition> is one of the names of MV_CONDITIONS below, <run> one of the names
# of MV_RUNS (scripts/validate_datasets.R): airway, pasilla, gse60450cov,
# gse60450basal. One process handles one (condition, run) pair so that the grid
# can be run in parallel; results are written per pair under validation/results/
# with the condition name as filename suffix.
#
# The engine itself is untouched: only mv_run_multiverse() arguments vary.

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
})
MV_REPO <- "/mnt/g/マイドライブ/toll_rnaseq"
source(file.path(MV_REPO, "scripts", "validate_common.R"))

MV_CONDITIONS <- list(
  # Experiment A: stricter call rule at the original cost (B = 20, target 0.10).
  strict_callrule = list(padj_cutoff = 0.01, lfc_cutoff = 2, direction_cutoff = 0.90,
                         B = 20L, target_efdr = 0.10),
  # Experiment B: default call rule, stricter eFDR target. mv_run_multiverse()
  # refuses targets below 0.10 unless B >= 100.
  B100_target001  = list(padj_cutoff = 0.05, lfc_cutoff = 1, direction_cutoff = 0.90,
                         B = 100L, target_efdr = 0.01),
  # Experiment C: both at once.
  strict_B100_target001 = list(padj_cutoff = 0.01, lfc_cutoff = 2, direction_cutoff = 0.90,
                               B = 100L, target_efdr = 0.01)
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: validate_tau_conditions.R <condition> <run>", call. = FALSE)
cond_name <- args[1]; run_name <- args[2]
if (!cond_name %in% names(MV_CONDITIONS))
  stop(sprintf("unknown condition %s; expected one of %s", cond_name,
               paste(names(MV_CONDITIONS), collapse = ", ")), call. = FALSE)
if (!run_name %in% names(MV_RUNS))
  stop(sprintf("unknown run %s; expected one of %s", run_name,
               paste(names(MV_RUNS), collapse = ", ")), call. = FALSE)

cond <- MV_CONDITIONS[[cond_name]]
d <- MV_RUNS[[run_name]]()
if (!is.null(d$note)) cat(d$note, "\n")

mv_validate(dataset = d$dataset, counts = d$counts, metadata = d$metadata,
            condition_col = d$condition_col, ref_level = d$ref_level,
            test_level = d$test_level, covariate = d$covariate,
            symbols = d$symbols, markers = d$markers,
            B = cond$B, seed = MV_SEED,
            padj_cutoff = cond$padj_cutoff, lfc_cutoff = cond$lfc_cutoff,
            direction_cutoff = cond$direction_cutoff, target_efdr = cond$target_efdr,
            condition_label = cond_name, tag_suffix = cond_name)
