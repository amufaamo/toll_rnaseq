#!/usr/bin/env Rscript
# Public-dataset validation 3/3: GSE60450 (Fu et al. 2015, Nat Cell Biol 17:365).
#
#   conda run -n toll Rscript scripts/validate_gse60450.R
#
# Mouse mammary gland, basal and luminal cell populations at three
# developmental stages (virgin, 18.5 day pregnancy, 2 day lactation), two
# animals per group (12 samples). This is the dataset used by the limma/edgeR
# "RNA-seq analysis is easy as 1-2-3" workflow (Law et al. 2016, F1000Research).
#
# Counts: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz
# Sample annotation is parsed from the GEO SOFT record saved at
# validation/data/GSE60450_samples.soft, fetched from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450&targ=gsm&form=text&view=quick
#
# Primary contrast: lactation vs pregnancy, pooling both cell populations, with
# the cell population as covariate (4 vs 4, balanced, 16 specifications).
# Secondary contrast: the same comparison within the basal population only
# (2 vs 2, no covariate, 8 specifications), run to see how the calibration
# behaves at the smallest group size the engine accepts.

source(file.path("/mnt/g/マイドライブ/toll_rnaseq", "scripts", "validate_common.R"))
data_dir <- file.path(MV_REPO, "validation", "data")

raw <- read.delim(gzfile(file.path(data_dir, "GSE60450_Lactation-GenewiseCounts.txt.gz")),
                  check.names = FALSE)
counts <- as.matrix(raw[, -(1:2)])
rownames(counts) <- as.character(raw$EntrezGeneID)
storage.mode(counts) <- "integer"
# Column headers are "<sample name>_<flowcell>_<barcode>_<lane>_R1".
colnames(counts) <- sub("_.*$", "", colnames(counts))

soft <- readLines(file.path(data_dir, "GSE60450_samples.soft"))
sample_name <- sub("^!Sample_description = Sample name: ", "",
                   grep("^!Sample_description = Sample name: ", soft, value = TRUE))
source_name <- sub("^!Sample_source_name_ch1 = ", "",
                   grep("^!Sample_source_name_ch1 = ", soft, value = TRUE))
stopifnot(length(sample_name) == 12, length(source_name) == 12)
ann <- data.frame(sample = sample_name, source = source_name, stringsAsFactors = FALSE)
ann$cell_type <- ifelse(grepl("basal cells", ann$source), "basal", "luminal")
ann$stage <- ifelse(grepl("virgin", ann$source), "virgin",
                    ifelse(grepl("pregnancy", ann$source), "pregnant", "lactate"))
ann <- ann[match(colnames(counts), ann$sample), , drop = FALSE]
stopifnot(identical(ann$sample, colnames(counts)), !anyNA(ann$stage))
cat("GEO sample annotation parsed from SOFT:\n"); print(ann[, c("sample", "cell_type", "stage")], row.names = FALSE)

metadata <- data.frame(sample = ann$sample, stage = ann$stage,
                       cell_type = ann$cell_type, stringsAsFactors = FALSE)

symbols <- setNames(rep(NA_character_, nrow(counts)), rownames(counts))
if (requireNamespace("org.Mm.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
  m <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
                                              keys = rownames(counts), keytype = "ENTREZID",
                                              columns = "SYMBOL"))
  m <- m[!duplicated(m$ENTREZID), , drop = FALSE]
  symbols[m$ENTREZID] <- m$SYMBOL
}

# Milk-protein genes that switch on at lactation, plus Mcl1, the focus of
# Fu et al. 2015 ("EGF-mediated induction of Mcl-1 at the switch to lactation").
markers <- c("Csn2", "Csn1s2b", "Wap", "Lalba", "Mcl1")

keep_stage <- metadata$stage %in% c("pregnant", "lactate")
mv_validate(dataset = "GSE60450", counts = counts[, keep_stage, drop = FALSE],
            metadata = metadata[keep_stage, , drop = FALSE],
            condition_col = "stage", ref_level = "pregnant", test_level = "lactate",
            covariate = "cell_type", symbols = symbols, markers = markers)

keep_basal <- keep_stage & metadata$cell_type == "basal"
mv_validate(dataset = "GSE60450basal", counts = counts[, keep_basal, drop = FALSE],
            metadata = metadata[keep_basal, , drop = FALSE],
            condition_col = "stage", ref_level = "pregnant", test_level = "lactate",
            covariate = NULL, symbols = symbols, markers = markers)
