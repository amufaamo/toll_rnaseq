# Dataset loaders shared by the public-dataset validation scripts.
#
# Each loader returns the arguments mv_validate() (scripts/validate_common.R)
# and mv_run_observed() (R/multiverse_engine.R) need for one run, so that the
# per-dataset validation scripts, the call-rule / eFDR-target sweeps
# (scripts/validate_tau_conditions.R) and the specification-concordance
# analysis (scripts/analyze_spec_concordance.R) all read the same data the
# same way.
#
# Requires MV_REPO to be defined (scripts/validate_common.R sets it).

mv_dataset_airway <- function() {
  suppressPackageStartupMessages({
    library(airway)
    library(SummarizedExperiment)
  })
  data(airway, package = "airway", envir = environment())
  airway <- get("airway", envir = environment())
  counts <- SummarizedExperiment::assay(airway)
  cd <- as.data.frame(SummarizedExperiment::colData(airway))
  metadata <- data.frame(sample = rownames(cd),
                         dex = as.character(cd$dex),
                         cell = as.character(cd$cell),
                         stringsAsFactors = FALSE)
  symbols <- setNames(as.character(SummarizedExperiment::rowData(airway)$symbol),
                      rownames(airway))
  # Glucocorticoid-response genes reported as induced by dexamethasone in Himes
  # et al. 2014 (PLoS ONE 9:e99625); CRISPLD2 is that paper's focus gene.
  list(dataset = "airway", counts = counts, metadata = metadata,
       condition_col = "dex", ref_level = "untrt", test_level = "trt",
       covariate = "cell", symbols = symbols,
       markers = c("CRISPLD2", "DUSP1", "KLF15", "PER1", "TSC22D3", "FKBP5", "ZBTB16", "CCL2"),
       note = sprintf("airway package version: %s", as.character(packageVersion("airway"))))
}

mv_dataset_pasilla <- function() {
  data_dir <- file.path(MV_REPO, "validation", "data")
  counts <- as.matrix(read.delim(file.path(data_dir, "pasilla_gene_counts.tsv"), row.names = "gene_id"))
  storage.mode(counts) <- "integer"

  ann <- read.csv(file.path(data_dir, "pasilla_sample_annotation.csv"), stringsAsFactors = FALSE)
  ann$sample <- sub("fb$", "", ann$file)
  ann <- ann[match(colnames(counts), ann$sample), , drop = FALSE]
  stopifnot(identical(ann$sample, colnames(counts)))
  metadata <- data.frame(sample = ann$sample,
                         condition = ann$condition,
                         type = sub("-", "_", ann$type),
                         stringsAsFactors = FALSE)

  symbols <- setNames(rep(NA_character_, nrow(counts)), rownames(counts))
  if (requireNamespace("org.Dm.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
    m <- suppressMessages(AnnotationDbi::select(org.Dm.eg.db::org.Dm.eg.db,
                                                keys = rownames(counts), keytype = "FLYBASE",
                                                columns = "SYMBOL"))
    m <- m[!duplicated(m$FLYBASE), , drop = FALSE]
    symbols[m$FLYBASE] <- m$SYMBOL
  }

  # ps (pasilla) is the RNAi target of Brooks et al. 2011 and must be knocked down.
  list(dataset = "pasilla", counts = counts, metadata = metadata,
       condition_col = "condition", ref_level = "untreated", test_level = "treated",
       covariate = "type", symbols = symbols, markers = "ps", note = NULL)
}

.mv_gse60450_raw <- function() {
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
  list(counts = counts, metadata = metadata, symbols = symbols, ann = ann)
}

# Milk-protein genes that switch on at lactation, plus Mcl1, the focus of
# Fu et al. 2015 ("EGF-mediated induction of Mcl-1 at the switch to lactation").
MV_GSE60450_MARKERS <- c("Csn2", "Csn1s2b", "Wap", "Lalba", "Mcl1")

mv_dataset_gse60450 <- function() {
  x <- .mv_gse60450_raw()
  keep <- x$metadata$stage %in% c("pregnant", "lactate")
  list(dataset = "GSE60450", counts = x$counts[, keep, drop = FALSE],
       metadata = x$metadata[keep, , drop = FALSE],
       condition_col = "stage", ref_level = "pregnant", test_level = "lactate",
       covariate = "cell_type", symbols = x$symbols, markers = MV_GSE60450_MARKERS,
       note = NULL, annotation = x$ann[, c("sample", "cell_type", "stage")])
}

mv_dataset_gse60450_basal <- function() {
  x <- .mv_gse60450_raw()
  keep <- x$metadata$stage %in% c("pregnant", "lactate") & x$metadata$cell_type == "basal"
  list(dataset = "GSE60450basal", counts = x$counts[, keep, drop = FALSE],
       metadata = x$metadata[keep, , drop = FALSE],
       condition_col = "stage", ref_level = "pregnant", test_level = "lactate",
       covariate = NULL, symbols = x$symbols, markers = MV_GSE60450_MARKERS,
       note = NULL)
}

# The four runs used throughout the validation, keyed by the short name the
# command-line drivers accept.
MV_RUNS <- list(
  airway         = mv_dataset_airway,
  pasilla        = mv_dataset_pasilla,
  gse60450cov    = mv_dataset_gse60450,
  gse60450basal  = mv_dataset_gse60450_basal
)
