#!/usr/bin/env Rscript
# Public-dataset validation 2/3: pasilla (Brooks et al. 2011, GSE18508).
#
#   conda run -n toll Rscript scripts/validate_pasilla.R
#
# Drosophila S2 cells, RNAi knockdown of pasilla vs untreated (3 treated,
# 4 untreated). The counts and sample annotation are the files shipped in the
# Bioconductor experiment package `pasilla`
# (inst/extdata/pasilla_gene_counts.tsv, inst/extdata/pasilla_sample_annotation.csv).
# The package itself does not install in this environment because its DEXSeq
# dependency needs system libraries, so the two files were taken from the
# package source tarball
# https://bioconductor.org/packages/3.22/data/experiment/src/contrib/pasilla_1.38.0.tar.gz
# and are stored under validation/data/.
# Contrast: condition treated vs untreated, with library type as covariate.

source(file.path("/mnt/g/マイドライブ/toll_rnaseq", "scripts", "validate_common.R"))
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

# org.Dm.eg.db maps FlyBase gene ids to symbols.
symbols <- setNames(rep(NA_character_, nrow(counts)), rownames(counts))
if (requireNamespace("org.Dm.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
  m <- suppressMessages(AnnotationDbi::select(org.Dm.eg.db::org.Dm.eg.db,
                                              keys = rownames(counts), keytype = "FLYBASE",
                                              columns = "SYMBOL"))
  m <- m[!duplicated(m$FLYBASE), , drop = FALSE]
  symbols[m$FLYBASE] <- m$SYMBOL
}

# ps (pasilla) is the RNAi target of Brooks et al. 2011 and must be knocked down.
markers <- c("ps")

mv_validate(dataset = "pasilla", counts = counts, metadata = metadata,
            condition_col = "condition", ref_level = "untreated", test_level = "treated",
            covariate = "type", symbols = symbols, markers = markers)
