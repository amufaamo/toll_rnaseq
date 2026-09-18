#!/usr/bin/env Rscript
# Public-dataset validation 1/3: airway (Himes et al. 2014, GSE52778).
#
#   conda run -n toll Rscript scripts/validate_airway.R
#
# Human airway smooth muscle cells, four donor cell lines each treated with
# dexamethasone or untreated (paired design, 8 samples). Data come from the
# Bioconductor experiment package `airway` (BiocManager::install("airway")).
# Contrast: dex trt vs untrt, with the donor cell line as covariate (16 specs).

suppressPackageStartupMessages({
  library(airway)
  library(SummarizedExperiment)
})
source(file.path("/mnt/g/マイドライブ/toll_rnaseq", "scripts", "validate_common.R"))

data(airway)
cat(sprintf("airway package version: %s\n", as.character(packageVersion("airway"))))
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
markers <- c("CRISPLD2", "DUSP1", "KLF15", "PER1", "TSC22D3", "FKBP5", "ZBTB16", "CCL2")

mv_validate(dataset = "airway", counts = counts, metadata = metadata,
            condition_col = "dex", ref_level = "untrt", test_level = "trt",
            covariate = "cell", symbols = symbols, markers = markers)
