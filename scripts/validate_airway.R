#!/usr/bin/env Rscript
# Public-dataset validation 1/3: airway (Himes et al. 2014, GSE52778).
#
#   conda run -n toll Rscript scripts/validate_airway.R
#
# Human airway smooth muscle cells, four donor cell lines each treated with
# dexamethasone or untreated (paired design, 8 samples). Data come from the
# Bioconductor experiment package `airway` (BiocManager::install("airway")).
# Contrast: dex trt vs untrt, with the donor cell line as covariate (16 specs).
# The loader lives in scripts/validate_datasets.R.

source(file.path("/mnt/g/マイドライブ/toll_rnaseq", "scripts", "validate_common.R"))

d <- mv_dataset_airway()
cat(d$note, "\n")
mv_validate(dataset = d$dataset, counts = d$counts, metadata = d$metadata,
            condition_col = d$condition_col, ref_level = d$ref_level,
            test_level = d$test_level, covariate = d$covariate,
            symbols = d$symbols, markers = d$markers)
