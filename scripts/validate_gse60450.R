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
# The loaders live in scripts/validate_datasets.R.

source(file.path("/mnt/g/マイドライブ/toll_rnaseq", "scripts", "validate_common.R"))

d <- mv_dataset_gse60450()
cat("GEO sample annotation parsed from SOFT:\n"); print(d$annotation, row.names = FALSE)
mv_validate(dataset = d$dataset, counts = d$counts, metadata = d$metadata,
            condition_col = d$condition_col, ref_level = d$ref_level,
            test_level = d$test_level, covariate = d$covariate,
            symbols = d$symbols, markers = d$markers)

b <- mv_dataset_gse60450_basal()
mv_validate(dataset = b$dataset, counts = b$counts, metadata = b$metadata,
            condition_col = b$condition_col, ref_level = b$ref_level,
            test_level = b$test_level, covariate = b$covariate,
            symbols = b$symbols, markers = b$markers)
