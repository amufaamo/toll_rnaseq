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
# The loader lives in scripts/validate_datasets.R.

source(file.path("/mnt/g/マイドライブ/toll_rnaseq", "scripts", "validate_common.R"))

d <- mv_dataset_pasilla()
mv_validate(dataset = d$dataset, counts = d$counts, metadata = d$metadata,
            condition_col = d$condition_col, ref_level = d$ref_level,
            test_level = d$test_level, covariate = d$covariate,
            symbols = d$symbols, markers = d$markers)
