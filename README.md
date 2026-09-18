# MultiverseDEG v4.0

MultiverseDEG is a [golem](https://thinkr-open.github.io/golem/)-structured R package
providing a Shiny GUI for **downstream** bulk RNA-seq analysis: from a count matrix to
publication-ready figures, without writing R code.

> **Scope.** v4 starts from a **count matrix**. It does *not* perform upstream processing
> (FASTQ QC, alignment, or quantification) — there is no fastp/STAR/Salmon/featureCounts
> step in this package. Generate counts with your own pipeline first, then upload the
> matrix here. (The unrelated legacy v3 application in `app.R` / `after_count/` /
> `before_count/` is retained in the repository for reference only and is not part of
> this package.)

## Analysis steps

The UI is a `bslib` (Bootstrap 5) `page_navbar` whose tabs follow the pipeline order.
Each tab stays locked with an empty-state card until its prerequisite step is done.

| # | Tab | Module | What it does |
|---|---|---|---|
| 1 | Upload | `mod_upload.R` | Count matrix (CSV/TSV/RDS) + sample metadata; format auto-detection |
| 2 | QC & Filter | `mod_qc.R` | Low-count gene filtering, library-size distribution |
| 3 | EDA | `mod_eda.R` | PCA, scree plot, sample-distance heatmap |
| 4 | DEG | `mod_deg.R` | DESeq2 Wald test, apeglm LFC shrinkage, optional ComBat-seq batch correction |
| 5 | Multiverse DEG | `mod_deg_multiverse.R` | Stability-ranked DEGs across a grid of defensible specifications, with bootstrap-calibrated effective FDR |
| 6 | Enrichment | `mod_enrichment.R` | GSEA (fgsea + MSigDB), GO/KEGG over-representation (clusterProfiler) |
| 7 | GSVA | `mod_gsva.R` | Single-sample pathway scoring |
| 8 | UpSet | `mod_deg_multi.R` | Intersection of DEG sets across biological contrasts |
| 9 | Time-series | `mod_timeseries.R` | maSigPro temporal expression modelling |
| 10 | Report | `mod_report.R` | Reproducible R script + HTML summary of every logged parameter |

Two cross-cutting components are not tabs:

- `mod_export.R` — a shared "Journal-Ready Export" modal (Nature/Cell/Science/Custom
  presets; PDF-cairo, SVG, TIFF 300 dpi, PNG 300 dpi; font and size controls).
- `class_AppState.R` — an R6 object holding all application state in `reactiveVal`
  wrappers, so state can be driven and inspected outside a Shiny session.

### Multiverse DEG

Rather than reporting one workflow, this tab runs every legal combination of filter ×
covariate × DE method for a single contrast and ranks genes by how consistently they are
called in the same direction. A parametric negative-binomial bootstrap from the fitted
condition-null model calibrates an effective FDR for the stability threshold. The
compute core (`R/multiverse_engine.R`) is Shiny-free and unit-tested; the design
rationale and its open questions are documented in
[`docs/design_multiverse_deg.md`](docs/design_multiverse_deg.md).

## Installation

Development happens in the conda environment `toll` (R 4.5.2).

```bash
conda activate toll
```

Hard dependencies are declared in `DESCRIPTION` under `Imports`
(`shiny`, `bslib`, `bsicons`, `golem`, `R6`, `DT`, `ggplot2`, `future`, `promises`,
`config`). Analysis packages are `Suggests`: every module checks with
`requireNamespace()` and reports a clear message instead of failing, so the app runs
with only the subset you need.

```r
# Imports
install.packages(c("shiny", "bslib", "bsicons", "golem", "R6", "DT",
                   "ggplot2", "future", "promises", "config"))

# Suggests (install what the tabs you use require)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2", "edgeR", "apeglm", "SummarizedExperiment",   # DEG + Multiverse DEG
  "sva",                                                  # ComBat-seq batch correction
  "fgsea", "msigdbr", "clusterProfiler", "enrichplot",    # Enrichment
  "org.Hs.eg.db", "org.Mm.eg.db",                         # Gene ID mapping
  "GSVA", "ComplexUpset", "maSigPro",                     # GSVA / UpSet / Time-series
  "ComplexHeatmap", "EnhancedVolcano", "ggrepel"          # Figures
))
```

## Running the app

From a clone of this repository:

```r
devtools::load_all(".")
run_app()
```

Or, once the package is installed:

```r
MultiverseDEG::run_app()
```

`run_app()` passes `...` through to `golem::with_golem_options()`.

## Input format

**Count matrix** — CSV, TSV, or RDS. Rows are genes, columns are samples, the first
column holds gene IDs, and a header row is required. Values must be raw integer counts
(not TPM/FPKM/CPM): DESeq2, edgeR, and the multiverse engine all reject non-integer
input.

**Sample metadata** — CSV or TSV. The **first column** holds sample names that match the
count matrix column names; remaining columns are experimental factors (condition, batch,
timepoint, …). You can also auto-generate a metadata skeleton from the count matrix
column names and edit the condition column.

`cd4_cd14.tsv` in the repository root is a human CD4 T-cell vs CD14 monocyte count
matrix (~26,000 genes) useful for exercising every tab.

## Development

```r
devtools::load_all(".")     # load the package
devtools::document()        # regenerate NAMESPACE and man/ from roxygen comments
devtools::test()            # run the testthat suite
```

Tests use testthat edition 3 and live in `tests/testthat/`. Tests that fit real models
skip themselves when `DESeq2`/`edgeR` are unavailable.

Because most analysis packages are `Suggests`, run `R CMD check` with suggested-package
enforcement disabled unless you have all of them installed. `_R_CHECK_CRAN_INCOMING_` is
also disabled: this package is not destined for CRAN, and that check only reports the
`.9000` development version suffix.

```bash
_R_CHECK_FORCE_SUGGESTS_=false _R_CHECK_CRAN_INCOMING_=false \
  Rscript -e 'rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))'
```

`.github/workflows/ci.yml` runs the same check on pushes and pull requests to `main` and
`v4.0-dev`.

The check is clean apart from two expected NOTEs:

- *"Namespace in Imports field not imported from: `future`"* — `future` is never called
  directly, but `promises::future_promise()` needs a `future` backend at runtime and
  `promises` only suggests it, so the dependency must be declared here.
- *"unable to verify current time"* — emitted when the check host cannot reach the
  network time service; it says nothing about the package.

## Repository layout

```
toll_rnaseq/
├── DESCRIPTION, NAMESPACE      # package metadata
├── R/                          # app shell, R6 state, modules, multiverse engine
├── inst/                       # golem-config.yml, app/www assets
├── man/                        # generated documentation
├── tests/testthat/             # test suite
├── docs/                       # design notes
├── cd4_cd14.tsv                # example count matrix
└── app.R, after_count/, before_count/, R_v3_backup/   # legacy v3, not packaged
```

## License

MIT — see [`LICENSE.md`](LICENSE.md).
