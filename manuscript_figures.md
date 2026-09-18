# MultiverseDEG: Manuscript Figures and Tables

Components prepared for the Software article submitted to *Briefings in Bioinformatics*
(manuscript: `manuscript_draft.md`; legends: `figure_legends.md`).

Figure sources in this repository:

| Item | Source | Rendered artefact |
| :--- | :--- | :--- |
| Figure 1 | `manuscript_fig1.mmd` (embedded verbatim below) | render with Mermaid |
| Figure 2 | `manuscript_fig2.mmd` (embedded verbatim below) | render with Mermaid |
| Figure 3 | `scripts/generate_fig3_specification_curve.R` | `fig3_specification_curve.pdf`, `fig3_specification_curve.png` |

## Figure 1: Architecture of MultiverseDEG v4.0

```mermaid
graph TD
    classDef input fill:#ecdbba,stroke:#e65100,stroke-width:2px,color:#333;
    classDef module fill:#2d4263,stroke:#0f3460,stroke-width:1.5px,color:#fff;
    classDef mv fill:#d55e00,stroke:#7f2704,stroke-width:3px,color:#fff;
    classDef state fill:#e8f4f8,stroke:#0f3460,stroke-width:2px,color:#333;
    classDef async fill:#ffb400,stroke:#c84b31,stroke-width:2px,color:#333;
    classDef output fill:#ecdbba,stroke:#e65100,stroke-width:2px,color:#333;

    IN["Gene-by-sample integer count matrix<br/>+ sample metadata table"]:::input --> M1

    subgraph APP["MultiverseDEG v4.0 - golem R package, Shiny with bslib v5 UI"]
        direction TB
        M1["1 Upload<br/>format validation"]:::module
        M2["2 QC and Filter<br/>expression thresholds"]:::module
        M3["3 EDA<br/>PCA, sample clustering"]:::module
        M4["4 DEG<br/>DESeq2 Wald + apeglm, edgeR"]:::module
        M5["5 Multiverse DEG<br/>specification grid, direction-consistent stability,<br/>calibrated eFDR"]:::mv
        M6["6 Enrichment<br/>clusterProfiler over-representation"]:::module
        M7["7 GSVA"]:::module
        M8["8 UpSet<br/>DEG sets across biological contrasts"]:::module
        M9["9 Time-series<br/>maSigPro"]:::module
        M10["10 Report"]:::module

        M1 --> M2 --> M3 --> M4 --> M5
        M5 -->|"multiverse:TEST_vs_REF, result type multiverse_stability"| M6
        M6 --> M7 --> M8 --> M9 --> M10

        STATE[("AppState (R6)<br/>uploaded data, per-module results,<br/>step-status map, parameter log")]:::state
        M1 <--> STATE
        M5 <--> STATE
        M10 <--> STATE

        ASYNC["ExtendedTask with future / promises<br/>bounded batches in background R processes"]:::async
        M4 -.-> ASYNC
        M5 -.-> ASYNC
        M6 -.-> ASYNC
        M7 -.-> ASYNC
        M9 -.-> ASYNC
    end

    M3 --> EXPORT
    M5 --> EXPORT
    M6 --> EXPORT
    EXPORT["Journal-ready export<br/>PDF via cairo_pdf and SVG<br/>Nature / Cell / Science page-width presets"]:::output
    M10 --> HTML["HTML report"]:::output
    STATE --> SCRIPT["Reproducible commented R script"]:::output
```

## Figure 2: The multiverse DEG workflow

```mermaid
graph TD
    classDef obs fill:#2d4263,stroke:#0f3460,stroke-width:1.5px,color:#fff;
    classDef null fill:#0072b2,stroke:#023858,stroke-width:1.5px,color:#fff;
    classDef calib fill:#d55e00,stroke:#7f2704,stroke-width:2px,color:#fff;
    classDef output fill:#ecdbba,stroke:#e65100,stroke-width:2px,color:#333;

    A["One fixed biological contrast<br/>raw integer counts, metadata"]:::obs

    A --> B["Build the specification grid<br/>filter x method x shrinkage x covariate<br/>family-legal combinations only"]:::obs
    B --> B2["Prune rank-deficient designs, confounded covariates<br/>and filters leaving too few genes; record the reason"]:::obs
    B2 --> C["Fit every retained specification<br/>DESeq2 Wald with and without apeglm,<br/>edgeR QL F-test, edgeR LRT (TMM)<br/>store log2FC, p and BH-adjusted p"]:::obs
    C --> D["Per-gene stability T = fraction of specifications<br/>calling the gene at padj &lt;= q and |log2FC| &gt;= a<br/>direction consistency C &gt;= 0.90 required"]:::obs

    A --> E["Fit the full negative-binomial model<br/>condition + covariates; assert the reconstructed<br/>mean equals the DESeq2 mu assay"]:::null
    E --> F["Null means: set only the condition coefficient to zero<br/>keep intercept, covariate coefficients<br/>and count-scale normalisation factors"]:::null
    F --> G["B parametric NB bootstrap replicates<br/>every filter and specification re-run from scratch<br/>retain only null discovery counts R_b(tau)"]:::null

    D --> H
    G --> H["eFDR(tau) = [(1 + sum_b R_b(tau)) / (B + 1)] / max(R_obs(tau), 1)<br/>evaluated at every attainable tau = k/m,<br/>made monotone by a cumulative minimum in ascending tau"]:::calib
    H --> I["Select the smallest positive tau meeting the eFDR target<br/>tau = 0 excluded; re-thresholding refits nothing"]:::calib

    I --> J["Stability-ranked gene table<br/>with per-gene eFDR"]:::output
    I --> K["Per-gene specification curve<br/>(Figure 3)"]:::output
    J --> L["Hand-off to the Enrichment and UpSet modules<br/>key multiverse:TEST_vs_REF<br/>result type multiverse_stability, pvalue = NA"]:::output
```

## Figure 3: Per-gene specification curve

Vector artefact: `fig3_specification_curve.pdf` (`cairo_pdf`, 183 x 130 mm, *Nature*
double-column width); raster preview: `fig3_specification_curve.png` (600 dpi).

The figure is produced by `scripts/generate_fig3_specification_curve.R`, which runs the
real engine (`R/multiverse_engine.R`) and the real plotting function
(`R/multiverse_plots.R::mv_plot_specification_curve()`) on the repository's synthetic
fixture (`tests/testthat/helper-multiverse.R::make_mv_counts()`), then writes the figure:

```
conda run -n toll Rscript scripts/generate_fig3_specification_curve.R
```

Run configuration and the summary statistics obtained (not hand-edited):

| Quantity | Value |
| :--- | :--- |
| Fixture | 400 genes, 40 true DE genes (log2FC +/- 1.5), n = 6 per group, seed 20260918 |
| Grid | 8 specifications (2 filters x {DESeq2 Wald, DESeq2 Wald + apeglm, edgeR QL/TMM, edgeR LRT/TMM}) |
| Bootstrap | B = 10, engine seed 1 |
| Call rule | adjusted *p* <= 0.05, \|log2FC\| >= 1, direction consistency >= 0.90, target eFDR 0.10 |
| Selected threshold | tau = 0.125 |
| Calls at the selected threshold | 33 |
| Realised false discovery proportion | 0.000 |
| True DE genes recovered | 33 of 40 |
| Gene plotted | g10 (stability 0.750, direction consistency 1.000, eFDR 0.0085; a true DE gene) |

Gene selection rule, fixed before inspection: among the genes called at the selected
threshold, the highest-ranked gene on which the multiverse is *not* unanimous
(stability < 1), ordered by eFDR, then stability, then gene identifier. A gene called in
every specification yields a degenerate curve that shows no disagreement.

**Resolved 2026-09-19.** An earlier revision of `manuscript_draft.md` reported 22 calls, a
realised false discovery proportion of 0.045 and 21 of 40 true DE genes recovered for a
run described with these same nominal fixture parameters (400 genes, 40 true DE, n = 6/
group, B = 10). That number was never backed by a script or fixture committed to this
repository: it traces to an ad hoc run from an earlier work session whose exact generator
was not saved. A 39-combination sweep over `make_mv_counts()` data seeds and engine seeds
found no seed reproducing it; the manuscript's numbers were consistent with a weaker true
effect size (log2FC of about 1.0) that does not exist in the repository. Rather than add an
unrecorded fixture to match an unreproducible number, `manuscript_draft.md`,
`cover_letter.md` and this file were updated to the number actually and reproducibly
obtained from the committed code and fixture, recorded in the table above.

## Table 1: Comparison with graphical RNA-seq analysis tools

| Feature | **MultiverseDEG v4.0** | **iDEP** [2] | **Galaxy** [1] | **DEBrowser** [3] |
| :--- | :---: | :---: | :---: | :---: |
| **Analysis entry point** | Count matrix + metadata | Count matrix | Raw reads or count matrix (tool collection) | Count matrix |
| **Deployment** | Local or self-hosted server (R package / Docker image) | Web application | Shared or cloud infrastructure | Not assessed |
| **Reproducible analysis script generation** | Yes | Yes | Not assessed | Not assessed |
| **Explicit multi-specification (multiverse) robustness reporting with a calibrated error rate** | **Yes** | No | No | No |

Notes, so the table is not read as claiming more than was checked:

- Comparators are restricted to the graphical tools cited in the manuscript. The last row
  reflects the prior-art comparison in the Introduction and in *Relation to existing
  tools*: these platforms document one selected analysis path and do not report the
  dependence of the result on that selection.
- Reproducible script generation is **not** a differentiator. iDEP already records the
  selected parameters and emits R / R Markdown code, and the manuscript makes no priority
  claim for this feature.
- Asynchronous computation, colourblind-safe defaults, journal page-width export presets,
  time-series support and input size limits are deliberately **not** tabulated: they are
  implemented in MultiverseDEG, but their presence or absence in the comparator tools was
  not verified, and an unverified cross in a comparison table is an unsupported claim.
- Rows describing an upstream FASTQ-to-counts pipeline and *de novo* assembly, present in
  earlier drafts, are removed. v4.0 has no upstream processing; it starts from a count
  matrix.

## Table 2: Measured single-fit cost

Simulated count matrix of 20,000 genes x 12 samples; elapsed time for one model fit.

| Operation | Elapsed time |
| :--- | :---: |
| DESeq2 (Wald) | 4.35 s |
| edgeR quasi-likelihood F-test | 2.54 s |

## Table 3: Cost of a calibrated multiverse run

Derived from Table 2 at a mixed mean of about 2.5 s per fit. A calibrated run is
approximately (*B* + 1) x \|*S*\| model fits plus the null fit and simulation.

| Configuration | Specifications \|*S*\| | Bootstrap *B* | Model fits | Approximate serial CPU time |
| :--- | :---: | :---: | :---: | :---: |
| Implemented grid, no covariate | 8 | 20 | 168 | ~7 CPU-minutes |
| Implemented grid, one covariate | 16 | 20 | 336 | ~14 CPU-minutes |
| Full design | 80 | 100 | 8,080 | ~5.6 CPU-hours |

These figures are hardware- and data-dependent and are reported to make the cost structure
explicit, not as a benchmark against other software. Multiverse analysis with
simulation-based calibration is one to three orders of magnitude more expensive than a
single pipeline, which is why the batched asynchronous scheduler and the no-refit
re-thresholding design are load-bearing.

## References cited in the tables

[1] Afgan, E., et al. (2018). *Nucleic Acids Research*, 46(W1), W537-W544.
[2] Ge, S. X., Son, E. W., & Yao, R. (2018). *BMC Bioinformatics*, 19(1), 534.
[3] Kucukural, A., et al. (2019). *BMC Genomics*, 20(1), 6.
