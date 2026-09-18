# Multiverse DEG analysis — design for MultiverseDEG v4

## 1. Motivation and positioning

Bulk RNA-seq DEG results depend on defensible choices: low-count filtering, normalisation, model family, nuisance-covariate adjustment, and calling thresholds. This module makes that choice set explicit, runs its legal members for one biological contrast, and returns a stability-ranked DEG list: a gene is useful downstream when it persists across the declared analysis multiverse. It follows multiverse analysis (Steegen et al., 2016, doi:10.1177/1745691616658637), specification curves (Simonsohn, Simmons & Nelson, 2020, doi:10.1038/s41562-020-0912-z), and vibration of effects (Patel, Burford & Ioannidis, 2015, doi:10.1016/j.jclinepi.2015.05.029).

This differs from iDEP and similar tools. iDEP records and reproduces one selected workflow, including selected parameters and generated R/R Markdown (Ge, Son & Yao, 2018, doi:10.1186/s12859-018-2486-6). Multiverse DEG asks instead: how much does this conclusion depend on which defensible workflow was selected? It neither treats paths as independent replications nor uses the union of all paths as discoveries. It also adds conservative, dataset-specific empirical-null calibration for stability.

It differs from the existing mod_deg_multi.R UpSet module. UpSet compares DEG sets from different biological contrasts. Multiverse DEG fixes one contrast and varies analytic specifications. A multiverse result can later be one named deg_results entry for comparison with biological contrasts in UpSet, but its name starts multiverse: so it cannot be mistaken for a biological contrast.

## 2. The multiverse: dimensions and legal specifications

The engine starts from integer raw counts, not counts_filtered(), because filtering is varied. It aligns samples to metadata as mod_deg.R does.

| Dimension | Default candidates | Why defensible |
|---|---|---|
| Filter | count >=10 in >=2 samples; >=10 in the smaller group size; >=20 in the smaller group size; edgeR filterByExpr using the design | Transparent abundance/prevalence screens and edgeR's design-aware filter are standard bulk-RNA-seq choices. |
| Normalisation | DESeq2 median-of-ratios; edgeR TMM/RLE/upper-quartile; limma-voom with edgeR TMM factors | Established composition-bias normalisations, each consumed only by its native family. |
| DE method | DESeq2 Wald; DESeq2 LRT; edgeR QL F-test; edgeR LRT; limma-voom moderated t | Mature negative-binomial and precision-weighted linear-model approaches. limma-voom is full-grid only, deferred from the Phase-2 MVP. |
| LFC shrinkage | DESeq2 Wald: apeglm on/off; all other arms: not applicable | apeglm shrinks a DESeq2 coefficient, not edgeR/limma estimates. |
| Covariate | condition-only; condition plus one user-nominated metadata covariate | A measured non-confounded nuisance covariate belongs in the design matrix rather than altered counts. |
| Optional ComBat-seq | Full-version, pre-filter sensitivity arm only when requireNamespace("sva", quietly=TRUE) | Never an MVP requirement; sva is absent in the verified environment. |

IHW is not an MVP axis: it is unavailable here and would make per-specification adjusted p-values less comparable. A later BH-versus-IHW extension is possible only behind requireNamespace("IHW", quietly=TRUE).

### Cardinality and cap

The full default is

\[
4\ {\rm filters}\times2\ {\rm covariate\ choices}\times
[2\ {\rm DESeq2\ Wald\ shrinkage\ states}+1\ {\rm DESeq2\ LRT}+
(3\ {\rm edgeR\ normalisations}\times2\ {\rm edgeR\ tests})+1\ {\rm voom}]
=80\ {\rm specifications}.
\]

This is a product over legal family-specific choices, not a false Cartesian product. The UI displays its count, warns above 100, and refuses more than 200 specifications in the MVP. A later option may sample uniformly within filter × covariate × method-family strata, retain inclusion probabilities, and use inverse-probability weights. This is unbiased for the declared full-grid average conditional on fixed selections, but increases variance and is not default.

The Phase-2 MVP preset is two filters (10 in the smaller group and filterByExpr) × four method paths (DESeq2 Wald unshrunk; DESeq2 Wald apeglm if available; edgeR QL/TMM; edgeR LRT/TMM) × covariate arms. The covariate axis is conditional, not a fixed 2: with no covariate column selected there is nothing to toggle, so `mv_build_specifications()` emits a single covariate arm (`NA`), giving **8 paths**. Selecting a covariate column adds a second arm (with/without that covariate in the design), giving **16 paths** (rank-deficient covariate/condition combinations are pruned from that 16, per the table below). If apeglm is absent, the shrinkage arm is visibly omitted and both counts drop by one method path (8→6, 16→12).

Padj and absolute-LFC rules are not multiplicative DE-fit axes. Every fit stores LFC, p-value, and BH padj. The supported no-refit call-rule grid is padj 0.01/0.05/0.10 × absolute LFC 0/0.5/1; bootstrap curves are stored for every selected pair. An arbitrary display-only threshold receives a not-bootstrap-calibrated badge.

| Family | Legal combination | Prune |
|---|---|---|
| DESeq2 Wald | native median-of-ratios; Wald; apeglm on/off | Never pass TMM/RLE/UQ factors alongside DESeq2 size factors (double normalisation); apeglm only for a unique Wald coefficient. |
| DESeq2 LRT | native normalisation; full ~batch+condition vs reduced ~batch, or full ~condition vs reduced ~1 | No apeglm. LFC is unshrunk fitted condition coefficient; p-value is LRT. |
| edgeR QL/LRT | TMM/RLE/UQ plus QL F-test or LRT | No apeglm; never provide pre-normalised/non-integer data to DGEList. |
| limma-voom | edgeR TMM factors, then voom/lmFit/eBayes | No apeglm; only TMM in default full grid to avoid unbalanced expansion. |
| Any | condition-only or full-rank batch+condition | Exclude one-level/missing/confounded covariates, filters leaving <2 samples/group, and failed fits. |

All filters are rerun per specification and bootstrap. A filtered-out gene is not imputed.

## 3. Stability

For legal specification set \(S\), \(m=|S|\), call rule \((q,a)\), and gene \(g\),

\[
I^+_{gs}=1\{g\ {\rm tested\ in}\ s,\ padj_{gs}\le q,\ LFC_{gs}\ge a\},\quad
I^-_{gs}=1\{g\ {\rm tested\ in}\ s,\ padj_{gs}\le q,\ LFC_{gs}\le-a\}.
\]

The equal-weight MVP estimator is

\[
stab^+_g=m^{-1}\sum_{s\in S}I^+_{gs},\quad
stab^-_g=m^{-1}\sum_{s\in S}I^-_{gs},\quad
T_g=\max(stab^+_g,stab^-_g).
\]

Direction is the maximizing sign and direction consistency is \(C_g=T_g/(stab^+_g+stab^-_g)\) when nonzero. A default call needs \(T_g\ge\tau\) and \(C_g\ge0.90\). A gene significant-up in half and significant-down in half has T=0.5, C=0.5, and is not stable.

Equal weights are appropriate because the selected grid is the target population, not because paths are independent. It can over-weight a family with many variants; show family counts and keep the default balanced. Cluster/inverse-family weighting is deferred. Filtered-out, untestable, or failed-fit genes receive zero indicators with denominator still m: missing is not significant. This is conservative for low-expression genes.

## 4. Calibrated effective FDR

### Adopted algorithm: full-grid parametric NB bootstrap

The result is called eFDR, not a formal FDR guarantee. Fit the **full** DESeq2 negative-binomial model (~batch+condition, or ~condition) first. This gives size factors and trended/MAP dispersions \(\hat\alpha_g\) without absorbing true condition effects into residual variance. Construct the null mean from the same fit by setting only the condition coefficient to zero, while retaining intercept, covariate coefficients, and count-scale normalisation factors: \(\hat\mu^0_{gj}=NF_{gj}2^{X_j\hat\beta^0_g}\). DESeq2 coefficients are verified log2-scale and its mu assay/count normalisation factors are count-scale. The engine reconstructs the full mean before zeroing condition and asserts equality with assays(dds)[["mu"]] to floating-point tolerance; a version change cannot silently mis-scale simulation. This replaces the earlier reduced-null fit because omitting condition inflates dispersion for real DE genes and biases expected null discoveries downward.

For b = 1,...,B:

1. Simulate integer counts \(Y^{(b)}_{gj}\sim NB(\hat\mu^0_{gj},size=1/\hat\alpha_g)\) with R rnbinom using mu and size. Preserve metadata, sample order, covariates, and condition labels.
2. Rerun every filter and legal specification. Do not reuse observed filter masks, normalisation factors, or p-values.
3. Compute \(T^{(b)}_g,C^{(b)}_g\) for each supported (q,a), then \(R_b(\tau)=\sum_g1\{T^{(b)}_g\ge\tau,C^{(b)}_g\ge0.90\}\).
4. Compute observed \(R_{obs}(\tau)\) identically. Store only bootstrap discovery counts, not bootstrap gene-by-specification arrays.

The displayed curve is

\[
\widehat{eFDR}(\tau)=\min\left[1,
\frac{(1+\sum_{b=1}^{B}R_b(\tau))/(B+1)}
{\max\{R_{obs}(\tau),1\}}\right].
\]

For zero observed calls eFDR is NA. The +1/B+1 smoothing is conservative for Monte-Carlo estimation. The curve is evaluated at every attainable stability value k/m, not a coarse grid. A gene with stability T is called at every threshold tau <= T, so its eFDR is min{raw eFDR(tau): tau <= T}. Therefore raw values are made monotone with a cumulative minimum in **ascending** tau order (with unavailable zero-call thresholds excluded), not by sweeping from stringent high tau downward. Each gene receives eFDR_g = eFDR(T_g), preserving the stability ranking. Tau zero is excluded from automatic selection because it imposes no cross-specification support. The UI shows B, expected null discoveries, and an uncertainty interval, then chooses the smallest positive tau meeting the user target unless a stricter tau is selected.

MVP allows B=20 only with an exploratory warning and no target below 0.10. B=100 is usable default (targets 0.10/0.20); B=500 is report-grade before offering 0.05.

### Alternatives

1. Condition-label permutation preserves observed margins but is sensitivity-only, offered if each group has >=5 samples, exchangeability holds (stratify within batch), and no confounding. At 3 vs 3 there are 20 labelled assignments but only 10 unique unsigned partitions, too coarse. At 5 vs 5 there are 126 unique unsigned partitions, hence the minimum offering rule. Real signals also mean permutation is not a complete-null generator.
2. Storey-style pi0 adjustment is fast but not used for calls: its target after nonlinear correlated stability selection is unclear. It may become descriptive only.

The NB bootstrap is adopted because it creates a complete null while retaining estimated means, dispersion, library sizes, and measured covariates, and reruns the whole decision rule.

**要検討 / Open question — formal guarantee.** Fitted-parameter reuse and correlated genes/specifications preclude exact finite-sample FDR control. Conservative fallback: call it eFDR, show uncertainty, and require simulation validation.

**要検討 / Open question — batch.** Perfect batch-condition confounding cannot be repaired; stop. Strong unmeasured batch can misspecify the null, so flag results exploratory.

**要検討 / Open question — n, zeros, outliers.** Require >=3/group and warn below five. Zero inflation/outlier domination can violate the NB null; no zero-inflated model is promised in MVP.

## 5. Specification curve

The primary curve is per gene. Columns are specifications ordered by that gene's unshrunk log2FC (ties by specification id). The upper panel plots signed -log10(padj), positive for positive LFC and negative for negative LFC, rather than LFC confidence intervals: it gives a common test-strength scale across DESeq2, edgeR, and limma, whose standard errors are not directly comparable. Use Okabe-Ito vermilion #D55E00 for significant/up, blue #0072B2 for significant/down, and grey #999999 otherwise.

The lower matrix aligns columns and uses rows for filter, normalisation, method/test, shrinkage, and covariate. Categories use Okabe-Ito blue, orange #E69F00, bluish green #009E73, purple #CC79A7, and grey, with labels/legend so colour is never the only encoding. A global companion summary shows stability distributions (median, 10th/90th percentiles) by family/filter; it does not claim a biologically meaningful median LFC.

Verified installed packages include patchwork 1.3.2, cowplot 1.2.0, gridExtra 2.3, but none is in DESCRIPTION. Use only ggplot2 4.0.3: one long data frame, shared integer spec_rank, and facet_grid(panel ~ ., scales="free_y", space="free_y"), so panels align with no dependency addition.

Literal in-app caption:

> Each column is one reasonable way to analyse the same samples. The top dots show how strongly this gene changes: dots above zero favour the test group and dots below zero favour the reference group. The tiles below say what changed. A gene is more robust when most columns point in the same direction and pass the threshold; more columns alone do not make it reliable.

## 6. Data model and downstream hand-off

Add reactiveVal fields, initialised and reset together:

- deg_multiverse: NULL or immutable run list containing run_id, contrast, specifications, stats, tested, stability, efdr_curve, bootstrap_counts, and metadata (seed, package versions, call rules, B, timings, warnings).
- deg_multiverse_progress: list(stage, completed, total, message), UI-only and updated in the main process.

Stats is a named [gene, specification, metric] double array with mandatory log2FoldChange, pvalue, padj and optional stat/baseMean. Tested is a logical gene-by-specification matrix. Stability has gene, positive/negative stability, maximum stability, direction consistency/direction, selected eFDR, call status. Bootstrap_counts is indexed by bootstrap/call rule/tau; efdr_curve is its aggregate.

20,000 genes × 200 specs × three doubles × 8 bytes = 96,000,000 bytes (91.6 MiB), plus about 4 MiB for tested; stat/baseMean add ~61 MiB. Approximately 110–170 MiB is acceptable but shown before launch. Never retain bootstrap per-gene arrays.

On selection, add one data frame under multiverse:treated_vs_control in state$deg_results(). It has gene, baseMean, log2FoldChange, lfcSE where comparable, stat, pvalue, padj. LFC is median tested-specification LFC in selected direction; pvalue is NA because there is no valid pooled raw p-value; padj is that gene's monotone eFDR_g evaluated at its own stability; stat is signed stability. Sort eFDR then stability. This matches actual Enrichment usage: ORA reads gene, padj, abs(log2FoldChange); GSEA accepts stat. GSVA reads counts, not deg_results; Report reads padj/LFC. The enrichment fallback must not use all-NA pvalue.

Log multiverse_deg contrast/levels; stable specification-table hash; requested/actual/skipped paths; filters; covariate/rank check; call-rule grid; consistency/tau/eFDR selection; bootstrap model/B/seed/package versions; timings; call count; warnings.

**要検討 / Open question — consumer semantics.** Existing modules call padj a per-gene BH q-value, but stability eFDR is threshold-level. MVP uses the compatibility field only with result_type="multiverse_stability" metadata and an explanatory downstream banner; a result-class contract is future work.

## 7. File layout and integration

Phase 2 will add:

- R/mod_deg_multiverse.R: thin UI/server, validation, scheduler, outputs/state hand-off.
- R/multiverse_engine.R: Shiny-free functions validate_multiverse_inputs(counts, metadata, contrast, covariate, options), build_multiverse_specifications(options, availability), run_multiverse_specification(counts, metadata, specification, contrast), run_multiverse_observed(...), compute_multiverse_stability(stats, tested, call_rule), fit_multiverse_null(counts, metadata, covariate), simulate_multiverse_null(null_fit, seed), bootstrap_multiverse_efdr(...), make_multiverse_handoff(run, selection), apply_multiverse_run(state, run, selection).
- R/multiverse_plots.R: plot_specification_curve(run, gene, call_rule), plot_global_multiverse_summary(run, call_rule).
- tests/testthat.R and tests/testthat/test-multiverse-*.R.

Insert “5 Multiverse DEG” after “4 DEG” in R/app_ui.R, renumber later tabs 6–10, and register mod_deg_multiverse_server("deg_multiverse", state=state) immediately after mod_deg_server. Add deg_multiverse="pending" to both AppState initialize and reset, and to app_server.R steps. The known reset omission of upset/timeseries remains otherwise untouched.

app_server.R renders badge_step outputs but .nav_label has no matching uiOutput placeholders. Phase 2 must add placeholders for all labels, including badge_deg_multiverse, or deliberately defer badges; adding only the step makes it invisible. Success sets deg_multiverse and compatibility deg to done, so unchanged Enrichment and UpSet gates consume the new deg_results entry.

## 8. Cost, scheduling, caching

Fits are approximately (B+1)|S| plus null fitting/simulation. A measured local toll benchmark on simulated 20,000 × 12 counts: DESeq2 4.35 s, edgeR QL 2.54 s, limma-voom 0.60 s elapsed. This is hardware/data dependent. Full 80 paths/B=100 is 8,080 fits, about 5.6 CPU-hours at a 2.5-s mixed mean. At that same ~2.5-s mixed mean: MVP without a covariate (8 paths)/B=20 is 168 fits (~7 serial CPU-minutes), MVP with a covariate (16 paths)/B=20 is 336 fits (~14 serial CPU-minutes); B=100 scales those to ~35 and ~45 minutes respectively. The synthetic validation run reported in §9/the manuscript used the no-covariate, 8-path grid.

Use the existing ExtendedTask plus promises future_promise idiom per bounded batch, not one opaque multiverse promise. A batch returns plain partial arrays/counts; its main-session callback merges it, updates deg_multiverse_progress, and queues the next batch. Suggested units: four observed specs, then one bootstrap replicate × four specs. Workers never call reactive APIs.

Do not set a private future plan in this module. Honour the app plan; app startup may set one shared future multisession plan via explicit MultiverseDEG.future_workers (default max(1, availableCores()-1)). Limit in-flight batches to that count and disable package-internal parallelism, preventing conflicts with existing async modules.

Cache a run hash of counts/metadata identity, contrast, options, package versions, seed. Changing tau, a supported call rule, table sort, or selected gene does no DE fit: re-threshold observed stats and look up bootstrap curves. Changing filter/method/covariate/B/seed creates a new run.

## 9. Phased implementation

### MVP (Phase 2)

Implement the 16-path grid; DESeq2 Wald/apeglm optional; edgeR TMM QL/LRT; stability/direction consistency; NB bootstrap B=20 exploratory/B=100 default; supported call-rule curves; per-gene faceted curve; batching/progress; CSV; deg_results hand-off; logging; and tests. Add edgeR to Suggests. Do not depend on sva/IHW.

Defer limma-voom, edgeR RLE/UQ, ComBat-seq/IHW, random sampling/weights, permutations, disk cache/resume, B=500 default, and result-class refactoring.

### Full

Add deferred paths/features, uncertainty bands, permutation sensitivity, family-balanced weights, global curve, disk-backed resume, and simulation studies for unmeasured batch/zero inflation/method imbalance.

## 10. Test plan

Create testthat edition-3 scaffolding. Generator: 1,000 genes; 6/group; log-normal baseline means (median 100); library multipliers 0.8–1.2; NB dispersion 0.15; 60 true up and 60 true down genes at log2FC +1.5/-1.5; 880 exact-null genes; seed 20260918. Test grid: two filters × condition-only DESeq2 Wald unshrunk and edgeR QL/TMM (apeglm disabled).

1. Invalid/double-normalised specs and rank-deficient covariates are pruned with reasons; unavailable packages skip only their arms.
2. Toy stats give stability 1 for a gene significant/up in every spec and 0 for one in none; mixed signs fail 0.90 consistency; calls are monotone as tau increases.
3. Generator, bootstrap counts, eFDR curve, and future sampling are deterministic under fixed seed.
4. Known-truth calibration: B=50, target eFDR 0.10; realised FDP among calls using 880 known nulls must be <=0.20 (nominal +0.10 absolute). This tolerates one finite dataset, 50 bootstraps, correlated genes, and discrete calls; it is a regression check, not proof.
5. Pure-null version (zero true DE, three fixed seeds, B=5 in the fast fixture): zero calls for a majority of seeds, never more than two calls, selected-tau eFDR is at least 0.10, and mean stability <0.05. These bounds tolerate rare discrete null excursions while detecting systematic selection; the B=50 known-truth validation remains the slow test.
6. Hand-calculated tests cover R_b, +1/B+1, cap at one, zero-observed.
7. Instantiate AppState$new() outside a Shiny server, apply engine result, and verify reactive fields, parameter log, named deg_results contract.
8. Plot functions return ggplot objects under ggplot2 4.0.3, share spec_rank, and hand-off has required columns, finite stat, sorted eFDR/stability, no pooled p-value.

**要検討 / Open question — test runtime.** B=50 may be slow for normal CI. Conservative fallback: B=10 fixture per commit plus required B=50 statistical-validation CI/nightly test; never silently reduce the report-grade test.

## Verified runtime interfaces and references

Verified in toll: R 4.5.2, DESeq2 1.50.2, edgeR 4.8.2, limma 3.66.0, ggplot2 4.0.3.

- DESeq2 DESeq: object, test Wald/LRT, fitType parametric/local/mean/glmGamPoi, sfType ratio/poscounts/iterate, full design, reduced.
- DESeq2 results: object, contrast/name, lfcThreshold, pAdjustMethod BH, filterFun, test Wald/LRT.
- DESeq2 lfcShrink: dds, coef/contrast/res, type apeglm/ashr/normal, lfcThreshold, apeMethod nbinomCR.
- edgeR filterByExpr.DGEList: y, design/group/lib.size; calcNormFactors.DGEList: object, method TMM/TMMwsp/RLE/upperquartile/none; estimateDisp.DGEList: y, design, trend/tagwise/robust; glmQLFit.DGEList: y, design, dispersion, abundance.trend, robust; glmQLFTest: glmfit, coef/contrast, poisson.bound; glmLRT: glmfit, coef/contrast.
- limma voom: counts, design, lib.size, normalize.method, plot; eBayes: fit, robust.

Installed but intentionally not adopted: patchwork 1.3.2, cowplot 1.2.0, gridExtra 2.3. sva and IHW are unavailable and optional.

References: Steegen et al. (2016), doi:10.1177/1745691616658637; Simonsohn et al. (2020), doi:10.1038/s41562-020-0912-z; Patel et al. (2015), doi:10.1016/j.jclinepi.2015.05.029; Love, Huber & Anders (2014), doi:10.1186/s13059-014-0550-8; Ge, Son & Yao (2018), doi:10.1186/s12859-018-2486-6.
