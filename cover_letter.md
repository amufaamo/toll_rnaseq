[Date]

Editor-in-Chief
Briefings in Bioinformatics

Dear Editor,

We are pleased to submit our manuscript titled "MultiverseDEG: multiverse differential expression analysis with calibrated stability in a graphical interface" for consideration as a Software article in Briefings in Bioinformatics.

Differential expression analysis of bulk RNA-seq requires a series of defensible but arbitrary choices — low-count filtering, normalisation, model family and test, log fold-change shrinkage, and nuisance-covariate adjustment — and the reported gene list depends on which combination the analyst happened to select. Existing graphical interfaces address the reproducibility of a selected workflow; iDEP, for example, already records the chosen parameters and emits the corresponding R code. What they do not address is that a perfectly reproducible analysis can still be a fragile one: one path is documented, and the other equally legal paths remain invisible.

MultiverseDEG v4.0 addresses that gap. Its central component is a multiverse DEG module that fixes one biological contrast, enumerates the legal analytic specifications, runs all of them, and ranks genes by *stability* — the fraction of specifications in which a gene is called, subject to a required direction consistency of 0.90. Because a stability threshold is not itself an error rate, the module calibrates stability against a dataset-specific parametric negative-binomial bootstrap null constructed from the full fitted model with the condition coefficient set to zero, and reports an effective FDR (eFDR) for each stability threshold. On a synthetic benchmark (400 genes, 40 true DE genes, n = 6 per group, B = 10) the eFDR curve was monotone non-increasing in the threshold, and the automatically selected threshold τ = 0.125 returned 22 calls with a realised false discovery proportion of 0.045 against a 0.10 target, recovering 21 of the 40 true DE genes.

We believe this is of interest to the readership for two reasons. First, the multiverse and specification-curve framework, and the related notion of vibration of effects in biomedical data, has so far been applied to RNA-seq mainly as a critique; here it is delivered as an operational analysis producing a ranked gene list with an attached error-rate estimate, a per-gene specification curve, and a hand-off into the downstream enrichment modules under an explicit result type. Second, it is available to users who do not write code, which is where undisclosed analytic flexibility is hardest to audit.

We wish to be explicit about scope. The package analyses bulk RNA-seq starting from a count matrix; alignment and quantification from FASTQ are not part of v4.0, and no end-to-end claim is made. The calibration reported here is demonstrated on synthetic data, and validation on public datasets is stated in the manuscript as the next step rather than as a completed result. The full design enumerates 80 legal specifications, of which an 8-path subset (16 with a nominated covariate) is implemented in the current release; this is stated in both the Abstract and the Implementation section.

The software is open-source (MIT) and freely available. All computation runs locally or on a user-controlled server, so count matrices and metadata need not be uploaded to a third-party host.
- **GitHub**: https://github.com/amufaamo/toll_rnaseq
- **Docker**: ghcr.io/amufaamo/easy-rna-seq:dev

We confirm that this manuscript has not been published elsewhere and is not under consideration by another journal. All authors have approved the manuscript and agree with its submission to Briefings in Bioinformatics.

Thank you for your time and consideration.

Sincerely,

Masakazu Hasegawa
[Institution]
amufaamo@gmail.com
