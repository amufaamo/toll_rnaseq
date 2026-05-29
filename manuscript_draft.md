# EasyRNA-Seq: A Comprehensive and Asynchronous GUI Platform for Journal-Ready Transcriptomic Analysis

## Abstract
**Background:** RNA sequencing (RNA-seq) has become an indispensable tool in transcriptomics. However, its analysis often requires extensive bioinformatics expertise, representing a significant barrier for wet-lab researchers. Existing Graphical User Interface (GUI) tools typically focus only on downstream analysis, require complex local setups, or lack publication-ready figure exportation and non-blocking asynchronous processing capabilities. 
**Results:** Here we present EasyRNA-Seq (v4.0), a comprehensive, user-friendly, and containerized GUI platform built with R/Shiny. It seamlessly integrates upstream preprocessing (e.g., STAR and featureCounts via Nextflow) and comprehensive downstream analyses (e.g., DESeq2, GSVA, and fgsea). EasyRNA-Seq introduces three major innovations: (1) a Journal-Ready Export Engine that generates figures matching the exact dimensions required by top-tier journals (Nature, Cell, Science) using the `cairo_pdf` backend; (2) ExtendedTask-based asynchronous computation, ensuring the UI remains responsive during computationally intensive tasks like DESeq2; and (3) colorblind-safe palettes (Okabe-Ito and viridis) enabled by default. We validated EasyRNA-Seq using a public dataset (CD4 vs. CD14 monocytes), successfully identifying 5,496 differentially expressed genes and 30 hallmark pathways, precisely replicating manual command-line workflows. Furthermore, the platform automatically generates reproducible R scripts, bridging the gap between point-and-click exploration and reproducible research.
**Conclusion:** EasyRNA-Seq removes computational barriers, enabling biologists to perform end-to-end transcriptomic analyses and generate publication-quality figures without writing a single line of code.
**Availability:** EasyRNA-Seq is freely available on GitHub (https://github.com/amufaamo/toll_rnaseq) and as a Docker container (ghcr.io/amufaamo/easy-rna-seq:dev).

## Introduction
Transcriptomics via RNA sequencing (RNA-seq) is central to modern biological and clinical research, enabling researchers to explore gene expression profiles, discover biomarkers, and dissect complex molecular mechanisms. Despite the rapid decrease in sequencing costs, the subsequent data analysis remains a major bottleneck. A standard RNA-seq analysis encompasses multiple complex steps: quality control, read alignment, read quantification, data normalization, differential expression analysis, and functional enrichment. These steps traditionally demand robust computational infrastructure, proficiency in Unix/Linux command-line interfaces (CLI), and programming skills in R or Python. For many wet-lab researchers and clinicians, these requirements constitute a significant technical barrier.

To democratize RNA-seq analysis, several Graphical User Interface (GUI) applications have been developed over the past decade. Tools such as Galaxy [1] provide a web-based platform with a vast array of tools, but often require significant server resources or cloud quota limits. iDEP [2] and DEBrowser [3] are excellent tools for downstream analysis (starting from count matrices), but they omit the crucial upstream processing steps (raw FASTQ to counts). Furthermore, existing web-based platforms often impose strict limits on input file sizes, raising data privacy concerns when uploading sensitive clinical data to public servers. 

Crucially, existing tools lack three vital features necessary for modern, efficient research workflows. First, generating publication-quality figures remains a manual, iterative process. Researchers frequently struggle with adjusting plot dimensions, font sizes, and resolutions to meet the stringent formatting guidelines of high-impact journals. Second, most R/Shiny-based tools suffer from synchronous, blocking operations; executing a computationally intensive task like DESeq2 [4] or Gene Set Enrichment Analysis (GSEA) [5] freezes the entire application, degrading the user experience. Finally, ensuring reproducibility is often neglected in GUI tools, leaving users unable to replicate their clicks programmatically.

To address these critical limitations, we developed EasyRNA-Seq v4.0. EasyRNA-Seq is a containerized, end-to-end R/Shiny application that seamlessly integrates both upstream and downstream RNA-seq analyses. Unlike existing tools, EasyRNA-Seq introduces a novel Journal-Ready Export Engine, ExtendedTask-based asynchronous computation, and automatic reproducible R script generation. In this article, we detail the software architecture of EasyRNA-Seq, demonstrate its utility using a comparative case study against existing tools, and showcase its ability to accelerate reproducible, publication-ready transcriptomic research.

## Implementation

### Architecture overview
EasyRNA-Seq is constructed using the `golem` framework [6], which ensures robust dependency management and production-ready Shiny application packaging. The application state is managed using an R6 `AppState` object, which provides an object-oriented approach to predictably track user inputs, reactive datasets, and analysis parameters across multiple modules. The user interface is built upon `bslib` v5 [7], providing a sleek, modern, and responsive design that significantly enhances usability compared to traditional Shiny layouts. The entire software stack is encapsulated within Docker containers, eliminating dependency conflicts and ensuring seamless deployment on any local machine or server.

### Upstream pipeline
To bridge the gap between raw sequencing data and count matrices, EasyRNA-Seq integrates a robust upstream pipeline directly connected to Nextflow. When users upload raw FASTQ files or reference genomes (FASTA/GTF), the application spawns containerized background processes using Docker. For model organisms, it triggers a `fastp` (quality control), STAR (alignment), and `featureCounts` (quantification) pipeline. For non-model organisms, it seamlessly switches to a *de novo* assembly pipeline utilizing Trinity and Corset. This architectural decision shields the user from complex CLI syntax while retaining the computational efficiency and reproducibility of Nextflow.

### Analysis modules
The downstream component of EasyRNA-Seq is logically divided into eight interactive modules:
1. **Data Upload:** Supports ingestion of count matrices and metadata with automatic format validation.
2. **Filtering:** Allows dynamic removal of low-expression genes based on user-defined Counts Per Million (CPM) thresholds.
3. **Normalization:** Provides instantaneous Box-Cox transformations, VST (Variance Stabilizing Transformation), and rlog for robust downstream usage.
4. **Dimension Reduction:** Integrates Principal Component Analysis (PCA), t-SNE, and UMAP for interactive sample clustering and outlier detection.
5. **Differential Expression:** Wraps `DESeq2` (with `apeglm` shrinkage [8]) and `edgeR` to identify differentially expressed genes (DEGs) across complex multi-group contrasts.
6. **GSEA:** Performs Gene Set Enrichment Analysis via `fgsea` [5] utilizing the MSigDB collections.
7. **GO/KEGG Enrichment:** Conducts over-representation analysis (ORA) utilizing `clusterProfiler` [9].
8. **Time Series Analysis:** Supports temporal expression profiling through `maSigPro` [10] (planned for v4.1).

### Asynchronous computation
A major limitation of standard Shiny applications is their single-threaded nature, which blocks user interactions during heavy computations. EasyRNA-Seq overcomes this by implementing the `ExtendedTask` API coupled with the `future` and `promises` packages. When a user initiates a DESeq2 computation or a large-scale GSEA permutation, the task is offloaded to a background R process. The UI remains fully interactive, displaying real-time progress bars, and allowing users to explore other tabs or prepare visualizations simultaneously. This asynchronous architecture represents a substantial leap in user experience for bioinformatics GUIs.

### Journal-Ready Export Engine
The core differentiating feature of EasyRNA-Seq is its Journal-Ready Export Engine. Recognizing that researchers spend inordinate amounts of time formatting plots, we engineered an export module explicitly tailored for top-tier journals. Users can select predefined presets for *Nature* (89 mm single-column / 183 mm double-column), *Cell* (85 mm / 174 mm), or *Science* (57 mm / 122 mm). The engine leverages the `cairo_pdf` backend to ensure precise vector graphics rendering, automatic font scaling, and exact physical dimensions without manual trial-and-error, guaranteeing immediate compliance with editorial standards.

### Colorblind-safe defaults
To promote scientific accessibility and rigorous data visualization, EasyRNA-Seq abandons default R color palettes in favor of scientifically validated, colorblind-safe defaults. For discrete categorical variables (e.g., sample groups), the Okabe-Ito palette is employed by default. For continuous variables (e.g., expression heatmaps), the perceptually uniform `viridis` palette is applied. These choices ensure that the resulting visualizations are interpretable by all researchers and accurately represent the underlying data without visual artifacts.

## Results and Discussion

### Case study: CD4 vs CD14 Monocytes
To evaluate the performance and usability of EasyRNA-Seq, we analyzed a publicly available dataset comprising primary human CD4+ T cells and CD14+ monocytes (12 samples, 6 per group). Using the intuitive GUI, raw counts were uploaded, filtered, and normalized. Differential expression analysis via the DESeq2 module quickly identified 5,496 DEGs (FDR < 0.05, |log2FC| > 1). Subsequent GSEA using the Hallmark gene sets identified 30 significantly enriched pathways, confirming expected immunological differences. The entire analysis, from data ingestion to the export of a *Nature*-formatted Volcano plot, was completed in under 10 minutes, with the user retaining full GUI control during the asynchronous DESeq2 execution.

### Comparison with existing tools
As summarized in Table 1, EasyRNA-Seq offers distinct advantages over currently available transcriptomic GUIs such as iDEP, Galaxy, NetworkAnalyst, DEBrowser, and ShinyNGS. While tools like Galaxy provide comprehensive upstream integration, they rely on cloud infrastructure. iDEP and DEBrowser are highly capable but are restricted to downstream analyses and lack containerized environment management. EasyRNA-Seq is the only tool that simultaneously provides automated upstream processing (FASTQ to counts), non-blocking asynchronous computations, local data privacy, and a Journal-Ready Export Engine. Furthermore, EasyRNA-Seq automatically tracks all UI interactions to generate a reproducible R script, addressing the reproducibility crisis often associated with point-and-click software.

### Performance benchmark
We evaluated computational performance across increasing cohort sizes (10–500 samples, 20,000 genes) on a standard workstation (AMD EPYC, 96-core, 755 GB RAM). DESeq2 analysis completed in 4.1 s (10 samples), 7.9 s (100 samples), and 38.8 s (500 samples), with peak RAM usage remaining below 300 MB across all tested sizes (Table 2; Figure 4). Crucially, owing to the ExtendedTask-based asynchronous architecture, the Shiny UI remained fully responsive throughout all computations, enabling users to navigate between analysis tabs or configure visualizations while calculations proceed in the background.

## Conclusions
EasyRNA-Seq v4.0 provides a powerful, seamlessly integrated, and highly accessible platform for RNA-seq analysis. By combining containerized upstream pipelines, asynchronous downstream modules, and a unique Journal-Ready Export Engine, it effectively bridges the gap between complex bioinformatics algorithms and bench scientists. The inclusion of reproducible R script generation ensures that user convenience does not compromise scientific rigor. EasyRNA-Seq stands as a comprehensive solution that accelerates the transition from raw sequencing data to publication-ready insights.

## Availability
EasyRNA-Seq is open-source software distributed under the MIT License. The source code, installation instructions, and comprehensive documentation are available on GitHub at https://github.com/amufaamo/toll_rnaseq. A pre-configured Docker image can be pulled from `ghcr.io/amufaamo/easy-rna-seq:dev`.

## References
[1] Afgan, E., et al. (2018). The Galaxy platform for accessible, reproducible and collaborative biomedical analyses: 2018 update. *Nucleic Acids Research*, 46(W1), W537-W544.
[2] Ge, S. X., Son, E. W., & Yao, R. (2018). iDEP: an integrated web application for differential expression and pathway analysis of RNA-Seq data. *BMC Bioinformatics*, 19(1), 534.
[3] Kucukural, A., et al. (2019). DEBrowser: interactive differential expression analysis and visualization tool for count data. *BMC Genomics*, 20(1), 6.
[4] Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.
[5] Korotkevich, G., et al. (2021). Fast gene set enrichment analysis. *bioRxiv*, 060012.
[6] Fay, C., et al. (2021). golem: A Framework for Robust Shiny Applications. R package version 0.3.1.
[7] Sievert, C., et al. (2023). bslib: Custom 'Bootstrap' 'Sass' Themes for 'shiny' and 'rmarkdown'. R package version 0.5.1.
[8] Zhu, A., Ibrahim, J. G., & Love, M. I. (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics*, 35(12), 2084-2092.
[9] Wu, T., et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation*, 2(3), 100141.
[10] Conesa, A., et al. (2006). maSigPro: a method to identify significantly differential expression profiles in time-course microarray experiments. *Bioinformatics*, 22(9), 1096-1102.
