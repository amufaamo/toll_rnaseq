# EasyRNA-Seq: Manuscript Figures and Tables

The following components are prepared for the "Applications Note" publication.

## Figure 1: System Architecture

```mermaid
graph TD
    classDef mainApp fill:#2d4263,stroke:#c84b31,stroke-width:2px,color:#fff;
    classDef docker fill:#e8f4f8,stroke:#0f3460,stroke-width:2px,color:#333;
    classDef data fill:#ecdbba,stroke:#e65100,stroke-width:2px,color:#333;

    A(Raw Reads<br>FASTQ)---|> B
    Ref(Reference Data<br>FASTA/GTF)---|> B

    A:::data; Ref:::data;

    subgraph EasyRNA_Seq_Platform [EasyRNA-Seq Shiny GUI Platform]
        B(Upstream Preprocessing<br>before_count):::mainApp
        D(Downstream Analysis<br>after_count):::mainApp
    end

    %% Docker Processes mapping
    B -.->|Spawns Containers| C1[fastp<br>QC & Trimming]:::docker
    B -.->|Spawns Containers| C2[STAR/featureCounts<br>Salmon Quant]:::docker
    B -.->|Spawns Containers| C3[Trinity & Corset<br>De novo Assembly]:::docker
    B -.->|Spawns Containers| C4[BUSCO & MultiQC<br>Evaluation]:::docker

    C1 & C2 & C3 & C4 -->|Generates| H[counts_matrix.csv]:::data
    H -->|Input| D
    
    D --> E1[Interactive Filtering]
    D --> E2[Dimension Reduction & Clustering<br>PCA, t-SNE, UMAP]
    D --> E3[Differential Expression<br>edgeR, DESeq2]
    D --> E4[Functional Enrichment<br>GSEA, GO, KEGG]
    D --> E5[Time Series Analysis<br>maSigPro]
    
    E1 & E2 & E3 & E4 & E5 ===> F[Publication-Ready<br>Figures & Tables]:::data
```

## Table 1: Comparison of Existing Transcriptomic GUI Tools

| Feature | **EasyRNA-Seq** | **iDEP** | **Galaxy** | **NetworkAnalyst** | **DEBrowser** | **ShinyNGS** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Execution Environment** | Local / Self-hosted Server | Cloud / Web-based | Web-based | Cloud / Web-based | Local / R Package | Local / Web-based |
| **Upstream Pipeline (FASTQ to Counts)** | ✅ (Docker-automated) | ❌ | ✅ | ❌ | ❌ | ❌ |
| **De novo Assembly (Non-model organisms)** | ✅ (Trinity + Corset) | ❌ | ⚠️ (Manual module setup) | ❌ | ❌ | ❌ |
| **Local Data Privacy (No Cloud Upload)** | ✅ | ❌ | ❌ | ❌ | ✅ (If local) | ✅ (If local) |
| **Tool Dependency Management** | Zero Setup (Containerized) | Zero Setup | Maintained by Cloud | Zero Setup | Manual (R Packages) | Manual (R Packages) |
| **Time Series Support** | ✅ (maSigPro) | ❌ | ❌ | ❌ | ❌ | ❌ |
| **Multi-group Complex Contrasts** | ✅ | ✅ | ✅ | ⚠️ | ✅ | ✅ |
| **Input File Size Limits** | Unlimited (Hard Drive capacity) | Restricted (< 100MB) | Quota-dependent | Moderate | RAM-dependent | RAM-dependent |
| **Journal-Ready Export (Nature/Cell/Science presets)** | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **Async computation (non-blocking UI)** | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **Colorblind-safe defaults** | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **Reproducible R script output** | ✅ | ⚠️ | ❌ | ❌ | ❌ | ❌ |

## Table 2: Supported Upstream Containers & Versions

| Tool | Version | Purpose |
| :--- | :---: | :--- |
| **fastp** | 0.23.4 | Read trimming and Quality Control |
| **STAR** | 2.7.11b | Reference-based Genome Alignment |
| **featureCounts** (subread) | 2.1.1 | Gene-level Read Summarization |
| **Salmon** | 1.10.0 | Transcript-level Quantification |
| **Trinity** | 2.15.1 | De novo Transcriptome Assembly |
| **Corset** | 1.09 | Unannotated/De novo Transcript Clustering |
| **BUSCO** | 5.5.0 | Assembly Completeness Assessment |
| **MultiQC** | 1.14 | Aggregate Report Generation |

## Table 2: Performance Benchmark (DESeq2 + GSEA, 20,000 genes)

| Samples | DESeq2 (s) | fgsea (s) | Total (s) | Peak RAM (MB) |
| :---: | :---: | :---: | :---: | :---: |
| 10 | 4.1 | 2.0 | 6.1 | 25 |
| 50 | 8.0 | 0.1 | 8.2 | 41 |
| 100 | 7.8 | 0.1 | 7.9 | 50 |
| 200 | 18.5 | 0.1 | 18.6 | 99 |
| 500 | 38.8 | 0.1 | 38.9 | 298 |

*Measured on AMD EPYC 96-core / 755 GB RAM. ExtendedTask ensures non-blocking UI throughout.*
