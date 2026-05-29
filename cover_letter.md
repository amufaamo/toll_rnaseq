[Date]

Editor-in-Chief
Briefings in Bioinformatics

Dear Editor,

We are pleased to submit our manuscript titled "EasyRNA-Seq v4.0" for consideration as a Software article in Briefings in Bioinformatics.

EasyRNA-Seq v4.0 represents a major update to our widely used RNA-seq analysis platform. In this version, we introduce three key novelties that significantly advance RNA-seq data analysis workflows:

1. **Journal-Ready Export Engine**: We have implemented an advanced plotting system that includes native presets for top-tier journals (Nature, Cell, Science) and fully supports `cairo_pdf`. This ensures users can generate publication-quality figures straight out of the box without requiring manual post-processing.
2. **ExtendedTask Asynchronous Processing**: By overhauling the underlying computational architecture, we incorporated a non-blocking UI for compute-intensive tasks. This profoundly improves user experience, system responsiveness, and scalability during large dataset analyses.
3. **Nextflow Direct Integration**: The tool now features a seamless end-to-end integration with Nextflow, fully automating the pipeline from raw FASTQ files to the final manuscript-ready figures.

Compared to existing tools such as iDEP, Galaxy, and DEBrowser, EasyRNA-Seq v4.0 uniquely bridges the gap between scalable pipeline execution and publication-oriented visualizations. None of our competitors currently offer this combination of asynchronous UI responsiveness, direct Nextflow integration, and immediate journal-ready plotting. 

The software is open-source and freely available. Source code and a containerized Docker image to ensure strict reproducibility can be found at:
- **GitHub**: https://github.com/amufaamo/toll_rnaseq
- **Docker**: ghcr.io/amufaamo/easy-rna-seq:dev

We confirm that this manuscript has not been published elsewhere and is not under consideration by another journal. All authors have approved the manuscript and agree with its submission to Briefings in Bioinformatics.

Thank you for your time and consideration.

Sincerely,

Masakazu Hasegawa
[Institution]
amufaamo@gmail.com
