# Briefings in Bioinformatics - Software Article Submission Checklist

## 1. Manuscript Structure & Mandatory Sections
- [ ] **Title Page**: Includes Title, Author list, Affiliations, and Corresponding author information.
- [ ] **Abstract**: Unstructured, typically 200–250 words.
- [ ] **Key Points**: 3 to 5 bullet points summarizing the main findings/impact (max ~100 characters per point).
- [ ] **Keywords**: 8-10 words (e.g., RNA-seq, GUI, Shiny, bioinformatics tools, reproducibility).
- [ ] **Introduction**: Context and significance of the software.
- [ ] **Methods / Implementation**: Detailed explanation of architecture (ExtendedTask, Nextflow integration, Journal-Ready Export Engine).
- [ ] **Results and Discussion**: Benchmarks or case studies demonstrating performance over tools like iDEP, Galaxy, and DEBrowser.
- [ ] **Conclusion**.
- [ ] **Data Availability Statement**: Statement clarifying how users can access the tool and data.
- [ ] **Code Availability**: GitHub and Docker links clearly provided.
- [ ] **Author Contributions**.
- [ ] **Conflict of Interest Statement**.
- [ ] **Funding Statement**.
- [ ] **References**.

## 2. Formatting & Word Limits
- [ ] **Word Count**: Verify length. Briefings in Bioinformatics generally targets 2,000 to 7,000 words depending on article category; keep concise and impactful.
- [ ] **Line Numbering**: Continuous line numbers included.
- [ ] **Spacing**: Double-spaced throughout.

## 3. Figures & Tables Requirements
- [ ] **Resolution**: Minimum 300 dpi for halftones/color, 600 dpi for line art. (Leverage `cairo_pdf` engine for vector-perfect exports).
- [ ] **Format**: EPS, TIFF, or high-resolution PDF.
- [ ] **File Size**: Typically keep individual figure files reasonable; check ScholarOne upload limits.
- [ ] **Captions**: Included sequentially at the end of the document or separately.

## 4. Software & Reproducibility Requirements
- [ ] **GitHub Repository**: Set to Public. Include a clear `README.md`.
- [ ] **Zenodo DOI**: Generated for the specific release (v4.0) to ensure long-term availability.
- [ ] **Docker Image**: Documented instructions on how to pull and run the environment.
- [ ] **Test Dataset**: Small demo dataset provided or easily downloaded for reviewers to test.

## 5. Submission System (ScholarOne Manuscripts)
- [ ] **ORCID iDs**: Collected and verified for all authors (especially the corresponding author).
- [ ] **Cover Letter**: Ready for upload.
- [ ] **Main Manuscript Document**: Uploaded in Word or PDF format.
- [ ] **Figures**: Uploaded individually in required high-resolution formats.
- [ ] **Supplementary Materials**: Uploaded (e.g., extended benchmark tables, additional plots).
