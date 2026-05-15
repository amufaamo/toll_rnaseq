FROM rocker/shiny-verse:latest

# Install system dependencies for Bioconductor packages
RUN apt-get update && apt-get install -y \
    libglpk-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpng-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('shiny', 'plotly', 'DT', 'shinycssloaders', 'pheatmap', 'ggplot2', 'tibble', 'writexl', 'dplyr', 'tidyr', 'dendextend', 'ggdendro', 'shinyjs', 'RColorBrewer', 'igraph', 'Rtsne', 'umap', 'msigdbr', 'ggrepel', 'rhandsontable'), repos='https://cran.rstudio.com/')"

# Install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(c('edgeR', 'DESeq2', 'AnnotationDbi', 'org.Hs.eg.db', 'org.Mm.eg.db', 'limma', 'fgsea', 'clusterProfiler', 'enrichplot', 'maSigPro'), ask=FALSE)"

# Create app directory
WORKDIR /srv/shiny-server/
COPY . /srv/shiny-server/

# Expose port
EXPOSE 3838

# Run app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 3838)"]
