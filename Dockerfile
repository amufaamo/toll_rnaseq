# Stage 1: builder
FROM rocker/r-ver:4.5.2 AS builder

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libuv1-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages. `remotes` (not the much heavier `devtools`, whose
# usethis/pkgdown/roxygen2/testthat/profvis dependency chain is unnecessary
# just to install a local package and has repeatedly failed to resolve in CI)
# installs the local package plus its declared Imports/Suggests below.
RUN R -e "install.packages(c('golem', 'shiny', 'bslib', 'bsicons', 'R6', 'DT', 'ggplot2', 'future', 'promises', 'qs2', 'svglite', 'remotes', 'msigdbr', 'fgsea'), repos='https://cloud.r-project.org/')"

# Install Bioconductor packages (Bioc 3.22)
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/'); BiocManager::install(version='3.22', ask=FALSE); BiocManager::install(c('DESeq2', 'apeglm', 'clusterProfiler', 'enrichplot', 'org.Hs.eg.db', 'org.Mm.eg.db', 'sva'), ask=FALSE)"

# Build and install the app
WORKDIR /app
COPY . .
RUN R -e "remotes::install_local(dependencies = FALSE, build = FALSE, upgrade = 'never')"

# Stage 2: runtime
FROM rocker/r-ver:4.5.2

# Install runtime system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy installed packages from builder
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

# Copy app files
WORKDIR /app
COPY --from=builder /app /app

EXPOSE 3838

# Start the application
CMD ["R", "-e", "options(shiny.port = 3838, shiny.host = '0.0.0.0'); MultiverseDEG::run_app()"]
