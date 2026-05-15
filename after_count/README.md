# Toll RNA-seq Analysis App

This Shiny application performs differential expression analysis (DEG) and visualization.

## How to Run (Recommended: Docker)

Using Docker ensures the application runs with the correct R version and package dependencies, regardless of your OS (Windows, Mac, Linux).

### Prerequisites
1.  **Install Docker Desktop**:
    *   **Mac**: Update/Install via [Docker Desktop for Mac](https://docs.docker.com/desktop/install/mac-install/).
    *   **Windows**: [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/).
    *   Once installed, **open the Docker Desktop application** and wait until the status bar icon stops animating (it must be running).

### Steps
1.  Open your terminal or command prompt.
2.  Check if docker is running:
    ```bash
    docker --version
    ```
    (If this says "command not found", check if Docker Desktop is installed and running).
3.  Navigate to this directory.
4.  Build and run the container:
    ```bash
    docker-compose up --build
    ```
    *(The first build may take 15-30 minutes to compile all R packages. Subsequent runs will be instant.)*

4.  Open your browser and go to:
    [http://localhost:3838](http://localhost:3838)

5.  To stop the app, press `Ctrl+C` in the terminal.

---

## How to Run (Local RStudio)

If you prefer to run it without Docker, you must ensure all dependencies are installed.

1.  Open `app.R` in RStudio.
2.  Install required packages:
    ```r
    install.packages(c("shiny", "plotly", "DT", "shinycssloaders", "pheatmap", "ggplot2", "tibble", "writexl", "dplyr", "tidyr"))
    
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(c("edgeR", "DESeq2", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "limma"))
    ```
3.  Click **Run App**.

---

## How to Develop (Modify Code)

Since the project folder is mounted into the container (via `docker-compose.yml`), you can edit files on your host machine and see changes immediately.

1.  **Edit Code**: Open `app.R` or `R/` files in your favorite editor (VS Code, RStudio, etc.).
2.  **Apply Changes**:
    *   **UI Changes**: Often just reloading the browser page is enough.
    *   **Server Logic**: You may need to restart the container (`Ctrl+C` then `docker-compose up`) to reload the R process.
3.  **Add New Packages**:
    *   If you need a new R package, add it to the `Dockerfile` (in the `RUN R -e ...` section).
    *   Then rebuild the container:
        ```bash
        docker-compose up --build
        ```

