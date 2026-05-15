# EasyRNA-Seq Integrated App (Upstream + Downstream)

options(shiny.maxRequestSize = 10000 * 1024^2) # 10GB Upload Limit

# --- Libraries ---
library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(processx)
library(readr)
library(shinycssloaders)
library(shinyjs)
library(DT)
library(plotly)
library(dplyr)
library(purrr)
library(data.table)
library(tibble)
library(edgeR)
library(ggplot2)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(Rtsne)
library(umap)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(ggrepel)
library(maSigPro)
library(rhandsontable)
library(org.Mm.eg.db) # Default

# --- Source Dependencies ---
source("before_count/EasyRNASeq_Preprocessor/docker_wrapper.R")
source("before_count/EasyRNASeq_Preprocessor/R/utils.R")

required_after <- c(
  "after_count/R/module_data_upload_metadata_new.R",
  "after_count/R/module_filtering.R",
  "after_count/R/module_processing.R",
  "after_count/R/module_dimension_reduction.R",
  "after_count/R/module_deg_analysis.R",
  "after_count/R/module_gsea.R",
  "after_count/R/module_go_enrichment_integrated.R",
  "after_count/R/module_timeseries_analysis.R"
)
for (m in required_after) source(m)

# --- UI Definition ---
ui <- navbarPage(
    title = "EasyRNA-Seq Integrated Platform v2.0",
    id = "main_nav",
    windowTitle = "EasyRNA-Seq",
    header = tags$style(HTML("
      /* Custom CSS */
      .navbar-brand { font-weight: bold; font-size: 22px; }
    ")),
    shinyjs::useShinyjs(),

    navbarMenu("Upstream Analysis (Read Count)", icon = icon("cogs"),
        # Tab 1: Reference Preparation
        tabPanel(
            "1. Reference Prep",
            sidebarLayout(
                sidebarPanel(
                    h4("Create Index"),
                    
                    radioButtons("index_type", "Index Type:", choices = c("Salmon (Transcriptome)" = "salmon", "STAR (Genome)" = "star"), inline = TRUE),
                    numericInput("ref_threads", "Threads", value = 4, min = 1, max = 32),
                    
                    radioButtons("ref_mode", "Source:", 
                                 choices = c("Local File" = "local", "URL Download" = "url", "Upload (Drag&Drop)" = "upload"), 
                                 inline = TRUE),
                    
                    # --- Local Mode ---
                    conditionalPanel(
                        condition = "input.ref_mode == 'local'",
                        shinyDirButton("ref_fasta_dir", "Select Fasta Folder", "Select folder containing Reference Fasta"),
                        textOutput("ref_fasta_path_display", inline = TRUE),
                        br(), br(),
                        h5("Select Reference Fasta (Transcriptome)"),
                        selectInput("ref_fasta_file", "Fasta File", choices = NULL),
                        uiOutput("ref_warning_ui"),
                        h5("Select GTF/GFF Annotation (Optional for Index, Used later)"),
                        shinyDirButton("ref_gtf_dir", "Select GTF Folder", "Select folder containing GTF"),
                        selectInput("ref_gtf_file", "GTF File", choices = NULL)
                    ),
    
                    # --- URL Mode ---
                    conditionalPanel(
                        condition = "input.ref_mode == 'url'",
                        textInput("url_fasta", "Fasta URL", placeholder = "https://example.com/transcriptome.fa.gz"),
                        textInput("url_gtf", "GTF URL (Optional)", placeholder = "https://example.com/annotation.gtf.gz"),
                        h5("Download Destination"),
                        shinyDirButton("download_dest_dir", "Select Folder", "Select folder to save downloaded/uploaded files"),
                        textOutput("download_dest_display"),
                        br(), br()
                    ),
    
                    # --- Upload Mode ---
                    conditionalPanel(
                        condition = "input.ref_mode == 'upload'",
                        fileInput("upload_fasta", "Upload Fasta (Transcriptome)", accept = c(".fa", ".fasta", ".fna", ".gz")),
                        fileInput("upload_gtf", "Upload GTF (Optional)", accept = c(".gtf", ".gff", ".gff3", ".gz")),
                        p("Files will be saved to the destination folder below:"),
                        shinyDirButton("upload_dest_dir", "Select Destination Folder", "Folder to save uploaded files"),
                        textOutput("upload_dest_display"),
                        br(), br()
                    ),
                    
                    hr(),
                    textInput("index_output_name", "Index Folder Name", value = "index_out"),
                    actionButton("btn_make_index", "Process & Create Index", class = "btn-primary", icon = icon("dna")),
                    br(), br(),
                    uiOutput("index_status_ui"),
                    hr(),
                    h5("Docker Status:"),
                    textOutput("docker_status")
                ),
                mainPanel(
                    h4("Console Output"),
                    verbatimTextOutput("log_ref")
                )
            )
        ),
    
        # Tab 2: QC & Quantification
        tabPanel(
            "2. QC & Quant",
            sidebarLayout(
                sidebarPanel(
                    h4("Analyze Samples"),
    
                    # Inputs
                    radioButtons("quant_method", "Quantification Method:", choices = c("Salmon", "STAR + featureCounts"), inline = TRUE),
                    h5("1. FASTQ Directory"),
                    shinyDirButton("fastq_dir", "Select Data Folder", "Select folder with FASTQ files"),
                    textOutput("fastq_dir_display"),
                    h5("2. Index Directory"),
                    shinyDirButton("index_dir_select", "Select Index Folder", "Select Index Folder"),
                    textOutput("index_dir_display"),
                    h5("3. Annotation (GTF)"),
                    textOutput("gtf_path_display_t2"),
                    p("Use GTF selected in Tab 1, or enter path:"),
                    textInput("manual_gtf_path", "Absolute Path to GTF", placeholder = "/path/to/genes.gtf"),
                    fileInput("upload_gtf_t2", "Or Upload GTF (Drag & Drop)", accept = c(".gtf", ".gff", ".gff3", ".gz")),
                    h5("4. Settings"),
                    numericInput("threads", "Threads", value = 4, min = 1, max = 32),
                    checkboxInput("t2_run_multiqc", "Run MultiQC Report", value = TRUE),
                    h5("5. Output Directory"),
                    shinyDirButton("output_dir_btn", "Select Output Folder", "Select folder for results"),
                    textOutput("output_dir_display"),
                    hr(),
                    h5("Detected Samples"),
                    tableOutput("sample_table"),
                    actionButton("btn_run_quant", "Start Analysis", class = "btn-success", icon = icon("play")),
                    br(), br(),
                    uiOutput("quant_status_ui")
                ),
                mainPanel(
                    h4("Progress"),
                    progressBar(id = "quant_progress", value = 0, display_pct = TRUE),
                    textOutput("current_activity"),
                    h4("Log Output"),
                    verbatimTextOutput("log_quant"),
                    h4("Results"),
                    uiOutput("results_ui")
                )
            )
        ),
    
        # Tab 3: Assembly
        tabPanel(
            "3. De novo Assembly",
            sidebarLayout(
                sidebarPanel(
                    h4(icon("dna"), "De novo Pipeline (Trinity + Salmon + Corset)"),
                    div(
                        class = "alert alert-warning",
                        icon("exclamation-triangle"), " Warning: This pipeline requires significant resources (32GB+ RAM recommended)."
                    ),
                    h5("1. FASTQ Directory"),
                    shinyDirButton("trinity_fastq_dir", "Select FASTQ Folder", "Select folder"),
                    textOutput("trinity_dir_display"),
                    
                    h5("2. Output Directory"),
                    shinyDirButton("trinity_out_dir_btn", "Select Output Folder", "Select output folder"),
                    textOutput("trinity_out_dir_display"),
    
                    h5("3. Settings"),
                    numericInput("trinity_mem", "Max Memory (GB)", value = 30, min = 4),
                    numericInput("trinity_cpu", "CPU Cores", value = 8, min = 1),
                    
                    checkboxInput("run_busco", "Run BUSCO Assessment", value = TRUE),
                    conditionalPanel(
                        condition = "input.run_busco == true",
                        selectInput("busco_lineage", "BUSCO Lineage", 
                                    choices = c("eukaryota_odb10", "vertebrata_odb10", "mammalia_odb10", "fungi_odb10", "bacteria_odb10", "metazoa_odb10", "archaea_odb10"),
                                    selected = "eukaryota_odb10")
                    ),
                    checkboxInput("t3_run_multiqc", "Run MultiQC Report", value = TRUE),
                    hr(),
                    h5("Detected Samples"),
                    tableOutput("trinity_sample_table"),
                    actionButton("btn_run_trinity_pipeline", "Start De novo Pipeline", class = "btn-danger", icon = icon("cogs")),
                    br(), br(),
                    uiOutput("trinity_status_ui")
                ),
                mainPanel(
                    h4("Progress"),
                    progressBar(id = "trinity_progress", value = 0, display_pct = TRUE),
                    textOutput("trinity_current_activity"),
                    h4("Log Output"),
                    verbatimTextOutput("log_trinity")
                )
            )
        )
    ),

    navbarMenu("Downstream Analysis (DEG & Vis)", icon = icon("chart-bar"),
        tabPanel("1. Data Input", dataUploadMetadataUI("dataTab")),
        tabPanel("2. Filtering", filteringUI("filterTab")),
        tabPanel("3. Processing", processingUI("procTab")),
        tabPanel("4. Dimension Reduction", dimensionReductionUI("dimRedTab")),
        tabPanel("5. DEG Analysis", degAnalysisUI("degTab")),
        tabPanel("6. GSEA", gseaUI("gseaTab")),
        tabPanel("7. GO Enrichment", goEnrichmentIntegratedUI("go_module")),
        tabPanel("8. Time Series", timeseriesAnalysisUI("timeseriesTab"))
    )
)

# --- Server Logic ---
server <- function(input, output, session) {

    # ==========================================
    # UPSTREAM SERVER LOGIC (before_count)
    # ==========================================
    volumes <- getVolumes()() # Get volumes once
    vals <- reactiveValues(
        docker_ok = FALSE,
        ref_dir = NULL,
        gtf_dir = NULL,
        fastq_dir = NULL,
        index_dir = NULL,
        trinity_dir = NULL,
        trinity_out_dir = NULL,

        # Logic for jobs
        job_queue = list(),
        job_current_index = 0,
        job_active = FALSE,
        process = NULL,
        log = c(),

        # Context
        current_task = "",
        download_dir = NULL,
        upload_dir = NULL,
        quant_out_dir = NULL
    )

    # Check startup
    observe({
        vals$docker_ok <- check_docker_status()
    })

    output$docker_status <- renderText({
        if (vals$docker_ok) "🟢 Running" else "🔴 Not Detected"
    })

    append_log <- function(msg, target = "all") {
        ts <- format(Sys.time(), "[%H:%M:%S] ")
        entry <- paste0(ts, msg)
        vals$log <- c(vals$log, entry)
    }

    # --- Tab 1: Ref Prep ---
    shinyDirChoose(input, "ref_fasta_dir", roots = volumes)
    shinyDirChoose(input, "ref_gtf_dir", roots = volumes)
    shinyDirChoose(input, "download_dest_dir", roots = volumes)
    shinyDirChoose(input, "upload_dest_dir", roots = volumes)

    observeEvent(input$ref_fasta_dir, {
        path <- parseDirPath(volumes, input$ref_fasta_dir)
        if (length(path) > 0) {
            vals$ref_dir <- path
            files <- list.files(path, pattern = "\\.(fa|fasta|fna)(\\.gz)?$")
            updateSelectInput(session, "ref_fasta_file", choices = files)
        }
    })

    observeEvent(input$ref_gtf_dir, {
        path <- parseDirPath(volumes, input$ref_gtf_dir)
        if (length(path) > 0) {
            vals$gtf_dir <- path
            files <- list.files(path, pattern = "\\.(gtf|gff|gff3)(\\.gz)?$")
            updateSelectInput(session, "ref_gtf_file", choices = files)
        }
    })

    observeEvent(input$download_dest_dir, {
        path <- parseDirPath(volumes, input$download_dest_dir)
        if (length(path) > 0) vals$download_dir <- path
    })
    
    observeEvent(input$upload_dest_dir, {
        path <- parseDirPath(volumes, input$upload_dest_dir)
        if (length(path) > 0) vals$upload_dir <- path
    })

    output$ref_fasta_path_display <- renderText({ vals$ref_dir })
    output$download_dest_display <- renderText({ vals$download_dir })
    output$upload_dest_display <- renderText({ vals$upload_dir })

    output$ref_warning_ui <- renderUI({
        req(input$ref_fasta_file)
        f <- tolower(input$ref_fasta_file)
        is_genome <- grepl("genomic", f) || (grepl("\\.fna", f) && !grepl("cdna", f) && !grepl("transcript", f) && !grepl("rna", f))

        if (is_genome) {
            div(class = "alert alert-danger", HTML("<b>Warning:</b> The selected file name suggests it is a <b>Genomic</b> FASTA.<br>Salmon requires a <b>Transcriptome (cDNA)</b> FASTA file for accurate quantification."))
        } else {
            NULL
        }
    })
    
    output$index_status_ui <- renderUI({
        if (vals$job_active) {
            div(class = "alert alert-info", icon("spinner", class = "fa-spin"), paste(" Processing:", vals$current_task))
        } else {
            div(class = "alert alert-success", icon("check"), " Ready / Idle")
        }
    })

    observeEvent(input$btn_make_index, {
        if (!vals$docker_ok) {
            append_log("Docker not running!")
            return()
        }
        new_queue <- list()
        final_fasta_path <- ""
        final_gtf_path <- "" 
        work_dir <- ""
        
        if (input$ref_mode == "local") {
            req(input$ref_fasta_file, vals$ref_dir)
            final_fasta_path <- file.path(vals$ref_dir, input$ref_fasta_file)
            work_dir <- vals$ref_dir
            if (!is.null(vals$gtf_dir) && input$ref_gtf_file != "") {
                final_gtf_path <- file.path(vals$gtf_dir, input$ref_gtf_file)
            }
        } else if (input$ref_mode == "url") {
            req(input$url_fasta, vals$download_dir)
            work_dir <- vals$download_dir
            f_url <- input$url_fasta
            f_name <- basename(f_url)
            if (grepl("\\?", f_name)) f_name <- strsplit(f_name, "\\?")[[1]][1]
            if (f_name == "") f_name <- "reference.fasta"
            final_fasta_path <- file.path(work_dir, f_name)
            new_queue[[length(new_queue) + 1]] <- list(type = "download", url = f_url, dest = final_fasta_path, desc = paste("Downloading FASTA:", f_name))
            
            if (input$url_gtf != "") {
                g_url <- input$url_gtf
                g_name <- basename(g_url)
                if (grepl("\\?", g_name)) g_name <- strsplit(g_name, "\\?")[[1]][1]
                if (g_name == "") g_name <- "annotation.gtf"
                final_gtf_path <- file.path(work_dir, g_name)
                new_queue[[length(new_queue) + 1]] <- list(type = "download", url = g_url, dest = final_gtf_path, desc = paste("Downloading GTF:", g_name))
            }
        } else if (input$ref_mode == "upload") {
            req(input$upload_fasta, vals$upload_dir)
            work_dir <- vals$upload_dir
            f_name <- input$upload_fasta$name
            final_fasta_path <- file.path(work_dir, f_name)
            append_log(paste("Copying uploaded FASTA to:", final_fasta_path))
            file.copy(input$upload_fasta$datapath, final_fasta_path, overwrite = TRUE)
            
            if (!is.null(input$upload_gtf)) {
                g_name <- input$upload_gtf$name
                final_gtf_path <- file.path(work_dir, g_name)
                append_log(paste("Copying uploaded GTF to:", final_gtf_path))
                file.copy(input$upload_gtf$datapath, final_gtf_path, overwrite = TRUE)
            }
        }
        
        out_dir_path <- file.path(work_dir, input$index_output_name)
        if (input$index_type == "salmon") {
             cmd_list <- run_salmon_index(final_fasta_path, out_dir_path, threads = input$ref_threads)
             desc_text <- "Building Salmon Index"
        } else {
             cmd_list <- run_star_index(final_fasta_path, final_gtf_path, out_dir_path, threads = input$ref_threads)
             desc_text <- "Building STAR Index"
        }
        
        new_queue[[length(new_queue) + 1]] <- list(type = "index", cmd = cmd_list, desc = desc_text)
        
        vals$job_queue <- new_queue
        vals$job_current_index <- 1
        vals$job_active <- TRUE
        vals$log <- c() 
        start_next_job()
    })

    # --- Tab 2: QC & Quant ---
    shinyDirChoose(input, "fastq_dir", roots = volumes)
    shinyDirChoose(input, "index_dir_select", roots = volumes)
    shinyDirChoose(input, "output_dir_btn", roots = volumes)

    observeEvent(input$fastq_dir, {
        path <- parseDirPath(volumes, input$fastq_dir)
        if (length(path) > 0) vals$fastq_dir <- path
    })
    observeEvent(input$index_dir_select, {
        path <- parseDirPath(volumes, input$index_dir_select)
        if (length(path) > 0) vals$index_dir <- path
    })
    observeEvent(input$output_dir_btn, {
        path <- parseDirPath(volumes, input$output_dir_btn)
        if (length(path) > 0) vals$quant_out_dir <- path
    })

    output$fastq_dir_display <- renderText({ vals$fastq_dir })
    output$index_dir_display <- renderText({ vals$index_dir })
    output$output_dir_display <- renderText({ vals$quant_out_dir })

    output$gtf_path_display_t2 <- renderText({
        if (!is.null(vals$gtf_dir) && input$ref_gtf_file != "") {
            file.path(vals$gtf_dir, input$ref_gtf_file)
        } else {
            "None selected in Tab 1"
        }
    })

    samples_reactive <- reactive({
        req(vals$fastq_dir)
        files <- list.files(vals$fastq_dir, pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE)
        if (length(files) == 0) return(NULL)
        df <- data.frame(path = files, stringsAsFactors = FALSE)
        df$is_r2 <- grepl("(_R2|_2)\\.(fastq|fq)", basename(files))
        r1_files <- df$path[!df$is_r2]
        samples <- list()
        for (r1 in r1_files) {
            base <- basename(r1)
            s_name <- gsub("(_R1|_1)\\.(fastq|fq)(\\.gz)?$", "", base)
            r2_cand <- files[basename(files) != base & grepl(s_name, basename(files))]
            r2 <- if (length(r2_cand) > 0) r2_cand[1] else NULL
            samples[[length(samples) + 1]] <- list(name = s_name, r1 = r1, r2 = r2)
        }
        return(samples)
    })

    output$sample_table <- renderTable({
        s <- samples_reactive()
        if (is.null(s)) return(data.frame(Status = "No FASTQ found"))
        do.call(rbind, lapply(s, function(x) data.frame(Sample = x$name, Type = if (is.null(x$r2)) "Single" else "Paired")))
    })
    
    output$quant_status_ui <- renderUI({
        if (vals$job_active && !is.null(vals$current_task) && !grepl("(Trinity|BUSCO|Corset)", vals$current_task)) {
            div(class = "alert alert-info", icon("spinner", class = "fa-spin"), paste(" Running:", vals$current_task))
        } else { NULL }
    })

    observeEvent(input$btn_run_quant, {
        if (!vals$docker_ok) { showNotification("Docker is not running!", type = "error"); return() }
        if (is.null(vals$fastq_dir)) { showNotification("Please select FASTQ Directory.", type = "error"); return() }
        if (is.null(vals$index_dir)) { showNotification("Please select Salmon Index Directory.", type = "error"); return() }
        if (is.null(vals$quant_out_dir)) { showNotification("Please select Output Directory.", type = "error"); return() }
        
        samples <- samples_reactive()
        if (is.null(samples)) { showNotification("No FASTQ files found in selected directory.", type = "error"); return() }

        gtf_path <- ""
        if (!is.null(input$upload_gtf_t2)) {
            gtf_path <- input$upload_gtf_t2$datapath
            append_log(paste("Using uploaded GTF (Tab 2):", input$upload_gtf_t2$name))
        } else if (input$manual_gtf_path != "") {
            gtf_path <- input$manual_gtf_path
        } else if (!is.null(vals$gtf_dir) && input$ref_gtf_file != "") {
             gtf_path <- file.path(vals$gtf_dir, input$ref_gtf_file)
        }
        
        if (gtf_path == "" || !file.exists(gtf_path)) { showNotification("GTF file required for aggregation!", type = "error"); return() }

        output_base <- vals$quant_out_dir
        if (!dir.exists(output_base)) dir.create(output_base, recursive = TRUE)

        queue <- list()
        quant_files <- c() 

        for (s in samples) {
            s_out <- file.path(output_base, s$name)
            fastp_res <- run_fastp(s$r1, s$r2, s_out, s$name)
            queue[[length(queue) + 1]] <- list(type = "fastp", cmd = fastp_res$cmd, desc = paste("QC (fastp):", s$name))

            clean_r1 <- file.path(s_out, fastp_res$clean_r1)
            clean_r2 <- if (!is.null(fastp_res$clean_r2)) file.path(s_out, fastp_res$clean_r2) else NULL

            if (input$quant_method == "Salmon") {
                salmon_wrapper_cmd <- run_salmon_quant(index_dir = vals$index_dir, r1 = clean_r1, r2 = clean_r2, output_dir = file.path(s_out, "salmon_quant"), threads = input$threads)
                queue[[length(queue) + 1]] <- list(type = "salmon", cmd = salmon_wrapper_cmd, desc = paste("Quant (Salmon):", s$name))
                quant_files[s$name] <- file.path(s_out, "salmon_quant", "quant.sf")
            } else {
                star_wrapper_cmd <- run_star_align(index_dir = vals$index_dir, r1 = clean_r1, r2 = clean_r2, gtf = gtf_path, output_dir = file.path(s_out, "star_align"), sample_name = s$name, threads = input$threads)
                queue[[length(queue) + 1]] <- list(type = "star", cmd = star_wrapper_cmd, desc = paste("Align (STAR):", s$name))
                bam_file <- file.path(s_out, "star_align", paste0(s$name, "_Aligned.sortedByCoord.out.bam"))
                fc_wrapper_cmd <- run_featurecounts(bam = bam_file, gtf = gtf_path, output_dir = file.path(s_out, "featureCounts"), sample_name = s$name, is_paired = !is.null(clean_r2), threads = input$threads)
                queue[[length(queue) + 1]] <- list(type = "featurecounts", cmd = fc_wrapper_cmd, desc = paste("Quant (featureCounts):", s$name))
                quant_files[s$name] <- file.path(s_out, "featureCounts", paste0(s$name, "_featurecounts.txt"))
            }
        }

        agg_type <- if (input$quant_method == "Salmon") "aggregation_salmon" else "aggregation_featurecounts"
        queue[[length(queue) + 1]] <- list(type = agg_type, quant_files = quant_files, gtf_path = gtf_path, out_csv = file.path(output_base, "counts_matrix.csv"), out_len = file.path(output_base, "gene_lengths.csv"), desc = "Aggregating counts")

        if (input$t2_run_multiqc) {
            mqc_out <- vals$quant_out_dir
            mqc_cmd <- run_multiqc(target_dir = vals$quant_out_dir, output_dir = mqc_out)
            queue[[length(queue) + 1]] <- list(type = "multiqc", cmd = mqc_cmd, desc = "Generating MultiQC Report")
        }

        vals$job_queue <- queue
        vals$job_current_index <- 1
        vals$job_active <- TRUE
        vals$log <- c()
        start_next_job()
    })

    # --- Tab 3: Trinity Pipeline ---
    shinyDirChoose(input, "trinity_fastq_dir", roots = volumes)
    shinyDirChoose(input, "trinity_out_dir_btn", roots = volumes)
    
    observeEvent(input$trinity_fastq_dir, {
        path <- parseDirPath(volumes, input$trinity_fastq_dir)
        if (length(path) > 0) vals$trinity_dir <- path
    })
    observeEvent(input$trinity_out_dir_btn, {
        path <- parseDirPath(volumes, input$trinity_out_dir_btn)
        if (length(path) > 0) vals$trinity_out_dir <- path
    })
    
    output$trinity_dir_display <- renderText({ vals$trinity_dir })
    output$trinity_out_dir_display <- renderText({ vals$trinity_out_dir })

    trinity_samples_reactive <- reactive({
        req(vals$trinity_dir)
        files <- list.files(vals$trinity_dir, pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE)
        if (length(files) == 0) return(NULL)
        df <- data.frame(path = files, stringsAsFactors = FALSE)
        df$is_r2 <- grepl("(_R2|_2)\\.(fastq|fq)", basename(files))
        r1_files <- df$path[!df$is_r2]
        samples <- list()
        for (r1 in r1_files) {
            base <- basename(r1)
            s_name <- gsub("(_R1|_1)\\.(fastq|fq)(\\.gz)?$", "", base)
            r2_cand <- files[basename(files) != base & grepl(s_name, basename(files))]
            r2 <- if (length(r2_cand) > 0) r2_cand[1] else NULL
            samples[[length(samples) + 1]] <- list(name = s_name, r1 = r1, r2 = r2)
        }
        return(samples)
    })

    output$trinity_sample_table <- renderTable({
        s <- trinity_samples_reactive()
        if (is.null(s)) return(data.frame(Status = "No FASTQ found"))
        do.call(rbind, lapply(s, function(x) data.frame(Sample = x$name, Type = if (is.null(x$r2)) "Single" else "Paired")))
    })

    output$trinity_status_ui <- renderUI({
        if (vals$job_active && !is.null(vals$current_task) && grepl("(fastp|Trinity|BUSCO|Salmon|Corset)", vals$current_task)) {
            div(class = "alert alert-info", icon("spinner", class = "fa-spin"), paste(" Running:", vals$current_task))
        } else { NULL }
    })

    observeEvent(input$btn_run_trinity_pipeline, {
        if (!vals$docker_ok) { showNotification("Docker not running!", type="error"); return() }
        req(vals$trinity_dir, vals$trinity_out_dir)
        samples <- trinity_samples_reactive()
        req(samples)

        output_base <- vals$trinity_out_dir
        if (!dir.exists(output_base)) dir.create(output_base, recursive = TRUE)

        queue <- list()
        clean_r1s <- c()
        clean_r2s <- c()

        for (s in samples) {
            s_out <- file.path(output_base, "fastp_qc", s$name)
            fastp_res <- run_fastp(s$r1, s$r2, s_out, s$name)
            queue[[length(queue) + 1]] <- list(type = "fastp", cmd = fastp_res$cmd, desc = paste("QC (fastp):", s$name))
            clean_r1s <- c(clean_r1s, file.path(s_out, fastp_res$clean_r1))
            if (!is.null(fastp_res$clean_r2)) clean_r2s <- c(clean_r2s, file.path(s_out, fastp_res$clean_r2))
        }

        trinity_out_dir <- file.path(output_base, "trinity_out")
        r2_list_arg <- if (length(clean_r2s) > 0) clean_r2s else NULL
        trinity_cmd <- run_trinity(r1_list = clean_r1s, r2_list = r2_list_arg, output_dir = trinity_out_dir, max_memory = paste0(input$trinity_mem, "G"), cpu = input$trinity_cpu)
        queue[[length(queue) + 1]] <- list(type = "trinity", cmd = trinity_cmd, desc = "Trinity Assembly")
        trinity_fasta <- file.path(trinity_out_dir, "Trinity.fasta")

        if (input$run_busco) {
            busco_cmd <- run_busco(fasta = trinity_fasta, lineage = input$busco_lineage, output_dir = file.path(output_base, "busco"), threads = input$trinity_cpu)
            queue[[length(queue) + 1]] <- list(type = "busco", cmd = busco_cmd, desc = "BUSCO Quality Check")
        }

        salmon_idx_dir <- file.path(output_base, "salmon_index")
        salmon_idx_cmd <- run_salmon_index(input_fasta = trinity_fasta, output_idx_dir = salmon_idx_dir, threads = input$trinity_cpu)
        queue[[length(queue) + 1]] <- list(type = "salmon_index", cmd = salmon_idx_cmd, desc = "Salmon Indexing (Trinity FASTA)")

        eq_files <- c()
        for (i in seq_along(samples)) {
            s <- samples[[i]]
            s_out <- file.path(output_base, "salmon_quant", s$name)
            sq_cmd <- run_salmon_quant(index_dir = salmon_idx_dir, r1 = clean_r1s[i], r2 = if(length(clean_r2s) >= i) clean_r2s[i] else NULL, output_dir = s_out, threads = input$trinity_cpu)
            idx <- which(sq_cmd == "quant")
            if (length(idx) > 0) sq_cmd <- append(sq_cmd, "--dumpEq", after = idx)
            queue[[length(queue) + 1]] <- list(type = "salmon_quant", cmd = sq_cmd, desc = paste("Salmon Quant:", s$name))
            eq_files <- c(eq_files, file.path(s_out, "aux_info", "eq_classes.txt"))
        }

        corset_out <- file.path(output_base, "corset")
        corset_cmd <- run_corset(eq_classes_files = eq_files, output_dir = corset_out, threads = input$trinity_cpu)
        queue[[length(queue) + 1]] <- list(type = "corset", cmd = corset_cmd, desc = "Corset Clustering")
        
        queue[[length(queue) + 1]] <- list(
            type = "corset_postprocess",
            corset_counts = file.path(corset_out, "counts.txt"),
            corset_clusters = file.path(corset_out, "clusters.txt"),
            out_csv = file.path(output_base, "counts_matrix.csv"),
            samples = sapply(samples, function(x) x$name),
            desc = "Formatting Count Matrix"
        )
        
        if (input$t3_run_multiqc) {
            mqc_out <- vals$trinity_out_dir
            mqc_cmd <- run_multiqc(target_dir = vals$trinity_out_dir, output_dir = mqc_out)
            queue[[length(queue) + 1]] <- list(type = "multiqc", cmd = mqc_cmd, desc = "Generating MultiQC Report")
        }

        vals$job_queue <- queue
        vals$job_current_index <- 1
        vals$job_active <- TRUE
        vals$log <- c()
        start_next_job()
    })

    # --- Job Execution Engine ---
    start_next_job <- function() {
        if (vals$job_current_index > length(vals$job_queue)) {
            vals$job_active <- FALSE
            append_log("All jobs completed successfully.")
            return()
        }

        job <- vals$job_queue[[vals$job_current_index]]
        vals$current_task <- job$desc
        append_log(paste("Starting:", job$desc))

        if (job$type %in% c("aggregation_salmon", "aggregation_featurecounts", "corset_postprocess")) {
            tryCatch({
                if (job$type == "aggregation_salmon") {
                    aggregate_salmon_counts(job$quant_files, job$gtf_path, job$out_csv, job$out_len)
                } else if (job$type == "aggregation_featurecounts") {
                    aggregate_featurecounts(job$quant_files, job$out_csv, job$out_len)
                } else if (job$type == "corset_postprocess") {
                    if (file.exists(job$corset_counts)) {
                       counts_df <- read.delim(job$corset_counts, header=FALSE, sep="\t", stringsAsFactors=FALSE)
                       colnames(counts_df) <- c("Geneid", job$samples)
                       write.csv(counts_df, job$out_csv, row.names=FALSE)
                    } else {
                       append_log("Error: Corset counts.txt not found.")
                    }
                }
                append_log("Post-processing finished.")
                updateProgressBar(session, "trinity_progress", value = (vals$job_current_index / length(vals$job_queue)) * 100)
                updateProgressBar(session, "quant_progress", value = (vals$job_current_index / length(vals$job_queue)) * 100)
                vals$job_current_index <- vals$job_current_index + 1
                start_next_job()
            }, error = function(e) {
                append_log(paste("Error in post-processing:", e$message))
                vals$job_active <- FALSE
            })
        } else if (job$type == "download") {
            tryCatch({
                ret <- download.file(job$url, job$dest, method = "auto", mode = "wb")
                if (ret != 0) stop(paste("Download returned status", ret))
                if (!file.exists(job$dest)) stop("File not found after download")
                append_log(paste("Downloaded:", job$dest))
                vals$job_current_index <- vals$job_current_index + 1
                start_next_job()
            }, error = function(e) {
                append_log(paste("Download Error:", e$message))
                vals$job_active <- FALSE
            })
        } else {
            d_cmd <- "docker"
            d_args <- job$cmd
            vals$process <- processx::process$new(command = d_cmd, args = d_args, stdout = "|", stderr = "|", cleanup = FALSE)
        }
    }

    observe({
        req(vals$job_active)
        invalidateLater(500)
        if (!is.null(vals$process) && vals$process$is_alive()) {
            out <- tryCatch(vals$process$read_output_lines(), error=function(e) character(0))
            err <- tryCatch(vals$process$read_error_lines(), error=function(e) character(0))
            if (length(out) > 0) append_log(tail(out, 1))
            if (length(err) > 0) append_log(tail(err, 1))
        } else if (!is.null(vals$process) && !vals$process$is_alive()) {
            status <- vals$process$get_exit_status()
            if (status == 0) {
                append_log("Command finished successfully.")
                updateProgressBar(session, "trinity_progress", value = (vals$job_current_index / length(vals$job_queue)) * 100)
                updateProgressBar(session, "quant_progress", value = (vals$job_current_index / length(vals$job_queue)) * 100)
                vals$process <- NULL
                vals$job_current_index <- vals$job_current_index + 1
                start_next_job()
            } else {
                append_log(paste("Process failed with status", status))
                vals$job_active <- FALSE
                vals$process <- NULL
            }
        }
    })

    output$log_ref <- renderText({ paste(vals$log, collapse = "\n") })
    output$log_quant <- renderText({ paste(vals$log, collapse = "\n") })
    output$log_trinity <- renderText({ paste(vals$log, collapse = "\n") })
    output$current_activity <- renderText({ if (vals$job_active) paste("Running:", vals$current_task) else "Idle" })
    output$trinity_current_activity <- renderText({ if (vals$job_active) paste("Running:", vals$current_task) else "Idle" })

    # ==========================================
    # DOWNSTREAM SERVER LOGIC (after_count)
    # ==========================================
    rv <- reactiveValues(
        merged_data = NULL,
        sample_metadata = NULL,
        file_info = NULL,
        gene_lengths = NULL,
        filtered_keep = NULL,
        background_genes_original = NULL,
        deg_results = NULL,
        selected_species = NULL,
        current_gene_id_type = "ENTREZID"
    )

    dataUploadMetadataServer("dataTab", rv)
    filteringServer("filterTab", rv)

    observe({
        req(rv$merged_data, rv$filtered_keep)
        gene_ids_from_merged_data <- rv$merged_data$Geneid
        if (!is.null(gene_ids_from_merged_data) && length(gene_ids_from_merged_data) == length(rv$filtered_keep)) {
            rv$background_genes_original <- gene_ids_from_merged_data[rv$filtered_keep]
        } else {
            rv$background_genes_original <- NULL
        }
    })

    processingServer("procTab", rv)
    dimensionReductionServer("dimRedTab", rv)
    degAnalysisServer("degTab", rv)
    gseaServer("gseaTab", rv)
    goEnrichmentIntegratedServer("go_module", deg_results_reactive=reactive(rv$deg_results), background_genes_reactive=reactive(rv$background_genes_original), selected_species_code_reactive=reactive(rv$selected_species))
    timeseriesAnalysisServer("timeseriesTab", rv)

}

shinyApp(ui, server)