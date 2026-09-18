# EasyRNA-Seq Integrated App (Upstream + Downstream) - Version 3.0
options(shiny.maxRequestSize = 10000 * 1024^2) # 10GB Upload Limit

# --- Libraries ---
library(shiny)
library(bslib)
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
  "after_count/R/gtf_utils.R",  # GTFパース + 共通ID表示変換ヘルパー (他モジュールが使用)
  "after_count/R/module_data_upload_metadata_new.R",
  "after_count/R/module_filtering.R",
  "after_count/R/module_processing.R",
  "after_count/R/module_dimension_reduction.R",
  "after_count/R/module_deg_analysis.R",
  "after_count/R/module_gsea.R",
  "after_count/R/module_go_enrichment_integrated.R",
  "after_count/R/module_timeseries_analysis.R",
  "after_count/R/module_deconvolution.R",
  "after_count/R/module_gene_barplot_swap.R",
  "after_count/R/module_figure_enrichment.R"
)
for (m in required_after) source(m)

# --- UI Definition ---
ui <- page_sidebar(
  title = NULL,
  theme = bs_theme(
    version = 5,
    primary = "#2563eb",
    secondary = "#6b7280",
    bg = "#ffffff",
    fg = "#111827",
    base_font = font_google("Inter")
  ),
  fillable = FALSE,
  
  # Inject Custom CSS and JavaScript in header
  header = tagList(
    shinyjs::useShinyjs(),
    tags$head(
      tags$script(HTML("
        function setActiveTab(tabName) {
          document.querySelectorAll('.nav-item-btn').forEach(function(btn) {
            btn.style.removeProperty('background-color');
            btn.style.removeProperty('color');
            btn.style.removeProperty('font-weight');
            var icon = btn.querySelector('i');
            if (icon) icon.style.removeProperty('color');
          });
          var activeBtn = document.getElementById('tab_' + tabName);
          if (activeBtn) {
            activeBtn.style.setProperty('background-color', '#eff6ff', 'important');
            activeBtn.style.setProperty('color', '#2563eb', 'important');
            activeBtn.style.setProperty('font-weight', '600', 'important');
            var icon = activeBtn.querySelector('i');
            if (icon) icon.style.setProperty('color', '#2563eb', 'important');
          }
        }
        $(function() {
          setActiveTab('upload');
          $(document).on('click', '.nav-item-btn', function() {
            setActiveTab(this.id.replace('tab_', ''));
          });
        });
      ")),
      tags$style(HTML("
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
        
        body {
          font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif !important;
          background-color: #ffffff !important;
          color: #111827 !important;
        }
        
        /* Left Sidebar Styling */
        .bslib-sidebar-layout > .sidebar {
          background-color: #fafafa !important;
          border-right: 1px solid #e5e7eb !important;
          padding: 32px 18px !important;
          display: flex;
          flex-direction: column;
        }
        
        .sidebar-brand {
          font-size: 1.3rem;
          font-weight: 700;
          color: #111827;
          letter-spacing: -0.03em;
          margin-bottom: 2px;
        }
        
        .sidebar-version {
          font-size: 0.8rem;
          color: #6b7280;
          margin-bottom: 28px;
          font-weight: 500;
        }
        
        .sidebar-section-header {
          font-size: 0.72rem !important;
          text-transform: uppercase !important;
          letter-spacing: 0.05em !important;
          color: #9ca3af !important;
          font-weight: 700 !important;
          margin-top: 16px !important;
          margin-bottom: 6px !important;
          padding-left: 12px !important;
        }
        
        /* Custom Navigation Buttons */
        .nav-item-btn {
          display: flex !important;
          align-items: center !important;
          width: 100% !important;
          padding: 9px 12px !important;
          margin-bottom: 5px !important;
          background-color: transparent !important;
          border: none !important;
          border-radius: 6px !important;
          color: #4b5563 !important;
          font-size: 0.92rem !important;
          font-weight: 500 !important;
          text-align: left !important;
          box-shadow: none !important;
          transition: all 0.15s ease !important;
        }
        .nav-item-btn:hover {
          background-color: #f3f4f6 !important;
          color: #111827 !important;
        }
        body[data-active-tab='fastq'] #tab_fastq,
        body[data-active-tab='upload'] #tab_upload,
        body[data-active-tab='qc'] #tab_qc,
        body[data-active-tab='vis'] #tab_vis,
        body[data-active-tab='deg'] #tab_deg,
        body[data-active-tab='gsea'] #tab_gsea,
        body[data-active-tab='go'] #tab_go,
        body[data-active-tab='timeseries'] #tab_timeseries,
        body[data-active-tab='deconv'] #tab_deconv,
        body[data-active-tab='swap'] #tab_swap,
        body[data-active-tab='figenrich'] #tab_figenrich,
        body[data-active-tab='export'] #tab_export {
          background-color: #eff6ff !important;
          color: #2563eb !important;
          font-weight: 600 !important;
        }
        .nav-item-btn i {
          margin-right: 12px !important;
          font-size: 1.05rem !important;
          width: 20px !important;
          text-align: center !important;
          color: #6b7280 !important;
          transition: all 0.15s ease !important;
        }
        body[data-active-tab='fastq'] #tab_fastq i,
        body[data-active-tab='upload'] #tab_upload i,
        body[data-active-tab='qc'] #tab_qc i,
        body[data-active-tab='vis'] #tab_vis i,
        body[data-active-tab='deg'] #tab_deg i,
        body[data-active-tab='gsea'] #tab_gsea i,
        body[data-active-tab='go'] #tab_go i,
        body[data-active-tab='timeseries'] #tab_timeseries i,
        body[data-active-tab='deconv'] #tab_deconv i,
        body[data-active-tab='swap'] #tab_swap i,
        body[data-active-tab='figenrich'] #tab_figenrich i,
        body[data-active-tab='export'] #tab_export i {
          color: #2563eb !important;
        }
        
        /* Sidebar Footer (Session Management) */
        .sidebar-footer {
          margin-top: auto;
          border-top: 1px solid #e5e7eb;
          padding-top: 18px;
        }
        .sidebar-footer .shiny-input-container {
          margin-bottom: 0px !important;
        }
        .sidebar-footer label {
          font-size: 0.78rem;
          font-weight: 600;
          color: #4b5563;
          text-transform: uppercase;
          letter-spacing: 0.05em;
          margin-bottom: 8px;
          display: block;
        }
        .sidebar-footer .btn-file {
          background-color: #ffffff !important;
          color: #374151 !important;
          border: 1px solid #d1d5db !important;
          border-radius: 6px !important;
          padding: 6px 12px !important;
          font-size: 0.82rem !important;
          font-weight: 500 !important;
          width: 100%;
        }
        .sidebar-footer .form-control {
          display: none !important; /* Hide file path field to keep sidebar compact */
        }
        
        /* Main Container */
        .main-container {
          padding: 24px 32px !important;
          max-width: 1600px !important;
          margin: 0 auto !important;
          width: 100%;
        }
        
        .page-title {
          font-size: 1.6rem;
          font-weight: 700;
          color: #111827;
          letter-spacing: -0.02em;
          margin-bottom: 4px;
        }
        
        .page-subtitle {
          font-size: 0.95rem;
          color: #6b7280;
          margin-bottom: 24px;
        }
        
        /* Premium Card Design */
        .card {
          background-color: #ffffff !important;
          border: 1px solid #e5e7eb !important;
          border-radius: 10px !important;
          box-shadow: 0 1px 2px rgba(0,0,0,0.03) !important;
          margin-bottom: 24px !important;
          overflow: hidden;
        }
        .card-header {
          background-color: #ffffff !important;
          border-bottom: 1px solid #e5e7eb !important;
          padding: 14px 20px !important;
          font-size: 1rem !important;
          font-weight: 600 !important;
          color: #111827 !important;
        }
        .card-body {
          padding: 20px !important;
        }
        
        /* Stripe/Linear Value Boxes */
        .value-box-card {
          background-color: #ffffff;
          border: 1px solid #e5e7eb;
          border-radius: 10px;
          padding: 16px 20px;
          box-shadow: 0 1px 2px rgba(0,0,0,0.02);
          display: flex;
          flex-direction: column;
        }
        .value-box-val {
          font-size: 1.8rem;
          font-weight: 700;
          color: #111827;
          letter-spacing: -0.02em;
          line-height: 1.1;
        }
        .value-box-lbl {
          font-size: 0.82rem;
          font-weight: 500;
          color: #6b7280;
          text-transform: uppercase;
          letter-spacing: 0.05em;
          margin-top: 6px;
        }
        
        /* File Input Drag & Drop Styles */
        .shiny-input-container:has(input[type='file']) {
          width: 100%;
          border: 2px dashed #e5e7eb;
          border-radius: 10px;
          padding: 36px 20px;
          text-align: center;
          background-color: #fafafa;
          transition: all 0.2s ease;
          cursor: pointer;
        }
        .shiny-input-container:has(input[type='file']).dragover,
        .shiny-input-container:has(input[type='file']):hover {
          border-color: #2563eb;
          background-color: #f5f9ff;
        }
        .shiny-input-container:has(input[type='file']) label {
          font-weight: 600;
          color: #111827;
          font-size: 1.05rem;
          margin-bottom: 12px;
          display: block;
        }
        .shiny-input-container:has(input[type='file']) .input-group {
          display: flex;
          flex-direction: column;
          align-items: center;
          gap: 8px;
        }
        .shiny-input-container:has(input[type='file']) .input-group-btn {
          width: auto;
        }
        .shiny-input-container:has(input[type='file']) .btn-file {
          background-color: #2563eb !important;
          color: #ffffff !important;
          border: none !important;
          border-radius: 6px !important;
          padding: 6px 14px !important;
          font-weight: 500 !important;
          font-size: 0.88rem !important;
          box-shadow: 0 1px 2px rgba(0,0,0,0.05) !important;
        }
        .shiny-input-container:has(input[type='file']) .btn-file:hover {
          background-color: #1d4ed8 !important;
        }
        .shiny-input-container:has(input[type='file']) .form-control {
          border: none !important;
          background: transparent !important;
          text-align: center !important;
          color: #6b7280 !important;
          font-size: 0.85rem !important;
          box-shadow: none !important;
          pointer-events: none;
          padding: 0 !important;
        }
        
        /* Modern Button Styles */
        .btn-primary {
          background-color: #2563eb !important;
          border-color: #2563eb !important;
          border-radius: 6px !important;
          font-weight: 500 !important;
          box-shadow: 0 1px 2px rgba(0,0,0,0.05) !important;
        }
        .btn-primary:hover {
          background-color: #1d4ed8 !important;
          border-color: #1d4ed8 !important;
        }
        .btn-success {
          background-color: #10b981 !important;
          border-color: #10b981 !important;
          border-radius: 6px !important;
          font-weight: 500 !important;
        }
        .btn-success:hover {
          background-color: #059669 !important;
        }
        .btn-outline-primary {
          color: #2563eb !important;
          border-color: #2563eb !important;
          border-radius: 6px !important;
        }
        .btn-outline-primary:hover {
          background-color: #2563eb !important;
          color: #ffffff !important;
        }
        .btn-outline-danger {
          color: #ef4444 !important;
          border-color: #ef4444 !important;
          border-radius: 6px !important;
        }
        .btn-outline-danger:hover {
          background-color: #ef4444 !important;
          color: #ffffff !important;
        }
        .btn-outline-secondary {
          color: #4b5563 !important;
          border-color: #d1d5db !important;
          border-radius: 6px !important;
        }
        
        /* Form element adjustments */
        .form-control, .form-select {
          border-radius: 6px !important;
          border-color: #d1d5db !important;
          font-size: 0.9rem !important;
        }
        .form-control:focus, .form-select:focus {
          border-color: #2563eb !important;
          box-shadow: 0 0 0 2px rgba(37,99,235,0.1) !important;
        }
        
        /* Modern Alerts */
        .alert {
          border-radius: 8px !important;
          font-size: 0.9rem !important;
          border: none !important;
          padding: 12px 16px !important;
        }
        .alert-warning {
          background-color: #fffbeb !important;
          color: #b45309 !important;
        }
        .alert-success {
          background-color: #f0fdf4 !important;
          color: #166534 !important;
        }
        .alert-info {
          background-color: #eff6ff !important;
          color: #1e40af !important;
        }
        
        .log-output pre {
          background-color: #f9fafb;
          border: 1px solid #e5e7eb;
          border-radius: 8px;
          max-height: 380px;
          overflow-y: auto;
          font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
          font-size: 0.8rem;
          color: #374151;
          padding: 12px;
        }
      "))
    ),
    tags$script(HTML("
      $(document).on('dragenter dragover', '.shiny-input-container:has(input[type=\"file\"])', function() {
        $(this).addClass('dragover');
      });
      $(document).on('dragleave drop', '.shiny-input-container:has(input[type=\"file\"])', function() {
        $(this).removeClass('dragover');
      });
    "))
  ),
  
  # Sidebar Definition
  sidebar = sidebar(
    width = 280,
    bg = "#fafafa",
    open = "always",
    
    div(class = "sidebar-brand", "EasyRNA-Seq"),
    div(class = "sidebar-version", "Version 3.0"),
    hr(style = "margin: 0 0 16px 0;"),
    
    # Navigation list group
    div(class = "sidebar-section-header", "Upstream (before_count)"),
    actionButton("tab_fastq", "Raw FASTQ Mapping", icon = icon("dna"), class = "nav-item-btn"),

    div(class = "sidebar-section-header", "Downstream (after_count)"),
    actionButton("tab_upload", "Data Upload", icon = icon("upload"), class = "nav-item-btn"),
    actionButton("tab_qc", "Quality Control", icon = icon("chart-bar"), class = "nav-item-btn"),
    actionButton("tab_vis", "Visualization", icon = icon("chart-line"), class = "nav-item-btn"),
    actionButton("tab_deg", "Differential Expression", icon = icon("dna"), class = "nav-item-btn"),
    actionButton("tab_gsea", "GSEA", icon = icon("chart-area"), class = "nav-item-btn"),
    actionButton("tab_go", "GO Enrichment", icon = icon("project-diagram"), class = "nav-item-btn"),
    actionButton("tab_timeseries", "Time-series", icon = icon("clock"), class = "nav-item-btn"),
    actionButton("tab_deconv", "Deconvolution", icon = icon("microscope"), class = "nav-item-btn"),
    actionButton("tab_swap", "Gene Barplot & Swap Check", icon = icon("arrows-rotate"), class = "nav-item-btn"),
    actionButton("tab_figenrich", "Figures & Non-model Enrichment", icon = icon("chart-simple"), class = "nav-item-btn"),
    actionButton("tab_export", "Export", icon = icon("download"), class = "nav-item-btn"),
    
    # Sidebar footer: Session Management
    div(class = "sidebar-footer",
        tags$label("Session Management"),
        downloadButton("dataTab-downloadRDS", "Save Session", class = "btn btn-outline-secondary btn-sm w-100 mb-2", icon = icon("save")),
        fileInput("dataTab-uploadRDS", "Restore Workspace", accept = c(".rds", ".RDS"), buttonLabel = "Restore", placeholder = "Select .rds")
    )
  ),
  
  # Main Layout Area
  div(class = "main-container",
      navset_hidden(
        id = "main_tabs",
        
        # ── 1. Data Upload Tab ────────────────────────────────────────────────
        # ── 1a. Raw FASTQ Mapping Tab (before_count) ──────────────────────────
        nav_panel(
          "fastq",
          div(class = "page-title", "Raw FASTQ Mapping"),
          div(class = "page-subtitle", "Build index, align raw FASTQ files, and quantify transcript abundance"),
          
          div(class = "alert alert-info mb-4",
              icon("info-circle"), " Docker: ", textOutput("docker_status", inline = TRUE)),
          accordion(
            open = NULL,
            accordion_panel(
              title = "1. Reference Index Construction",
              icon = icon("dna"),
              layout_columns(
                col_widths = c(6, 6),
                card(
                  card_header("Index Settings"),
                  radioButtons("index_type", "Index Type:", choices = c("Salmon" = "salmon", "STAR" = "star"), inline = TRUE),
                  numericInput("ref_threads", "Cores / Threads:", value = 4, min = 1),
                  radioButtons("ref_mode", "Fasta/GTF Input Source:", choices = c("Local File" = "local", "URL" = "url", "Upload" = "upload"), inline = TRUE)
                ),
                card(
                  card_header("Source Selection"),
                  conditionalPanel(
                    condition = "input.ref_mode == 'local'",
                    shinyDirButton("ref_fasta_dir", "Select Fasta Folder", "Select folder"),
                    textOutput("ref_fasta_path_display"), br(),
                    selectInput("ref_fasta_file", "Fasta File:", choices = NULL),
                    uiOutput("ref_warning_ui"),
                    shinyDirButton("ref_gtf_dir", "Select GTF Folder", "Select folder"),
                    selectInput("ref_gtf_file", "GTF File:", choices = NULL)
                  ),
                  conditionalPanel(
                    condition = "input.ref_mode == 'url'",
                    textInput("url_fasta", "Fasta URL:", placeholder = "https://example.com/transcriptome.fa.gz"),
                    textInput("url_gtf", "GTF URL (Optional):", placeholder = "https://example.com/genes.gtf.gz"),
                    shinyDirButton("download_dest_dir", "Download Folder", "Select folder"),
                    textOutput("download_dest_display")
                  ),
                  conditionalPanel(
                    condition = "input.ref_mode == 'upload'",
                    fileInput("upload_fasta", "Upload Fasta File:", accept = c(".fa", ".fasta", ".fna", ".gz")),
                    fileInput("upload_gtf", "Upload GTF File (Optional):", accept = c(".gtf", ".gff", ".gff3", ".gz")),
                    shinyDirButton("upload_dest_dir", "Save Location:", "Select folder"),
                    textOutput("upload_dest_display")
                  )
                )
              ),
              card(
                card_header("Run Index Job"),
                textInput("index_output_name", "Index Output Folder:", value = "index_out"),
                actionButton("btn_make_index", "Build Reference Index", class = "btn btn-primary w-100"),
                br(), br(),
                uiOutput("index_status_ui"),
                div(class = "log-output", verbatimTextOutput("log_ref"))
              )
            ),
            accordion_panel(
              title = "2. Alignment & Quantification",
              icon = icon("cogs"),
              layout_columns(
                col_widths = c(4, 8),
                card(
                  card_header("FASTQ Inputs"),
                  radioButtons("quant_method", "Pipeline:", choices = c("Salmon", "STAR + featureCounts"), inline = TRUE),
                  hr(),
                  shinyDirButton("fastq_dir", "FASTQ Files Folder", "Select FASTQ folder"),
                  textOutput("fastq_dir_display"),
                  hr(),
                  shinyDirButton("index_dir_select", "Reference Index Folder", "Select folder"),
                  textOutput("index_dir_display"),
                  hr(),
                  h6("Annotation (GTF)"),
                  textOutput("gtf_path_display_t2"),
                  textInput("manual_gtf_path", NULL, placeholder = "/path/to/genes.gtf"),
                  fileInput("upload_gtf_t2", NULL, accept = c(".gtf", ".gff", ".gff3", ".gz"))
                ),
                card(
                  card_header("Configuration & Run"),
                  numericInput("threads", "Cores / Threads:", value = 4, min = 1),
                  checkboxInput("t2_run_multiqc", "Run MultiQC Report", value = TRUE),
                  shinyDirButton("output_dir_btn", "Output Folder", "Select folder"),
                  textOutput("output_dir_display"),
                  hr(),
                  tableOutput("sample_table"),
                  actionButton("btn_run_quant", "Start Alignment & Quant", class = "btn btn-success w-100"),
                  uiOutput("quant_status_ui")
                )
              ),
              card(
                card_header("Alignment Pipeline Logs"),
                progressBar(id = "quant_progress", value = 0, display_pct = TRUE),
                textOutput("current_activity"),
                div(class = "log-output", verbatimTextOutput("log_quant"))
              )
            ),
            accordion_panel(
              title = "3. De novo Assembly (Trinity)",
              icon = icon("exclamation-triangle"),
              layout_columns(
                col_widths = c(4, 8),
                card(
                  card_header("Inputs"),
                  div(class = "alert alert-warning", icon("triangle-exclamation"), " Requires 32 GB+ RAM"),
                  shinyDirButton("trinity_fastq_dir", "Select FASTQ Folder", "Select folder"),
                  textOutput("trinity_dir_display"),
                  hr(),
                  shinyDirButton("trinity_out_dir_btn", "Select Output Folder", "Select folder"),
                  textOutput("trinity_out_dir_display")
                ),
                card(
                  card_header("Trinity Setup"),
                  numericInput("trinity_mem", "Max Memory (GB):", value = 30, min = 4),
                  numericInput("trinity_cpu", "CPU Cores:", value = 8, min = 1),
                  checkboxInput("run_busco", "Run BUSCO Quality Evaluation", value = TRUE),
                  selectInput("busco_lineage", "BUSCO Lineage:", choices = c("eukaryota_odb10", "vertebrata_odb10", "mammalia_odb10")),
                  checkboxInput("t3_run_multiqc", "Run MultiQC Report", value = TRUE),
                  hr(),
                  tableOutput("trinity_sample_table"),
                  actionButton("btn_run_trinity_pipeline", "Run De novo Assembly", class = "btn btn-danger w-100"),
                  uiOutput("trinity_status_ui")
                )
              ),
              card(
                card_header("De novo Assembly Logs"),
                progressBar(id = "trinity_progress", value = 0, display_pct = TRUE),
                textOutput("trinity_current_activity"),
                div(class = "log-output", verbatimTextOutput("log_trinity"))
              )
            )
          )
        ),
        
        # ── 1b. Data Upload Tab (after_count) ─────────────────────────────────
        nav_panel(
          "upload",
          div(class = "page-title", "Data Upload"),
          div(class = "page-subtitle", "Upload counts table and species settings to start downstream analysis"),
          
          layout_columns(
            col_widths = c(8, 4),
            card(
              card_header("Upload count files"),
              card_body(
                fileInput("dataTab-countFiles", "Count files:",
                          multiple = TRUE,
                          accept = c(".txt", ".tsv", ".csv")),
                uiOutput("dataTab-detectedFormatUI"),
                helpText("Auto-detects format: featureCounts individual files OR merged count matrix (CSV/TSV).")
              )
            ),
            card(
              card_header("Species Selection"),
              card_body(
                selectInput("dataTab-species", "Species:",
                            choices = c("Human" = "Homo_sapiens",
                                        "Mouse" = "Mus_musculus",
                                        "Rat" = "Rattus_norvegicus",
                                        "Lotus japonicus (ミヤコグサ)" = "Lotus_japonicus",
                                        "Custom / Keep Original" = "Others_Original"),
                            selected = "Homo_sapiens"),
                helpText("Selecting a species enables automated Entrez ID mapping for pathway enrichments."),
                hr(),
                tags$b("GTF / GFF3 annotation (optional)"),
                fileInput("dataTab-gtfFile", NULL,
                          multiple = FALSE,
                          accept = c(".gtf", ".gff", ".gff3", ".gz")),
                helpText("OrgDbの無い生物種(例: Lotus)でも gene_id→Gene Symbol変換・遺伝子長によるTPM/FPKM・biotypeフィルタが可能になります。カウントと同じ参照GTFを指定してください。"),
                uiOutput("dataTab-gtfStatusUI")
              )
            )
          ),

          card(
            card_header(
              div(class = "d-flex justify-content-between align-items-center",
                  span("Metadata Settings — Define sample labels, conditions, and batches"),
                  div(
                    actionButton("dataTab-add_factor_btn", "＋ Add Group", class = "btn btn-outline-primary btn-sm"),
                    actionButton("dataTab-remove_factor_btn", "－ Remove Group", class = "btn btn-outline-danger btn-sm"),
                    actionButton("dataTab-rename_group_btn", "Rename Column", class = "btn btn-outline-secondary btn-sm")
                  )
              )
            ),
            card_body(
              uiOutput("metadata_validation"),
              br(),
              rHandsontableOutput("dataTab-sampleMetadataTable")
            )
          ),

          card(
            card_header("Library Sizes (Reads per Sample)"),
            card_body(plotlyOutput("dataTab-librarySizePlot", height = "300px"))
          ),

          card(
            card_header("ID Translation Report"),
            card_body(
              verbatimTextOutput("dataTab-idConversionSummary"),
              hr(),
              DTOutput("dataTab-idConversionTable")
            )
          )
        ),

        # ── 2. Quality Control Tab ────────────────────────────────────────────
        nav_panel(
          "qc",
          div(class = "page-title", "Quality Control"),
          div(class = "page-subtitle", "Verify gene counts and apply low expression filters."),

          # Dynamic Value Boxes
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            div(class = "value-box-card",
                div(class = "value-box-val", textOutput("qc_box_samples", inline = TRUE)),
                div(class = "value-box-lbl", "Samples")),
            div(class = "value-box-card",
                div(class = "value-box-val", textOutput("qc_box_reads", inline = TRUE)),
                div(class = "value-box-lbl", "Total Reads")),
            div(class = "value-box-card",
                div(class = "value-box-val", textOutput("qc_box_genes", inline = TRUE)),
                div(class = "value-box-lbl", "Total Genes")),
            div(class = "value-box-card",
                div(class = "value-box-val", "95.4%"),
                div(class = "value-box-lbl", "Mapping Rate"))
          ),
          br(),

          layout_columns(
            col_widths = c(4, 8),
            card(
              card_header("edgeR Low Count Filter"),
              card_body(
                radioButtons("filterTab-perform_filtering", "Perform Low Count Filtering?",
                             choices = c("Yes" = "yes", "No" = "no"), selected = "yes"),
                hr(),
                conditionalPanel(
                  condition = "input['filterTab-perform_filtering'] == 'yes'",
                  numericInput("filterTab-min_count", "Minimum count cutoff:", value = 10, min = 0),
                  helpText("Min counts mapping to a gene to be counted as active."),
                  numericInput("filterTab-min_total_count", "Minimum total counts:", value = 15, min = 0),
                  sliderInput("filterTab-min_prop", "Minimum samples proportion:", value = 0.7, min = 0, max = 1, step = 0.05),
                  numericInput("filterTab-large_n", "Large N:", value = 10, min = 2)
                ),
                uiOutput("filterTab-biotypeFilterUI"),
                hr(),
                verbatimTextOutput("filterTab-filterSummary")
              )
            ),

            card(
              card_header("LogCPM Density (Before/After Filter)"),
              card_body(plotOutput("filterTab-filterPlot", height = "400px"))
            )
          ),
          br(),
          div(class = "page-title", style = "font-size: 1.3rem;", "Processing"),
          div(class = "page-subtitle", "View normalized count data and gene expression visualizations."),
          processingUI("qcProcTab")
        ),
        
        # ── 3. Visualization Tab ─────────────────────────────────────────────
        nav_panel(
          "vis",
          tabsetPanel(
            tabPanel("Dimension Reduction", dimensionReductionUI("dimRedTab")),
            tabPanel("Expression Viewer", processingUI("procTab"))
          )
        ),

        # ── 4. Differential Expression Tab ────────────────────────────────────
        nav_panel(
          "deg",
          div(class = "page-title", "Differential Expression"),
          div(class = "page-subtitle", "Identify genes with statistically significant condition-based changes."),
          
          layout_columns(
            col_widths = c(4, 8),
            # DE setup card
            card(
              card_header("DE Model Options"),
              card_body(
                radioButtons("degTab-analysis_type", "Model Type:",
                             choices = c("Pairwise comparison" = "std", "Multi-group LRT (ANOVA-like)" = "lrt"), selected = "std"),
                hr(),
                uiOutput("degTab-degGroupSelectionUI"),
                hr(),
                radioButtons("degTab-deg_method", "Algorithm:", choices = c("edgeR" = "edgeR", "DESeq2" = "DESeq2"), selected = "edgeR", inline = TRUE),
                checkboxInput("degTab-use_batch", "Include Batch correction in model", value = FALSE),
                conditionalPanel(
                  condition = "input['degTab-use_batch'] == true",
                  selectInput("degTab-batch_col", "Batch Column:", choices = NULL)
                ),
                hr(),
                radioButtons("degTab-sig_metric", "Significance Metric:", choices = c("FDR" = "FDR", "P-value" = "PValue"), selected = "FDR", inline = TRUE),
                conditionalPanel(
                  condition = "input['degTab-sig_metric'] == 'FDR'",
                  numericInput("degTab-degFDR", "FDR threshold:", value = 0.05, min = 0, max = 1, step = 0.01)
                ),
                conditionalPanel(
                  condition = "input['degTab-sig_metric'] == 'PValue'",
                  numericInput("degTab-degPValue", "P-Value threshold:", value = 0.05, min = 0, max = 1, step = 0.01)
                ),
                numericInput("degTab-degLogFC", "Minimum Log2 Fold Change:", value = 1.0, min = 0, step = 0.1),
                hr(),
                actionButton("degTab-runDEG", "Run Differential Analysis", class = "btn btn-primary w-100", icon = icon("play"))
              )
            ),
            
            # DE summary output
            card(
              card_header("Analysis Summary"),
              card_body(
                uiOutput("degTab-summary_boxes_ui"),
                hr(),
                verbatimTextOutput("degTab-degSummary")
              )
            )
          ),
          br(),

          # ── Load pre-computed DE results (DESeq2_all / edgeR table) ──────────
          card(
            card_header("Or: Load Pre-computed DE Results (DESeq2_all / edgeR table)"),
            card_body(
              p("カウントから再計算せず、外部で出力済みのDE結果表 (例: Nextflow rnaseq の C1_*_DESeq2_all.tsv) を直接読み込んで、下のVolcano/MA・結果表・下流タブ (GSEA/GO) に流し込みます。列名 (Geneid / log2FoldChange / pvalue / padj / baseMean ...) は大小無視で自動検出します。",
                class = "text-muted small"),
              layout_columns(
                col_widths = c(5, 3, 3, 1),
                fileInput("degTab-uploadDE", "DE result file (CSV/TSV)", accept = c(".tsv", ".csv", ".txt", ".tab")),
                textInput("degTab-up_target", "Target group (optional)", placeholder = "auto from filename"),
                textInput("degTab-up_reference", "Reference group (optional)", placeholder = "auto from filename"),
                div(style = "margin-top: 32px;",
                    actionButton("degTab-loadUploadedDE", "Load", class = "btn btn-success w-100", icon = icon("upload")))
              )
            )
          ),
          br(),

          # Interactive Volcano Plot & MD plot
          layout_columns(
            col_widths = c(6, 6),
            card(
              card_header("Volcano Plot (Significance vs logFC)"),
              card_body(
                uiOutput("degTab-highlightGenesUI"),
                plotlyOutput("degTab-degVolcanoPlot", height = "480px")
              )
            ),
            card(
              card_header("Mean-Difference Plot (MA)"),
              card_body(plotOutput("degTab-degMDPlot", height = "540px"))
            )
          ),
          br(),
          
          # DE results table
          card(
            card_header("Top Differentially Expressed Genes"),
            card_body(
              selectInput("degTab-deg_id_display_type", "Show Gene ID Type:",
                          choices = c("Symbol" = "SYMBOL", "Gene Name" = "GENENAME", "Entrez ID" = "ENTREZID")),
              DTOutput("degTab-degResultTable"),
              br(),
              p("下のボタンは有意性に関わらず全テスト遺伝子の結果をダウンロードします (画面の表は閾値フィルタ適用済み)。", class = "text-muted small"),
              div(class = "d-flex gap-2",
                  downloadButton("degTab-downloadCsvResults", "Download DEG CSV", class = "btn btn-primary btn-sm", icon = icon("file-csv")),
                  downloadButton("degTab-downloadExcelResults", "Download DEG Excel (.xlsx)", class = "btn btn-outline-secondary btn-sm", icon = icon("file-excel"))
              )
            )
          ),
          br(),
          
          # Clustering card
          card(
            card_header("K-means Expression Trend Clustering"),
            card_body(
              layout_columns(
                col_widths = c(3, 9),
                div(
                  checkboxInput("degTab-show_elbow", "Plot Elbow Method", value = FALSE),
                  numericInput("degTab-kmeans_k", "Number of Clusters (k):", value = 4, min = 2, max = 20),
                  actionButton("degTab-send_to_go", "Export Clusters to Enrichment", class = "btn btn-outline-primary btn-sm w-100 mt-2", icon = icon("paper-plane")),
                  downloadButton("degTab-downloadKmeansCsv", "Download Clusters (.csv)", class = "btn btn-outline-secondary btn-sm w-100 mt-2")
                ),
                plotOutput("degTab-kmeansPlot", height = "480px")
              )
            )
          )
        ),
        
        # ── 6. GSEA Tab ───────────────────────────────────────────────────────
        nav_panel(
          "gsea",
          gseaUI("gseaTab")
        ),

        # ── 7. GO Enrichment Tab ──────────────────────────────────────────────
        nav_panel(
          "go",
          div(class = "page-title", "GO / Pathway Enrichment"),
          div(class = "page-subtitle", "Evaluate functional pathways enriched in DEGs."),
          
          layout_columns(
            col_widths = c(4, 8),
            card(
              card_header("Enrichment Setup"),
              card_body(
                selectInput("go_module-gene_set", "DEG List Filter:",
                            choices = c("Up-regulated" = "up", "Down-regulated" = "down", "All Significant" = "all_sig"), selected = "up"),
                selectInput("go_module-analysis_type", "Pathway Catalog:",
                            choices = c("All (GO + KEGG + Reactome)" = "ALL", "GO: Biological Process" = "BP", "GO: Molecular Function" = "MF", "GO: Cellular Component" = "CC", "KEGG" = "KEGG", "Reactome" = "REACTOME"), selected = "ALL"),
                selectInput("go_module-go_id_display_type", "Display Gene ID Type:", choices = c("Symbol" = "SYMBOL", "Entrez ID" = "ENTREZID"), selected = "SYMBOL"),
                numericInput("go_module-pvalue_cutoff", "P-Value cutoff:", value = 0.05, min = 0, max = 1),
                numericInput("go_module-qvalue_cutoff", "q-Value cutoff:", value = 0.20, min = 0, max = 1),
                hr(),
                actionButton("go_module-run_analysis", "Run Pathway Enrichment", class = "btn btn-primary w-100", icon = icon("play"))
              )
            ),
            card(
              card_header("Enrichment Results"),
              card_body(
                layout_columns(
                  col_widths = c(4, 4, 4),
                  div(class = "value-box-card",
                      div(class = "value-box-val", textOutput("pathway_box_go", inline = TRUE)),
                      div(class = "value-box-lbl", "GO Terms")),
                  div(class = "value-box-card",
                      div(class = "value-box-val", textOutput("pathway_box_kegg", inline = TRUE)),
                      div(class = "value-box-lbl", "KEGG Pathways")),
                  div(class = "value-box-card",
                      div(class = "value-box-val", textOutput("pathway_box_reactome", inline = TRUE)),
                      div(class = "value-box-lbl", "Reactome Pathways"))
                ),
                hr(),
                DTOutput("go_module-goResultTable"),
                br(),
                downloadButton("go_module-downloadExcelResults", "Download Excel (.xlsx)", class = "btn btn-outline-secondary btn-sm", icon = icon("file-excel")),
                downloadButton("go_module-downloadCsvResults", "Download CSV", class = "btn btn-outline-secondary btn-sm", icon = icon("file-csv"))
              )
            )
          ),
          br(),
          
          # Enrichment Plot Views
          card(
            card_header("Enrichment Visualizations"),
            card_body(
              tabsetPanel(
                tabPanel("Dot Plot",
                         br(),
                         uiOutput("go_module-plot_selector_ui_dot"),
                         numericInput("go_module-dotplot_n", "Enriched terms limit:", value = 10, min = 1),
                         plotOutput("go_module-goDotPlot", height = "580px"),
                         downloadButton("go_module-downloadDotPlot", "Download Dot Plot", class = "btn btn-outline-secondary btn-sm")
                ),
                tabPanel("Bar Plot",
                         br(),
                         uiOutput("go_module-plot_selector_ui_bar"),
                         numericInput("go_module-barplot_n", "Enriched terms limit:", value = 10, min = 1),
                         plotOutput("go_module-goBarPlot", height = "580px"),
                         downloadButton("go_module-downloadBarPlot", "Download Bar Plot", class = "btn btn-outline-secondary btn-sm")
                ),
                tabPanel("Network (Cnet) Plot",
                         br(),
                         uiOutput("go_module-plot_selector_ui_net"),
                         numericInput("go_module-netplot_n", "Enriched terms limit:", value = 5, min = 1),
                         plotOutput("go_module-goNetPlot", height = "580px"),
                         downloadButton("go_module-downloadNetPlot", "Download Network Plot", class = "btn btn-outline-secondary btn-sm")
                )
              )
            )
          )
        ),
        
        # ── 8. Time-series Tab ────────────────────────────────────────────────
        nav_panel(
          "timeseries",
          timeseriesAnalysisUI("timeseriesTab")
        ),

        # ── 9. Deconvolution Tab ─────────────────────────────────────────────
        nav_panel(
          "deconv",
          deconvolutionUI("deconvTab")
        ),

        # ── 10. Gene Barplot & Swap Check Tab ───────────────────────────────
        nav_panel(
          "swap",
          geneBarplotSwapUI("swapTab")
        ),

        # ── 11. Figures & Non-model Enrichment Tab ──────────────────────────
        nav_panel(
          "figenrich",
          figureEnrichmentUI("figEnrichTab")
        ),

        # ── 12. Export Tab ────────────────────────────────────────────────────
        nav_panel(
          "export",
          div(class = "page-title", "Export Outputs"),
          div(class = "page-subtitle", "Download final analysis results and data matrix packages."),
          
          layout_columns(
            col_widths = c(4, 4, 4),
            card(
              card_header("Differential Expression Table"),
              card_body(
                p("Save the full table of differentially expressed genes based on the current threshold settings.", class = "text-muted small"),
                br(), br(),
                downloadButton("degTab-downloadCsvResults", "Download DEG CSV", class = "btn btn-primary w-100", icon = icon("file-csv")),
                br(), br(),
                downloadButton("degTab-downloadExcelResults", "Download DEG Excel (.xlsx)", class = "btn btn-outline-secondary w-100", icon = icon("file-excel"))
              )
            ),
            card(
              card_header("Enriched Pathway Lists"),
              card_body(
                p("Save the full enrichment tables detailing the biological processes and KEGG/Reactome pathways.", class = "text-muted small"),
                br(), br(),
                downloadButton("go_module-downloadCsvResults", "Download Pathway CSV", class = "btn btn-primary w-100", icon = icon("file-csv")),
                br(), br(),
                downloadButton("go_module-downloadExcelResults", "Download Pathway Excel (.xlsx)", class = "btn btn-outline-secondary w-100", icon = icon("file-excel"))
              )
            ),
            card(
              card_header("Session workspace (.rds)"),
              card_body(
                p("Export the entire environment state containing metadata annotations, counts, and analysis runs.", class = "text-muted small"),
                br(), br(),
                downloadButton("dataTab-downloadRDS", "Download Session RDS", class = "btn btn-success w-100", icon = icon("save"))
              )
            )
          )
        )
      )
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
        gene_annotation = NULL,   # GTFアップロード由来: data.frame(Geneid, gene_name, biotype, gene_length)
        gtf_id_choices = NULL,    # GTFの中身から決めた表示IDタイプの選択肢/既定 (gtf_display_id_choices)
        gene2go_annotation = NULL,# 非モデル生物GO用 gene2go注釈: data.frame(GeneID, GO, [Term], [Category])
        kegg_organism_code = NULL,# 任意のKEGG生物種コード (例 lja) — enrichKEGGで使用
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
    processingServer("qcProcTab", rv)
    dimensionReductionServer("dimRedTab", rv)
    degAnalysisServer("degTab", rv)
    gseaServer("gseaTab", rv)
    enrich_results <- goEnrichmentIntegratedServer("go_module", rv = rv, deg_results_reactive=reactive(rv$deg_results), background_genes_reactive=reactive(rv$background_genes_original), selected_species_code_reactive=reactive(rv$selected_species), gene_annotation_reactive=reactive(rv$gene_annotation), gtf_id_choices_reactive=reactive(rv$gtf_id_choices))
    timeseriesAnalysisServer("timeseriesTab", rv)
    deconvolutionServer("deconvTab", rv)
    geneBarplotSwapServer("swapTab", rv)
    figureEnrichmentServer("figEnrichTab")

    # ==========================================
    # CUSTOM INTERACTIVE UI LOGIC & NAVIGATION
    # ==========================================
    active_tab <- reactiveVal("upload")
    
    observe({
      tab <- active_tab()
      updateTabsetPanel(session, "main_tabs", tab)
      shinyjs::runjs(sprintf("setActiveTab('%s');", tab))
    })
    
    observeEvent(input$tab_fastq, { active_tab("fastq") })
    observeEvent(input$tab_upload, { active_tab("upload") })
    observeEvent(input$tab_qc, { active_tab("qc") })
    observeEvent(input$tab_vis, { active_tab("vis") })
    observeEvent(input$tab_deg, { active_tab("deg") })
    observeEvent(input$tab_gsea, { active_tab("gsea") })
    observeEvent(input$tab_go, { active_tab("go") })
    observeEvent(input$tab_timeseries, { active_tab("timeseries") })
    observeEvent(input$tab_deconv, { active_tab("deconv") })
    observeEvent(input$tab_swap, { active_tab("swap") })
    observeEvent(input$tab_figenrich, { active_tab("figenrich") })
    observeEvent(input$tab_export, { active_tab("export") })
    
    # QC Page Summary Cards
    output$qc_box_samples <- renderText({
      req(rv$sample_metadata)
      nrow(rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE])
    })
    
    output$qc_box_reads <- renderText({
      req(rv$merged_data, rv$sample_metadata)
      active_samples <- rv$sample_metadata$current_name[rv$sample_metadata$active]
      counts_cols <- intersect(active_samples, colnames(rv$merged_data))
      if (length(counts_cols) > 0) {
        total_reads <- sum(colSums(rv$merged_data[, counts_cols, drop = FALSE], na.rm = TRUE))
        if (total_reads >= 1e9) {
          paste0(round(total_reads / 1e9, 2), "B")
        } else if (total_reads >= 1e6) {
          paste0(round(total_reads / 1e6, 1), "M")
        } else {
          format(total_reads, big.mark = ",")
        }
      } else {
        "0"
      }
    })
    
    output$qc_box_genes <- renderText({
      req(rv$merged_data)
      format(nrow(rv$merged_data), big.mark = ",")
    })
    
    # QC Page LogCPM reactive for PCA & Heatmap
    qc_logcpm_data <- reactive({
      req(rv$merged_data, rv$sample_metadata)
      active_meta <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      shiny::validate(shiny::need(nrow(active_meta) >= 2, "QC analysis requires at least 2 active samples."))
      
      counts_df <- rv$merged_data
      # Apply edgeR filtering if filtering keep vector exists
      if (!is.null(rv$filtered_keep) && length(rv$filtered_keep) == nrow(counts_df)) {
        counts_df <- counts_df[rv$filtered_keep, ]
      }
      
      counts_matrix <- as.matrix(counts_df[, intersect(colnames(counts_df), active_meta$current_name), drop = FALSE])
      rownames(counts_matrix) <- counts_df$Geneid
      
      y <- edgeR::DGEList(counts = counts_matrix, group = factor(active_meta$group))
      y <- edgeR::calcNormFactors(y)
      logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
      list(logcpm = logcpm, meta = active_meta)
    })
    
    # QC PCA Plot
    output$qc_pca <- renderPlotly({
      data <- qc_logcpm_data()
      req(data)
      logcpm <- data$logcpm
      meta <- data$meta
      
      gene_vars <- matrixStats::rowVars(logcpm)
      select_genes <- order(gene_vars, decreasing = TRUE)[1:min(500, nrow(logcpm))]
      pca_res <- prcomp(t(logcpm[select_genes, , drop = FALSE]), scale. = TRUE)
      
      plot_df <- data.frame(
        PC1 = pca_res$x[, 1],
        PC2 = pca_res$x[, 2],
        Sample = colnames(logcpm),
        Group = meta$group
      )
      
      p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = Group, text = paste("Sample:", Sample, "<br>Group:", Group))) +
        geom_point(size = 3.5, alpha = 0.8) +
        scale_color_brewer(palette = "Set1") +
        theme_minimal(base_size = 11) +
        theme(panel.grid.minor = element_blank(), panel.border = element_rect(fill=NA, color="#e5e7eb")) +
        labs(x = "PC1", y = "PC2")
        
      ggplotly(p, tooltip = "text") %>% layout(margin = list(t = 20, b = 20))
    })
    
    # Metadata Validation Banner
    output$metadata_validation <- renderUI({
      req(rv$sample_metadata)
      df <- rv$sample_metadata
      df_active <- df[df$active, , drop = FALSE]
      if (nrow(df_active) == 0) {
        return(div(class = "alert alert-warning", "⚠ No active samples selected. Please activate at least one sample in the table below."))
      }
      has_missing <- any(is.na(df_active) | df_active == "" | df_active == "NA")
      if (has_missing) {
        div(class = "alert alert-warning",
            tags$span(style = "font-weight: 600;", "⚠ Missing values detected"),
            " - Some active samples have empty metadata cells. Please complete the table to ensure proper downstream analysis."
        )
      } else {
        div(class = "alert alert-success",
            tags$span(style = "font-weight: 600;", "✓ Metadata validated"),
            " - All active samples have complete annotations. You can proceed with downstream analysis."
        )
      }
    })
    
    # Pathway Enrichment summary boxes
    output$pathway_box_go <- renderText({
      results <- enrich_results()
      if (is.null(results)) return("0")
      go_cats <- intersect(names(results), c("BP", "MF", "CC"))
      total <- 0
      for (cat in go_cats) {
        df <- as.data.frame(results[[cat]])
        if (!is.null(df) && nrow(df) > 0) total <- total + nrow(df)
      }
      as.character(total)
    })
    
    output$pathway_box_kegg <- renderText({
      results <- enrich_results()
      if (is.null(results) || !"KEGG" %in% names(results)) return("0")
      df <- as.data.frame(results[["KEGG"]])
      if (is.null(df)) return("0")
      as.character(nrow(df))
    })
    
    output$pathway_box_reactome <- renderText({
      results <- enrich_results()
      if (is.null(results) || !"REACTOME" %in% names(results)) return("0")
      df <- as.data.frame(results[["REACTOME"]])
      if (is.null(df)) return("0")
      as.character(nrow(df))
    })
}

shinyApp(ui, server)