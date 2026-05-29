#' Report Module — Reproducible R script + HTML summary report
#' @import shiny
#' @importFrom bslib layout_sidebar sidebar card card_header card_body value_box
#' @importFrom bsicons bs_icon
#' @noRd
mod_report_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 350,
      card(
        card_header(bsicons::bs_icon("file-earmark-text"), " Report Options"),
        card_body(
          checkboxGroupInput(
            ns("sections"), "Include sections",
            choices  = c("QC & Filter" = "qc", "EDA" = "eda",
                         "DEG"         = "deg", "Enrichment" = "enrichment"),
            selected = c("qc", "eda", "deg")
          ),
          tags$hr(),
          h6("Figure Export Defaults"),
          selectInput(ns("fig_preset"), "Journal preset",
                      choices  = c("Nature 1-col (89mm)", "Nature 2-col (183mm)",
                                   "Cell 1-col (85mm)", "Science 1-col (57mm)", "Custom"),
                      selected = "Nature 1-col (89mm)"),
          selectInput(ns("fig_font"), "Font",
                      choices  = c("Helvetica", "Arial", "Times New Roman"),
                      selected = "Helvetica"),
          tags$hr(),
          downloadButton(ns("dl_rscript"), "R Script",
                         class = "btn-outline-secondary w-100 mb-2"),
          downloadButton(ns("dl_report"),  "HTML Report",
                         class = "btn-primary w-100")
        )
      )
    ),
    uiOutput(ns("main_content"))
  )
}

#' @noRd
mod_report_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$main_content <- renderUI({
      if (!state$step_ready("upload")) {
        return(ui_empty_state("lock", "Report Locked",
                              "Complete at least the Upload step."))
      }
      tagList(
        # ── Pipeline summary cards ───────────────────────────────────────────
        fluidRow(
          column(2, .status_box("Upload",     state$pipeline_status()$upload)),
          column(2, .status_box("QC",         state$pipeline_status()$qc)),
          column(2, .status_box("EDA",        state$pipeline_status()$eda)),
          column(2, .status_box("DEG",        state$pipeline_status()$deg)),
          column(2, .status_box("Enrichment", state$pipeline_status()$enrichment)),
          column(2, .status_box("Report",     "pending"))
        ),
        # ── Parameter log ────────────────────────────────────────────────────
        card(
          card_header("Analysis Parameters (for reproducibility)"),
          card_body(
            verbatimTextOutput(ns("params_summary"))
          )
        ),
        # ── DEG summary ──────────────────────────────────────────────────────
        if (state$step_ready("deg") && length(state$deg_results()) > 0) {
          results  <- state$deg_results()
          contrast <- names(results)[length(results)]
          res      <- results[[length(results)]]
          card(
            card_header(paste("DEG Summary:", contrast)),
            card_body(
              fluidRow(
                column(4, value_box("Sig. DEG",
                  sum(!is.na(res$padj) & res$padj < 0.05),
                  bsicons::bs_icon("asterisk"), theme = "primary")),
                column(4, value_box("Up",
                  sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1),
                  bsicons::bs_icon("arrow-up"), theme = "danger")),
                column(4, value_box("Down",
                  sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1),
                  bsicons::bs_icon("arrow-down"), theme = "info"))
              )
            )
          )
        }
      )
    })

    output$params_summary <- renderPrint({
      log <- state$params_log()
      if (length(log) == 0) {
        cat("No parameters logged yet. Run analyses first.\n")
        return(invisible(NULL))
      }
      cat("=== EasyRNA-Seq Analysis Parameters ===\n")
      cat(sprintf("Generated: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      prev_step <- ""
      for (nm in names(log)) {
        e <- log[[nm]]
        if (e$step != prev_step) {
          cat(sprintf("\n[%s]\n", toupper(e$step)))
          prev_step <- e$step
        }
        cat(sprintf("  %-25s = %s\n", e$key, paste(e$value, collapse = ", ")))
      }
    })

    # ── R Script download ─────────────────────────────────────────────────────
    output$dl_rscript <- downloadHandler(
      filename = "EasyRNASeq_analysis.R",
      content  = function(file) {
        log <- state$params_log()
        preset_dims <- list(
          "Nature 1-col (89mm)"  = c(89,  89),
          "Nature 2-col (183mm)" = c(183, 120),
          "Cell 1-col (85mm)"    = c(85,  85),
          "Science 1-col (57mm)" = c(57,  57),
          "Custom"               = c(120, 90)
        )
        dims <- preset_dims[[input$fig_preset]] %||% c(89, 89)

        lines <- c(
          "# ============================================================",
          "# EasyRNA-Seq v4.0 — Reproducible R Script",
          paste0("# Generated : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
          paste0("# Journal   : ", input$fig_preset),
          "# ============================================================",
          "",
          "library(DESeq2)",
          "library(ggplot2)",
          if (isTRUE("apeglm" %in% rownames(installed.packages()))) "library(apeglm)",
          "",
          "# ── Load data ──────────────────────────────────────────────────",
          '# counts  <- read.csv("counts.csv", row.names=1, check.names=FALSE)',
          '# metadata <- read.csv("metadata.csv")',
          "",
          "# ── Parameters ─────────────────────────────────────────────────"
        )

        for (nm in names(log)) {
          e <- log[[nm]]
          val_str <- if (is.character(e$value)) paste0('"', e$value, '"') else e$value
          lines <- c(lines,
            sprintf("# [%s] %s = %s", e$step, e$key,
                    paste(val_str, collapse = ", ")))
        }

        # DESeq2 template
        deg_params <- log[grep("^deg\\.", names(log))]
        cond  <- deg_params[["deg.condition_col"]]$value %||% "condition"
        ref   <- deg_params[["deg.ref_level"]]$value     %||% "control"
        test  <- deg_params[["deg.test_level"]]$value    %||% "treated"

        lines <- c(lines, "",
          "# ── DESeq2 ─────────────────────────────────────────────────────",
          sprintf('meta$%s <- relevel(factor(meta$%s), ref = "%s")', cond, cond, ref),
          sprintf('dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta, design=~%s)', cond),
          "dds <- DESeq(dds)",
          sprintf('res <- lfcShrink(dds, coef="%s_%s_vs_%s", type="apeglm")', cond, test, ref),
          "res_df <- as.data.frame(res)",
          "",
          "# ── Volcano plot ────────────────────────────────────────────────",
          "res_df$sig <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,",
          "                     ifelse(res_df$log2FoldChange > 0, 'Up', 'Down'), 'NS')",
          "p_volcano <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color=sig)) +",
          "  geom_point(alpha=0.5, size=1.2) +",
          '  scale_color_manual(values=c(Up="#D55E00", Down="#0072B2", NS="#AAAAAA")) +',
          "  theme_minimal(base_size=7)",
          "",
          sprintf('ggsave("volcano.pdf", plot=p_volcano, width=%s, height=%s, units="mm",', dims[1], dims[2]),
          '       device=cairo_pdf)'
        )

        writeLines(lines, file)
      }
    )

    # ── HTML Report download ──────────────────────────────────────────────────
    output$dl_report <- downloadHandler(
      filename = paste0("EasyRNASeq_report_",
                        format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
      content  = function(file) {
        # Build inline HTML report (no Quarto/rmarkdown dependency)
        log    <- state$params_log()
        counts <- state$counts()
        filt   <- state$counts_filtered()
        results <- state$deg_results()

        sections <- c(
          "<html><head>",
          "<meta charset='UTF-8'>",
          "<style>",
          "body{font-family:Helvetica,Arial,sans-serif;max-width:900px;margin:2em auto;font-size:13px}",
          "h1{color:#2C3E50}h2{color:#18BC9C;border-bottom:1px solid #eee;padding-bottom:.3em}",
          "table{border-collapse:collapse;width:100%}",
          "th,td{padding:4px 8px;border:1px solid #ddd;text-align:left}",
          "th{background:#f5f5f5}",
          "pre{background:#f8f8f8;padding:1em;border-radius:4px;overflow:auto}",
          ".badge{display:inline-block;padding:2px 6px;border-radius:3px;font-size:11px}",
          ".done{background:#18BC9C;color:white}.pending{background:#ccc;color:#555}",
          "</style></head><body>",
          "<h1>EasyRNA-Seq v4.0 — Analysis Report</h1>",
          sprintf("<p>Generated: %s</p>", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
          "<hr>"
        )

        # Pipeline status
        status <- state$pipeline_status()
        sections <- c(sections, "<h2>Pipeline Status</h2><table><tr>")
        for (step in names(status)) {
          cls <- if (status[[step]] == "done") "done" else "pending"
          sections <- c(sections,
            sprintf("<td><span class='badge %s'>%s</span><br>%s</td>",
                    cls, status[[step]], step))
        }
        sections <- c(sections, "</tr></table>")

        # Data summary
        if (!is.null(counts)) {
          sections <- c(sections,
            "<h2>Data Summary</h2><table>",
            sprintf("<tr><th>Genes (raw)</th><td>%s</td></tr>",    format(nrow(counts), big.mark=",")),
            sprintf("<tr><th>Samples</th><td>%s</td></tr>",         ncol(counts)),
            if (!is.null(filt))
              sprintf("<tr><th>Genes (filtered)</th><td>%s</td></tr>", format(nrow(filt), big.mark=",")),
            "</table>")
        }

        # Parameters
        if (length(log) > 0) {
          sections <- c(sections, "<h2>Analysis Parameters</h2><pre>")
          prev_step <- ""
          for (nm in names(log)) {
            e <- log[[nm]]
            if (e$step != prev_step) {
              sections <- c(sections, sprintf("\n[%s]", toupper(e$step)))
              prev_step <- e$step
            }
            sections <- c(sections,
              sprintf("  %-25s = %s", e$key, paste(e$value, collapse=", ")))
          }
          sections <- c(sections, "</pre>")
        }

        # DEG table (top 50)
        if (length(results) > 0) {
          key <- names(results)[length(results)]
          res <- results[[length(results)]]
          sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
          sig <- head(sig[order(sig$padj), ], 50)
          sections <- c(sections,
            sprintf("<h2>Top DEGs: %s (FDR &lt; 0.05, top 50)</h2>", key),
            "<table><tr><th>Gene</th><th>log2FC</th><th>FDR</th></tr>")
          for (i in seq_len(nrow(sig))) {
            gene <- if ("gene" %in% colnames(sig)) sig$gene[i] else rownames(sig)[i]
            lfc  <- round(sig$log2FoldChange[i], 3)
            padj <- format(sig$padj[i], digits=3, scientific=TRUE)
            col  <- if (!is.na(lfc) && lfc > 0) "#FFD7CC" else "#CCE5FF"
            sections <- c(sections,
              sprintf("<tr><td>%s</td><td style='background:%s'>%s</td><td>%s</td></tr>",
                      gene, col, lfc, padj))
          }
          sections <- c(sections, "</table>")
        }

        sections <- c(sections, "</body></html>")
        writeLines(sections, file)
      }
    )
  })
}

# ── Helper: pipeline status badge ─────────────────────────────────────────────
.status_box <- function(label, status) {
  theme <- if (status == "done") "success" else "light"
  icon_name <- if (status == "done") "check-circle-fill" else "circle"
  value_box(label, status, bsicons::bs_icon(icon_name), theme = theme)
}
