# ★★★ Load required packages at the beginning of the module ★★★
library(shiny)
library(shinyjs) # Using shinyjs for conditional UI display is often more robust
library(edgeR)
library(purrr) # to use map
library(dplyr) # to use bind_rows, mutate
library(tidyr) # to use pivot_longer
library(tibble)# to use rownames_to_column
library(ggplot2)
library(shinycssloaders)
library(data.table) # to use fread, setnames etc.

filteringUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 4,
      # Widget to enable/disable filtering
      radioButtons(ns("perform_filtering"), "Perform Filtering?",
                   choices = c("Yes" = "yes", "No" = "no"),
                   selected = "yes"),
      hr(),
      # Conditional panel that only shows when filtering is enabled
      shinyjs::useShinyjs(), # Set up shinyjs
      div(id = ns("filtering_params"),
          h4("Filtering Parameters"),
          helpText("Filter genes based on expression levels using edgeR's filterByExpr function. The group factor is automatically considered."),
          numericInput(ns("min_count"), "Minimum Count (min.count)", value = 10, min = 0, step = 1),
          helpText("Minimum count required for a gene to be considered expressed."),
          numericInput(ns("min_total_count"), "Minimum Total Count (min.total.count)", value = 15, min = 0, step = 1),
          helpText("Minimum total count across all samples for a gene to be kept."),
          numericInput(ns("large_n"), "Large N (large.n)", value = 10, min = 2, step = 1),
          helpText("Number of samples where min.prop is applied."),
          sliderInput(ns("min_prop"), "Minimum Proportion (min.prop)", value = 0.7, min = 0, max = 1, step = 0.05),
          helpText("Minimum proportion of samples that must have a count of at least min.count.")
      )
    ),
    mainPanel(
      width = 8,
      h4("Filtering Summary"),
      withSpinner(verbatimTextOutput(ns("filterSummary")), type = 6),
      hr(),
      h4("LogCPM Distribution (Before vs. After Filtering)"),
      helpText("This plot shows how low-expressed genes are removed by filtering."),
      withSpinner(plotOutput(ns("filterPlot"), height="500px"), type = 6)
    )
  )
}

filteringServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    
    # Show/hide the parameter panel based on the radio button selection
    observe({
      shinyjs::toggle(id = "filtering_params", condition = input$perform_filtering == "yes")
    })
    
    filtering_results <- reactive({
      req(rv$merged_data, rv$sample_metadata, input$perform_filtering)
      
      # ★★★ 修正箇所: アクティブなサンプルのみを抽出 ★★★
      active_samples_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      validate(need(nrow(active_samples_metadata) > 0, "解析対象のサンプルが選択されていません。データ入力タブで少なくとも1つのサンプルをチェックしてください。"))
      active_sample_names <- active_samples_metadata$current_name
      
      counts_df <- rv$merged_data
      # Geneid列とアクティブなサンプル列のみを保持
      counts_df_active <- counts_df[, c("Geneid", intersect(colnames(counts_df), active_sample_names)), drop = FALSE]
      
      count_matrix <- counts_df_active[,-1, drop=FALSE]
      rownames(count_matrix) <- counts_df_active$Geneid
      
      # メタデータもアクティブなものにフィルタリング
      samples_metadata <- active_samples_metadata
      
      sample_order <- colnames(count_matrix)
      samples_metadata_ordered <- samples_metadata[match(sample_order, samples_metadata$current_name),]
      group <- factor(samples_metadata_ordered$group)
      
      y <- DGEList(counts = count_matrix, group = group)
      n_genes_before <- nrow(y)
      validate(need(n_genes_before > 0, "Count data is not available."))
      
      if (input$perform_filtering == "no") {
        # If filtering is disabled, keep all genes
        keep <- rep(TRUE, n_genes_before)
      } else {
        # If filtering is enabled, use filterByExpr
        keep <- tryCatch({
          filterByExpr(y, group = group, min.count = input$min_count, min.total.count = input$min_total_count, large.n = input$large_n, min.prop = input$min_prop)
        }, error = function(e) {
          showNotification(paste("filterByExpr error:", e$message), type="error")
          NULL
        })
        validate(need(!is.null(keep), "An error occurred during filtering calculation."))
      }
      
      n_genes_after <- sum(keep)
      kept_genes_ratio <- round(mean(keep) * 100, 1)
      
      # このyはフィルタリング前のものなので、プロット用のy_unfilteredとして返す
      # keepベクトルを適用した後のデータは後続のプロットで作成
      list(
        n_genes_before = n_genes_before,
        n_genes_after = n_genes_after,
        kept_genes_ratio = kept_genes_ratio,
        keep_vector = keep,
        y_unfiltered = y
      )
    })
    
    observe({
      results <- filtering_results()
      if(is.list(results) && !is.null(results$keep_vector)){
        rv$filtered_keep <- results$keep_vector
        message("Filtering results (keep vector) have been stored in rv.")
      } else {
        rv$filtered_keep <- NULL
        message("Filtering results have been reset.")
      }
    })
    
    output$filterSummary <- renderPrint({
      results <- filtering_results()
      req(results)
      
      cat("--- Filtering Results (Active Samples Only) ---\n") # タイトル修正
      if(input$perform_filtering == "no"){
        cat("Filtering is currently disabled.\n\n")
        cat("Total number of genes:", results$n_genes_before, "\n")
      } else {
        if(is.na(results$n_genes_after)) {
          cat("Error or calculation not performed.\n")
        } else {
          cat("Total genes before filtering:", results$n_genes_before, "\n")
          cat("Genes after filtering:", results$n_genes_after, "\n")
          cat("Percentage of genes kept:", results$kept_genes_ratio, "%\n")
          cat("\nParameters used:\n")
          cat("  min.count =", input$min_count, "\n")
          cat("  min.total.count =", input$min_total_count, "\n")
          cat("  large.n =", input$large_n, "\n")
          cat("  min.prop =", input$min_prop, "\n")
          cat("  (Group information was automatically considered)\n")
        }
      }
      cat("\nNote: This is a preview.\nThe filtering result from the last preview will be applied in the subsequent tabs.\n")
    })
    
    output$filterPlot <- renderPlot({
      results <- filtering_results()
      # ★★★ 修正箇所: y_unfiltered を使用 ★★★
      req(results, results$keep_vector, results$y_unfiltered)
      
      y_unfiltered <- results$y_unfiltered
      keep <- results$keep_vector
      logcpm_unfiltered <- edgeR::cpm(y_unfiltered, log = TRUE, prior.count = 2)
      
      validate(need(length(keep) == nrow(logcpm_unfiltered), "Mismatch between filtering results and data."))
      
      df_unfiltered <- as.data.frame(logcpm_unfiltered) %>%
        tibble::rownames_to_column("Geneid") %>%
        tidyr::pivot_longer(-Geneid, names_to = "Sample", values_to = "logCPM") %>%
        mutate(Status = "Before Filtering")
      
      if (input$perform_filtering == "yes" && sum(keep) > 0) {
        logcpm_filtered <- logcpm_unfiltered[keep, , drop = FALSE]
        
        df_filtered <- as.data.frame(logcpm_filtered) %>%
          tibble::rownames_to_column("Geneid") %>%
          tidyr::pivot_longer(-Geneid, names_to = "Sample", values_to = "logCPM") %>%
          mutate(Status = "After Filtering")
        
        plot_df <- bind_rows(df_unfiltered, df_filtered)
        plot_df$Status <- factor(plot_df$Status, levels = c("Before Filtering", "After Filtering"))
        
        # Plot with both distributions
        ggplot(plot_df, aes(x = logCPM, color = Status)) +
          geom_density(alpha = 0.5, linewidth = 1) +
          scale_color_manual(values = c("Before Filtering" = "grey50", "After Filtering" = "steelblue")) +
          labs(title = "Comparison of LogCPM Distributions (Active Samples)",
               x = "Log2(CPM + 2)",
               y = "Density",
               color = "Filtering Status") +
          theme_bw(base_size = 14) +
          theme(legend.position = "top")
        
      } else {
        # Plot only the 'Before' distribution
        ggplot(df_unfiltered, aes(x = logCPM)) +
          geom_density(alpha = 0.5, linewidth = 1, color = "grey50", fill = "grey50") +
          labs(title = "LogCPM Distribution (Active Samples)",
               subtitle = "Filtering is disabled",
               x = "Log2(CPM + 2)",
               y = "Density") +
          theme_bw(base_size = 14)
      }
    })
    
  })
}