# R/module_dimension_reduction.R
# (Dendrogramを静的プロットとして描画し、ラベル回転を確実にする修正)

# --- 必要なライブラリ ---
library(shiny)
library(shinycssloaders)
library(edgeR)
library(ggplot2)
library(plotly)
library(dendextend)
library(ggdendro)
library(RColorBrewer)
library(pheatmap)
library(igraph)
library(Rtsne)
library(umap)
library(matrixStats)
library(dplyr)

dimensionReductionUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h4("サンプル間類似性の可視化"),
    helpText("サンプル間の全体的な類似性を、選択した手法で可視化します。フィルタリング後のデータに基づき計算されます。"),
    hr(),
    fluidRow(
      column(4, selectInput(ns("dim_reduction_method"), "可視化手法:", choices = c("MDS (edgeR)", "PCA", "t-SNE", "UMAP", "階層的クラスタリング (デンドログラム)" = "Dendrogram", "相関ヒートマップ" = "Heatmap", "ネットワークグラフ" = "Network"), selected = "MDS (edgeR)")),
      column(8,
             conditionalPanel( condition = paste0("input['", ns("dim_reduction_method"), "'] == 'MDS (edgeR)'"), numericInput(ns("mds_top_genes"), "使用遺伝子数(Top N):", 500, min=50, step=50) ),
             conditionalPanel( condition = paste0("input['", ns("dim_reduction_method"), "'] == 'PCA'"), numericInput(ns("pca_top_genes"), "使用遺伝子数(Top N変動):", 500, min=50, step=50) ),
             conditionalPanel( condition = paste0("input['", ns("dim_reduction_method"), "'] == 't-SNE'"), fluidRow(column(6, numericInput(ns("tsne_top_genes"), "使用遺伝子数(Top N変動):", 500, min=50, step=50)), column(6, numericInput(ns("tsne_perplexity"), "Perplexity:", 10, min=2, step=1))), helpText("Perplexityは サンプル数(n)に対し 2 ~ floor((n-1)/3) の範囲で設定。") ),
             conditionalPanel( condition = paste0("input['", ns("dim_reduction_method"), "'] == 'UMAP'"), fluidRow(column(6, numericInput(ns("umap_top_genes"), "使用遺伝子数(Top N変動):", 500, min=50, step=50)), column(6, numericInput(ns("umap_n_neighbors"), "近傍数(n_neighbors):", 15, min=2, step=1))) ),
             conditionalPanel(
               condition = paste0("input['", ns("dim_reduction_method"), "'] == 'Dendrogram'"),
               fluidRow(
                 column(6, selectInput(ns("dist_method"), "距離計算 (サンプル間):", choices = c("euclidean", "maximum", "manhattan", "canberra", "minkowski"), selected = "euclidean")),
                 column(6, selectInput(ns("hclust_method"), "クラスタリング:", choices = c("ward.D2", "ward.D", "complete", "average"), selected = "ward.D2"))
               ),
               helpText("logCPM値に基づいてサンプル間の距離を計算し、クラスタリングします。"),
               hr(),
               downloadButton(ns("downloadPlot"), "プロットをダウンロード (.png)", icon = icon("download"))
             ),
             conditionalPanel(
               condition = paste0("input['", ns("dim_reduction_method"), "'] == 'Heatmap'"),
               fluidRow(
                 column(4, selectInput(ns("heatmap_dist_method"), "距離計算 (クラスタリング用):", choices = c("euclidean", "maximum", "manhattan", "canberra", "minkowski"), selected = "euclidean")),
                 column(4, selectInput(ns("heatmap_hclust_method"), "クラスタリング:", choices = c("ward.D2", "ward.D", "complete", "average"), selected = "ward.D2")),
                 column(4, selectInput(ns("cor_method"), "相関係数 (色表示用):", choices = c("pearson", "spearman"), selected = "pearson"))
               ),
               checkboxInput(ns("heatmap_show_values"), "相関係数値を表示", value = FALSE),
               checkboxInput(ns("heatmap_cluster_cols"), "列(サンプル)をクラスタリング", value = TRUE),
               checkboxInput(ns("heatmap_cluster_rows"), "行(サンプル)をクラスタリング", value = TRUE)
             ),
             conditionalPanel( condition = paste0("input['", ns("dim_reduction_method"), "'] == 'Network'"), fluidRow(column(4, selectInput(ns("network_cor_method"), "相関係数:", choices = c("pearson", "spearman"), selected = "pearson")), column(4, sliderInput(ns("network_cor_threshold"), "相関閾値(エッジ):", min = 0, max = 1, value = 0.8, step = 0.05)), column(4, selectInput(ns("network_layout"), "レイアウト:", choices = c("nicely", "Fruchterman-Reingold"="fr", "Circle"="circle", "Kamada-Kawai"="kk"), selected = "nicely")) ) )
      )
    ),
    hr(),
    withSpinner(uiOutput(ns("plotArea")), type = 6)
  )
}

dimensionReductionServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    prepared_data <- reactive({
      req(rv$merged_data, rv$sample_metadata, rv$filtered_keep)
      
      active_samples_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      validate(need(nrow(active_samples_metadata) > 0, "解析対象のサンプルが選択されていません。データ入力タブで少なくとも1つのサンプルをチェックしてください。"))
      active_sample_names <- active_samples_metadata$current_name
      
      keep_vector <- rv$filtered_keep
      counts_df <- rv$merged_data
      validate(need(length(keep_vector) == nrow(counts_df), "フィルタリング情報とカウントデータの行数が一致しません。フィルタリングタブを再実行してください。"))
      
      counts_df_active <- counts_df[, c("Geneid", intersect(colnames(counts_df), active_sample_names)), drop = FALSE]
      
      filtered_counts_df <- counts_df_active[keep_vector, ]
      validate(need(nrow(filtered_counts_df) > 0, "フィルタリングの結果、遺伝子が残りませんでした。"))
      
      count_matrix <- filtered_counts_df[,-1, drop=FALSE]
      rownames(count_matrix) <- filtered_counts_df$Geneid
      
      sample_order <- colnames(count_matrix)
      samples_metadata_ordered <- active_samples_metadata[match(sample_order, active_samples_metadata$current_name),]
      group <- factor(samples_metadata_ordered$group)
      
      y <- DGEList(counts = count_matrix, group = group)
      y <- calcNormFactors(y)
      logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
      validate(need(nrow(logcpm) >= 2, "logCPM計算の対象となる遺伝子が少なすぎます。"))
      
      list(y = y, logcpm = logcpm, group = group, samples_metadata = samples_metadata_ordered)
    })
    
    plot_object_reactive <- reactiveVal(NULL)
    
    dimension_reduction_results <- reactive({
      prep_data <- prepared_data(); req(prep_data)
      selected_method <- input$dim_reduction_method
      y <- prep_data$y; logcpm <- prep_data$logcpm; group <- prep_data$group
      samples_metadata <- prep_data$samples_metadata; n_samples <- ncol(y)
      plot_object <- NULL; plot_type <- "other"; title_suffix <- ""
      
      results <- tryCatch({
        # ★★★ ここでプロットのタイプを決定 ★★★
        plot_type <- if (selected_method == "Dendrogram") "static" else "plotly"
        
        if (selected_method == "Dendrogram") {
          req(input$dist_method, input$hclust_method)
          validate(need(n_samples >= 3, "クラスタリングには最低3サンプル必要です。"))
          data_for_dist <- t(logcpm)
          dist_matrix <- dist(data_for_dist, method = input$dist_method)
          hclust_obj <- hclust(dist_matrix, method = input$hclust_method)
          
          dend <- as.dendrogram(hclust_obj)
          ddata <- dendro_data(dend, type = "rectangle")
          
          label_df <- ddata$labels %>%
            mutate(group = factor(samples_metadata$group[match(label, samples_metadata$current_name)]))
          
          unique_groups <- levels(label_df$group)
          num_colors <- length(unique_groups)
          color_palette <- RColorBrewer::brewer.pal(max(3, min(num_colors, 9)), "Set1")[1:num_colors]
          names(color_palette) <- unique_groups
          
          y_range <- range(ddata$segments$y, na.rm = TRUE)
          y_min_for_text <- - (y_range[2] * 0.1)
          
          p <- ggplot() +
            geom_segment(data = segment(ddata), aes(x = x, y = y, xend = xend, yend = yend)) +
            geom_text(data = label_df,
                      aes(x = x, y = y_min_for_text, label = label, color = group),
                      angle = 90,
                      hjust = 1,
                      vjust = 0.5,
                      size = 3.5) +
            scale_color_manual(name = "Group", values = color_palette) +
            scale_y_continuous(expand = expansion(mult = c(0.25, 0.1))) +
            labs(y = paste("Height (", input$dist_method, "distance )"), x = "") +
            theme_bw(base_size = 12) +
            theme(
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              legend.position = "top",
              plot.margin = margin(t = 20, r = 20, b = 40, l = 20),
              clip = "off"
            )
          title_suffix <- paste("(Dist:", input$dist_method, ", Linkage:", input$hclust_method, ")")
          plot_object <- p
          plot_object_reactive(p)
          
        } else if (selected_method %in% c("MDS (edgeR)", "PCA", "t-SNE", "UMAP")) {
          plot_df <- NULL; x_lab <- "Dim 1"; y_lab <- "Dim 2";
          if (selected_method == "MDS (edgeR)") {
            req(input$mds_top_genes); validate(need(n_samples >= 3, "MDSには最低3サンプル必要です。"));
            mds_results <- plotMDS(y, top = input$mds_top_genes, plot = FALSE);
            plot_df <- data.frame( Dim1 = mds_results$x, Dim2 = mds_results$y, Sample = colnames(y), Group = group, stringsAsFactors = FALSE );
            var_explained <- mds_results$var.explained; x_lab <- paste0("Leading logFC dim 1 (", round(var_explained[1]*100, 1), "%)"); y_lab <- paste0("Leading logFC dim 2 (", round(var_explained[2]*100, 1), "%)");
            title_suffix <- paste("(Top", input$mds_top_genes, "genes)")
          } else {
            top_genes_n <- switch(selected_method, "PCA"=input$pca_top_genes, "t-SNE"=input$tsne_top_genes, "UMAP"=input$umap_top_genes, 500);
            req(top_genes_n); gene_vars <- matrixStats::rowVars(logcpm); n_select <- min(top_genes_n, nrow(logcpm));
            validate(need(n_select >= 2, "変動上位遺伝子が少なすぎます。")); select_genes <- order(gene_vars, decreasing = TRUE)[1:n_select];
            data_for_dim_red <- t(logcpm[select_genes, , drop = FALSE]); title_suffix <- paste("(Top", n_select, "variable genes)");
            if (selected_method == "PCA") {
              validate(need(n_samples >= 3, "PCAには最低3サンプル必要です。")); pca_results <- prcomp(data_for_dim_red, scale. = TRUE); pca_summary <- summary(pca_results);
              var_explained <- pca_summary$importance[2, 1:2]; plot_df <- data.frame( Dim1 = pca_results$x[, 1], Dim2 = pca_results$x[, 2], Sample = rownames(pca_results$x), Group = group[match(rownames(pca_results$x), colnames(y))], stringsAsFactors = FALSE );
              x_lab <- paste0("PC 1 (", round(var_explained[1]*100, 1), "%)"); y_lab <- paste0("PC 2 (", round(var_explained[2]*100, 1), "%)")
            } else if (selected_method == "t-SNE") {
              min_samples_tsne <- 7; validate(need(n_samples >= min_samples_tsne, paste0("t-SNEには最低",min_samples_tsne,"サンプル必要です。")));
              req(input$tsne_perplexity); input_perplexity <- input$tsne_perplexity; max_valid_perplexity <- floor((n_samples-1) / 3); min_valid_perplexity <- 2;
              validate(need(input_perplexity >= min_valid_perplexity && input_perplexity <= max_valid_perplexity, paste0("Perplexityは",min_valid_perplexity,"-",max_valid_perplexity,"の範囲で設定してください。")));
              tsne_results <- Rtsne::Rtsne(data_for_dim_red, dims = 2, perplexity = input_perplexity, verbose = FALSE, check_duplicates = FALSE, pca = FALSE);
              plot_df <- data.frame( Dim1 = tsne_results$Y[, 1], Dim2 = tsne_results$Y[, 2], Sample = rownames(data_for_dim_red), Group = group[match(rownames(data_for_dim_red), colnames(y))], stringsAsFactors = FALSE );
              title_suffix <- paste(title_suffix, ", Perp:", input$tsne_perplexity)
            } else if (selected_method == "UMAP") {
              min_samples_umap <- 10; validate(need(n_samples >= min_samples_umap, paste0("UMAPには最低",min_samples_umap,"サンプルを推奨します。")));
              req(input$umap_n_neighbors); umap_config <- umap::umap.defaults; umap_config$n_neighbors <- input$umap_n_neighbors; umap_config$n_neighbors <- min(umap_config$n_neighbors, n_samples - 1);
              validate(need(umap_config$n_neighbors >= 2, "近傍数(n_neighbors)が2未満です。")); umap_results <- umap::umap(data_for_dim_red, config = umap_config);
              plot_df <- data.frame( Dim1 = umap_results$layout[, 1], Dim2 = umap_results$layout[, 2], Sample = rownames(data_for_dim_red), Group = group[match(rownames(data_for_dim_red), colnames(y))], stringsAsFactors = FALSE );
              title_suffix <- paste(title_suffix, ", Neigh:", umap_config$n_neighbors)
            }
          }
          p <- ggplot(plot_df, aes(x = Dim1, y = Dim2, color = Group, text = paste("Sample:", Sample, "<br>Group:", Group))) +
            geom_point(size = 3, alpha = 0.8) +
            labs(x = x_lab, y = y_lab) + theme_bw(base_size = 12) +
            theme(legend.title = element_text(face="bold")) +
            geom_text(aes(label=Sample), size=3, vjust=-0.5, hjust=0.5, check_overlap = TRUE)
          plot_object <- p
          plot_object_reactive(p)
        } else if (selected_method == "Heatmap") {
          plot_type <- "static"
          req(input$heatmap_dist_method, input$heatmap_hclust_method, input$cor_method)
          validate(need(n_samples >= 2, "相関計算には最低2サンプル必要です。"))
          cor_matrix <- cor(logcpm, method = input$cor_method)
          annotation_col <- data.frame(Group = group, row.names = colnames(logcpm))
          plot_object <- list( matrix = cor_matrix, annotation_col = annotation_col, clustering_method = input$heatmap_hclust_method, cluster_cols = input$heatmap_cluster_cols, cluster_rows = input$heatmap_cluster_rows, display_numbers = input$heatmap_show_values, dist_method = input$heatmap_dist_method )
          title_suffix <- paste("(Cor:", input$cor_method, ", Dist:", input$heatmap_dist_method, ", Linkage:", input$heatmap_hclust_method, ")")
        } else if (selected_method == "Network") {
          plot_type <- "static"
          req(input$network_cor_threshold, input$network_layout, input$network_cor_method)
          validate(need(n_samples >= 2, "相関計算には最低2サンプル必要です。"))
          cor_matrix <- cor(logcpm, method = input$network_cor_method)
          adj_matrix <- (abs(cor_matrix) >= input$network_cor_threshold) & (cor_matrix != 1)
          graph_obj <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
          validate(need(igraph::gsize(graph_obj) > 0, "指定した相関閾値ではサンプル間の繋がり（エッジ）がありません。"))
          igraph::V(graph_obj)$group <- group[match(igraph::V(graph_obj)$name, colnames(logcpm))]
          layout_func <- switch(input$network_layout, "fr"=igraph::layout_with_fr, "circle"=igraph::layout_in_circle, "kk"=igraph::layout_with_kk, igraph::layout_nicely)
          layout_coords <- layout_func(graph_obj)
          plot_object <- list(graph = graph_obj, layout = layout_coords)
          title_suffix <- paste("(Cor:", input$network_cor_method, ">", input$network_cor_threshold, ", Layout:", input$network_layout, ")")
        }
        list(plot_object = plot_object, plot_type = plot_type, method = selected_method, title_suffix = title_suffix)
      }, error = function(e) { validate(paste(selected_method, "計算中にエラー:", e$message)) })
      results
    })
    
    output$plotArea <- renderUI({
      results <- dimension_reduction_results(); req(results)
      ns <- session$ns
      if (results$plot_type == "plotly") {
        plotlyOutput(ns("interactivePlot"), height = "600px")
      } else if (results$plot_type == "static") {
        plotOutput(ns("staticPlot"), height = "600px")
      } else {
        tags$p("プロットタイプを判別できませんでした。")
      }
    })
    
    output$interactivePlot <- renderPlotly({
      results <- dimension_reduction_results()
      req(results, results$plot_type == "plotly", inherits(results$plot_object, "ggplot"))
      plot_title <- paste(results$method, "Plot for Active Samples", results$title_suffix)
      p <- results$plot_object
      ggplotly(p, tooltip = "text") %>%
        layout(title = list(text = plot_title, x = 0.05),
               legend = list(title = list(text = '<b> Group </b>'), orientation = "h", y = -0.2, xanchor="center", x=0.5))
    })
    
    output$staticPlot <- renderPlot({
      results <- dimension_reduction_results()
      req(results, results$plot_type == "static")
      plot_title <- paste(results$method, "Plot for Active Samples", results$title_suffix)
      plot_obj_data <- results$plot_object
      
      if (results$method == "Heatmap" || results$method == "Network") {
        if (results$method == "Heatmap") {
          req(is.list(plot_obj_data), all(c("matrix", "annotation_col") %in% names(plot_obj_data)))
          if(anyNA(plot_obj_data$matrix)) { validate("ヒートマップ描画エラー: 相関行列にNAが含まれています。") }
          color_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(100)
          pheatmap::pheatmap( plot_obj_data$matrix, annotation_col = plot_obj_data$annotation_col,
                              clustering_distance_rows = if(plot_obj_data$cluster_rows) plot_obj_data$dist_method else NA,
                              clustering_distance_cols = if(plot_obj_data$cluster_cols) plot_obj_data$dist_method else NA,
                              clustering_method = plot_obj_data$clustering_method, cluster_rows = plot_obj_data$cluster_rows,
                              cluster_cols = plot_obj_data$cluster_cols, show_rownames = (ncol(plot_obj_data$matrix) <= 40),
                              show_colnames = (ncol(plot_obj_data$matrix) <= 40), display_numbers = plot_obj_data$display_numbers,
                              number_format = "%.2f", fontsize = 8, main = plot_title, color = color_palette, na_col = "grey50" )
        } else if (results$method == "Network") {
          req(is.list(plot_obj_data), all(c("graph", "layout") %in% names(plot_obj_data)))
          g <- plot_obj_data$graph; l <- plot_obj_data$layout
          unique_groups <- unique(igraph::V(g)$group)
          if (length(unique_groups) > 0) {
            num_colors <- length(unique_groups)
            if (num_colors <= 9) { group_colors <- RColorBrewer::brewer.pal(max(3, num_colors), "Set1")[1:num_colors]
            } else { group_colors <- rainbow(num_colors) }
            names(group_colors) <- unique_groups
            vertex_colors <- group_colors[as.character(igraph::V(g)$group)]
          } else { vertex_colors <- "skyblue"; group_colors <- NULL }
          plot(g, layout = l, vertex.color = vertex_colors, vertex.label = igraph::V(g)$name, vertex.size = 10, vertex.label.cex = 0.8, edge.color = "grey50", main = plot_title)
          if (!is.null(group_colors) && length(unique_groups) > 0) {
            legend("topright", legend=unique_groups, fill=group_colors, border=NA, cex=0.8, title="Group")
          }
        }
      } else if (results$method == "Dendrogram") {
        # ★★★ Dendrogramはggplotオブジェクトを直接描画 ★★★
        req(inherits(plot_obj_data, "ggplot"))
        print(plot_obj_data + ggtitle(plot_title))
      }
    })
    
    output$downloadPlot <- downloadHandler(
      filename = function() {
        results <- dimension_reduction_results(); req(results)
        paste0(results$method, "_plot_", Sys.Date(), ".png")
      },
      content = function(file) {
        p_to_save <- plot_object_reactive(); req(p_to_save, inherits(p_to_save, "ggplot"))
        ggsave(file, plot = p_to_save, device = "png", width = 10, height = 8, dpi = 300)
      }
    )
    
  })
}