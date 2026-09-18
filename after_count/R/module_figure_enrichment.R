# R/module_figure_enrichment.R
# タブ10: 「アップロード作図 & 非モデル生物エンリッチメント」
# 既存モジュール（msigdbr/OrgDbベース）には手を加えず、自己完結で以下を提供:
#   (1) DEG表 → Volcano / MA プロット
#   (2) GO/KEGG ORA結果 → dotplot / barplot
#   (3) GSEA結果 → NESバープロット
#   (4) カウント行列 → Heatmap（選択遺伝子）
#   (5) 非モデル生物: DEG表 + eggNOG注釈 → GO/KEGG ORA+GSEA（clusterProfiler term2gene）
#
# 入力ファイルは Nextflow rnaseq パイプライン（--run_func_enrich）の出力をそのまま受け付ける:
#   DEG:      Geneid, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
#   ORA:      ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count
#   GSEA:     ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, ..., core_enrichment
#   eggNOG:   *.emapper.annotations（#query, ..., GOs(10), KEGG_Pathway(13)）
# CSV/TSV 自動判別。先頭の # コメント行は無視。

library(GO.db) # (5)のONTOLOGY/TERM取得で未qualifiedのGO.dbオブジェクトを参照するため必須

# ---- 共通ヘルパ ----------------------------------------------------------
.fe_read_table <- function(path) {
  # 先頭の # コメント行を除去してから fread（CSV/TSV自動判別）
  lines <- readLines(path, warn = FALSE)
  data_lines <- lines[!grepl("^#", lines)]
  data_lines <- data_lines[nzchar(trimws(data_lines))]
  if (length(data_lines) == 0) return(NULL)
  tmp <- tempfile()
  writeLines(data_lines, tmp)
  on.exit(unlink(tmp))
  df <- tryCatch(
    as.data.frame(data.table::fread(tmp, sep = "auto", header = TRUE,
                                    check.names = FALSE)),
    error = function(e) NULL
  )
  df
}

.fe_find_col <- function(df, candidates) {
  if (is.null(df)) return(NULL)
  cn <- colnames(df)
  for (c in candidates) {
    hit <- cn[tolower(cn) == tolower(c)]
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

.fe_split_terms <- function(v, g) {
  keep <- !is.na(v) & v != "-" & v != ""
  v <- v[keep]; g <- g[keep]
  if (!length(v)) return(data.frame(term = character(), gene = character()))
  do.call(rbind, Map(function(t, gene) {
    tt <- strsplit(t, ",")[[1]]
    if (!length(tt)) return(NULL)
    data.frame(term = tt, gene = gene, stringsAsFactors = FALSE)
  }, v, g))
}

# ============================ UI =========================================
figureEnrichmentUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h4("📊 アップロード作図 & 非モデル生物エンリッチメント"),
    helpText("Nextflow RNA-seqパイプライン等の中間ファイル（DEG / GO・KEGG / GSEA / カウント）を",
             "アップロードして図を生成します。OrgDbの無い非モデル生物は eggNOG注釈から",
             "GO/KEGGエンリッチメントを直接実行できます。"),
    tabsetPanel(
      id = ns("subtabs"),

      # (1) Volcano / MA -------------------------------------------------
      tabPanel("① Volcano / MA",
        sidebarLayout(
          sidebarPanel(width = 4,
            fileInput(ns("deg_file"), "DEG表 (.tsv/.csv)", accept = c(".tsv",".csv",".txt")),
            radioButtons(ns("deg_plot_type"), "プロット種別:",
                         c("Volcano" = "volcano", "MA" = "ma"), inline = TRUE),
            uiOutput(ns("deg_col_ui")),
            numericInput(ns("deg_padj_th"), "padj 閾値:", value = 0.05, min = 0, step = 0.01),
            numericInput(ns("deg_lfc_th"), "|log2FC| 閾値:", value = 1.0, min = 0, step = 0.5),
            numericInput(ns("deg_label_n"), "ラベル表示する上位遺伝子数:", value = 10, min = 0),
            numericInput(ns("deg_base_size"), "文字サイズ:", value = 14, min = 6),
            downloadButton(ns("dl_deg_pdf"), "PDF保存", icon = icon("file-pdf"))
          ),
          mainPanel(width = 8,
            withSpinner(plotOutput(ns("deg_plot"), height = "560px"), type = 6),
            br(), DTOutput(ns("deg_sig_table"))
          )
        )
      ),

      # (2) ORA dotplot/barplot -----------------------------------------
      tabPanel("② GO/KEGG (ORA)",
        sidebarLayout(
          sidebarPanel(width = 4,
            fileInput(ns("ora_file"), "ORA結果 (.tsv/.csv)", accept = c(".tsv",".csv",".txt")),
            radioButtons(ns("ora_plot_type"), "プロット種別:",
                         c("dotplot" = "dot", "barplot" = "bar"), inline = TRUE),
            numericInput(ns("ora_show_n"), "表示項目数:", value = 20, min = 1),
            selectInput(ns("ora_order"), "並び替え:",
                        c("p.adjust" = "padj", "Count" = "count")),
            numericInput(ns("ora_base_size"), "文字サイズ:", value = 13, min = 6),
            downloadButton(ns("dl_ora_pdf"), "PDF保存", icon = icon("file-pdf"))
          ),
          mainPanel(width = 8,
            withSpinner(plotOutput(ns("ora_plot"), height = "600px"), type = 6),
            br(), DTOutput(ns("ora_table"))
          )
        )
      ),

      # (3) GSEA NES bar -------------------------------------------------
      tabPanel("③ GSEA",
        sidebarLayout(
          sidebarPanel(width = 4,
            fileInput(ns("gsea_file"), "GSEA結果 (.tsv/.csv)", accept = c(".tsv",".csv",".txt")),
            numericInput(ns("gsea_show_n"), "表示項目数（|NES|上位）:", value = 20, min = 1),
            numericInput(ns("gsea_padj_th"), "p.adjust 閾値:", value = 0.05, min = 0, step = 0.01),
            numericInput(ns("gsea_base_size"), "文字サイズ:", value = 13, min = 6),
            downloadButton(ns("dl_gsea_pdf"), "PDF保存", icon = icon("file-pdf"))
          ),
          mainPanel(width = 8,
            withSpinner(plotOutput(ns("gsea_plot"), height = "600px"), type = 6),
            br(), DTOutput(ns("gsea_table"))
          )
        )
      ),

      # (4) Heatmap ------------------------------------------------------
      tabPanel("④ Heatmap",
        sidebarLayout(
          sidebarPanel(width = 4,
            fileInput(ns("hm_file"), "カウント/発現行列 (.csv/.tsv) 1列目=遺伝子ID",
                      accept = c(".tsv",".csv",".txt")),
            textAreaInput(ns("hm_genes"), "対象遺伝子（空欄=分散上位を自動選択）:",
                          placeholder = "1行1遺伝子, または カンマ区切り", height = "100px"),
            numericInput(ns("hm_top_n"), "自動選択時の分散上位数:", value = 50, min = 2),
            checkboxInput(ns("hm_scale"), "行方向にZスコア化 (scale)", value = TRUE),
            checkboxInput(ns("hm_log"), "log2(x+1) 変換", value = TRUE),
            numericInput(ns("hm_base_size"), "文字サイズ:", value = 10, min = 5),
            downloadButton(ns("dl_hm_pdf"), "PDF保存", icon = icon("file-pdf"))
          ),
          mainPanel(width = 8,
            withSpinner(plotOutput(ns("hm_plot"), height = "640px"), type = 6)
          )
        )
      ),

      # (5) 非モデル生物 eggNOG エンリッチメント -------------------------
      tabPanel("⑤ 非モデル生物 (eggNOG)",
        sidebarLayout(
          sidebarPanel(width = 4,
            helpText("OrgDbの無い生物向け。DEG表とeggNOG-mapper注釈から",
                     "GO/KEGGのORA+GSEAをclusterProfilerで実行します。"),
            fileInput(ns("egg_deg"), "DEG表 (.tsv/.csv)", accept = c(".tsv",".csv",".txt")),
            fileInput(ns("egg_ann"), "eggNOG注釈 (*.emapper.annotations)",
                      accept = c(".annotations",".tsv",".txt")),
            fileInput(ns("egg_q2g"), "query2gene.tsv（任意: query≠遺伝子IDの場合）",
                      accept = c(".tsv",".txt")),
            selectInput(ns("egg_db"), "対象:",
                        c("GO Biological Process" = "BP", "GO Molecular Function" = "MF",
                          "GO Cellular Component" = "CC", "KEGG Pathway" = "KEGG")),
            radioButtons(ns("egg_method"), "解析:",
                         c("ORA (DEG過剰代表)" = "ora", "GSEA (ランキング)" = "gsea"), inline = TRUE),
            numericInput(ns("egg_padj_th"), "DEG: padj 閾値:", value = 0.05, min = 0, step = 0.01),
            numericInput(ns("egg_lfc_th"), "DEG: |log2FC| 閾値:", value = 1.0, min = 0, step = 0.5),
            numericInput(ns("egg_show_n"), "プロット表示項目数:", value = 20, min = 1),
            actionButton(ns("egg_run"), "実行", icon = icon("play")),
            hr(),
            downloadButton(ns("dl_egg_tsv"), "結果テーブル (.tsv)", icon = icon("download")),
            downloadButton(ns("dl_egg_pdf"), "プロットPDF", icon = icon("file-pdf"))
          ),
          mainPanel(width = 8,
            withSpinner(plotOutput(ns("egg_plot"), height = "600px"), type = 6),
            br(), withSpinner(DTOutput(ns("egg_table")), type = 6)
          )
        )
      )
    )
  )
}

# ============================ Server =====================================
figureEnrichmentServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---------- (1) DEG: Volcano / MA -------------------------------------
    deg_data <- reactive({
      req(input$deg_file)
      df <- .fe_read_table(input$deg_file$datapath)
      validate(need(!is.null(df) && ncol(df) >= 2, "DEG表を読み込めませんでした。"))
      df
    })

    output$deg_col_ui <- renderUI({
      df <- deg_data(); cn <- colnames(df)
      id_def   <- .fe_find_col(df, c("Geneid","gene","ID","GeneSymbol")) %||% cn[1]
      lfc_def  <- .fe_find_col(df, c("log2FoldChange","logFC","log2FC")) %||% cn[1]
      padj_def <- .fe_find_col(df, c("padj","FDR","adj.P.Val","p.adjust","qvalue")) %||% cn[1]
      base_def <- .fe_find_col(df, c("baseMean","logCPM","AveExpr","basemean"))
      tagList(
        selectInput(ns("deg_id_col"),   "遺伝子ID列:", cn, selected = id_def),
        selectInput(ns("deg_lfc_col"),  "log2FC列:",  cn, selected = lfc_def),
        selectInput(ns("deg_padj_col"), "padj列:",    cn, selected = padj_def),
        if (input$deg_plot_type == "ma")
          selectInput(ns("deg_base_col"), "発現量列(MA用):", cn,
                      selected = base_def %||% cn[1])
      )
    })

    deg_plot_obj <- reactive({
      df <- deg_data()
      req(input$deg_id_col, input$deg_lfc_col, input$deg_padj_col)
      d <- data.frame(
        id   = as.character(df[[input$deg_id_col]]),
        lfc  = suppressWarnings(as.numeric(df[[input$deg_lfc_col]])),
        padj = suppressWarnings(as.numeric(df[[input$deg_padj_col]])),
        stringsAsFactors = FALSE
      )
      d <- d[is.finite(d$lfc) & !is.na(d$padj), ]
      validate(need(nrow(d) > 0, "有効な行がありません（列指定を確認）。"))
      d$sig <- ifelse(d$padj < input$deg_padj_th & abs(d$lfc) > input$deg_lfc_th,
                      ifelse(d$lfc > 0, "Up", "Down"), "NS")
      d$sig <- factor(d$sig, levels = c("Up","Down","NS"))
      cols <- c(Up = "#e41a1c", Down = "#377eb8", NS = "grey70")
      if (input$deg_plot_type == "volcano") {
        d$neglog <- -log10(pmax(d$padj, .Machine$double.xmin))
        p <- ggplot(d, aes(lfc, neglog, color = sig)) +
          geom_point(alpha = 0.6, size = 1.4) +
          geom_vline(xintercept = c(-1,1) * input$deg_lfc_th, linetype = 2, color = "grey50") +
          geom_hline(yintercept = -log10(input$deg_padj_th), linetype = 2, color = "grey50") +
          labs(x = "log2 Fold Change", y = "-log10(padj)", color = "")
      } else {
        req(input$deg_base_col)
        d$expr <- suppressWarnings(as.numeric(df[[input$deg_base_col]][match(d$id, as.character(df[[input$deg_id_col]]))]))
        d <- d[is.finite(d$expr), ]
        d$xexpr <- log10(pmax(d$expr, 1))
        p <- ggplot(d, aes(xexpr, lfc, color = sig)) +
          geom_point(alpha = 0.6, size = 1.4) +
          geom_hline(yintercept = c(-1,1) * input$deg_lfc_th, linetype = 2, color = "grey50") +
          labs(x = "log10(mean expression)", y = "log2 Fold Change", color = "")
      }
      p <- p + scale_color_manual(values = cols) + theme_bw(base_size = input$deg_base_size)
      if (input$deg_label_n > 0) {
        top <- d[order(d$padj), ]; top <- head(top[top$sig != "NS", ], input$deg_label_n)
        if (nrow(top) > 0) {
          yv <- if (input$deg_plot_type == "volcano") "neglog" else "lfc"
          xv <- if (input$deg_plot_type == "volcano") "lfc" else "xexpr"
          p <- p + ggrepel::geom_text_repel(
            data = top, aes_string(x = xv, y = yv, label = "id"),
            size = input$deg_base_size * 0.22, max.overlaps = 50, show.legend = FALSE)
        }
      }
      p
    })

    output$deg_plot <- renderPlot({ print(deg_plot_obj()) })
    output$deg_sig_table <- renderDT({
      df <- deg_data()
      req(input$deg_id_col, input$deg_lfc_col, input$deg_padj_col)
      d <- df[, c(input$deg_id_col, input$deg_lfc_col, input$deg_padj_col)]
      colnames(d) <- c("Geneid","log2FC","padj")
      d <- d[is.finite(suppressWarnings(as.numeric(d$log2FC))) & !is.na(as.numeric(d$padj)), ]
      d <- d[as.numeric(d$padj) < input$deg_padj_th & abs(as.numeric(d$log2FC)) > input$deg_lfc_th, ]
      datatable(d[order(as.numeric(d$padj)), ], rownames = FALSE,
                options = list(pageLength = 10, scrollX = TRUE))
    })
    output$dl_deg_pdf <- downloadHandler(
      filename = function() paste0("DEG_", input$deg_plot_type, "_", Sys.Date(), ".pdf"),
      content = function(file) { pdf(file, width = 8, height = 7); print(deg_plot_obj()); dev.off() }
    )

    # ---------- (2) ORA dotplot / barplot --------------------------------
    ora_data <- reactive({
      req(input$ora_file)
      df <- .fe_read_table(input$ora_file$datapath)
      validate(need(!is.null(df), "ORA結果を読み込めませんでした。"))
      idc <- .fe_find_col(df, c("ID","term")); dsc <- .fe_find_col(df, c("Description","name"))
      pac <- .fe_find_col(df, c("p.adjust","padj","qvalue","FDR"))
      cnt <- .fe_find_col(df, c("Count","count"))
      gr  <- .fe_find_col(df, c("GeneRatio"))
      validate(need(!is.null(pac), "p.adjust列が見つかりません。"))
      out <- data.frame(
        ID = if (!is.null(idc)) as.character(df[[idc]]) else seq_len(nrow(df)),
        Description = if (!is.null(dsc)) as.character(df[[dsc]]) else
                      if (!is.null(idc)) as.character(df[[idc]]) else "",
        p.adjust = suppressWarnings(as.numeric(df[[pac]])),
        Count = if (!is.null(cnt)) suppressWarnings(as.numeric(df[[cnt]])) else NA_real_,
        stringsAsFactors = FALSE
      )
      if (!is.null(gr)) {
        rr <- sapply(strsplit(as.character(df[[gr]]), "/"),
                     function(x) { x <- as.numeric(x); if (length(x)==2 && x[2]>0) x[1]/x[2] else NA })
        out$GeneRatio <- rr
        if (all(is.na(out$Count))) out$Count <- sapply(strsplit(as.character(df[[gr]]), "/"),
                                                        function(x) as.numeric(x[1]))
      }
      out[is.finite(out$p.adjust), ]
    })

    ora_plot_obj <- reactive({
      d <- ora_data()
      ord <- if (input$ora_order == "count" && !all(is.na(d$Count))) order(-d$Count) else order(d$p.adjust)
      d <- d[ord, ]; d <- head(d, input$ora_show_n)
      d$Description <- factor(d$Description, levels = rev(d$Description))
      xval <- if (!is.null(d$GeneRatio) && any(is.finite(d$GeneRatio))) "GeneRatio" else "Count"
      if (input$ora_plot_type == "dot") {
        p <- ggplot(d, aes_string(x = xval, y = "Description",
                                  size = "Count", color = "p.adjust")) +
          geom_point() + scale_color_gradient(low = "#e41a1c", high = "#377eb8") +
          labs(y = NULL)
      } else {
        p <- ggplot(d, aes_string(x = "Count", y = "Description", fill = "p.adjust")) +
          geom_col() + scale_fill_gradient(low = "#e41a1c", high = "#377eb8") +
          labs(y = NULL)
      }
      p + theme_bw(base_size = input$ora_base_size)
    })
    output$ora_plot <- renderPlot({ print(ora_plot_obj()) })
    output$ora_table <- renderDT({ datatable(ora_data(), rownames = FALSE,
                                              options = list(pageLength = 10, scrollX = TRUE)) })
    output$dl_ora_pdf <- downloadHandler(
      filename = function() paste0("ORA_", input$ora_plot_type, "_", Sys.Date(), ".pdf"),
      content = function(file) { pdf(file, width = 9, height = 7); print(ora_plot_obj()); dev.off() }
    )

    # ---------- (3) GSEA NES bar -----------------------------------------
    gsea_data <- reactive({
      req(input$gsea_file)
      df <- .fe_read_table(input$gsea_file$datapath)
      validate(need(!is.null(df), "GSEA結果を読み込めませんでした。"))
      idc <- .fe_find_col(df, c("ID","pathway","term"))
      dsc <- .fe_find_col(df, c("Description","name"))
      nesc<- .fe_find_col(df, c("NES")); pac <- .fe_find_col(df, c("p.adjust","padj","qvalue"))
      validate(need(!is.null(nesc) && !is.null(pac), "NES列またはp.adjust列が見つかりません。"))
      data.frame(
        Description = if (!is.null(dsc)) as.character(df[[dsc]]) else as.character(df[[idc]]),
        NES = suppressWarnings(as.numeric(df[[nesc]])),
        p.adjust = suppressWarnings(as.numeric(df[[pac]])),
        stringsAsFactors = FALSE
      )
    })
    gsea_plot_obj <- reactive({
      d <- gsea_data()
      d <- d[is.finite(d$NES) & d$p.adjust < input$gsea_padj_th, ]
      validate(need(nrow(d) > 0, "閾値を満たす項目がありません。"))
      d <- d[order(-abs(d$NES)), ]; d <- head(d, input$gsea_show_n)
      d$dir <- ifelse(d$NES > 0, "Activated", "Suppressed")
      d$Description <- factor(d$Description, levels = d$Description[order(d$NES)])
      ggplot(d, aes(x = NES, y = Description, fill = dir)) +
        geom_col() +
        scale_fill_manual(values = c(Activated = "#e41a1c", Suppressed = "#377eb8")) +
        labs(y = NULL, fill = "") + theme_bw(base_size = input$gsea_base_size)
    })
    output$gsea_plot <- renderPlot({ print(gsea_plot_obj()) })
    output$gsea_table <- renderDT({ datatable(gsea_data(), rownames = FALSE,
                                               options = list(pageLength = 10, scrollX = TRUE)) })
    output$dl_gsea_pdf <- downloadHandler(
      filename = function() paste0("GSEA_NESbar_", Sys.Date(), ".pdf"),
      content = function(file) { pdf(file, width = 9, height = 7); print(gsea_plot_obj()); dev.off() }
    )

    # ---------- (4) Heatmap ----------------------------------------------
    hm_matrix <- reactive({
      req(input$hm_file)
      df <- .fe_read_table(input$hm_file$datapath)
      validate(need(!is.null(df) && ncol(df) >= 3, "発現行列を読み込めませんでした。"))
      rn <- as.character(df[[1]]); m <- df[, -1, drop = FALSE]
      num <- vapply(m, function(x) is.numeric(x) || !any(is.na(suppressWarnings(as.numeric(x)))), logical(1))
      m <- m[, num, drop = FALSE]
      m <- as.matrix(sapply(m, function(x) suppressWarnings(as.numeric(x))))
      rownames(m) <- rn
      m[is.finite(rowSums(m)), , drop = FALSE]
    })
    hm_plot_fun <- reactive({
      m <- hm_matrix()
      genes_in <- trimws(unlist(strsplit(input$hm_genes, "[,\n\r]+")))
      genes_in <- genes_in[nzchar(genes_in)]
      if (length(genes_in) > 0) {
        sel <- intersect(genes_in, rownames(m))
        validate(need(length(sel) >= 2, "指定遺伝子が行列内に2つ以上見つかりません。"))
        m <- m[sel, , drop = FALSE]
      } else {
        v <- matrixStats::rowVars(m)
        m <- m[order(-v)[seq_len(min(input$hm_top_n, nrow(m)))], , drop = FALSE]
      }
      if (input$hm_log) m <- log2(m + 1)
      scale_arg <- if (input$hm_scale) "row" else "none"
      function() pheatmap::pheatmap(m, scale = scale_arg, fontsize = input$hm_base_size,
                                    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                                    silent = FALSE)
    })
    output$hm_plot <- renderPlot({ hm_plot_fun()() })
    output$dl_hm_pdf <- downloadHandler(
      filename = function() paste0("Heatmap_", Sys.Date(), ".pdf"),
      content = function(file) { pdf(file, width = 8, height = 9); print(hm_plot_fun()()); dev.off() }
    )

    # ---------- (5) 非モデル生物 eggNOG エンリッチメント -----------------
    egg_results <- eventReactive(input$egg_run, {
      req(input$egg_deg, input$egg_ann)
      deg <- .fe_read_table(input$egg_deg$datapath)
      validate(need(!is.null(deg), "DEG表を読み込めませんでした。"))
      idc  <- .fe_find_col(deg, c("Geneid","gene","ID"))
      lfcc <- .fe_find_col(deg, c("log2FoldChange","logFC","log2FC"))
      pac  <- .fe_find_col(deg, c("padj","FDR","p.adjust"))
      stc  <- .fe_find_col(deg, c("stat","t","logFC_x_signed_p"))
      validate(need(!is.null(idc), "DEG表に遺伝子ID列がありません。"))

      # # で始まる行（##コメント + #queryヘッダ）を全て除去し、列名は手動付与
      em <- read.delim(input$egg_ann$datapath, comment.char = "#",
                       header = FALSE, quote = "", stringsAsFactors = FALSE)
      # 列名（emapper v2標準）
      cn <- c("query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl",
              "COG_category","Description","Preferred_name","GOs","EC","KEGG_ko",
              "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE",
              "KEGG_TC","CAZy","BiGG_Reaction","PFAMs")
      colnames(em) <- cn[seq_len(ncol(em))]

      # query -> gene マッピング
      if (!is.null(input$egg_q2g)) {
        q2g <- as.data.frame(data.table::fread(input$egg_q2g$datapath, header = FALSE))
        colnames(q2g)[1:2] <- c("query","gene")
        em <- merge(em, q2g[,1:2], by = "query")
      } else {
        em$gene <- em$query
      }

      # term2gene 構築
      if (input$egg_db == "KEGG") {
        t2g <- unique(.fe_split_terms(em$KEGG_Pathway, em$gene))
        t2g <- t2g[grepl("^ko[0-9]{5}$", t2g$term), ]
        t2n <- tryCatch({
          km <- clusterProfiler:::kegg_list("pathway", "ko")
          data.frame(term = sub("^path:", "", km$from), name = km$to)
        }, error = function(e) NULL)
      } else {
        t2g <- unique(.fe_split_terms(em$GOs, em$gene))
        t2g <- t2g[grepl("^GO:", t2g$term), ]
        ids <- unique(t2g$term)
        gn <- suppressMessages(AnnotationDbi::select(GO.db, keys = ids,
                  columns = c("TERM","ONTOLOGY"), keytype = "GOID"))
        gn <- gn[gn$ONTOLOGY == input$egg_db, ]
        t2g <- t2g[t2g$term %in% gn$GOID, ]
        t2n <- data.frame(term = gn$GOID, name = gn$TERM)
      }
      validate(need(nrow(t2g) > 0, "注釈から有効な term2gene を構築できませんでした。"))
      universe <- unique(em$gene)

      if (input$egg_method == "ora") {
        sig <- deg[!is.na(suppressWarnings(as.numeric(deg[[pac]]))) &
                   as.numeric(deg[[pac]]) < input$egg_padj_th &
                   abs(as.numeric(deg[[lfcc]])) > input$egg_lfc_th, ]
        gene <- intersect(as.character(sig[[idc]]), universe)
        validate(need(length(gene) > 0, "閾値を満たすDEGが背景に含まれません。"))
        res <- clusterProfiler::enricher(gene, TERM2GENE = t2g, TERM2NAME = t2n,
                  universe = universe, pvalueCutoff = 1, qvalueCutoff = 1)
      } else {
        rk_col <- stc %||% lfcc
        validate(need(!is.null(rk_col), "ランキングに使える列(stat/logFC)がありません。"))
        d2 <- deg[!is.na(suppressWarnings(as.numeric(deg[[rk_col]]))), ]
        gl <- setNames(as.numeric(d2[[rk_col]]), as.character(d2[[idc]]))
        gl <- sort(gl[names(gl) %in% universe], decreasing = TRUE)
        set.seed(1234)
        res <- clusterProfiler::GSEA(gl, TERM2GENE = t2g, TERM2NAME = t2n,
                  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, eps = 0, seed = TRUE)
      }
      validate(need(!is.null(res) && nrow(as.data.frame(res)) > 0, "有意な結果が得られませんでした。"))
      list(obj = res, method = input$egg_method)
    })

    output$egg_table <- renderDT({
      r <- egg_results(); datatable(as.data.frame(r$obj), rownames = FALSE,
                                    options = list(pageLength = 10, scrollX = TRUE))
    })
    egg_plot_obj <- reactive({
      r <- egg_results(); n <- input$egg_show_n
      if (r$method == "ora") {
        enrichplot::dotplot(r$obj, showCategory = n) + ggtitle("eggNOG ORA")
      } else {
        df <- as.data.frame(r$obj)
        df <- df[order(-abs(df$NES)), ]; df <- head(df, n)
        df$Description <- factor(df$Description, levels = df$Description[order(df$NES)])
        df$dir <- ifelse(df$NES > 0, "Activated", "Suppressed")
        ggplot(df, aes(x = NES, y = Description, fill = dir)) + geom_col() +
          scale_fill_manual(values = c(Activated = "#e41a1c", Suppressed = "#377eb8")) +
          labs(y = NULL, fill = "") + theme_bw(base_size = 13) + ggtitle("eggNOG GSEA")
      }
    })
    output$egg_plot <- renderPlot({ print(egg_plot_obj()) })
    output$dl_egg_tsv <- downloadHandler(
      filename = function() paste0("eggNOG_", input$egg_db, "_", input$egg_method, "_", Sys.Date(), ".tsv"),
      content = function(file) {
        r <- egg_results(); write.table(as.data.frame(r$obj), file, sep = "\t",
                                         quote = FALSE, row.names = FALSE)
      }
    )
    output$dl_egg_pdf <- downloadHandler(
      filename = function() paste0("eggNOG_", input$egg_db, "_", input$egg_method, "_", Sys.Date(), ".pdf"),
      content = function(file) { pdf(file, width = 9, height = 7); print(egg_plot_obj()); dev.off() }
    )

    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
  })
}
