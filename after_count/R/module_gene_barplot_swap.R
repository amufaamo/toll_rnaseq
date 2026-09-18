# R/module_gene_barplot_swap.R
# Tab 9: 鍵遺伝子 Barplot + サンプルスワップ検証 + AI解釈エクスポート
#
# 追加背景: 星野先生 opal RNA-seq で判明した2ニーズ (2026-04-10)
#   1. バッチ別・条件別の鍵遺伝子 barplot
#   2. サンプルスワップ疑いの相関検証
#   3. 解析結果をAIに送れるプロンプト生成

library(shiny)
library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(shinycssloaders)
library(DT)
library(AnnotationDbi)
library(plotly)

# ============================================================
# UI
# ============================================================

geneBarplotSwapUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h4("Tab 9: 鍵遺伝子 Barplot・スワップ検証・AI解釈"),
    helpText(
      "特定遺伝子の発現量をバッチ・条件別に比較します。",
      "サンプルスワップ（実験ミスによるラベル間違い）の検証と、",
      "解析結果をAIに送るためのプロンプト生成も行えます。"
    ),
    hr(),

    tabsetPanel(

      # -------- Tab A: 鍵遺伝子 barplot --------
      tabPanel("鍵遺伝子 Barplot",
        br(),
        fluidRow(
          column(5,
            textAreaInput(ns("gene_list"),
              label = "遺伝子名（1行1遺伝子、Symbol または Ensembl ID）:",
              value = "Txnip\nIfna4\nIfnb1",
              rows  = 6,
              placeholder = "Txnip\nIfna4\n..."
            )
          ),
          column(4,
            selectInput(ns("color_by"), "色分け:",
                        choices = c("batch", "group", "stim"), selected = "batch"),
            checkboxInput(ns("free_y"), "Y軸をスケール別々にする", value = TRUE),
            numericInput(ns("ncol"), "列数:", value = 2, min = 1, max = 4),
            checkboxInput(ns("show_logcpm"), "log2CPMで表示（デフォルトはCPM）", value = FALSE)
          ),
          column(3,
            br(),
            actionButton(ns("run_barplot"), "Barplot を描画",
                         class = "btn-primary btn-lg", width = "100%"),
            br(), br(),
            helpText("遺伝子名が見つからない場合は、",
                     "Ensembl IDで試してください。",
                     "例: ENSMUSG00000038393 (Txnip)")
          )
        ),
        withSpinner(plotOutput(ns("barplot"), height = "550px")),
        fluidRow(
          column(3, downloadButton(ns("dl_bar_pdf"), "PDF")),
          column(3, downloadButton(ns("dl_bar_png"), "PNG"))
        ),
        hr(),
        h5("D19刺激応答（log2FC: D19 / medium）バッチ比較"),
        helpText(
          "D19あり/なしの比をバッチごとに比較します。",
          "Batch Iの応答が弱い場合は除外を検討してください。",
          "※メタデータに 'stim' 列（値: D19, med）が必要です。"
        ),
        withSpinner(plotOutput(ns("response_plot"), height = "350px")),
        fluidRow(
          column(3, downloadButton(ns("dl_resp_pdf"), "PDF")),
          column(3, downloadButton(ns("dl_resp_csv"), "CSV"))
        ),
        br(),
        withSpinner(DTOutput(ns("response_table")))
      ),

      # -------- Tab B: スワップ検証 --------
      tabPanel("サンプルスワップ検証",
        br(),
        helpText(
          "実験ミスによるサンプルラベルの入れ替わりを確認します。",
          "疑いのあるサンプルを選択し、他サンプルとの相関を比較してください。",
          "スワップがある場合、本来の条件より別条件との相関が高くなります。"
        ),
        fluidRow(
          column(3,
            selectInput(ns("swap_target"),
                        "疑わしいサンプル:",
                        choices = NULL)
          ),
          column(3,
            selectInput(ns("expected_group"),
                        "本来あるべき条件グループ:",
                        choices = NULL)
          ),
          column(3,
            selectInput(ns("suspected_group"),
                        "スワップ先の疑いがある条件グループ:",
                        choices = NULL)
          ),
          column(3,
            br(),
            actionButton(ns("run_swap"), "検証を実行",
                         class = "btn-warning btn-lg", width = "100%")
          )
        ),
        withSpinner(plotOutput(ns("swap_barplot"), height = "380px")),
        verbatimTextOutput(ns("swap_result_text")),
        fluidRow(
          column(3, downloadButton(ns("dl_swap_pdf"), "barplot PDF")),
          column(3, downloadButton(ns("dl_swap_csv"), "相関テーブル CSV"))
        ),
        hr(),
        h5("全サンプル 相関ヒートマップ"),
        withSpinner(plotOutput(ns("cor_heatmap"), height = "580px")),
        fluidRow(
          column(3, downloadButton(ns("dl_heatmap_pdf"), "ヒートマップ PDF"))
        ),
        br(),
        withSpinner(DTOutput(ns("swap_table")))
      ),

      # -------- Tab C: AI解釈エクスポート --------
      tabPanel("🤖 AIに解釈を依頼",
        br(),
        helpText(
          "解析結果を Claude などのAIに送るためのプロンプトを自動生成します。",
          "プロンプトをコピーして、AIのチャット画面に貼り付けてください。"
        ),
        hr(),
        h5("1. 解析コンテキストを入力してください"),
        fluidRow(
          column(6,
            textInput(ns("ai_experiment"),
                      "実験の概要:",
                      value = "グルコース/アロース条件下でのD19刺激に対するマウス免疫細胞の転写応答"),
            textAreaInput(ns("ai_question"),
                          "AIへの質問（自由記述）:",
                          value = "以下の結果を生物学的に解釈し、論文の考察に使える説明を日本語で書いてください。",
                          rows = 3)
          ),
          column(6,
            checkboxGroupInput(ns("ai_include"),
                               "プロンプトに含める内容:",
                               choices  = c(
                                 "サンプル情報・実験デザイン" = "design",
                                 "スワップ検証結果"          = "swap",
                                 "D19応答（鍵遺伝子）"       = "response"
                               ),
                               selected = c("design", "swap", "response"))
          )
        ),
        fluidRow(
          column(3,
            actionButton(ns("gen_prompt"), "プロンプト生成",
                         class = "btn-success btn-lg", width = "100%")
          )
        ),
        br(),
        h5("2. 生成されたプロンプト（コピーしてAIに貼り付けてください）"),
        withSpinner(verbatimTextOutput(ns("ai_prompt"))),
        fluidRow(
          column(4,
            downloadButton(ns("dl_prompt_txt"), "プロンプトをテキストファイルに保存")
          )
        ),
        br(),
        div(
          style = "background-color: #f0f8ff; padding: 12px; border-radius: 6px;",
          h6("使い方のヒント"),
          tags$ul(
            tags$li("Claude Code や ChatGPT に貼り付けてください"),
            tags$li("日本語での返答を希望する場合は質問欄に「日本語で」と記入してください"),
            tags$li("図のPDFをAIに送りたい場合は各タブのダウンロードボタンを使ってください")
          )
        )
      )
    )
  )
}

# ============================================================
# Server
# ============================================================

geneBarplotSwapServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- 共通: 正規化データを準備（dimension_reductionと同じパターン）----
    prepared_data <- reactive({
      req(rv$merged_data, rv$sample_metadata, rv$filtered_keep)

      active_meta <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      validate(need(nrow(active_meta) > 0,
                    "解析対象のサンプルが選択されていません。Tab 1 でサンプルを確認してください。"))

      active_names <- active_meta$current_name
      counts_df    <- rv$merged_data
      keep_vec     <- rv$filtered_keep

      validate(need(length(keep_vec) == nrow(counts_df),
                    "フィルタリング情報が古いです。Tab 2（フィルタリング）を再実行してください。"))

      counts_active <- counts_df[, c("Geneid", intersect(colnames(counts_df), active_names)), drop=FALSE]
      counts_filt   <- counts_active[keep_vec, ]

      count_mat <- counts_filt[, -1, drop=FALSE]
      rownames(count_mat) <- counts_filt$Geneid

      sample_order <- colnames(count_mat)
      meta_ord     <- active_meta[match(sample_order, active_meta$current_name), ]
      group        <- factor(meta_ord$group)

      y      <- DGEList(counts = count_mat, group = group)
      y      <- calcNormFactors(y)
      logcpm <- edgeR::cpm(y, log = TRUE,  prior.count = 2)
      cpm    <- edgeR::cpm(y, log = FALSE, prior.count = 1)

      list(logcpm = logcpm, cpm = cpm, meta = meta_ord, group = group)
    })

    # ---- セレクトボックスの更新 ----
    observe({
      req(rv$sample_metadata)
      meta    <- rv$sample_metadata[rv$sample_metadata$active, , drop=FALSE]
      samples <- meta$current_name
      groups  <- unique(meta$group)

      updateSelectInput(session, "swap_target",   choices = samples, selected = samples[1])
      updateSelectInput(session, "expected_group",  choices = groups,  selected = groups[1])
      updateSelectInput(session, "suspected_group", choices = groups,
                        selected = if (length(groups) > 1) groups[2] else groups[1])

      # color_by の選択肢をメタデータの列名に合わせる
      meta_cols <- setdiff(colnames(meta), c("id","current_name","active","time"))
      updateSelectInput(session, "color_by", choices = meta_cols, selected = "batch")
    })

    # ---- 遺伝子名 → Ensembl ID 解決 ----
    resolve_genes <- reactive({
      req(prepared_data(), input$gene_list)
      mat   <- prepared_data()$logcpm
      genes <- trimws(strsplit(input$gene_list, "\n")[[1]])
      genes <- genes[nchar(genes) > 0]

      ensembl_re <- "^ENS[A-Z]+G[0-9]+"
      found <- list()
      for (g in genes) {
        if (grepl(ensembl_re, g) && g %in% rownames(mat)) {
          found[[g]] <- g
        } else {
          # Symbol → Ensembl 変換を試みる（OrgDbが利用可能な場合）
          if (!is.null(rv$selected_species)) {
            orgdb_name <- c(
              "Homo_sapiens" = "org.Hs.eg.db", "Mus_musculus" = "org.Mm.eg.db",
              "Rattus_norvegicus" = "org.Rn.eg.db"
            )[[rv$selected_species]]
            if (!is.null(orgdb_name) && requireNamespace(orgdb_name, quietly=TRUE)) {
              orgdb <- get(orgdb_name)
              ens_ids <- tryCatch(
                suppressMessages(mapIds(orgdb, keys=g, column="ENSEMBL",
                                        keytype="SYMBOL", multiVals="first")),
                error = function(e) NA_character_
              )
              ens_ids <- ens_ids[!is.na(ens_ids) & ens_ids %in% rownames(mat)]
              if (length(ens_ids) > 0) {
                found[[g]] <- ens_ids[1]
                next
              }
            }
          }
          # フォールバック: symbolそのままで検索
          if (g %in% rownames(mat)) found[[g]] <- g
        }
      }
      found
    })

    # ---- Barplot 描画 ----
    barplot_react <- eventReactive(input$run_barplot, {
      req(prepared_data(), resolve_genes())
      pd   <- prepared_data()
      mat  <- if (input$show_logcpm) pd$logcpm else pd$cpm
      gmap <- resolve_genes()
      validate(need(length(gmap) > 0,
                    "指定した遺伝子がデータに見つかりませんでした。遺伝子名またはEnsembl IDを確認してください。"))

      mat_sub <- mat[unlist(gmap), , drop=FALSE]
      df <- mat_sub %>%
        as.data.frame() %>%
        rownames_to_column("ensembl") %>%
        mutate(gene = names(gmap)[match(ensembl, unlist(gmap))]) %>%
        pivot_longer(-c(ensembl, gene), names_to="sample", values_to="value") %>%
        left_join(pd$meta %>% rename(sample=current_name), by="sample") %>%
        mutate(sample = factor(sample, levels=pd$meta$current_name))

      color_var <- input$color_by
      y_lab <- if (input$show_logcpm) "log2CPM" else "CPM (normalized)"
      scale_y  <- if (input$free_y) "free_y" else "fixed"

      p <- ggplot(df, aes(x=sample, y=value, fill=.data[[color_var]])) +
        geom_bar(stat="identity") +
        facet_wrap(~gene, scales=scale_y, ncol=input$ncol) +
        theme_bw(base_size=11) +
        theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
              strip.text=element_text(face="bold"),
              legend.position="top") +
        labs(title="Key Gene Expression by Sample",
             subtitle=paste0("edgeR TMM-normalized ", y_lab),
             x=NULL, y=y_lab, fill=color_var)
      p
    })

    output$barplot <- renderPlot({ req(barplot_react()); barplot_react() })

    # ---- D19応答 barplot ----
    response_data_react <- eventReactive(input$run_barplot, {
      req(prepared_data(), resolve_genes())
      pd   <- prepared_data()
      mat  <- pd$cpm
      gmap <- resolve_genes()
      if (length(gmap)==0) return(NULL)

      meta <- pd$meta
      if (!all(c("stim","batch") %in% colnames(meta))) return(NULL)

      mat_sub <- mat[unlist(gmap), , drop=FALSE]
      df <- mat_sub %>%
        as.data.frame() %>%
        rownames_to_column("ensembl") %>%
        mutate(gene=names(gmap)[match(ensembl, unlist(gmap))]) %>%
        pivot_longer(-c(ensembl, gene), names_to="sample", values_to="cpm") %>%
        left_join(meta %>% rename(sample=current_name), by="sample")

      group_vars <- c("gene","batch","stim")
      if ("glucose" %in% colnames(meta)) group_vars <- c(group_vars,"glucose")

      resp <- df %>%
        group_by(across(all_of(group_vars))) %>%
        summarise(mean_cpm=mean(cpm, na.rm=TRUE), .groups="drop") %>%
        pivot_wider(names_from=stim, values_from=mean_cpm) %>%
        mutate(
          log2FC  = log2((.data[["D19"]]+0.1) / (.data[["med"]]+0.1)),
          x_label = if ("glucose" %in% colnames(.)) paste0(batch,"_",glucose) else as.character(batch)
        )
      resp
    })

    output$response_plot <- renderPlot({
      req(response_data_react())
      resp <- response_data_react()
      ggplot(resp, aes(x=x_label, y=log2FC, fill=batch)) +
        geom_bar(stat="identity") +
        geom_hline(yintercept=0, linetype="dashed") +
        facet_wrap(~gene, scales="free_y") +
        theme_bw(base_size=11) +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        labs(title="D19 Stimulation Response per Batch",
             x="Batch", y="log2FC (D19 / medium)", fill="Batch")
    })

    output$response_table <- renderDT({
      req(response_data_react())
      datatable(response_data_react() %>% mutate(across(where(is.numeric), ~round(.x,3))),
                options=list(pageLength=10, scrollX=TRUE), rownames=FALSE)
    })

    # ---- スワップ検証 ----
    swap_result <- eventReactive(input$run_swap, {
      req(prepared_data(), input$swap_target, input$expected_group, input$suspected_group)
      pd     <- prepared_data()
      mat    <- pd$logcpm
      target <- input$swap_target
      validate(need(target %in% colnames(mat), paste("サンプル", target, "が見つかりません。")))

      cor_mat <- cor(mat, method="pearson")
      cor_vec <- sort(cor_mat[target, ], decreasing=TRUE)

      meta <- pd$meta %>% rename(sample=current_name)
      group_cors <- data.frame(sample=names(cor_vec), correlation=as.numeric(cor_vec)) %>%
        filter(sample != target) %>%
        left_join(meta, by="sample") %>%
        arrange(desc(correlation))

      exp_cor  <- mean(group_cors$correlation[group_cors$group==input$expected_group],  na.rm=TRUE)
      susp_cor <- mean(group_cors$correlation[group_cors$group==input$suspected_group], na.rm=TRUE)

      list(target=target, cor_mat=cor_mat, group_cors=group_cors,
           exp_cor=exp_cor, susp_cor=susp_cor,
           exp_group=input$expected_group, susp_group=input$suspected_group,
           is_swap = susp_cor > exp_cor, meta=pd$meta)
    })

    output$swap_barplot <- renderPlot({
      req(swap_result())
      res <- swap_result()
      df  <- res$group_cors %>%
        mutate(
          highlight = case_when(
            group == res$exp_group  ~ paste0(res$exp_group,  " (expected)"),
            group == res$susp_group ~ paste0(res$susp_group, " (swap hypothesis)"),
            TRUE ~ "Other"),
          sample = factor(sample, levels=rev(res$group_cors$sample))
        )
      colors <- setNames(c("#D55E00","#0072B2","#999999"),
                         c(paste0(res$exp_group," (expected)"),
                           paste0(res$susp_group," (swap hypothesis)"),
                           "Other"))
      ggplot(df, aes(x=sample, y=correlation, fill=highlight)) +
        geom_bar(stat="identity") + coord_flip() +
        geom_hline(yintercept=res$exp_cor,  linetype="dashed", color="#D55E00", linewidth=1) +
        geom_hline(yintercept=res$susp_cor, linetype="dotted", color="#0072B2", linewidth=1) +
        scale_fill_manual(values=colors, drop=FALSE) +
        theme_bw(base_size=12) +
        labs(title=paste("Sample Swap Test:", res$target),
             subtitle=sprintf("r(vs %s)=%.3f  |  r(vs %s)=%.3f",
                              res$exp_group, res$exp_cor, res$susp_group, res$susp_cor),
             x=NULL, y="Pearson correlation (log2CPM)", fill="Condition") +
        ylim(0.8, 1.0)
    })

    output$swap_result_text <- renderPrint({
      req(swap_result())
      res <- swap_result()
      cat("=== Sample Swap Test Result ===\n")
      cat(sprintf("対象サンプル: %s\n", res$target))
      cat(sprintf("期待グループ   (%s): avg r = %.4f\n", res$exp_group,  res$exp_cor))
      cat(sprintf("スワップ疑い   (%s): avg r = %.4f\n", res$susp_group, res$susp_cor))
      if (res$is_swap) {
        cat("\n>>> 警告: スワップ先疑いのグループとの相関がより高い\n")
        cat("    サンプルスワップの可能性を支持する結果です。\n")
        cat("    実験担当者に分注手順の確認を依頼してください。\n")
      } else {
        cat("\n>>> 期待グループとの相関がより高い。\n")
        cat("    サンプルスワップの可能性は低いです。\n")
      }
    })

    output$cor_heatmap <- renderPlot({
      req(swap_result())
      res  <- swap_result()
      meta <- res$meta
      ann_col <- meta[, intersect(c("batch","glucose","stim","group"), colnames(meta)), drop=FALSE]
      rownames(ann_col) <- meta$current_name
      pheatmap(res$cor_mat,
               annotation_col=ann_col, annotation_row=ann_col,
               color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
               breaks=seq(0.85, 1.0, length.out=101),
               display_numbers=round(res$cor_mat,3),
               number_format="%.3f", fontsize_number=7,
               main=paste("Sample Correlation Heatmap →", res$target, "swap check"),
               cluster_rows=TRUE, cluster_cols=TRUE)
    })

    output$swap_table <- renderDT({
      req(swap_result())
      datatable(swap_result()$group_cors %>% mutate(correlation=round(correlation,4)),
                options=list(pageLength=15, scrollX=TRUE), rownames=FALSE)
    })

    # ============================================================
    # AI解釈プロンプト生成
    # ============================================================

    ai_prompt_text <- eventReactive(input$gen_prompt, {
      req(rv$sample_metadata)
      meta <- rv$sample_metadata[rv$sample_metadata$active, , drop=FALSE]

      lines <- c()
      lines <- c(lines,
        "# RNA-seq 解析結果 — AI解釈依頼",
        "",
        paste0("**実験概要**: ", input$ai_experiment),
        "",
        input$ai_question,
        "",
        "---"
      )

      # 実験デザイン情報
      if ("design" %in% input$ai_include) {
        lines <- c(lines,
          "",
          "## 実験デザイン",
          paste0("- サンプル数: ", nrow(meta), " サンプル"),
          paste0("- 生物種: ", rv$selected_species %||% "不明"),
          ""
        )
        if ("group" %in% colnames(meta)) {
          grp_tbl <- table(meta$group)
          lines <- c(lines, "**条件グループ:**")
          for (g in names(grp_tbl)) {
            lines <- c(lines, paste0("- ", g, ": ", grp_tbl[[g]], " サンプル"))
          }
        }
        if ("batch" %in% colnames(meta)) {
          batch_tbl <- table(meta$batch)
          lines <- c(lines, "", "**バッチ構成:**")
          for (b in names(batch_tbl)) {
            lines <- c(lines, paste0("- バッチ ", b, ": ", batch_tbl[[b]], " サンプル"))
          }
        }
      }

      # スワップ検証結果
      if ("swap" %in% input$ai_include && !is.null(swap_result()) &&
          tryCatch(!is.null(swap_result()$target), error=function(e) FALSE)) {
        res <- swap_result()
        lines <- c(lines,
          "",
          "## サンプルスワップ検証結果",
          paste0("- 対象サンプル: ", res$target),
          paste0("- 期待グループ (", res$exp_group, ") との平均相関: r = ", round(res$exp_cor, 4)),
          paste0("- スワップ疑い (", res$susp_group, ") との平均相関: r = ", round(res$susp_cor, 4)),
          paste0("- 判定: ", if (res$is_swap) "スワップの可能性あり（スワップ先疑いのほうが相関高い）" else "スワップの可能性は低い"),
          "",
          "**上位相関サンプル Top5:**"
        )
        top5 <- head(res$group_cors, 5)
        for (i in 1:nrow(top5)) {
          lines <- c(lines,
            sprintf("  %d. %s (group=%s, r=%.4f)",
                    i, top5$sample[i], top5$group[i], top5$correlation[i]))
        }
      }

      # D19応答
      if ("response" %in% input$ai_include &&
          !is.null(response_data_react()) &&
          tryCatch(nrow(response_data_react()) > 0, error=function(e) FALSE)) {
        resp <- response_data_react()
        lines <- c(lines,
          "",
          "## D19刺激応答（log2FC: D19 / medium）",
          "| 遺伝子 | バッチ | 糖種 | log2FC |",
          "|--------|--------|------|--------|"
        )
        for (i in 1:nrow(resp)) {
          r <- resp[i, ]
          sugar <- if ("glucose" %in% colnames(r)) r$glucose else "-"
          lines <- c(lines,
            sprintf("| %s | %s | %s | %.2f |",
                    r$gene, r$batch, sugar, r$log2FC))
        }
      }

      lines <- c(lines, "", "---", "",
                 "上記のRNA-seq解析結果について、生物学的な観点から解釈をお願いします。",
                 "特に以下の点についてコメントいただけると助かります：",
                 "1. 実験結果の妥当性（バッチ間の一貫性）",
                 "2. 異常サンプルが存在する場合の対処方針",
                 "3. 解析結果から得られる生物学的知見",
                 "4. 論文の考察に使えるポイント")

      paste(lines, collapse="\n")
    })

    output$ai_prompt <- renderText({
      if (input$gen_prompt == 0) {
        return("← 「プロンプト生成」ボタンを押すと、AIに送るためのテキストが生成されます。")
      }
      req(ai_prompt_text())
      ai_prompt_text()
    })

    # ============================================================
    # ダウンロード
    # ============================================================

    output$dl_bar_pdf <- downloadHandler(
      filename = "barplot_key_genes.pdf",
      content  = function(f) {
        req(barplot_react())
        ggsave(f, barplot_react(), width=12, height=8, device=cairo_pdf)
      }
    )
    output$dl_bar_png <- downloadHandler(
      filename = "barplot_key_genes.png",
      content  = function(f) {
        req(barplot_react())
        ggsave(f, barplot_react(), width=12, height=8, dpi=200)
      }
    )
    output$dl_resp_pdf <- downloadHandler(
      filename = "D19_response_by_batch.pdf",
      content  = function(f) {
        req(response_data_react())
        p <- ggplot(response_data_react(),
                    aes(x=x_label, y=log2FC, fill=batch)) +
          geom_bar(stat="identity") + geom_hline(yintercept=0, linetype="dashed") +
          facet_wrap(~gene, scales="free_y") + theme_bw(base_size=11)
        ggsave(f, p, width=10, height=6, device=cairo_pdf)
      }
    )
    output$dl_resp_csv <- downloadHandler(
      filename = "D19_response_by_batch.csv",
      content  = function(f) {
        req(response_data_react())
        write.csv(response_data_react(), f, row.names=FALSE)
      }
    )
    output$dl_swap_pdf <- downloadHandler(
      filename = "sample_swap_test.pdf",
      content  = function(f) {
        req(swap_result())
        # renderPlotと同じ内容を再現
        pdf(f, width=9, height=6)
        res <- swap_result()
        df  <- res$group_cors %>%
          mutate(
            highlight = case_when(
              group == res$exp_group  ~ paste0(res$exp_group,  " (expected)"),
              group == res$susp_group ~ paste0(res$susp_group, " (swap hypothesis)"),
              TRUE ~ "Other"),
            sample = factor(sample, levels=rev(res$group_cors$sample))
          )
        colors <- setNames(c("#D55E00","#0072B2","#999999"),
                           c(paste0(res$exp_group," (expected)"),
                             paste0(res$susp_group," (swap hypothesis)"),
                             "Other"))
        p <- ggplot(df, aes(x=sample, y=correlation, fill=highlight)) +
          geom_bar(stat="identity") + coord_flip() +
          geom_hline(yintercept=res$exp_cor,  linetype="dashed", color="#D55E00") +
          geom_hline(yintercept=res$susp_cor, linetype="dotted", color="#0072B2") +
          scale_fill_manual(values=colors, drop=FALSE) +
          theme_bw(base_size=12) +
          labs(title=paste("Sample Swap Test:", res$target),
               x=NULL, y="Pearson correlation", fill="Condition") +
          ylim(0.8, 1.0)
        print(p)
        dev.off()
      }
    )
    output$dl_swap_csv <- downloadHandler(
      filename = "sample_swap_correlations.csv",
      content  = function(f) {
        req(swap_result())
        write.csv(swap_result()$group_cors, f, row.names=FALSE)
      }
    )
    output$dl_heatmap_pdf <- downloadHandler(
      filename = "sample_correlation_heatmap.pdf",
      content  = function(f) {
        req(swap_result())
        res  <- swap_result()
        meta <- res$meta
        ann_col <- meta[, intersect(c("batch","glucose","stim","group"), colnames(meta)), drop=FALSE]
        rownames(ann_col) <- meta$current_name
        pdf(f, width=10, height=9)
        pheatmap(res$cor_mat,
                 annotation_col=ann_col, annotation_row=ann_col,
                 color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
                 breaks=seq(0.85,1.0,length.out=101),
                 display_numbers=round(res$cor_mat,3),
                 number_format="%.3f", fontsize_number=7,
                 cluster_rows=TRUE, cluster_cols=TRUE)
        dev.off()
      }
    )
    output$dl_prompt_txt <- downloadHandler(
      filename = "ai_prompt.txt",
      content  = function(f) {
        req(ai_prompt_text())
        writeLines(ai_prompt_text(), f)
      }
    )
  })
}

# NULL-coalescing operator (base R で未定義の場合のフォールバック)
`%||%` <- function(a, b) if (!is.null(a)) a else b
