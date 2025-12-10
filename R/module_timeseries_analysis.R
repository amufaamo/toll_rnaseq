# R/module_timeseries_analysis.R
# ★★★ maSigProとrhandsontableをロード ★★★
# (shiny::runApp()の前に 'rhandsontable' と 'maSigPro' が
#  インストール＆ロードされている必要があります)
# library(maSigPro)
# library(rhandsontable) 

# --- 1. UI定義 ---
timeseriesAnalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        width = 5, 
        h3("時系列解析 (maSigPro)"),
        p("maSigProを用いて、時系列データから発現パターンが有意に変化する遺伝子を検出します。"),
        
        # --- ★★★ 修正箇所: edesign (実験デザイン) の定義UI ★★★ ---
        hr(),
        h4("1. 実験デザイン(edesign)の定義"),
        helpText("「フィルタリング」タブで選択したアクティブなサンプルを読み込み、maSigPro用の実験デザイン（Time, Replicate）を定義してください。"),
        
        actionButton(ns("load_active_samples"), "アクティブなサンプルを読み込み/リセット", class = "btn-info"),
        br(), br(),
        
        helpText("テーブル内の Time と Replicate を直接編集してください。"),
        rHandsontableOutput(ns("edesign_table")),
        br(),
        
        # ★★★ 削除: グループ（群）を追加するUIを削除 ★★★
        # (単一系列の解析を前提とするため)
        
        # --- ★★★ 修正ここまで ★★★ ---
        
        hr(),
        h4("2. maSigPro パラメータ設定"),
        
        numericInput(ns("degree"), "回帰モデルの次数 (Degree)", value = 2, min = 1, max = 5),
        helpText("3点の時系列の場合、次数(Degree)は 1 (直線) または 2 (曲線/放物線) を推奨します。"),
        numericInput(ns("alpha"), "有意水準 (Alpha)", value = 0.05, min = 0.001, max = 0.1),
        numericInput(ns("rsq_cutoff"), "R-squared カットオフ", value = 0.7, min = 0, max = 1),
        
        actionButton(ns("run_analysis"), "解析実行", class = "btn-primary"),
        
        hr(),
        h4("3. 結果の可視化"),
        selectInput(ns("sig_gene_group"), "表示する遺伝子群", choices = NULL),
        
        # --- ★★★ 修正箇所: kのデフォルト値を変更 ★★★ ---
        numericInput(ns("k_clusters"), "クラスター数 (k)", value = 5, min = 2, max = 20),
        helpText("3点の時系列(例: 2d, 3d, 4d)の場合、k=5 (上昇/減少/山型/谷型/変動なし)あたりが解釈しやすいデフォルトです。"),
        # --- ★★★ 修正ここまで ★★★ ---
        
        actionButton(ns("plot_genes"), "選択した遺伝子群をプロット")
        
      ),
      mainPanel(
        width = 7, 
        tabsetPanel(
          tabPanel("解析結果サマリー", 
                   withSpinner(DT::DTOutput(ns("sig_genes_summary_table"))),
                   verbatimTextOutput(ns("maSigPro_summary"))
          ),
          tabPanel("遺伝子プロット (see.genes)", 
                   p("「結果の可視化」セクションで表示したい遺伝子群を選択し、「プロット」ボタンを押してください。"),
                   helpText("似た発現パターンを持つ遺伝子同士がクラスタリングされて表示されます。"),
                   withSpinner(plotOutput(ns("gene_plots"), height = "800px"))
          ),
          tabPanel("クラスタリングプロット (PlotGroups)",
                   helpText("上記で選択した遺伝子群について、各クラスターの平均プロファイルを表示します。"),
                   withSpinner(plotOutput(ns("cluster_plot"), height = "800px"))
          )
        )
      )
    )
  )
}

# --- 2. Serverロジック ---
timeseriesAnalysisServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    
    # --- リアクティブ値 ---
    results <- reactiveValues(
      sig_genes = NULL,
      design = NULL, # make.design.matrixの結果
      data = NULL # 正規化済みカウントデータ
    )
    
    # edesignテーブルを保持するReactiveVal
    edesign_data <- reactiveVal(NULL)
    
    # --- ★★★ 修正箇所: edesignテーブルのロジック (シンプル化) ★★★ ---
    
    # 1. アクティブなサンプルを読み込む
    observeEvent(input$load_active_samples, {
      req(rv$sample_metadata)
      
      active_meta <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      
      validate(need(nrow(active_meta) > 0, "アクティブなサンプルがありません。「データ入力」タブでサンプルを選択してください。"))
      
      # maSigProのedesignの基本形を作成
      # ★★★ 修正: 単一グループ（'Group'）を前提とし、'1' を入れておく ★★★
      df <- data.frame(
        SampleName = active_meta$current_name,
        Time = 0,
        Replicate = 1,
        Group = 1, # 単一グループを 'Group' とし、すべて '1' に
        stringsAsFactors = FALSE
      )
      
      edesign_data(df)
      showNotification("アクティブなサンプルを読み込みました。Time と Replicate を編集してください。", type = "message")
    })
    
    # 2. ★★★ 削除: 実験群の列を追加するロジックを削除 ★★★
    
    # 3. rhandsontableを描画する
    output$edesign_table <- renderRHandsontable({
      df <- edesign_data()
      req(df)
      
      rhandsontable(df, readOnly = FALSE, rowHeaders = NULL, stretchH = "all") %>%
        hot_col("SampleName", readOnly = TRUE) %>% # SampleName列は編集不可
        hot_col("Group", readOnly = TRUE) %>%      # ★★★ 追加: Group列も編集不可に
        hot_col(c("Time", "Replicate"), type = "numeric") # 数値入力を強制
    })
    
    # --- ★★★ 修正箇所: 解析実行ロジック (vars = 'all') ★★★ ---
    observeEvent(input$run_analysis, {
      
      # フィルタリング結果と、編集されたedesignテーブルを要求
      req(rv$filtered_keep, rv$merged_data, input$edesign_table)
      
      edesign_from_table <- hot_to_r(input$edesign_table)
      
      # バリデーション
      validate(
        need(all(c("Time", "Replicate", "Group") %in% colnames(edesign_from_table)), "Time, Replicate, Group列が必要です。"),
        need(is.numeric(edesign_from_table$Time), "Time列は数値である必要があります。"),
        need(is.numeric(edesign_from_table$Replicate), "Replicate列は数値である必要があります。")
      )
      
      shiny::withProgress(message = 'maSigPro解析を実行中...', value = 0, {
        
        # 1. データの準備
        incProgress(0.1, detail = "データを準備しています...")
        
        # (A) maSigPro用の 'edesign' を作成
        rownames(edesign_from_table) <- edesign_from_table$SampleName
        edesign <- edesign_from_table[, -which(colnames(edesign_from_table) == "SampleName"), drop = FALSE]
        
        # (B) フィルタリング済みの 'counts' データを準備
        counts_filtered <- rv$merged_data[rv$filtered_keep, , drop = FALSE]
        rownames(counts_filtered) <- counts_filtered$Geneid
        counts_filtered <- counts_filtered[, -1, drop = FALSE] # Geneid列を削除
        
        # (C) edesign と counts のサンプル順を一致させる ★最重要★
        sample_order <- rownames(edesign)
        
        validate(need(all(sample_order %in% colnames(counts_filtered)), "edesignのサンプル名がカウントデータに見つかりません。"))
        
        counts <- counts_filtered[, sample_order, drop = FALSE]
        
        # 2. maSigProの実行
        incProgress(0.3, detail = "正規化とデザイン行列の作成中...")
        
        # 正規化
        dge <- DGEList(counts = counts)
        dge <- calcNormFactors(dge)
        norm_counts <- cpm(dge, log = FALSE)
        
        # make.design.matrix (ユーザー定義のedesignを使用)
        design <- tryCatch({
          maSigPro::make.design.matrix(edesign, degree = input$degree)
        }, error = function(e) {
          showNotification(paste("make.design.matrixでエラー:", e$message), type = "error")
          NULL
        })
        validate(need(!is.null(design), "デザイン行列の作成に失敗しました。"))
        
        incProgress(0.5, detail = "p.vectorを実行中...")
        fit <- tryCatch({
          maSigPro::p.vector(norm_counts, design, Q = 1, MT.adjust = "BH", min.obs = 3)
        }, error = function(e) {
          showNotification(paste("p.vectorでエラー:", e$message), type = "error")
          NULL
        })
        validate(need(!is.null(fit), "p.vectorの実行に失敗しました。"))
        
        incProgress(0.7, detail = "T.fitを実行中...")
        tstep <- tryCatch({
          maSigPro::T.fit(fit, step.method = "backward", alfa = input$alpha)
        }, error = function(e) {
          showNotification(paste("T.fitでエラー:", e$message), type = "error")
          NULL
        })
        validate(need(!is.null(tstep), "T.fitの実行に失敗しました。"))
        
        
        incProgress(0.9, detail = "有意遺伝子を抽出中...")
        # ★★★ 修正: vars = "groups" から "all" に変更 ★★★
        # "all" を指定すると、時間経過(Time)によって有意に変動する
        # すべての遺伝子（単一系列解析の目的）を取得します。
        sigs <- tryCatch({
          maSigPro::get.siggenes(tstep, rsq = input$rsq_cutoff, vars = "all")
        }, error = function(e) {
          showNotification(paste("get.siggenesでエラー:", e$message), type = "error")
          NULL
        })
        validate(need(!is.null(sigs), "有意遺伝子の抽出に失敗しました。"))
        
        # 結果を保存
        results$sig_genes <- sigs
        results$design <- design 
        results$data <- norm_counts 
        
        # UIの更新
        if (!is.null(sigs$sig.genes) && length(names(sigs$sig.genes)) > 0) {
          updateSelectInput(session, "sig_gene_group", choices = names(sigs$sig.genes))
          showNotification("maSigPro解析が完了しました。", type = "message")
        } else {
          updateSelectInput(session, "sig_gene_group", choices = NULL)
          showNotification("解析は完了しましたが、有意な遺伝子群が見つかりませんでした。", type = "warning")
        }
        
        incProgress(1, detail = "完了")
      })
    })
    
    # --- 解析結果サマリー ---
    output$sig_genes_summary_table <- DT::renderDT({
      req(results$sig_genes)
      # $summaryはリストかもしれないので、data.frameに変換
      summary_df <- as.data.frame(do.call(rbind, results$sig_genes$summary))
      DT::datatable(summary_df, options = list(scrollX = TRUE, pageLength = 10), rownames = TRUE)
    })
    
    output$maSigPro_summary <- renderPrint({
      req(results$sig_genes)
      print(results$sig_genes$summary)
    })
    
    # --- 遺伝子プロット ---
    # k (クラスター数) はUIから取得 (input$k_clusters)
    observeEvent(input$plot_genes, {
      req(results$sig_genes, results$data, results$design, input$sig_gene_group)
      
      sig_genes_data <- results$sig_genes$sig.genes[[input$sig_gene_group]]
      
      validate(
        need(!is.null(sig_genes_data), "選択された群の有意遺伝子データがありません。"),
        need(nrow(sig_genes_data) > 0, "選択された群に有意遺伝子がありません。"),
        need(nrow(sig_genes_data) >= input$k_clusters, paste("クラスター数(k)は遺伝子数(", nrow(sig_genes_data), ")以下にしてください。"))
      )
      
      output$gene_plots <- renderPlot({
        maSigPro::see.genes(sig_genes_data, 
                            edesign = results$design$edesign, 
                            show.fit = TRUE, 
                            dis = results$design$dis,
                            cluster.data = 1, 
                            k = input$k_clusters, # ★UIからkを取得 (デフォルト5)
                            newX11 = FALSE,
                            main = paste("Gene Profiles:", input$sig_gene_group, "(k=", input$k_clusters, ")"))
      })
      
      output$cluster_plot <- renderPlot({
        see_genes_result <- maSigPro::see.genes(sig_genes_data, 
                                                edesign = results$design$edesign, 
                                                show.fit = TRUE, 
                                                dis = results$design$dis,
                                                cluster.data = 1, 
                                                k = input$k_clusters,
                                                newX11 = FALSE) # プロットはしない
        
        maSigPro::PlotGroups(results$data[rownames(sig_genes_data), ], 
                             edesign = results$design$edesign,
                             dis = results$design$dis,
                             groups = see_genes_result$cut, # ★see.genesのクラスタリング結果を使用
                             newX11 = FALSE,
                             main = paste("Cluster Averages:", input$sig_gene_group, "(k=", input$k_clusters, ")"))
      })
      
      showNotification("プロットを更新しました。", type = "message")
      
    })
    
  })
}