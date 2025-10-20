# R/module_timeseries_analysis.R

# --- 1. UI定義 ---
timeseriesAnalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        h3("時系列解析 (maSigPro)"),
        p("maSigProを用いて、時系列データから発現パターンが有意に変化する遺伝子を検出します。"),
        
        # --- 入力パラメータ ---
        selectInput(ns("time_col"), "時間を示す列を選択", choices = NULL),
        selectInput(ns("replicate_col"), "生物学的反復を示す列を選択", choices = NULL),
        selectInput(ns("group_cols"), "比較する実験群の列を選択 (1つ以上)", choices = NULL, multiple = TRUE),
        
        numericInput(ns("degree"), "回帰モデルの次数 (Degree)", value = 2, min = 1, max = 5),
        numericInput(ns("alpha"), "有意水準 (Alpha)", value = 0.05, min = 0.001, max = 0.1),
        numericInput(ns("rsq_cutoff"), "R-squared カットオフ", value = 0.7, min = 0, max = 1),
        
        actionButton(ns("run_analysis"), "解析実行", class = "btn-primary"),
        
        hr(),
        h4("結果の可視化"),
        # UI for selecting which significant genes to plot
        selectInput(ns("sig_gene_group"), "表示する遺伝子群", choices = NULL),
        actionButton(ns("plot_genes"), "選択した遺伝子群をプロット")
        
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("解析結果サマリー", 
                   withSpinner(DT::DTOutput(ns("sig_genes_summary_table"))),
                   verbatimTextOutput(ns("maSigPro_summary"))
          ),
          tabPanel("遺伝子プロット", 
                   p("「結果の可視化」セクションで表示したい遺伝子群を選択し、「プロット」ボタンを押してください。"),
                   withSpinner(plotOutput(ns("gene_plots"), height = "800px"))
          ),
          tabPanel("クラスタリングプロット",
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
      design = NULL,
      data = NULL
    )
    
    # --- メタデータ列をUIに反映 ---
    observe({
      req(rv$sample_metadata)
      choices <- colnames(rv$sample_metadata)
      updateSelectInput(session, "time_col", choices = choices, selected = choices[1])
      updateSelectInput(session, "replicate_col", choices = choices, selected = choices[2])
      updateSelectInput(session, "group_cols", choices = choices, selected = choices[3])
    })
    
    # --- 解析実行 ---
    observeEvent(input$run_analysis, {
      req(rv$filtered_keep, rv$merged_data, rv$sample_metadata,
          input$time_col, input$replicate_col, length(input$group_cols) > 0)
      
      shiny::withProgress(message = 'maSigPro解析を実行中...', value = 0, {
        
        # 1. データの準備
        incProgress(0.1, detail = "データを準備しています...")
        
        # フィルタリングされたカウントデータを取得
        counts <- rv$merged_data[rv$filtered_keep, -1] # Geneid列を除外
        rownames(counts) <- rv$merged_data$Geneid[rv$filtered_keep]
        
        # design matrixの作成
        meta <- rv$sample_metadata
        
        # maSigPro用のedesignを作成
        # 必要な列: Time, Replicate, そして実験群の列
        # ユーザーが選択した列名がmaSigProの要求する'Time', 'Replicate'と異なる場合があるので、ここで列名を変更
        edesign <- meta[, c(input$time_col, input$replicate_col, input$group_cols), drop = FALSE]
        colnames(edesign)[1] <- "Time"
        colnames(edesign)[2] <- "Replicate"
        
        # 時間の列は数値である必要がある
        if(!is.numeric(edesign$Time)) {
            shiny::showNotification("時間を示す列は数値である必要があります。", type = "error")
            return(NULL)
        }

        # 2. maSigProの実行
        incProgress(0.3, detail = "p.vectorを実行中...")
        
        # 正規化 (TMMなど) - この例ではedgeRのDGEListとcalcNormFactorsを流用
        dge <- DGEList(counts = counts)
        dge <- calcNormFactors(dge)
        norm_counts <- cpm(dge, log = FALSE)
        
        # make.design.matrix
        design <- maSigPro::make.design.matrix(edesign, degree = input$degree)
        
        # p.vector
        fit <- maSigPro::p.vector(norm_counts, design, Q = 1, MT.adjust = "BH", min.obs = 3)
        
        incProgress(0.6, detail = "T.fitを実行中...")
        # T.fit
        tstep <- maSigPro::T.fit(fit, step.method = "backward", alfa = input$alpha)
        
        incProgress(0.8, detail = "有意遺伝子を抽出中...")
        # get.siggenes
        sigs <- maSigPro::get.siggenes(tstep, rsq = input$rsq_cutoff, vars = "groups")
        
        # 結果を保存
        results$sig_genes <- sigs
        results$design <- design
        results$data <- norm_counts
        
        # UIの更新
        updateSelectInput(session, "sig_gene_group", choices = names(sigs$sig.genes))
        
        shiny::showNotification("maSigPro解析が完了しました。", type = "message")
        incProgress(1, detail = "完了")
      })
    })
    
    # --- 解析結果サマリー ---
    output$sig_genes_summary_table <- DT::renderDT({
      req(results$sig_genes)
      DT::datatable(results$sig_genes$summary, options = list(scrollX = TRUE))
    })
    
    output$maSigPro_summary <- renderPrint({
      req(results$sig_genes)
      results$sig_genes$summary
    })
    
    # --- 遺伝子プロット ---
    observeEvent(input$plot_genes, {
        req(results$sig_genes, results$data, results$design, input$sig_gene_group)
        
        output$gene_plots <- renderPlot({
            maSigPro::see.genes(results$sig_genes$sig.genes[[input$sig_gene_group]], 
                                edesign = results$design$edesign, 
                                show.fit = TRUE, 
                                dis = results$design$dis,
                                cluster.data = 1, 
                                k = 9, # k-meansのクラスタ数 (例)
                                newX11 = FALSE)
        })

        output$cluster_plot <- renderPlot({
            maSigPro::PlotGroups(results$data[results$sig_genes$sig.genes[[input$sig_gene_group]][,1],], 
                                 edesign = results$design$edesign, 
                                 show.fit = TRUE, 
                                 dis = results$design$dis,
                                 newX11 = FALSE)
        })
    })

  })
}
