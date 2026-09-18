# R/module_data_upload_metadata.R
# (エイリアス検索によるID変換率向上版)

library(shiny)
library(data.table)
library(purrr)
library(dplyr)
library(shinycssloaders)
library(plotly)
library(AnnotationDbi)
library(DT)
library(rhandsontable)

# --- Species UI Choices ---
# (変更なし)
species_choices_ui <- c(
  "Human (Homo sapiens)" = "Homo_sapiens",
  "Mouse (Mus musculus)" = "Mus_musculus",
  "Rat (Rattus norvegicus)" = "Rattus_norvegicus",
  "Fly (Drosophila melanogaster)" = "Drosophila_melanogaster",
  "Worm (Caenorhabditis elegans)" = "Caenorhabditis_elegans",
  "Zebrafish (Danio rerio)" = "Danio_rerio",
  "Yeast (Saccharomyces cerevisiae)" = "Saccharomyces_cerevisiae",
  "Arabidopsis (Arabidopsis thaliana)" = "Arabidopsis_thaliana",
  "Bovine (Bos taurus)" = "Bos_taurus",
  "Chicken (Gallus gallus)" = "Gallus_gallus",
  "Canine (Canis familiaris)" = "Canis_familiaris",
  "Rhesus (Macaca mulatta)" = "Macaca_mulatta",
  "Chimp (Pan troglodytes)" = "Pan_troglodytes",
  "Pig (Sus scrofa)" = "Sus_scrofa",
  "Xenopus (Xenopus laevis)" = "Xenopus_laevis",
  "Anopheles (Anopheles gambiae)" = "Anopheles_gambiae",
  "E. coli K12 (Escherichia coli K12)" = "Escherichia_coli_K12",
  "E. coli Sakai (Escherichia coli Sakai)" = "Escherichia_coli_Sakai",
  "Malaria (Plasmodium falciparum)" = "Plasmodium_falciparum",
  "Myxococcus xanthus DK 1622 (Myxococcus xanthus)" = "Myxococcus_xanthus_DK1622",
  "Others (Keep Original ID)" = "Others_Original"
)

# --- Species to OrgDb package mapping ---
# (変更なし)
orgdb_species_map <- c(
  "Homo_sapiens" = "org.Hs.eg.db",
  "Mus_musculus" = "org.Mm.eg.db",
  "Rattus_norvegicus" = "org.Rn.eg.db",
  "Drosophila_melanogaster" = "org.Dm.eg.db",
  "Caenorhabditis_elegans" = "org.Ce.eg.db",
  "Danio_rerio" = "org.Dr.eg.db",
  "Saccharomyces_cerevisiae" = "org.Sc.sgd.db",
  "Arabidopsis_thaliana" = "org.At.tair.db",
  "Bos_taurus" = "org.Bt.eg.db",
  "Gallus_gallus" = "org.Gg.eg.db",
  "Canis_familiaris" = "org.Cf.eg.db",
  "Macaca_mulatta" = "org.Mmu.eg.db",
  "Pan_troglodytes" = "org.Pt.eg.db",
  "Sus_scrofa" = "org.Ss.eg.db",
  "Xenopus_laevis" = "org.Xl.eg.db",
  "Anopheles_gambiae" = "org.Ag.eg.db",
  "Escherichia_coli_K12" = "org.EcK12.eg.db",
  "Escherichia_coli_Sakai" = "org.EcSakai.eg.db",
  "Plasmodium_falciparum" = "org.Pf.plasmo.db",
  "Myxococcus_xanthus_DK1622" = "org.Mxanthus.db"
)

# --- Helper functions (ID detection/Conversion) ---
detect_gene_id_type <- function(ids) {
  ids_clean <- na.omit(ids[ids != "" & !is.na(ids)])
  if (length(ids_clean) == 0) {
    return("UNKNOWN")
  }
  n_total <- length(ids_clean)
  n_sample <- min(n_total, 1000)
  ids_sample <- sample(ids_clean, n_sample)
  if (mean(grepl("^ENS[A-Z0-9]*[FPTG]\\d{10,}(\\.\\d+)?$", ids_sample, ignore.case = TRUE)) > 0.8) {
    return("ENSEMBL")
  }
  if (mean(grepl("^[NX][CMPWTZ]_[0-9]+(\\.\\d+)?$", ids_sample, ignore.case = TRUE)) > 0.8) {
    return("REFSEQ")
  }
  if (mean(grepl("^[0-9]+$", ids_sample)) > 0.9) {
    return("ENTREZID")
  }
  if (mean(grepl("^([A-Za-z][A-Za-z0-9-]*[A-Za-z0-9]|[A-Za-z])$", ids_sample)) > 0.7 && mean(grepl("^ENS", ids_sample, T)) < 0.1 && mean(grepl("^[NX][M_]", ids_sample, T)) < 0.1 && mean(grepl("^[0-9]+$", ids_sample)) < 0.1) {
    return("SYMBOL")
  }
  return("UNKNOWN")
}

infer_species_from_gene_ids <- function(ids) {
  ids_clean <- na.omit(as.character(ids[ids != "" & !is.na(ids)]))
  if (length(ids_clean) == 0) {
    return(NULL)
  }
  n_sample <- min(length(ids_clean), 1000)
  ids_sample <- sample(ids_clean, n_sample)
  species_patterns <- c(
    "Mus_musculus" = "^ENSMUSG",
    "Homo_sapiens" = "^ENSG",
    "Rattus_norvegicus" = "^ENSRNOG"
  )
  match_rates <- vapply(species_patterns, function(pattern) {
    mean(grepl(pattern, ids_sample, ignore.case = TRUE))
  }, numeric(1))
  best_match <- names(which.max(match_rates))
  if (length(best_match) == 1 && match_rates[[best_match]] > 0.8) {
    return(best_match)
  }
  NULL
}

# ★★★ 修正: エイリアス検索を追加した関数 ★★★
convert_or_pass_ids <- function(original_ids, detected_keytype, selected_species_code) {
  if (selected_species_code == "Others_Original") {
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(FALSE, length(original_ids)), failed_count = 0, final_id_type = detected_keytype %||% "ORIGINAL_UNCONVERTED", conversion_map = NULL))
  }
  if (selected_species_code == "Lotus_japonicus") {
    lotus_entrez <- gsub("^LOC", "", as.character(original_ids))
    debug_map <- data.frame(
      Original = as.character(original_ids),
      Converted = lotus_entrez,
      Status = "Success",
      stringsAsFactors = FALSE
    )
    return(list(processed_ids = lotus_entrez,
                failed_indices = rep(FALSE, length(original_ids)),
                failed_count = 0,
                final_id_type = "ENTREZID",
                conversion_map = debug_map))
  }
  if (detected_keytype == "ENTREZID") {
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(FALSE, length(original_ids)), failed_count = 0, final_id_type = "ENTREZID", conversion_map = NULL))
  }

  # OrgDbの準備
  orgdb_pkg_name <- orgdb_species_map[[selected_species_code]]
  if (is.null(orgdb_pkg_name) || !requireNamespace(orgdb_pkg_name, quietly = TRUE)) {
    warning(paste0("OrgDb package '", orgdb_pkg_name, "' not found."))
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(TRUE, length(original_ids)), failed_count = length(original_ids), final_id_type = detected_keytype %||% "UNKNOWN", conversion_map = NULL))
  }

  require(orgdb_pkg_name, character.only = TRUE, quietly = TRUE)
  org_db <- get(orgdb_pkg_name)
  keys_to_convert <- as.character(original_ids)
  if (detected_keytype == "ENSEMBL") {
    keys_to_convert <- gsub("\\..*$", "", keys_to_convert)
  }

  if (!detected_keytype %in% keytypes(org_db)) {
    warning(paste0("Keytype '", detected_keytype, "' not supported by '", orgdb_pkg_name, "'."))
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(TRUE, length(original_ids)), failed_count = length(original_ids), final_id_type = detected_keytype %||% "UNKNOWN", conversion_map = NULL))
  }

  message(paste0("[ID_CONV] 1st pass: Converting ", length(unique(keys_to_convert)), " unique '", detected_keytype, "' IDs..."))
  unique_input_keys <- unique(keys_to_convert)

  # 1. まず通常の変換（Symbol -> Entrez）
  entrez_map <- tryCatch(
    suppressMessages(mapIds(org_db, keys = unique_input_keys, column = "ENTREZID", keytype = detected_keytype, multiVals = "first")),
    error = function(e) {
      warning(e$message)
      NULL
    }
  )

  if (is.null(entrez_map)) {
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(TRUE, length(original_ids)), failed_count = length(original_ids), final_id_type = detected_keytype %||% "UNKNOWN", conversion_map = NULL))
  }

  # 2. 失敗したものを特定
  failed_keys_1st <- names(entrez_map)[is.na(entrez_map)]

  # 3. エイリアス（別名）で再挑戦 (Symbolの場合のみ有効)
  if (detected_keytype == "SYMBOL" && length(failed_keys_1st) > 0 && "ALIAS" %in% keytypes(org_db)) {
    message(paste0("[ID_CONV] 2nd pass: Trying ALIAS lookup for ", length(failed_keys_1st), " failed IDs..."))

    alias_map <- tryCatch(
      suppressMessages(mapIds(org_db, keys = failed_keys_1st, column = "ENTREZID", keytype = "ALIAS", multiVals = "first")),
      error = function(e) {
        NULL
      }
    )

    if (!is.null(alias_map)) {
      # 成功したものを元のマップに統合
      rescued_count <- sum(!is.na(alias_map))
      message(paste0("[ID_CONV] Rescued ", rescued_count, " IDs using ALIAS."))
      entrez_map[names(alias_map)] <- alias_map
    }
  }

  # 結果の整形
  converted_values <- entrez_map[keys_to_convert]
  failed_indices_logical <- is.na(converted_values)
  num_failed <- sum(failed_indices_logical)

  final_ids <- as.character(converted_values)
  # 失敗した箇所は元のIDに戻す
  final_ids[failed_indices_logical] <- as.character(original_ids[failed_indices_logical])

  # デバッグ用マップを作成 (Original -> Entrez)
  debug_map <- data.frame(
    Original = unique_input_keys,
    Converted = entrez_map[unique_input_keys],
    Status = ifelse(is.na(entrez_map[unique_input_keys]), "Failed", "Success"),
    stringsAsFactors = FALSE
  )

  return(list(
    processed_ids = final_ids,
    failed_indices = failed_indices_logical,
    failed_count = num_failed,
    final_id_type = "ENTREZID_MIXED",
    conversion_map = debug_map
  ))
}

# --- Helper function: Auto-detect file format ---
detect_file_format <- function(files_info) {
  if (is.null(files_info) || nrow(files_info) == 0) return(NULL)
  if (nrow(files_info) > 1) return("individual")

  fpath <- files_info$datapath[1]
  fname <- files_info$name[1]
  ext <- tolower(tools::file_ext(fname))
  sep_char <- if (ext == "csv") "," else "\t"

  lines <- tryCatch(readLines(fpath, n = 10, warn = FALSE), error = function(e) character(0))
  data_lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  if (length(data_lines) == 0) return("individual")

  cols <- trimws(strsplit(data_lines[1], sep_char)[[1]])
  featurecounts_markers <- c("Chr", "Start", "End", "Strand", "Length")
  if (any(featurecounts_markers %in% cols)) return("individual")
  return("merged")
}

# --- Helper function: Trim common strings ---
# (変更なし)
find_lcp <- function(strs) {
  if (length(strs) < 2) {
    return("")
  }
  char_lists <- strsplit(strs, "")
  min_len <- min(sapply(char_lists, length))
  if (min_len == 0) {
    return("")
  }
  lcp <- ""
  for (i in 1:min_len) {
    char_to_check <- char_lists[[1]][i]
    if (all(sapply(char_lists, function(x) x[i] == char_to_check))) {
      lcp <- paste0(lcp, char_to_check)
    } else {
      break
    }
  }
  return(lcp)
}
find_lcs <- function(strs) {
  if (length(strs) < 2) {
    return("")
  }
  rev_strs <- sapply(strs, function(x) intToUtf8(rev(utf8ToInt(x))))
  rev_lcs <- find_lcp(rev_strs)
  lcs <- intToUtf8(rev(utf8ToInt(rev_lcs)))
  return(lcs)
}
trim_common_strings <- function(strs) {
  if (length(strs) < 2) {
    return(strs)
  }
  lcp <- find_lcp(strs)
  lcs <- find_lcs(strs)
  trimmed_strs <- sub(paste0("^", lcp), "", strs)
  trimmed_strs <- sub(paste0(lcs, "$"), "", trimmed_strs)
  return(trimmed_strs)
}

dataUploadMetadataUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4(icon("upload"), "1. ファイルアップロード"),
      tags$div(class = "alert alert-info", style = "padding: 10px;",
               icon("info-circle"), " カウントデータをアップロードします。結合済みの行列(CSV/TSV)か、featureCountsが出力した個別ファイルかを選べます。"),
      radioButtons(ns("inputType"), "データの形式:",
        choices = c(
          "個別の featureCounts ファイル (複数可)" = "individual",
          "結合済みカウント行列 (CSV/TSV)" = "merged"
        ),
        selected = "individual"
      ),
      conditionalPanel(
        condition = paste0("input['", ns("inputType"), "'] == 'individual'"),
        fileInput(ns("featureCountsFiles"), "featureCounts出力ファイル", multiple = TRUE, accept = c(".txt", ".tsv")),
        helpText("各ファイルの1列目(Geneid)、6列目(Length)、7列目(Count)を結合します。")
      ),
      conditionalPanel(
        condition = paste0("input['", ns("inputType"), "'] == 'merged'"),
        fileInput(ns("mergedCountFile"), "カウント行列ファイル", multiple = FALSE, accept = c(".csv", ".txt", ".tsv")),
        helpText("行名がGene Symbolの場合、自動でEntrez IDに変換します。")
      ),
      hr(),
      h4(icon("paw"), "2. 生物種選択"),
      selectInput(ns("species"), label = "解析対象の生物種またはID処理方法:", choices = species_choices_ui, selected = "Homo_sapiens"),
      hr(),
      h4(icon("dna"), "2b. GTF/GFF アノテーション (任意)"),
      fileInput(ns("gtfFile"), "GTF / GFF3 ファイル", multiple = FALSE,
                accept = c(".gtf", ".gff", ".gff3", ".gz")),
      helpText(icon("info-circle"), " OrgDbの無い生物種でも gene_id→Gene Symbol 変換、遺伝子長によるTPM/FPKM、biotypeフィルタが可能になります。カウントファイルと同じ参照GTFを指定してください。"),
      uiOutput(ns("gtfStatusUI")),
      hr(),
      h4(icon("project-diagram"), "2c. カスタムアノテーション (任意)"),
      fileInput(ns("customAnnotationFile"), "eggNOG-mapper アノテーション (.emapper.annotations)", multiple = FALSE, accept = c(".annotations", ".txt", ".tsv")),
      helpText(icon("info-circle"), " 非モデル生物でGO / KEGG / GSEA解析を行う場合にアップロードしてください。"),
      uiOutput(ns("customAnnotationStatusUI")),
      hr(),
      h4(icon("sitemap"), "2d. 非モデル生物の GO 注釈 (OrgDbが無い生物用)"),
      helpText(icon("info-circle"), " OrgDbの無い生物でGOエンリッチメントを行うための注釈です。下のいずれかで設定できます（同梱済みの生物は不要）。"),
      # --- 方法A: Tax ID を入れるだけで自動取得 (コマンド不要) ---
      tags$b("方法A: NCBI Tax ID から自動取得"),
      div(class = "input-group",
        textInput(ns("taxIdInput"), NULL, value = "", placeholder = "例: 34305 (ミヤコグサ)"),
      ),
      actionButton(ns("fetchGene2go"), "この生物のGO注釈を取得", icon = icon("download"), class = "btn btn-primary btn-sm w-100"),
      helpText(icon("info-circle"), " NCBIのGO注釈を取得します。NCBI Taxonomyで種名検索→Tax ID。",
               tags$b("初回のみ約1.3GBのDB取得で数分かかります"), "（以降は即時／同梱生物は不要）。"),
      br(),
      # --- 方法B: gene2go ファイルを手動アップロード (上級者向け) ---
      tags$details(
        tags$summary("方法B: gene2go ファイルを手動アップロード (上級者向け)"),
        fileInput(ns("gene2goFile"), "gene2go ファイル (GeneID, GO, [Term], [Category] 列)", multiple = FALSE,
                  accept = c(".tsv", ".txt", ".csv", ".gz")),
        helpText(icon("info-circle"), " NCBI gene2go形式（タブ区切り; 列: GeneID, GO, 任意でTerm/Category=Process/Function/Component）。GeneIDはカウント行列の遺伝子IDと一致させてください。")
      ),
      uiOutput(ns("gene2goStatusUI")),
      textInput(ns("keggOrgCode"), "KEGG生物種コード (任意, 例: lja)", value = "", placeholder = "lja"),
      helpText(icon("info-circle"), " KEGGがサポートする生物種コードを指定すると、その生物でKEGGパスウェイ解析が可能になります（オンライン取得）。"),
      hr(),
      h4(icon("save"), "3. セッション管理"),
      downloadButton(ns("downloadRDS"), "セッション保存"),
      br(), br(),
      fileInput(ns("uploadRDS"), "セッション復元", accept = c(".rds", ".RDS"))
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel(
          "サンプル情報 & QC",
          h4("サンプル別 合計リード数"),
          uiOutput(ns("lib_size_color_ui")),
          withSpinner(plotlyOutput(ns("librarySizePlot")), type = 6),
          hr(),
          h4("サンプル情報入力"),
          actionButton(ns("add_factor_btn"), "＋ グループを追加 (Add Group)", icon = icon("plus")),
          actionButton(ns("remove_factor_btn"), "－ グループを削除 (Remove Group)", icon = icon("minus")),
          actionButton(ns("rename_group_btn"), "列名を変更 (Rename)", icon = icon("edit")),
          br(), br(),
          withSpinner(rHandsontableOutput(ns("sampleMetadataTable")), type = 6)
        ),
        tabPanel(
          "ID変換レポート (デバッグ用)",
          h4("Gene ID 変換状況"),
          p("ここで、どのような遺伝子IDが変換に成功/失敗したかを確認できます。失敗が多い場合、入力ファイルのID形式や生物種が合っているか確認してください。"),
          verbatimTextOutput(ns("idConversionSummary")),
          hr(),
          h5("変換詳細テーブル (検索可能)"),
          DTOutput(ns("idConversionTable"))
        )
      )
    )
  )
}

dataUploadMetadataServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    rv$current_gene_id_type <- reactiveVal("ENTREZID")
    rv$conversion_debug_info <- reactiveVal(NULL) # デバッグ情報用

    # --- GTF/GFF アノテーション ---
    gtf_raw <- reactiveVal(NULL)  # data.frame(gene_id, gene_name, biotype, gene_length)
    observeEvent(input$gtfFile, {
      req(input$gtfFile)
      withProgress(message = "GTF/GFFを解析中...", value = 0.5, {
        ann <- tryCatch(parse_gtf_annotation(input$gtfFile$datapath), error = function(e) {
          showNotification(paste("GTF解析エラー:", e$message), type = "error", duration = 15)
          NULL
        })
        gtf_raw(ann)
        if (!is.null(ann)) {
          # GTFの中身から「表示する遺伝子IDタイプ」の選択肢と既定値を決定し、
          # 各タブ(DEG / GO)のドロップダウンへ反映できるよう rv に格納する。
          rv$gtf_id_choices <- gtf_display_id_choices(ann)
          sel_label <- names(rv$gtf_id_choices$choices)[
            rv$gtf_id_choices$choices == rv$gtf_id_choices$selected][1]
          showNotification(paste0("GTF読み込み完了: ", nrow(ann),
                                  " 遺伝子のアノテーションを取得しました。表示IDタイプを「",
                                  sel_label, "」に自動設定しました。"),
                           type = "message", duration = 8)
        } else {
          rv$gtf_id_choices <- NULL
        }
      })
    })

    output$gtfStatusUI <- renderUI({
      ann <- gtf_raw()
      if (is.null(ann)) return(NULL)
      n_len <- sum(!is.na(ann$gene_length))
      n_bt  <- length(unique(ann$biotype[ann$biotype != "unknown"]))
      tags$div(class = "alert alert-success", style = "padding: 8px; font-size: 0.85rem;",
        icon("check-circle"),
        HTML(paste0(" GTF適用中: <b>", nrow(ann), "</b> 遺伝子 / 遺伝子長 <b>", n_len,
                    "</b> 件 / biotype <b>", n_bt, "</b> 種類"))
      )
    })

    # --- Custom Annotation (eggNOG-mapper) ---
    observeEvent(input$customAnnotationFile, {
      req(input$customAnnotationFile)
      withProgress(message = "カスタムアノテーションを解析中...", value = 0.5, {
        tryCatch({
          # Read file, skipping comments
          lines <- readLines(input$customAnnotationFile$datapath)
          data_lines <- lines[!grepl("^##", lines)]
          # Find the header line, it starts with #query
          header_idx <- which(grepl("^#query", data_lines))
          if (length(header_idx) > 0) {
            data_lines[header_idx] <- sub("^#", "", data_lines[header_idx])
          }
          
          df <- read.table(text = data_lines, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
          
          if (!"query" %in% colnames(df)) {
            stop("eggNOG-mapperのフォーマットとして認識できませんでした ('query' 列が見つかりません)。")
          }
          
          # Parse GO
          go_t2g <- NULL
          if ("GOs" %in% colnames(df)) {
            go_list <- lapply(1:nrow(df), function(i) {
              gos <- df$GOs[i]
              if (!is.na(gos) && gos != "" && gos != "-") {
                data.frame(term = trimws(strsplit(gos, ",")[[1]]), gene = df$query[i], stringsAsFactors = FALSE)
              } else {
                NULL
              }
            })
            go_t2g <- do.call(rbind, go_list)
            if (!is.null(go_t2g) && nrow(go_t2g) > 0) go_t2g <- unique(go_t2g)
          }
          
          # Parse KO
          ko_t2g <- NULL
          if ("KEGG_ko" %in% colnames(df)) {
            ko_list <- lapply(1:nrow(df), function(i) {
              kos <- df$KEGG_ko[i]
              if (!is.na(kos) && kos != "" && kos != "-") {
                kos_split <- strsplit(kos, ",")[[1]]
                kos_split <- gsub("^ko:", "", kos_split)
                data.frame(term = trimws(kos_split), gene = df$query[i], stringsAsFactors = FALSE)
              } else {
                NULL
              }
            })
            ko_t2g <- do.call(rbind, ko_list)
            if (!is.null(ko_t2g) && nrow(ko_t2g) > 0) ko_t2g <- unique(ko_t2g)
          }
          
          # Parse KEGG Pathway
          pathway_t2g <- NULL
          if ("KEGG_Pathway" %in% colnames(df)) {
            pathway_list <- lapply(1:nrow(df), function(i) {
              pws <- df$KEGG_Pathway[i]
              if (!is.na(pws) && pws != "" && pws != "-") {
                pws_split <- strsplit(pws, ",")[[1]]
                data.frame(term = trimws(pws_split), gene = df$query[i], stringsAsFactors = FALSE)
              } else {
                NULL
              }
            })
            pathway_t2g <- do.call(rbind, pathway_list)
            if (!is.null(pathway_t2g) && nrow(pathway_t2g) > 0) pathway_t2g <- unique(pathway_t2g)
          }
          
          rv$custom_annotations <- list(
            go_t2g = go_t2g,
            ko_t2g = ko_t2g,
            pathway_t2g = pathway_t2g
          )
          
          showNotification("カスタムアノテーションの読み込みに成功しました。", type = "message", duration = 6)
        }, error = function(e) {
          showNotification(paste("カスタムアノテーション読み込みエラー:", e$message), type = "error", duration = 15)
          rv$custom_annotations <- NULL
        })
      })
    })
    
    output$customAnnotationStatusUI <- renderUI({
      ca <- rv$custom_annotations
      if (is.null(ca)) return(NULL)
      go_n <- if(!is.null(ca$go_t2g)) nrow(ca$go_t2g) else 0
      ko_n <- if(!is.null(ca$ko_t2g)) nrow(ca$ko_t2g) else 0
      pw_n <- if(!is.null(ca$pathway_t2g)) nrow(ca$pathway_t2g) else 0
      
      tags$div(class = "alert alert-success", style = "padding: 8px; font-size: 0.85rem;",
        icon("check-circle"),
        HTML(paste0(" カスタムアノテーション適用中: <br>GO: <b>", go_n, "</b> 件, KO: <b>", ko_n, "</b> 件, Pathway: <b>", pw_n, "</b> 件"))
      )
    })

    # --- 方法A: Tax ID から NCBI gene2go を自動取得 (コマンド不要) ---
    observeEvent(input$fetchGene2go, {
      tax <- trimws(input$taxIdInput %||% "")
      if (!grepl("^[0-9]+$", tax)) {
        showNotification("Tax ID を数値で入力してください (例: ミヤコグサ=34305)。", type = "error", duration = 8)
        return(NULL)
      }
      withProgress(message = "GO注釈を取得中...", value = 0.05, {
        tryCatch({
          df <- fetch_gene2go_by_taxid(tax, progress = function(f, m) setProgress(value = f, message = m))
          rv$gene2go_annotation <- df
          showNotification(paste0("GO注釈取得完了 (Tax ID ", tax, "): ", nrow(df), " 注釈 / 遺伝子 ",
                                  length(unique(df$GeneID)), " 件。GOエンリッチメントが可能になりました。"),
                           type = "message", duration = 10)
        }, error = function(e) {
          showNotification(paste("GO注釈の取得に失敗:", conditionMessage(e)), type = "error", duration = 15)
        })
      })
    })

    # --- 方法B: gene2go アノテーション 手動アップロード ---
    # NCBI gene2go 形式のタブ区切り表 (列: GeneID, GO, 任意で Term, Category) を読み込み、
    # rv$gene2go_annotation に格納する。GOモジュールがこれを TERM2GENE として使用する。
    observeEvent(input$gene2goFile, {
      req(input$gene2goFile)
      withProgress(message = "gene2goを解析中...", value = 0.5, {
        tryCatch({
          path <- input$gene2goFile$datapath
          is_gz <- grepl("\\.gz$", input$gene2goFile$name, ignore.case = TRUE)
          con <- if (is_gz) gzfile(path) else path
          # 区切りはタブ優先、ダメならカンマで再試行
          df <- tryCatch(
            read.table(con, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, fill = TRUE, comment.char = ""),
            error = function(e) NULL
          )
          if (is.null(df) || ncol(df) < 2) {
            con2 <- if (is_gz) gzfile(path) else path
            df <- read.csv(con2, header = TRUE, stringsAsFactors = FALSE)
          }
          # 列名を正規化 (GeneID / GO の別名を許容)
          cn <- tolower(colnames(df))
          gid_i <- which(cn %in% c("geneid", "gene_id", "gene", "entrezid", "entrez"))[1]
          go_i  <- which(cn %in% c("go", "go_id", "goid", "term_id"))[1]
          if (is.na(gid_i) || is.na(go_i)) {
            stop("必須列が見つかりません。GeneID 列と GO 列を含むタブ区切り表が必要です。")
          }
          out <- data.frame(GeneID = as.character(df[[gid_i]]), GO = as.character(df[[go_i]]), stringsAsFactors = FALSE)
          term_i <- which(cn %in% c("term", "go_term", "name", "description"))[1]
          cat_i  <- which(cn %in% c("category", "ontology", "namespace", "aspect"))[1]
          if (!is.na(term_i)) out$Term <- as.character(df[[term_i]])
          if (!is.na(cat_i)) {
            # BP/MF/CC や biological_process 等を NCBI表記(Process/Function/Component)へ正規化
            raw <- tolower(as.character(df[[cat_i]]))
            out$Category <- ifelse(grepl("bp|process|biological", raw), "Process",
                            ifelse(grepl("mf|function|molecular", raw), "Function",
                            ifelse(grepl("cc|component|cellular", raw), "Component", as.character(df[[cat_i]]))))
          }
          out <- out[!is.na(out$GeneID) & nzchar(out$GeneID) & !is.na(out$GO) & nzchar(out$GO), , drop = FALSE]
          out <- unique(out)
          if (nrow(out) == 0) stop("有効な GeneID-GO ペアがありません。")
          rv$gene2go_annotation <- out
          showNotification(paste0("gene2go読み込み完了: ", nrow(out), " 注釈 / 遺伝子 ",
                                  length(unique(out$GeneID)), " 件。非モデル生物のGO解析が可能になりました。"),
                           type = "message", duration = 8)
        }, error = function(e) {
          showNotification(paste("gene2go読み込みエラー:", e$message), type = "error", duration = 15)
          rv$gene2go_annotation <- NULL
        })
      })
    })

    output$gene2goStatusUI <- renderUI({
      g <- rv$gene2go_annotation
      if (is.null(g)) return(NULL)
      has_cat <- "Category" %in% colnames(g)
      tags$div(class = "alert alert-success", style = "padding: 8px; font-size: 0.85rem;",
        icon("check-circle"),
        HTML(paste0(" gene2go適用中: <b>", nrow(g), "</b> 注釈 / 遺伝子 <b>",
                    length(unique(g$GeneID)), "</b> 件",
                    if (has_cat) " / ontology別(BP/MF/CC)対応" else " / ontology列なし(全GO一括)"))
      )
    })

    # KEGG生物種コード (任意): GOモジュールが enrichKEGG(organism=...) で使用
    observeEvent(input$keggOrgCode, {
      code <- trimws(input$keggOrgCode %||% "")
      rv$kegg_organism_code <- if (nzchar(code)) code else NULL
    }, ignoreInit = TRUE)


    # 元ID -> 内部Geneid の対応表 (カウント処理時に各経路で格納)
    id_map_store <- reactiveVal(NULL)  # data.frame(orig, Geneid)

    # id_map_store と gtf_raw から rv$gene_annotation を構築し、遺伝子長を上書きする。
    apply_gtf_annotation <- function() {
      gtf <- gtf_raw()
      map_df <- id_map_store()
      if (is.null(gtf) || is.null(map_df)) { rv$gene_annotation <- NULL; return(invisible(NULL)) }
      dt <- data.table::data.table(orig = as.character(map_df$orig), Geneid = as.character(map_df$Geneid))
      g  <- data.table::as.data.table(gtf)
      m  <- merge(dt, g, by.x = "orig", by.y = "gene_id", all.x = TRUE)
      ann <- m[, list(
        gene_name        = { v <- gene_name[!is.na(gene_name) & gene_name != ""]; if (length(v) > 0) v[1] else Geneid[1] },
        biotype          = { v <- biotype[!is.na(biotype)]; if (length(v) > 0) v[1] else "unknown" },
        gene_length      = { v <- gene_length[!is.na(gene_length)]; if (length(v) > 0) max(v) else NA_real_ },
        gene_description = { v <- gene_description[!is.na(gene_description) & gene_description != ""]; if (length(v) > 0) v[1] else Geneid[1] }
      ), by = Geneid]
      rv$gene_annotation <- as.data.frame(ann)
      # GTF由来の遺伝子長が得られた場合は rv$gene_lengths を上書き
      if (any(!is.na(ann$gene_length))) {
        rv$gene_lengths <- data.frame(Geneid = ann$Geneid, Length = ann$gene_length, stringsAsFactors = FALSE)
      }
    }

    # GTFが (カウントデータの後に) アップロード/変更された場合も再構築
    observeEvent(gtf_raw(), {
      if (!is.null(id_map_store())) apply_gtf_annotation()
    }, ignoreNULL = FALSE)

    data_processing_trigger <- reactive({
      list(countFiles = input$countFiles, species = input$species)
    })

    observeEvent(data_processing_trigger(), {
      trigger <- data_processing_trigger()
      files_info <- trigger$countFiles
      species_code <- trigger$species

      if (is.null(files_info)) return()

      inputType <- detect_file_format(files_info)
      if (is.null(inputType)) return()

      message("--- Data Processing Start ---")
      rv$merged_data <- NULL
      rv$sample_metadata <- NULL
      rv$gene_lengths <- NULL
      rv$filtered_keep <- NULL
      rv$deg_results <- NULL
      rv$background_genes_original <- NULL
      rv$conversion_debug_info(NULL)
      rv$selected_species <- species_code

      tryCatch(
        {
          # === A. 個別ファイル ===
          if (inputType == "individual") {
            file_paths <- files_info$datapath
            original_filenames <- files_info$name
            names_no_ext <- sub("\\.[^.]*$", "", original_filenames)
            trimmed_names <- trim_common_strings(names_no_ext)
            initial_sample_names <- make.unique(trimmed_names)

            raw_data_list <- map(file_paths, ~ {
              dt <- fread(.x, header = TRUE, sep = "\t", stringsAsFactors = FALSE, select = c(1, 6, 7), na.strings = c("NA", "NaN", ""))
              setnames(dt, c("OriginalGeneid", "Length", "Count"))
              dt
            })

            # IDクリーニング
            first_file_ids <- trimws(as.character(raw_data_list[[1]]$OriginalGeneid))
            first_file_ids <- gsub('^"|"$', "", first_file_ids) # 引用符除去

            detected_keytype <- detect_gene_id_type(first_file_ids)

            all_original <- map(raw_data_list, ~ {
              ids <- trimws(as.character(.x$OriginalGeneid))
              gsub('^"|"$', "", ids)
            })
            common_ids <- Reduce(intersect, all_original)
            if (length(common_ids) == 0) stop("共通のGeneIDが見つかりません。")
            inferred_species_code <- infer_species_from_gene_ids(common_ids)
            if (!is.null(inferred_species_code) && inferred_species_code != species_code && species_code != "Others_Original") {
              message(paste0("Detected Ensembl species '", inferred_species_code, "'. Overriding selected species '", species_code, "' for ID conversion."))
              species_code <- inferred_species_code
              rv$selected_species <- species_code
              updateSelectInput(session, "species", selected = species_code)
              showNotification("入力Gene IDから生物種を自動補正しました。", type = "message", duration = 6)
            }

            conversion_res <- convert_or_pass_ids(common_ids, detected_keytype, species_code)
            rv$current_gene_id_type(conversion_res$final_id_type)
            rv$conversion_debug_info(conversion_res$conversion_map)

            map_df <- data.frame(OriginalGeneid = common_ids, Geneid = conversion_res$processed_ids, stringsAsFactors = FALSE)
            id_map_store(data.frame(orig = map_df$OriginalGeneid, Geneid = map_df$Geneid, stringsAsFactors = FALSE))

            processed_list <- map2(raw_data_list, all_original, function(dt, clean_ids) {
              dt[, OriginalGeneid := clean_ids] # クリーニング済みIDで上書き
              dt <- dt[OriginalGeneid %in% map_df$OriginalGeneid, ]
              merged <- merge(dt, map_df, by = "OriginalGeneid")
              merged[, OriginalGeneid := NULL]
              merged
            })

            rv$gene_lengths <- as.data.frame(processed_list[[1]][, .(Geneid, Length)][!duplicated(Geneid)])

            # Aggregate counts by Geneid (dplyr ensures unique Geneids per sample)
            count_frames <- map2(processed_list, initial_sample_names, function(dt, sname) {
              df <- as.data.frame(dt)[, c("Geneid", "Count")]
              df %>%
                group_by(Geneid) %>%
                summarise(!!sname := sum(Count, na.rm = TRUE), .groups = "drop") %>%
                as.data.frame()
            })

            # Merge all samples (base merge avoids many-to-many issues)
            merged_df <- Reduce(function(a, b) merge(a, b, by = "Geneid", all = TRUE), count_frames)
            merged_df[is.na(merged_df)] <- 0
            rv$merged_data <- merged_df

            # === B. 結合済みマトリックス ===
          } else if (inputType == "merged") {
            infile <- list(datapath = files_info$datapath[1], name = files_info$name[1])
            ext <- tools::file_ext(infile$name)
            sep_char <- if (tolower(ext) %in% c("tsv", "txt")) "\t" else ","

            all_lines_df <- read.delim(infile$datapath, sep = sep_char, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, quote = "")
            if (nrow(all_lines_df) < 2) stop("ファイルに行が少なすぎます。")

            # ヘッダー処理 (1行目)
            raw_header <- as.character(all_lines_df[1, -1])
            raw_header <- raw_header[!is.na(raw_header) & raw_header != ""]

            # データ行探索
            is_valid_count_row <- function(row_vec) {
              if (length(row_vec) < 2) {
                return(FALSE)
              }
              vals <- as.character(row_vec[-1])
              if (any(is.na(vals)) || any(vals == "")) {
                return(FALSE)
              }
              vals_clean <- gsub(",", "", vals)
              nums <- suppressWarnings(as.numeric(vals_clean))
              return(!any(is.na(nums)))
            }
            data_start_idx <- -1
            for (i in 2:nrow(all_lines_df)) {
              if (is_valid_count_row(all_lines_df[i, ])) {
                data_start_idx <- i
                break
              }
            }
            if (data_start_idx == -1) stop("有効な数値データ行が見つかりませんでした。")

            candidate_data <- all_lines_df[data_start_idx:nrow(all_lines_df), ]
            keep_rows_idx <- apply(candidate_data, 1, is_valid_count_row)
            final_data <- candidate_data[keep_rows_idx, , drop = FALSE]
            if (nrow(final_data) == 0) stop("データ行が残りませんでした。")

            # ★★★ IDクリーニング強化 ★★★
            raw_ids <- as.character(final_data[, 1])
            clean_ids <- trimws(raw_ids) # 空白除去
            clean_ids <- gsub('^"|"$', "", clean_ids) # 引用符除去
            clean_ids <- gsub("^'|'$", "", clean_ids) # シングルクォート除去

            counts_part <- final_data[, -1, drop = FALSE]
            counts_mat <- apply(counts_part, 2, function(x) as.numeric(gsub(",", "", x)))

            if (length(clean_ids) != nrow(counts_mat)) stop("ID列とデータ列の行数が不一致。")

            clean_mat <- as.matrix(counts_mat)
            rownames(clean_mat) <- clean_ids

            if (ncol(clean_mat) > length(raw_header)) {
              extra_cols <- (length(raw_header) + 1):ncol(clean_mat)
              raw_header <- c(raw_header, paste0("Sample_", extra_cols))
            } else if (ncol(clean_mat) < length(raw_header)) {
              raw_header <- raw_header[1:ncol(clean_mat)]
            }
            colnames(clean_mat) <- raw_header
            initial_sample_names <- raw_header

            original_ids <- rownames(clean_mat)
            detected_keytype <- detect_gene_id_type(original_ids)
            inferred_species_code <- infer_species_from_gene_ids(original_ids)
            if (!is.null(inferred_species_code) && inferred_species_code != species_code && species_code != "Others_Original") {
              message(paste0("Detected Ensembl species '", inferred_species_code, "'. Overriding selected species '", species_code, "' for ID conversion."))
              species_code <- inferred_species_code
              rv$selected_species <- species_code
              updateSelectInput(session, "species", selected = species_code)
              showNotification("入力Gene IDから生物種を自動補正しました。", type = "message", duration = 6)
            }

            if (species_code != "Others_Original") {
              message(paste0("Matrix Input: Detected ID type '", detected_keytype, "'. Converting to Entrez..."))
              conversion_res <- convert_or_pass_ids(original_ids, detected_keytype, species_code)
              rv$current_gene_id_type(conversion_res$final_id_type)
              rv$conversion_debug_info(conversion_res$conversion_map)

              map_df <- data.frame(Original = original_ids, New = conversion_res$processed_ids, stringsAsFactors = FALSE)
              id_map_store(data.frame(orig = map_df$Original, Geneid = map_df$New, stringsAsFactors = FALSE))
              clean_mat_subset <- clean_mat[map_df$Original, , drop = FALSE]
              df_for_agg <- as.data.frame(clean_mat_subset)
              df_for_agg$Geneid <- map_df$New

              agg_df <- df_for_agg %>%
                group_by(Geneid) %>%
                summarise(across(everything(), sum)) %>%
                ungroup()
              rv$merged_data <- as.data.frame(agg_df)
            } else {
              rv$current_gene_id_type(detected_keytype %||% "Original")
              rv$conversion_debug_info(data.frame(Original = original_ids, Converted = "Skipped", Status = "Skipped"))
              df_res <- as.data.frame(clean_mat)
              df_res$Geneid <- rownames(clean_mat)
              df_res <- df_res[, c("Geneid", setdiff(colnames(df_res), "Geneid"))]
              rv$merged_data <- df_res
              # 変換なし: 元IDと内部Geneidは同一
              id_map_store(data.frame(orig = df_res$Geneid, Geneid = df_res$Geneid, stringsAsFactors = FALSE))
            }
            rv$gene_lengths <- data.frame(Geneid = rv$merged_data$Geneid, Length = 1, stringsAsFactors = FALSE)
            showNotification(paste0("読み込み完了。データ開始:", data_start_idx, "行目"), type = "message")
          }

          # GTFアノテーションが在れば内部Geneid基準で結合し、遺伝子長を上書き
          apply_gtf_annotation()

          # === 共通: サンプルメタデータ ===
          current_samples <- setdiff(colnames(rv$merged_data), "Geneid")
          rv$sample_metadata <- data.frame(
            id = paste0("sample_", 1:length(current_samples)),
            current_name = current_samples,
            group = rep("GroupA", length(current_samples)),
            active = rep(TRUE, length(current_samples)),
            stringsAsFactors = FALSE
          )
        },
        error = function(e) {
          showNotification(paste("エラー:", e$message), type = "error", duration = 15)
          rv$merged_data <- NULL
          rv$sample_metadata <- NULL
        }
      )
    })

    # --- 検出フォーマットバッジ ---
    output$detectedFormatUI <- renderUI({
      files_info <- input$countFiles
      if (is.null(files_info)) return(NULL)
      fmt <- detect_file_format(files_info)
      n <- nrow(files_info)
      if (fmt == "individual") {
        div(class = "alert alert-info mt-2 mb-0 py-2",
            tags$strong("Detected:"),
            paste0(" featureCounts individual files (", n, " file", if(n>1)"s" else "", ")"))
      } else {
        div(class = "alert alert-success mt-2 mb-0 py-2",
            tags$strong("Detected:"),
            paste0(" Merged count matrix — ", files_info$name[1]))
      }
    })

    # --- 列の追加・削除のハンドラ ---
    observeEvent(input$add_factor_btn, {
      req(rv$sample_metadata)
      df <- rv$sample_metadata

      base_names <- c("group2", "group3", "group4", "group5", "group6")
      existing_cols <- colnames(df)
      new_col <- NA
      for (nm in base_names) {
        if (!(nm %in% existing_cols)) {
          new_col <- nm
          break
        }
      }
      if (is.na(new_col)) {
        new_col <- paste0("group", length(existing_cols)) # fallback
      }

      df[[new_col]] <- rep("Condition", nrow(df))
      rv$sample_metadata <- df
    })

    observeEvent(input$remove_factor_btn, {
      req(rv$sample_metadata)
      df <- rv$sample_metadata
      cols <- colnames(df)
      factor_cols <- setdiff(cols, c("id", "current_name", "active", "group", "time"))
      if (length(factor_cols) > 0) {
        last_factor <- factor_cols[length(factor_cols)]
        df[[last_factor]] <- NULL
        rv$sample_metadata <- df
      } else {
        showNotification("削除できる追加グループがありません。", type = "warning", duration = 3)
      }
    })

    observeEvent(input$rename_group_btn, {
      req(rv$sample_metadata)
      meta <- rv$sample_metadata
      factor_cols <- setdiff(colnames(meta), c("id", "current_name", "active", "time"))

      showModal(modalDialog(
        title = "グループ列名の変更",
        lapply(factor_cols, function(fc) {
          textInput(ns(paste0("rename_input_", fc)), paste("現在の列名:", fc), value = fc)
        }),
        easyClose = TRUE,
        footer = tagList(
          modalButton("キャンセル"),
          actionButton(ns("confirm_rename_btn"), "変更を保存", class = "btn-primary")
        )
      ))
    })

    observeEvent(input$confirm_rename_btn, {
      req(rv$sample_metadata)
      meta <- rv$sample_metadata
      factor_cols <- setdiff(colnames(meta), c("id", "current_name", "active", "time"))

      new_names <- sapply(factor_cols, function(fc) {
        val <- input[[paste0("rename_input_", fc)]]
        if (is.null(val) || nchar(trimws(val)) == 0) {
          return(fc)
        }
        make.names(trimws(val)) # Use safe variable names
      })

      new_names <- make.unique(as.character(new_names))

      for (i in seq_along(factor_cols)) {
        colnames(meta)[names(meta) == factor_cols[i]] <- new_names[i]
      }

      rv$sample_metadata <- meta
      removeModal()
      showNotification("列名を変更しました。", type = "message")
    })

    # --- UI生成 ---
    output$sampleMetadataTable <- renderRHandsontable({
      if (is.null(rv$sample_metadata)) {
        return(NULL)
      }

      # UI表示用に列順序を調整
      # id列は内部管理用とし、編集不可にする
      # 存在する列だけを抽出して並べ替え
      all_cols <- colnames(rv$sample_metadata)
      base_cols <- c("active", "current_name", "group")
      factor_cols <- setdiff(all_cols, c(base_cols, "id", "time"))
      cols_to_show <- c(base_cols, factor_cols)

      # time列があれば末尾に追加
      if ("time" %in% all_cols) {
        cols_to_show <- c(base_cols, factor_cols, "time")
      }

      df_display <- rv$sample_metadata[, cols_to_show, drop = FALSE]

      rhp <- rhandsontable(df_display, rowHeaders = NULL, stretchH = "all") %>%
        hot_col("active", halign = "center", type = "checkbox") %>%
        hot_col("current_name", header = "Sample Name") %>%
        hot_col("group", header = "Group") %>%
        hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)

      if ("time" %in% all_cols) {
        rhp <- rhp %>% hot_col("time", header = "Time", type = "numeric")
      }

      for (col in factor_cols) {
        rhp <- rhp %>% hot_col(col, header = col)
      }

      return(rhp)
    })

    # --- サンプル情報更新監視 ---
    observeEvent(input$sampleMetadataTable, {
      req(input$sampleMetadataTable, rv$merged_data)
      new_df <- hot_to_r(input$sampleMetadataTable)

      # current_name の重複チェック
      if (any(duplicated(new_df$current_name))) {
        new_df$current_name <- make.unique(new_df$current_name)
        showNotification("Duplicate sample names detected and fixed.", type = "warning")
      }

      current_meta <- rv$sample_metadata

      # 非表示列 (id 等) を current_meta から復元
      hidden_cols <- setdiff(colnames(current_meta), colnames(new_df))
      for (col in hidden_cols) {
        new_df[[col]] <- current_meta[[col]]
      }

      # 変更検知 (列ごとに比較)
      is_changed <- FALSE
      shared_cols <- intersect(colnames(current_meta), colnames(new_df))
      for (col in shared_cols) {
        if (!identical(current_meta[[col]], new_df[[col]])) {
          is_changed <- TRUE
          break
        }
      }

      if (is_changed) {
        # merged_dataの列名更新 (current_nameが変わった場合)
        if (!identical(current_meta$current_name, new_df$current_name)) {
          old_names <- current_meta$current_name
          new_names <- new_df$current_name
          current_cols <- colnames(rv$merged_data)
          matched_idx <- match(old_names, current_cols)

          if (any(is.na(matched_idx))) {
            warning("[Metadata Update] Column name mismatch detected.")
          } else {
            colnames(rv$merged_data)[matched_idx] <- new_names
          }
        }
        # rv$sample_metadata の列順序に合わせて保存
        rv$sample_metadata <- new_df[, colnames(rv$sample_metadata)]
      }
    })

    # --- デバッグ情報出力 ---
    output$idConversionSummary <- renderPrint({
      req(rv$conversion_debug_info())
      df <- rv$conversion_debug_info()
      if (is.null(df) || nrow(df) == 0) {
        return("変換情報なし")
      }

      n_total <- nrow(df)
      n_success <- sum(df$Status == "Success")
      n_failed <- sum(df$Status == "Failed")

      cat("Total IDs:", n_total, "\n")
      cat("Success:", n_success, "(", round(n_success / n_total * 100, 1), "%)\n")
      cat("Failed:", n_failed, "(", round(n_failed / n_total * 100, 1), "%)\n\n")

      if (n_failed > 0) {
        cat("--- Failed IDs (Sample) ---\n")
        print(head(df[df$Status == "Failed", "Original"], 10))
        cat("\n(ヒント: これらのIDが正しいGene Symbolか、空白が含まれていないか確認してください。)\n")
      }
    })

    output$idConversionTable <- renderDT({
      req(rv$conversion_debug_info())
      datatable(rv$conversion_debug_info(), style = "bootstrap5", class = "table-hover table-sm", options = list(pageLength = 10, scrollX = TRUE), filter = "top")
    })

    # --- プロットなど ---
    output$lib_size_color_ui <- renderUI({
      req(rv$sample_metadata)
      meta <- rv$sample_metadata
      factor_cols <- setdiff(colnames(meta), c("id", "current_name", "active", "time"))
      if (length(factor_cols) == 0) {
        return(NULL)
      }
      selectInput(ns("lib_size_color"), "色分けするグループ:", choices = factor_cols, selected = factor_cols[1])
    })

    plot_data <- reactive({
      req(rv$merged_data, rv$sample_metadata)
      act_meta <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      shiny::validate(shiny::need(nrow(act_meta) > 0, "プロット対象なし"))
      cols <- intersect(act_meta$current_name, colnames(rv$merged_data))
      if (length(cols) == 0) {
        return(NULL)
      }
      num_dat <- rv$merged_data[, cols, drop = FALSE]
      if (!all(sapply(num_dat, is.numeric))) {
        return(NULL)
      }
      tot <- colSums(num_dat, na.rm = TRUE)

      color_col <- input$lib_size_color
      if (is.null(color_col) || !(color_col %in% colnames(act_meta))) {
        factor_cols <- setdiff(colnames(act_meta), c("id", "current_name", "active", "time"))
        color_col <- if (length(factor_cols) > 0) factor_cols[1] else NULL
      }

      df <- data.frame(current_name = names(tot), TotalReads = tot, stringsAsFactors = FALSE)
      if (!is.null(color_col)) {
        df <- left_join(df, act_meta[, c("current_name", color_col)], by = "current_name")
        df$color_var <- df[[color_col]]
      } else {
        df$color_var <- "All"
      }
      df$current_name <- factor(df$current_name, levels = act_meta$current_name)
      df
    })

    output$librarySizePlot <- renderPlotly({
      d <- plot_data()
      req(d, nrow(d) > 0)
      plot_ly(d, x = ~current_name, y = ~TotalReads, color = ~color_var, type = "bar") %>% layout(margin = list(b = 100))
    })
    output$downloadRDS <- downloadHandler(filename = function() {
      paste0("session_", Sys.Date(), ".rds")
    }, content = function(f) {
      saveRDS(reactiveValuesToList(rv), f)
    })
    observeEvent(input$uploadRDS, {
      req(input$uploadRDS)
      d <- readRDS(input$uploadRDS$datapath)
      for (n in names(d)) if (n %in% names(rv)) rv[[n]] <- d[[n]]
    })

    `%||%` <- function(a, b) if (!is.null(a)) a else b
  })
}
