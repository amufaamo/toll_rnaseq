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
  ids_clean <- na.omit(ids[ids != "" & !is.na(ids)]); if (length(ids_clean) == 0) return("UNKNOWN")
  n_total <- length(ids_clean); n_sample <- min(n_total, 1000); ids_sample <- sample(ids_clean, n_sample)
  if (mean(grepl("^ENS[A-Z0-9]*[FPTG]\\d{10,}(\\.\\d+)?$", ids_sample, ignore.case = TRUE)) > 0.8) return("ENSEMBL")
  if (mean(grepl("^[NX][CMPWTZ]_[0-9]+(\\.\\d+)?$", ids_sample, ignore.case = TRUE)) > 0.8) return("REFSEQ")
  if (mean(grepl("^[0-9]+$", ids_sample)) > 0.9) return("ENTREZID")
  if (mean(grepl("^([A-Za-z][A-Za-z0-9-]*[A-Za-z0-9]|[A-Za-z])$", ids_sample)) > 0.7 && mean(grepl("^ENS",ids_sample,T))<0.1 && mean(grepl("^[NX][M_]",ids_sample,T))<0.1 && mean(grepl("^[0-9]+$",ids_sample))<0.1) return("SYMBOL")
  return("UNKNOWN")
}

# ★★★ 修正: エイリアス検索を追加した関数 ★★★
convert_or_pass_ids <- function(original_ids, detected_keytype, selected_species_code) {
  if (selected_species_code == "Others_Original") {
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(FALSE, length(original_ids)), failed_count = 0, final_id_type = detected_keytype %||% "ORIGINAL_UNCONVERTED", conversion_map = NULL))
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
  if (detected_keytype == "ENSEMBL") { keys_to_convert <- gsub("\\..*$", "", keys_to_convert) }
  
  if (!detected_keytype %in% keytypes(org_db)) {
    warning(paste0("Keytype '", detected_keytype, "' not supported by '", orgdb_pkg_name, "'."))
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(TRUE, length(original_ids)), failed_count = length(original_ids), final_id_type = detected_keytype %||% "UNKNOWN", conversion_map = NULL))
  }
  
  message(paste0("[ID_CONV] 1st pass: Converting ", length(unique(keys_to_convert)), " unique '", detected_keytype, "' IDs..."))
  unique_input_keys <- unique(keys_to_convert)
  
  # 1. まず通常の変換（Symbol -> Entrez）
  entrez_map <- tryCatch(
    suppressMessages(mapIds(org_db, keys = unique_input_keys, column = "ENTREZID", keytype = detected_keytype, multiVals = "first")),
    error = function(e) { warning(e$message); NULL }
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
      error = function(e) { NULL }
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

# --- Helper function: Trim common strings ---
# (変更なし)
find_lcp <- function(strs) {
  if (length(strs) < 2) return(""); char_lists <- strsplit(strs, ""); min_len <- min(sapply(char_lists, length)); if (min_len == 0) return("")
  lcp <- ""; for (i in 1:min_len) { char_to_check <- char_lists[[1]][i]; if (all(sapply(char_lists, function(x) x[i] == char_to_check))) { lcp <- paste0(lcp, char_to_check) } else { break } }; return(lcp)
}
find_lcs <- function(strs) {
  if (length(strs) < 2) return(""); rev_strs <- sapply(strs, function(x) intToUtf8(rev(utf8ToInt(x)))); rev_lcs <- find_lcp(rev_strs); lcs <- intToUtf8(rev(utf8ToInt(rev_lcs))); return(lcs)
}
trim_common_strings <- function(strs) {
  if (length(strs) < 2) return(strs); lcp <- find_lcp(strs); lcs <- find_lcs(strs); trimmed_strs <- sub(paste0("^", lcp), "", strs); trimmed_strs <- sub(paste0(lcs, "$"), "", trimmed_strs); return(trimmed_strs)
}

dataUploadMetadataUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("1. ファイルアップロード"),
      radioButtons(ns("inputType"), "データの形式:",
                   choices = c("個別の featureCounts ファイル (複数可)" = "individual",
                               "結合済みカウント行列 (CSV/TSV)" = "merged"),
                   selected = "individual"),
      
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
      h4("2. 生物種選択"),
      selectInput(ns("species"), label = "解析対象の生物種またはID処理方法:", choices = species_choices_ui, selected = "Homo_sapiens"),
      hr(),
      h4("3. セッション管理"),
      downloadButton(ns("downloadRDS"), "セッション保存"),
      br(),br(),
      fileInput(ns("uploadRDS"), "セッション復元", accept = c(".rds", ".RDS"))
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("サンプル情報 & QC",
                 h4("サンプル別 合計リード数"),
                 withSpinner(plotlyOutput(ns("librarySizePlot")), type = 6),
                 hr(),
                 h4("サンプル情報入力"),
                 withSpinner(uiOutput(ns("sampleMetadataInput")), type = 6)
        ),
        tabPanel("ID変換レポート (デバッグ用)",
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
    
    data_processing_trigger <- reactive({
      list(inputType = input$inputType, individualFiles = input$featureCountsFiles, mergedFile = input$mergedCountFile, species = input$species)
    })
    
    observeEvent(data_processing_trigger(), {
      trigger <- data_processing_trigger()
      inputType <- trigger$inputType
      species_code <- trigger$species
      
      if (inputType == "individual" && is.null(trigger$individualFiles)) return()
      if (inputType == "merged" && is.null(trigger$mergedFile)) return()
      
      message("--- Data Processing Start ---")
      rv$merged_data <- NULL; rv$sample_metadata <- NULL; rv$gene_lengths <- NULL
      rv$filtered_keep <- NULL; rv$deg_results <- NULL; rv$background_genes_original <- NULL
      rv$conversion_debug_info(NULL)
      rv$selected_species <- species_code
      
      tryCatch({
        # === A. 個別ファイル ===
        if (inputType == "individual") {
          file_paths <- trigger$individualFiles$datapath
          original_filenames <- trigger$individualFiles$name
          names_no_ext <- sub("\\.[^.]*$", "", original_filenames)
          trimmed_names <- trim_common_strings(names_no_ext)
          initial_sample_names <- make.unique(trimmed_names)
          
          raw_data_list <- map(file_paths, ~{
            dt <- fread(.x, header = TRUE, sep = "\t", stringsAsFactors = FALSE, select = c(1, 6, 7), na.strings = c("NA", "NaN", ""))
            setnames(dt, c("OriginalGeneid", "Length", "Count"))
            dt
          })
          
          # IDクリーニング
          first_file_ids <- trimws(as.character(raw_data_list[[1]]$OriginalGeneid))
          first_file_ids <- gsub('^"|"$', '', first_file_ids) # 引用符除去
          
          detected_keytype <- detect_gene_id_type(first_file_ids)
          
          all_original <- map(raw_data_list, ~ {
            ids <- trimws(as.character(.x$OriginalGeneid))
            gsub('^"|"$', '', ids)
          })
          common_ids <- Reduce(intersect, all_original)
          if(length(common_ids) == 0) stop("共通のGeneIDが見つかりません。")
          
          conversion_res <- convert_or_pass_ids(common_ids, detected_keytype, species_code)
          rv$current_gene_id_type(conversion_res$final_id_type)
          rv$conversion_debug_info(conversion_res$conversion_map)
          
          map_df <- data.frame(OriginalGeneid = common_ids, Geneid = conversion_res$processed_ids, stringsAsFactors = FALSE)
          
          processed_list <- map2(raw_data_list, all_original, function(dt, clean_ids) {
            dt[, OriginalGeneid := clean_ids] # クリーニング済みIDで上書き
            dt <- dt[OriginalGeneid %in% map_df$OriginalGeneid, ]
            merged <- merge(dt, map_df, by="OriginalGeneid")
            merged[, OriginalGeneid := NULL]
            merged
          })
          
          rv$gene_lengths <- as.data.frame(processed_list[[1]][, .(Geneid, Length)][!duplicated(Geneid)])
          
          count_list <- map(processed_list, ~ .x[, .(Geneid, Count)])
          merged_dt <- count_list[[1]]
          setnames(merged_dt, "Count", initial_sample_names[1])
          
          if(length(count_list) > 1){
            merged_df <- as.data.frame(merged_dt)
            for(i in 2:length(count_list)){
              curr <- count_list[[i]]
              setnames(curr, "Count", initial_sample_names[i])
              merged_df <- full_join(merged_df, as.data.frame(curr), by="Geneid")
            }
            merged_dt <- as.data.table(merged_df)
          }
          for(col in setdiff(names(merged_dt), "Geneid")) set(merged_dt, which(is.na(merged_dt[[col]])), col, 0)
          
          if(any(duplicated(merged_dt$Geneid))){
            merged_dt <- merged_dt[, lapply(.SD, sum, na.rm=TRUE), by=Geneid, .SDcols=setdiff(names(merged_dt), "Geneid")]
          }
          rv$merged_data <- as.data.frame(merged_dt)
          
        # === B. 結合済みマトリックス ===
        } else if (inputType == "merged") {
          req(trigger$mergedFile)
          infile <- trigger$mergedFile
          ext <- tools::file_ext(infile$name)
          sep_char <- if(tolower(ext) %in% c("tsv", "txt")) "\t" else ","
          
          all_lines_df <- read.delim(infile$datapath, sep = sep_char, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, quote="")
          if(nrow(all_lines_df) < 2) stop("ファイルに行が少なすぎます。")
          
          # ヘッダー処理 (1行目)
          raw_header <- as.character(all_lines_df[1, -1])
          raw_header <- raw_header[!is.na(raw_header) & raw_header != ""]
          
          # データ行探索
          is_valid_count_row <- function(row_vec) {
             if(length(row_vec) < 2) return(FALSE)
             vals <- as.character(row_vec[-1])
             if(any(is.na(vals)) || any(vals == "")) return(FALSE)
             vals_clean <- gsub(",", "", vals)
             nums <- suppressWarnings(as.numeric(vals_clean))
             return(!any(is.na(nums)))
          }
          data_start_idx <- -1
          for(i in 2:nrow(all_lines_df)) {
             if(is_valid_count_row(all_lines_df[i, ])) { data_start_idx <- i; break }
          }
          if(data_start_idx == -1) stop("有効な数値データ行が見つかりませんでした。")
          
          candidate_data <- all_lines_df[data_start_idx:nrow(all_lines_df), ]
          keep_rows_idx <- apply(candidate_data, 1, is_valid_count_row)
          final_data <- candidate_data[keep_rows_idx, , drop = FALSE]
          if(nrow(final_data) == 0) stop("データ行が残りませんでした。")
          
          # ★★★ IDクリーニング強化 ★★★
          raw_ids <- as.character(final_data[, 1])
          clean_ids <- trimws(raw_ids)          # 空白除去
          clean_ids <- gsub('^"|"$', '', clean_ids) # 引用符除去
          clean_ids <- gsub("^'|'$", "", clean_ids) # シングルクォート除去
          
          counts_part <- final_data[, -1, drop = FALSE]
          counts_mat <- apply(counts_part, 2, function(x) as.numeric(gsub(",", "", x)))
          
          if(length(clean_ids) != nrow(counts_mat)) stop("ID列とデータ列の行数が不一致。")
          
          clean_mat <- as.matrix(counts_mat)
          rownames(clean_mat) <- clean_ids
          
          if(ncol(clean_mat) > length(raw_header)) {
             extra_cols <- (length(raw_header)+1):ncol(clean_mat)
             raw_header <- c(raw_header, paste0("Sample_", extra_cols))
          } else if(ncol(clean_mat) < length(raw_header)) {
             raw_header <- raw_header[1:ncol(clean_mat)]
          }
          colnames(clean_mat) <- raw_header
          initial_sample_names <- raw_header
          
          original_ids <- rownames(clean_mat)
          detected_keytype <- detect_gene_id_type(original_ids)
          
          if (species_code != "Others_Original") {
             message(paste0("Matrix Input: Detected ID type '", detected_keytype, "'. Converting to Entrez..."))
             conversion_res <- convert_or_pass_ids(original_ids, detected_keytype, species_code)
             rv$current_gene_id_type(conversion_res$final_id_type)
             rv$conversion_debug_info(conversion_res$conversion_map)
             
             map_df <- data.frame(Original = original_ids, New = conversion_res$processed_ids, stringsAsFactors = FALSE)
             clean_mat_subset <- clean_mat[map_df$Original, , drop=FALSE]
             df_for_agg <- as.data.frame(clean_mat_subset)
             df_for_agg$Geneid <- map_df$New
             
             agg_df <- df_for_agg %>% group_by(Geneid) %>% summarise(across(everything(), sum)) %>% ungroup()
             rv$merged_data <- as.data.frame(agg_df)
             
          } else {
            rv$current_gene_id_type(detected_keytype %||% "Original")
            rv$conversion_debug_info(data.frame(Original=original_ids, Converted="Skipped", Status="Skipped"))
            df_res <- as.data.frame(clean_mat)
            df_res$Geneid <- rownames(clean_mat)
            df_res <- df_res[, c("Geneid", setdiff(colnames(df_res), "Geneid"))]
            rv$merged_data <- df_res
          }
          rv$gene_lengths <- data.frame(Geneid = rv$merged_data$Geneid, Length = 1, stringsAsFactors = FALSE)
          showNotification(paste0("読み込み完了。データ開始:", data_start_idx, "行目"), type = "message")
        }
        
        # === 共通: サンプルメタデータ ===
        current_samples <- setdiff(colnames(rv$merged_data), "Geneid")
        rv$sample_metadata <- data.frame(
          id = paste0("sample_", 1:length(current_samples)),
          current_name = current_samples,
          group = rep("GroupA", length(current_samples)),
          time = rep(0, length(current_samples)),
          active = rep(TRUE, length(current_samples)),
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        showNotification(paste("エラー:", e$message), type = "error", duration = 15)
        rv$merged_data <- NULL; rv$sample_metadata <- NULL
      })
    })
    
    # --- UI生成 ---
    output$sampleMetadataInput <- renderUI({
      if (is.null(rv$sample_metadata)) { return(tags$p("ファイルをアップロードしてください。")) }
      input_tags <- map(1:nrow(rv$sample_metadata), ~{
        sample_info <- rv$sample_metadata[.x, ]
        tagList(fluidRow(
            column(2, checkboxInput(inputId = ns(paste0("active_sample_", sample_info$id)), label = "解析対象", value = sample_info$active)),
            column(4, textInput(inputId = ns(paste0("sample_name_", sample_info$id)), label = paste0("サンプル ", .x, " 名前"), value = sample_info$current_name)),
            column(3, textInput(inputId = ns(paste0("group_name_", sample_info$id)), label = "グループ", value = sample_info$group)),
            column(3, numericInput(inputId = ns(paste0("time_", sample_info$id)), label = "Time", value = sample_info$time %||% 0, min = 0, step = 1))
        ))
      })
      do.call(tagList, input_tags)
    })
    
    # --- サンプル情報更新監視 ---
    observe({
      all_inputs <- reactive({
        req(rv$sample_metadata)
        lapply(1:nrow(rv$sample_metadata), function(i) {
          sid <- rv$sample_metadata$id[i]
          list(input[[paste0("sample_name_", sid)]], input[[paste0("group_name_", sid)]], input[[paste0("time_", sid)]], input[[paste0("active_sample_", sid)]])
        })
      })
      debounced_inputs <- debounce(all_inputs, 500)
      observeEvent(debounced_inputs(), {
        req(rv$sample_metadata, rv$merged_data)
        n_s <- nrow(rv$sample_metadata); new_names <- character(n_s); new_grps <- character(n_s); new_times <- numeric(n_s); new_acts <- logical(n_s); valid <- TRUE
        isolate({
          vals <- reactiveValuesToList(input)
          for(i in 1:n_s){
            sid <- rv$sample_metadata$id[i]
            nm <- vals[[paste0("sample_name_", sid)]]; grp <- vals[[paste0("group_name_", sid)]]; tm <- vals[[paste0("time_", sid)]]; act <- vals[[paste0("active_sample_", sid)]]
            if(is.null(nm)||is.null(grp)||is.null(tm)||is.null(act)||!is.numeric(tm)) { valid <- FALSE; break }
            new_names[i] <- nm; new_grps[i] <- grp; new_times[i] <- tm; new_acts[i] <- act
          }
        })
        if(valid){
          if(!identical(new_names, rv$sample_metadata$current_name) || !identical(new_grps, rv$sample_metadata$group) || !identical(new_times, rv$sample_metadata$time) || !identical(new_acts, rv$sample_metadata$active)){
            unq_nms <- make.unique(new_names)
            if(any(unq_nms != new_names)){ showNotification("サンプル名重複修正", type="warning"); new_names <- unq_nms }
            if(!is.null(rv$merged_data) && (ncol(rv$merged_data)-1 == length(new_names))){
              rv$sample_metadata$current_name <- new_names; rv$sample_metadata$group <- new_grps; rv$sample_metadata$time <- new_times; rv$sample_metadata$active <- new_acts
              colnames(rv$merged_data) <- c("Geneid", new_names)
            }
          }
        }
      })
    })
    
    # --- デバッグ情報出力 ---
    output$idConversionSummary <- renderPrint({
      req(rv$conversion_debug_info)
      df <- rv$conversion_debug_info()
      if (is.null(df) || nrow(df) == 0) return("変換情報なし")
      
      n_total <- nrow(df)
      n_success <- sum(df$Status == "Success")
      n_failed <- sum(df$Status == "Failed")
      
      cat("Total IDs:", n_total, "\n")
      cat("Success:", n_success, "(", round(n_success/n_total*100, 1), "%)\n")
      cat("Failed:", n_failed, "(", round(n_failed/n_total*100, 1), "%)\n\n")
      
      if(n_failed > 0) {
        cat("--- Failed IDs (Sample) ---\n")
        print(head(df[df$Status == "Failed", "Original"], 10))
        cat("\n(ヒント: これらのIDが正しいGene Symbolか、空白が含まれていないか確認してください。)\n")
      }
    })
    
    output$idConversionTable <- renderDT({
      req(rv$conversion_debug_info)
      datatable(rv$conversion_debug_info(), options = list(pageLength = 10, scrollX = TRUE), filter="top")
    })
    
    # --- プロットなど ---
    plot_data <- reactive({
      req(rv$merged_data, rv$sample_metadata)
      act_meta <- rv$sample_metadata[rv$sample_metadata$active, , drop=FALSE]
      validate(need(nrow(act_meta) > 0, "プロット対象なし"))
      cols <- intersect(act_meta$current_name, colnames(rv$merged_data))
      if(length(cols) == 0) return(NULL)
      num_dat <- rv$merged_data[, cols, drop=FALSE]
      if(!all(sapply(num_dat, is.numeric))) return(NULL)
      tot <- colSums(num_dat, na.rm=TRUE)
      df <- data.frame(current_name = names(tot), TotalReads = tot, stringsAsFactors=FALSE)
      df <- left_join(df, act_meta[, c("current_name", "group")], by="current_name")
      df$current_name <- factor(df$current_name, levels = act_meta$current_name)
      df
    })
    output$librarySizePlot <- renderPlotly({
      d <- plot_data(); req(d, nrow(d)>0)
      plot_ly(d, x=~current_name, y=~TotalReads, color=~group, type='bar') %>% layout(margin=list(b=100))
    })
    output$downloadRDS <- downloadHandler(filename=function(){paste0("session_", Sys.Date(), ".rds")}, content=function(f){saveRDS(reactiveValuesToList(rv), f)})
    observeEvent(input$uploadRDS, { req(input$uploadRDS); d <- readRDS(input$uploadRDS$datapath); for(n in names(d)) if(n %in% names(rv)) rv[[n]] <- d[[n]] })
    
    `%||%` <- function(a, b) if (!is.null(a)) a else b
  })
}
