# R/module_data_upload_metadata.R
# (セッション管理、サンプル選択、Time入力、および結合済みマトリックスのアップロード[上下メタデータ除去強化版])

library(shiny)
library(data.table)
library(purrr)
library(dplyr)
library(shinycssloaders)
library(plotly)
library(AnnotationDbi)

# --- Species UI Choices ---
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

convert_or_pass_ids <- function(original_ids, detected_keytype, selected_species_code) {
  if (selected_species_code == "Others_Original") {
    message("[ID_PROC_UPLOAD] Species is 'Others_Original'. Using original IDs directly.")
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(FALSE, length(original_ids)), failed_count = 0, final_id_type = detected_keytype %||% "ORIGINAL_UNCONVERTED"))
  }
  if (detected_keytype == "ENTREZID") {
    message("[ID_CONV_UPLOAD] Gene IDs are already EntrezID.")
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(FALSE, length(original_ids)), failed_count = 0, final_id_type = "ENTREZID"))
  }
  if (detected_keytype == "UNKNOWN" || is.null(selected_species_code) || !nzchar(selected_species_code)) {
    warning("[ID_CONV_UPLOAD] Cannot convert to EntrezID: Input type is UNKNOWN or species not selected for conversion.")
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(TRUE, length(original_ids)), failed_count = length(original_ids), final_id_type = detected_keytype %||% "UNKNOWN"))
  }
  orgdb_pkg_name <- orgdb_species_map[[selected_species_code]]
  if (is.null(orgdb_pkg_name) || !nzchar(orgdb_pkg_name)) {
    warning(paste0("[ID_CONV_UPLOAD] No OrgDb package found for species '", selected_species_code, "'. Cannot convert to EntrezID."))
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(TRUE, length(original_ids)), failed_count = length(original_ids), final_id_type = detected_keytype %||% "UNKNOWN"))
  }
  if (!requireNamespace(orgdb_pkg_name, quietly = TRUE)) { stop(paste0("[ID_CONV_UPLOAD] OrgDb package '", orgdb_pkg_name, "' not installed. Please install it. Cannot convert to EntrezID.")) }
  require(orgdb_pkg_name, character.only = TRUE, quietly = TRUE); org_db <- get(orgdb_pkg_name)
  keys_to_convert <- as.character(original_ids); if (detected_keytype == "ENSEMBL") { keys_to_convert <- gsub("\\..*$", "", keys_to_convert) }
  
  if (!detected_keytype %in% keytypes(org_db)) { 
    warning(paste0("[ID_CONV_UPLOAD] Keytype '", detected_keytype, "' not supported by '", orgdb_pkg_name, "'. Cannot convert."))
    return(list(processed_ids = as.character(original_ids), failed_indices = rep(TRUE, length(original_ids)), failed_count = length(original_ids), final_id_type = detected_keytype %||% "UNKNOWN")) 
  }
  
  message(paste0("[ID_CONV_UPLOAD] Converting ", length(unique(keys_to_convert)), " unique '", detected_keytype, "' IDs to ENTREZID using '", orgdb_pkg_name, "'."))
  unique_input_keys <- unique(keys_to_convert)
  entrez_map <- tryCatch(suppressMessages(mapIds(org_db, keys = unique_input_keys, column = "ENTREZID", keytype = detected_keytype, multiVals = "first")), error = function(e) { warning(paste0("[ID_CONV_UPLOAD] Error during mapIds for EntrezID: ", e$message)); NULL })
  
  if (is.null(entrez_map)) { return(list(processed_ids = as.character(original_ids), failed_indices = rep(TRUE, length(original_ids)), failed_count = length(original_ids), final_id_type = detected_keytype %||% "UNKNOWN")) }
  
  converted_values <- entrez_map[keys_to_convert]
  failed_indices_logical <- is.na(converted_values)
  num_failed <- sum(failed_indices_logical)
  
  if (num_failed > 0) { message(paste0("[ID_CONV_UPLOAD] ", num_failed, " out of ", length(original_ids), " original IDs could not be converted to EntrezID from type '", detected_keytype, "'.")) }
  
  return(list(processed_ids = as.character(converted_values), failed_indices = failed_indices_logical, failed_count = num_failed, final_id_type = "ENTREZID"))
}

# --- Helper function: Trim common strings ---
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
      h4("1. ファイルアップロード new"),
      
      radioButtons(ns("inputType"), "データの形式:",
                   choices = c("個別の featureCounts ファイル (複数可)" = "individual",
                               "結合済みカウント行列 (CSV/TSV)" = "merged"),
                   selected = "individual"),
      
      # 個別ファイル用UI
      conditionalPanel(
        condition = paste0("input['", ns("inputType"), "'] == 'individual'"),
        fileInput(ns("featureCountsFiles"), "featureCounts出力ファイル (.txt, .tsv)", 
                  multiple = TRUE, accept = c("text/plain", ".txt", ".tsv"), 
                  buttonLabel = "複数ファイル選択...", placeholder = "ファイル未選択"),
        helpText("各ファイルの1列目(Geneid)、6列目(Length)、7列目(Count)を結合します。")
      ),
      
      # 結合済みファイル用UI
      conditionalPanel(
        condition = paste0("input['", ns("inputType"), "'] == 'merged'"),
        fileInput(ns("mergedCountFile"), "カウント行列ファイル (.csv, .txt, .tsv)",
                  multiple = FALSE, accept = c(".csv", ".txt", ".tsv"),
                  buttonLabel = "ファイル選択...", placeholder = "ファイル未選択"),
        helpText(HTML("<b>注意:</b> 行名がGene Symbolの場合、自動でEntrez IDに変換します。<br>ファイルの上部や下部にメタデータ(文字列のみの行)が含まれていても、自動的に除外して読み込みます。"))
      ),
      
      hr(),
      h4("2. 生物種選択"),
      selectInput(ns("species"), label = "解析対象の生物種またはID処理方法:", choices = species_choices_ui, selected = "Homo_sapiens"),
      hr(),
      h4("3. セッション管理"),
      p("現在の作業状態を保存、または以前の状態を復元します。"),
      downloadButton(ns("downloadRDS"), "現在のセッションを保存 (.rds)"),
      br(),br(),
      fileInput(ns("uploadRDS"), "セッションを復元 (.rds)", accept = c(".rds", ".RDS")),
      helpText("注意: セッションを復元すると、現在の作業内容は上書きされます。")
    ),
    mainPanel(
      width = 8,
      h4("サンプル別 合計リード数"),
      withSpinner(plotlyOutput(ns("librarySizePlot")), type = 6),
      hr(),
      h4("サンプル情報入力"),
      helpText("各サンプルのチェックボックス、名前、グループ、Time(数値)を編集してください。"),
      withSpinner(uiOutput(ns("sampleMetadataInput")), type = 6),
      helpText("注意: ファイルを再アップロードするか生物種を変更すると、入力情報はリセットされることがあります。")
    )
  )
}

dataUploadMetadataServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv$current_gene_id_type <- reactiveVal("ENTREZID")
    
    # データ処理のトリガー
    data_processing_trigger <- reactive({
      list(
        inputType = input$inputType,
        individualFiles = input$featureCountsFiles,
        mergedFile = input$mergedCountFile,
        species = input$species
      )
    })
    
    observeEvent(data_processing_trigger(), {
      trigger <- data_processing_trigger()
      inputType <- trigger$inputType
      species_code <- trigger$species
      
      if (inputType == "individual" && is.null(trigger$individualFiles)) return()
      if (inputType == "merged" && is.null(trigger$mergedFile)) return()
      
      message("--- Data Processing Start ---")
      message(paste("Input Type:", inputType))
      
      # リセット
      rv$merged_data <- NULL; rv$sample_metadata <- NULL; rv$gene_lengths <- NULL
      rv$filtered_keep <- NULL; rv$deg_results <- NULL; rv$background_genes_original <- NULL
      rv$selected_species <- species_code
      
      tryCatch({
        
        # === A. 個別 featureCounts ファイルの処理 ===
        if (inputType == "individual") {
          file_paths <- trigger$individualFiles$datapath
          original_filenames <- trigger$individualFiles$name
          names_no_ext <- sub("\\.[^.]*$", "", original_filenames)
          trimmed_names <- trim_common_strings(names_no_ext)
          initial_sample_names <- make.unique(trimmed_names)
          
          # ファイル読み込み
          raw_data_list <- map(file_paths, ~{ 
            # 必要な列がない場合のエラーメッセージを強化
            tryCatch({
              dt <- fread(.x, header = TRUE, sep = "\t", stringsAsFactors = FALSE, select = c(1, 6, 7), na.strings = c("NA", "NaN", ""))
              validate(need(ncol(dt) >= 3, paste(basename(.x),"に必要な列(Geneid, Length, Count)がありません。")))
              setnames(dt, c("OriginalGeneid", "Length", "Count"))
              dt
            }, error = function(e) {
              stop(paste0("ファイル '", basename(.x), "' の読み込みに失敗: ", e$message, 
                          "\n※結合済みマトリックスの場合は「データの形式」で「結合済みカウント行列」を選択してください。"))
            })
          })
          
          # ID変換の前処理
          first_file_ids <- raw_data_list[[1]]$OriginalGeneid
          detected_keytype <- detect_gene_id_type(first_file_ids)
          
          all_original <- map(raw_data_list, ~ unique(.x$OriginalGeneid))
          common_ids <- Reduce(intersect, all_original)
          if(length(common_ids) == 0) stop("共通のGeneIDが見つかりません。")
          
          conversion_res <- convert_or_pass_ids(common_ids, detected_keytype, species_code)
          rv$current_gene_id_type(conversion_res$final_id_type)
          
          valid_idx <- !conversion_res$failed_indices
          if(sum(valid_idx) == 0) stop("ID変換可能な遺伝子がありませんでした。")
          
          map_df <- data.frame(OriginalGeneid = common_ids[valid_idx], Geneid = conversion_res$processed_ids[valid_idx], stringsAsFactors = FALSE)
          map_df <- map_df[!duplicated(map_df$Geneid), ]
          
          # データ結合
          processed_list <- map(raw_data_list, ~{
            dt <- .x[OriginalGeneid %in% map_df$OriginalGeneid, ]
            merged <- merge(dt, map_df, by="OriginalGeneid")
            merged[, OriginalGeneid := NULL]
            merged
          })
          
          # Length情報の保存
          rv$gene_lengths <- as.data.frame(processed_list[[1]][, .(Geneid, Length)][!duplicated(Geneid)])
          
          # カウント行列作成
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
          
          # === B. 結合済みマトリックス (CSV/TSV) の処理 [強化版] ===
        } else if (inputType == "merged") {
          req(trigger$mergedFile)
          infile <- trigger$mergedFile
          
          ext <- tools::file_ext(infile$name)
          sep_char <- if(tolower(ext) %in% c("tsv", "txt")) "\t" else ","
          
          # ヘッダーなしで全て読み込み、自力で解析する（メタデータ対策）
          raw_df <- read.delim(infile$datapath, sep = sep_char, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
          
          # 1列目を除いたデータ部分を行列化
          if(ncol(raw_df) < 2) stop("ファイルに2列以上のデータが見つかりません。")
          data_vals_mat <- as.matrix(raw_df[, -1, drop = FALSE])
          
          # 各セルが数値として解釈可能かチェック
          is_num_mat <- suppressWarnings(matrix(!is.na(as.numeric(data_vals_mat)), nrow = nrow(data_vals_mat)))
          
          # 行ごとの数値セル率を計算。サンプル列が全て数値である行を「データ行」とみなす
          # (1行でもデータがあればOK)
          row_is_data <- rowSums(is_num_mat) == ncol(is_num_mat)
          
          data_indices <- which(row_is_data)
          if(length(data_indices) == 0) stop("有効な数値データ行が見つかりませんでした。")
          
          # データ行の開始位置の1つ上の行を「ヘッダー」とみなす
          first_data_idx <- data_indices[1]
          header_row_idx <- first_data_idx - 1
          
          if(header_row_idx < 1) {
            # ヘッダーがない場合(1行目からデータ)はV1, V2...を使うか、エラーにする
            # 通常はエラー推奨だが、ここでは仮の名前を付ける
            warning("ヘッダー行が見つかりません。1行目からデータとみなします。")
            sample_names <- paste0("Sample_", 1:ncol(data_vals_mat))
          } else {
            # ヘッダー行を取得
            sample_names <- as.character(raw_df[header_row_idx, -1])
          }
          
          # IDと数値データを抽出
          clean_ids <- as.character(raw_df[data_indices, 1])
          clean_counts <- data_vals_mat[data_indices, , drop = FALSE]
          
          # マトリックス構築
          clean_mat <- matrix(suppressWarnings(as.numeric(clean_counts)), nrow = nrow(clean_counts))
          rownames(clean_mat) <- clean_ids
          colnames(clean_mat) <- sample_names
          
          initial_sample_names <- sample_names
          
          # ID変換処理
          original_ids <- rownames(clean_mat)
          detected_keytype <- detect_gene_id_type(original_ids)
          
          if (species_code != "Others_Original") {
            message(paste0("Matrix Input: Detected ID type '", detected_keytype, "'. Converting to Entrez..."))
            conversion_res <- convert_or_pass_ids(original_ids, detected_keytype, species_code)
            rv$current_gene_id_type(conversion_res$final_id_type)
            
            map_df <- data.frame(Original = original_ids, New = conversion_res$processed_ids, stringsAsFactors = FALSE)
            map_df <- map_df[!is.na(map_df$New), ]
            
            if(nrow(map_df) == 0) stop("ID変換の結果、有効な遺伝子が残りませんでした。")
            
            clean_mat_subset <- clean_mat[map_df$Original, , drop=FALSE]
            df_for_agg <- as.data.frame(clean_mat_subset)
            df_for_agg$Geneid <- map_df$New
            
            agg_df <- df_for_agg %>%
              group_by(Geneid) %>%
              summarise(across(everything(), sum)) %>%
              ungroup()
            
            rv$merged_data <- as.data.frame(agg_df)
            
          } else {
            rv$current_gene_id_type(detected_keytype %||% "Original")
            df_res <- as.data.frame(clean_mat)
            df_res$Geneid <- rownames(clean_mat)
            df_res <- df_res[, c("Geneid", setdiff(colnames(df_res), "Geneid"))]
            rv$merged_data <- df_res
          }
          
          # ★★★ Gene Lengthのダミー生成 ★★★
          rv$gene_lengths <- data.frame(Geneid = rv$merged_data$Geneid, Length = 1, stringsAsFactors = FALSE)
          
          showNotification("結合済みマトリックスを読み込みました。ヘッダー以外のメタデータ行は自動的に除外されました。", type = "message", duration = 8)
        }
        
        # === 共通: サンプルメタデータの初期化 ===
        current_samples <- setdiff(colnames(rv$merged_data), "Geneid")
        
        rv$sample_metadata <- data.frame(
          id = paste0("sample_", 1:length(current_samples)),
          current_name = current_samples,
          group = rep("GroupA", length(current_samples)),
          time = rep(0, length(current_samples)),
          active = rep(TRUE, length(current_samples)),
          stringsAsFactors = FALSE
        )
        message("--- Sample metadata initialized ---")
        
      }, error = function(e) {
        error_msg <- paste("データ処理エラー:", e$message)
        showNotification(error_msg, type = "error", duration = 15)
        rv$merged_data <- NULL; rv$sample_metadata <- NULL
      })
    })
    
    # --- UI生成 ---
    output$sampleMetadataInput <- renderUI({
      if (is.null(rv$sample_metadata)) { return(tags$p("ファイルをアップロードしてください。")) }
      input_tags <- map(1:nrow(rv$sample_metadata), ~{
        sample_info <- rv$sample_metadata[.x, ]
        tagList(
          fluidRow(
            column(2, checkboxInput(inputId = ns(paste0("active_sample_", sample_info$id)), label = "解析対象", value = sample_info$active)),
            column(4, textInput(inputId = ns(paste0("sample_name_", sample_info$id)), label = paste0("サンプル ", .x, " 名前"), value = sample_info$current_name)),
            column(3, textInput(inputId = ns(paste0("group_name_", sample_info$id)), label = "グループ", value = sample_info$group)),
            column(3, numericInput(inputId = ns(paste0("time_", sample_info$id)), label = "Time", value = sample_info$time %||% 0, min = 0, step = 1))
          )
        )
      })
      do.call(tagList, input_tags)
    })
    
    # --- サンプル情報更新の監視 ---
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
        n_s <- nrow(rv$sample_metadata); new_names <- char(n_s); new_grps <- char(n_s); new_times <- num(n_s); new_acts <- logi(n_s); valid <- TRUE
        isolate({
          vals <- reactiveValuesToList(input)
          for(i in 1:n_s){
            sid <- rv$sample_metadata$id[i]
            nm <- vals[[paste0("sample_name_", sid)]]; grp <- vals[[paste0("group_name_", sid)]]; tm <- vals[[paste0("time_", sid)]]; act <- vals[[paste0("active_sample_", sid)]]
            if(is.null(nm) || is.null(grp) || is.null(tm) || is.null(act) || !is.numeric(tm)) { valid <- FALSE; break }
            new_names[i] <- nm; new_grps[i] <- grp; new_times[i] <- tm; new_acts[i] <- act
          }
        })
        if(valid){
          if(!identical(new_names, rv$sample_metadata$current_name) || !identical(new_grps, rv$sample_metadata$group) || !identical(new_times, rv$sample_metadata$time) || !identical(new_acts, rv$sample_metadata$active)){
            unq_nms <- make.unique(new_names)
            if(any(unq_nms != new_names)){ showNotification("サンプル名重複のため自動修正されました。", type="warning"); new_names <- unq_nms }
            if(!is.null(rv$merged_data) && (ncol(rv$merged_data)-1 == length(new_names))){
              rv$sample_metadata$current_name <- new_names; rv$sample_metadata$group <- new_grps; rv$sample_metadata$time <- new_times; rv$sample_metadata$active <- new_acts
              colnames(rv$merged_data) <- c("Geneid", new_names)
            }
          }
        }
      })
    })
    
    # --- プロット用データ ---
    plot_data <- reactive({
      req(rv$merged_data, rv$sample_metadata)
      act_meta <- rv$sample_metadata[rv$sample_metadata$active, , drop=FALSE]
      validate(need(nrow(act_meta) > 0, "プロット対象のサンプルがありません。"))
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
      plot_ly(d, x=~current_name, y=~TotalReads, color=~group, type='bar', text=~paste("Sample:", current_name, "<br>Reads:", format(TotalReads, big.mark=",")), hoverinfo='text') %>%
        layout(xaxis=list(title="", tickangle=-45), yaxis=list(title="Total Reads"), margin=list(b=100))
    })
    
    output$downloadRDS <- downloadHandler(filename=function(){paste0("session_", Sys.Date(), ".rds")}, content=function(f){saveRDS(reactiveValuesToList(rv), f)})
    observeEvent(input$uploadRDS, { req(input$uploadRDS); d <- readRDS(input$uploadRDS$datapath); for(n in names(d)) if(n %in% names(rv)) rv[[n]] <- d[[n]] })
    
    char <- character; num <- numeric; logi <- logical
    `%||%` <- function(a, b) if (!is.null(a)) a else b
  })
}