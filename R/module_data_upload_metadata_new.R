# (セッション管理機能とサンプル選択機能を統合)

library(shiny)
library(data.table) # fread, setnames, setDT, copy, set
library(purrr)      # map
library(dplyr)      # full_join, left_join
library(shinycssloaders) # withSpinner
library(plotly)     # plotlyOutput
library(AnnotationDbi) # mapIds, columns, keytypes

# --- Species UI Choices (UI表示用、値は orgdb_species_map のキーと一致) ---
species_choices_ui <- c( # ユーザーに表示する名前 = rv$selected_species に格納される値
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

# --- Species to OrgDb package mapping (グローバルまたはヘルパースクリプトに配置推奨) ---
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

# --- Helper function: Detect Gene ID type ---
detect_gene_id_type <- function(ids) {
  ids_clean <- na.omit(ids[ids != "" & !is.na(ids)])
  if (length(ids_clean) == 0) return("UNKNOWN")
  n_total <- length(ids_clean)
  n_sample <- min(n_total, 1000)
  ids_sample <- sample(ids_clean, n_sample)
  
  ensembl_pattern <- "^ENS[A-Z0-9]*[FPTG]\\d{10,}(\\.\\d+)?$"
  ensembl_match_rate <- mean(grepl(ensembl_pattern, ids_sample, ignore.case = TRUE))
  if (ensembl_match_rate > 0.8) return("ENSEMBL")
  
  refseq_pattern <- "^[NX][CMPWTZ]_[0-9]+(\\.\\d+)?$"
  refseq_match_rate <- mean(grepl(refseq_pattern, ids_sample, ignore.case = TRUE))
  if (refseq_match_rate > 0.8) return("REFSEQ")
  
  entrez_pattern <- "^[0-9]+$"
  entrez_match_rate <- mean(grepl(entrez_pattern, ids_sample))
  if (entrez_match_rate > 0.9) return("ENTREZID")
  
  symbol_pattern <- "^([A-Za-z][A-Za-z0-9-]*[A-Za-z0-9]|[A-Za-z])$"
  symbol_match_rate <- mean(grepl(symbol_pattern, ids_sample))
  if (symbol_match_rate > 0.7 && ensembl_match_rate < 0.1 && refseq_match_rate < 0.1 && entrez_match_rate < 0.1) {
    return("SYMBOL")
  }
  
  return("UNKNOWN")
}

# --- Helper function: Convert IDs to Entrez (or pass through if "Others_Original") ---
convert_or_pass_ids <- function(original_ids, detected_keytype, selected_species_code) {
  if (selected_species_code == "Others_Original") {
    message("[ID_PROC_UPLOAD] Species is 'Others_Original'. Using original IDs directly.")
    return(list(processed_ids = as.character(original_ids),
                failed_indices = rep(FALSE, length(original_ids)),
                failed_count = 0,
                final_id_type = detected_keytype %||% "ORIGINAL_UNCONVERTED"))
  }
  
  if (detected_keytype == "ENTREZID") {
    message("[ID_CONV_UPLOAD] Gene IDs are already EntrezID.")
    return(list(processed_ids = as.character(original_ids),
                failed_indices = rep(FALSE, length(original_ids)),
                failed_count = 0,
                final_id_type = "ENTREZID"))
  }
  if (detected_keytype == "UNKNOWN" || is.null(selected_species_code) || !nzchar(selected_species_code)) {
    warning("[ID_CONV_UPLOAD] Cannot convert to EntrezID: Input type is UNKNOWN or species not selected for conversion.")
    return(list(processed_ids = as.character(original_ids),
                failed_indices = rep(TRUE, length(original_ids)),
                failed_count = length(original_ids),
                final_id_type = detected_keytype %||% "UNKNOWN"))
  }
  
  orgdb_pkg_name <- orgdb_species_map[[selected_species_code]]
  if (is.null(orgdb_pkg_name) || !nzchar(orgdb_pkg_name)) {
    warning(paste0("[ID_CONV_UPLOAD] No OrgDb package found for species '", selected_species_code, "'. Cannot convert to EntrezID."))
    return(list(processed_ids = as.character(original_ids),
                failed_indices = rep(TRUE, length(original_ids)),
                failed_count = length(original_ids),
                final_id_type = detected_keytype %||% "UNKNOWN"))
  }
  if (!requireNamespace(orgdb_pkg_name, quietly = TRUE)) {
    stop(paste0("[ID_CONV_UPLOAD] OrgDb package '", orgdb_pkg_name, "' not installed. Please install it. Cannot convert to EntrezID."))
  }
  
  require(orgdb_pkg_name, character.only = TRUE, quietly = TRUE)
  org_db <- get(orgdb_pkg_name)
  
  keys_to_convert <- as.character(original_ids)
  if (detected_keytype == "ENSEMBL") {
    keys_to_convert <- gsub("\\..*$", "", keys_to_convert)
  }
  
  if (!"ENTREZID" %in% columns(org_db)) {
    warning(paste0("[ID_CONV_UPLOAD] 'ENTREZID' column not found in OrgDb '", orgdb_pkg_name, "'. Cannot convert."))
    return(list(processed_ids = as.character(original_ids),
                failed_indices = rep(TRUE, length(original_ids)),
                failed_count = length(original_ids),
                final_id_type = detected_keytype %||% "UNKNOWN"))
  }
  if (!detected_keytype %in% keytypes(org_db)) {
    warning(paste0("[ID_CONV_UPLOAD] Keytype '", detected_keytype, "' not supported by '", orgdb_pkg_name, "'. Cannot convert."))
    return(list(processed_ids = as.character(original_ids),
                failed_indices = rep(TRUE, length(original_ids)),
                failed_count = length(original_ids),
                final_id_type = detected_keytype %||% "UNKNOWN"))
  }
  
  message(paste0("[ID_CONV_UPLOAD] Converting ", length(unique(keys_to_convert)), " unique '", detected_keytype,
                 "' IDs to ENTREZID using '", orgdb_pkg_name, "'."))
  
  unique_input_keys <- unique(keys_to_convert)
  entrez_map <- tryCatch(
    suppressMessages(mapIds(org_db,
                            keys = unique_input_keys,
                            column = "ENTREZID",
                            keytype = detected_keytype,
                            multiVals = "first")),
    error = function(e) {
      warning(paste0("[ID_CONV_UPLOAD] Error during mapIds for EntrezID: ", e$message))
      NULL
    }
  )
  
  if (is.null(entrez_map)) {
    return(list(processed_ids = as.character(original_ids),
                failed_indices = rep(TRUE, length(original_ids)),
                failed_count = length(original_ids),
                final_id_type = detected_keytype %||% "UNKNOWN"))
  }
  
  converted_values <- entrez_map[keys_to_convert]
  
  failed_indices_logical <- is.na(converted_values)
  num_failed <- sum(failed_indices_logical)
  
  if (num_failed > 0) {
    message(paste0("[ID_CONV_UPLOAD] ", num_failed, " out of ", length(original_ids),
                   " original IDs could not be converted to EntrezID from type '", detected_keytype, "'."))
  }
  return(list(processed_ids = as.character(converted_values),
              failed_indices = failed_indices_logical,
              failed_count = num_failed,
              final_id_type = "ENTREZID"))
}

# --- Helper function: Trim common strings from sample names ---
find_lcp <- function(strs) { # Longest Common Prefix
  if (length(strs) < 2) return("")
  char_lists <- strsplit(strs, "")
  min_len <- min(sapply(char_lists, length))
  if (min_len == 0) return("")
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

find_lcs <- function(strs) { # Longest Common Suffix
  if (length(strs) < 2) return("")
  rev_strs <- sapply(strs, function(x) intToUtf8(rev(utf8ToInt(x))))
  rev_lcs <- find_lcp(rev_strs)
  lcs <- intToUtf8(rev(utf8ToInt(rev_lcs)))
  return(lcs)
}

trim_common_strings <- function(strs) {
  if (length(strs) < 2) return(strs)
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
      h4("1. ファイルアップロード"),
      fileInput(ns("featureCountsFiles"),
                "featureCounts出力ファイルを選択 (.txt, .tsv)",
                multiple = TRUE, accept = c("text/plain", ".txt", ".tsv"),
                buttonLabel = "ファイルを選択...", placeholder = "ファイルが選択されていません"),
      helpText(HTML("各ファイルの1列目(Geneid)、6列目(Length)、7列目(カウント)を読み込みます。<br><b>「生物種」で特定の種を選択した場合:</b> GeneidをEntrez IDに変換します。<br><b>「生物種」で「Others」を選択した場合:</b> Geneidを変換せず、元のIDのまま使用します。")),
      hr(),
      h4("2. 生物種選択"),
      selectInput(ns("species"),
                  label = "解析対象の生物種またはID処理方法を選択:",
                  choices = species_choices_ui,
                  selected = "Homo_sapiens"
      ),
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
      helpText("各サンプルのチェックボックス、名前、グループを編集してください。チェックを外すと解析対象から除外されます。"), ### ★ 変更 ★
      withSpinner(uiOutput(ns("sampleMetadataInput")), type = 6),
      helpText("注意: ファイルを再アップロードするか生物種を変更すると、入力したサンプル情報はリセットされることがあります。")
    )
  )
}

dataUploadMetadataServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv$current_gene_id_type <- reactiveVal("ENTREZID")
    
    data_processing_trigger <- reactive({
      list(files = input$featureCountsFiles, species = input$species)
    })
    
    observeEvent(data_processing_trigger(), {
      triggered_data <- data_processing_trigger()
      files_info <- triggered_data$files
      species_code <- triggered_data$species
      
      req(files_info)
      
      if (identical(files_info, rv$file_info) && identical(species_code, rv$selected_species)) {
        return()
      }
      
      rv$file_info <- files_info
      rv$selected_species <- species_code
      
      file_paths <- rv$file_info$datapath
      original_filenames <- rv$file_info$name
      
      names_no_ext <- sub("\\.[^.]*$", "", original_filenames)
      trimmed_names <- trim_common_strings(names_no_ext)
      initial_sample_names <- make.unique(trimmed_names)
      
      message("--- File Upload/Species Change Event Start ---")
      message(paste("Selected species/ID mode for processing:", rv$selected_species))
      
      rv$merged_data <- NULL; rv$sample_metadata <- NULL; rv$gene_lengths <- NULL
      rv$filtered_keep <- NULL; rv$deg_results <- NULL; rv$background_genes_original <- NULL
      
      tryCatch({
        raw_data_list <- map(file_paths, ~{
          dt <- fread(.x, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      select = c(1, 6, 7), na.strings = c("NA", "NaN", ""))
          validate(need(ncol(dt) >= 3, paste(basename(.x),"に必要な列がありません。")))
          setnames(dt, c("OriginalGeneid", "Length", "Count"))
          validate(need(is.numeric(dt$Length), paste(basename(.x),"のLength列が数値ではありません。")))
          validate(need(is.numeric(dt$Count), paste(basename(.x),"のCount列が数値ではありません。")))
          dt
        })
        message("--- ファイル読み込み完了 ---")
        validate(need(length(raw_data_list) > 0, "ファイルの読み込みに失敗しました。"))
        
        first_file_ids <- raw_data_list[[1]]$OriginalGeneid
        validate(need(length(first_file_ids) > 0, "最初のファイルにGeneIDが見つかりません。"))
        detected_original_keytype <- detect_gene_id_type(first_file_ids)
        
        if((is.null(detected_original_keytype) || detected_original_keytype == "UNKNOWN") && rv$selected_species != "Others_Original") {
          stop("アップロードされたファイルのGeneIDタイプを特定できませんでした。'Others'を選択してください。")
        }
        
        id_processing_message <- if (rv$selected_species == "Others_Original") {
          paste0("元のGeneIDタイプを '", detected_original_keytype %||% "UNKNOWN", "' と推定。このIDを直接使用します。")
        } else {
          paste0("元のGeneIDタイプを '", detected_original_keytype %||% "UNKNOWN", "' と推定。EntrezIDへ変換します...")
        }
        showNotification(id_processing_message, type="message", duration=7)
        
        all_original_geneids_list <- map(raw_data_list, ~ unique(.x$OriginalGeneid))
        common_original_geneids <- Reduce(intersect, all_original_geneids_list)
        
        if(length(common_original_geneids) == 0){
          stop("アップロードされたファイル間で共通のGeneIDが見つかりませんでした。")
        }
        
        conversion_result <- convert_or_pass_ids(
          original_ids = common_original_geneids,
          detected_keytype = detected_original_keytype,
          selected_species_code = rv$selected_species
        )
        
        rv$current_gene_id_type(conversion_result$final_id_type)
        message(paste0("--- Current internal GeneID type set to: ", rv$current_gene_id_type(), " ---"))
        
        valid_common_indices <- !conversion_result$failed_indices
        num_successfully_processed_common <- sum(valid_common_indices)
        
        if(conversion_result$failed_count > 0 && rv$selected_species != "Others_Original"){
          showNotification(
            paste0(conversion_result$failed_count, "個の共通GeneIDが変換できず、除外されます。"),
            type="warning", duration = 10
          )
        }
        validate(need(num_successfully_processed_common > 0,
                      "共通のGeneIDを処理できませんでした。生物種やID形式を確認してください。"))
        
        common_original_ids_valid <- common_original_geneids[valid_common_indices]
        common_processed_ids_valid <- conversion_result$processed_ids[valid_common_indices]
        
        processed_to_original_map_df <- data.frame(
          OriginalGeneid = common_original_ids_valid,
          ProcessedIdCol = common_processed_ids_valid,
          stringsAsFactors = FALSE
        )
        processed_to_original_map_df <- processed_to_original_map_df[!duplicated(processed_to_original_map_df$ProcessedIdCol), ]
        
        processed_raw_data_list <- map(raw_data_list, ~{
          dt_temp <- .x
          dt_temp_filtered <- dt_temp[OriginalGeneid %in% processed_to_original_map_df$OriginalGeneid, ]
          dt_temp_merged_processed <- merge(dt_temp_filtered, processed_to_original_map_df, by="OriginalGeneid", all.x=FALSE, all.y=FALSE)
          dt_temp_merged_processed[, OriginalGeneid := NULL]
          setnames(dt_temp_merged_processed, "ProcessedIdCol", "Geneid")
          return(dt_temp_merged_processed)
        })
        
        first_processed_file <- processed_raw_data_list[[1]]
        rv_gene_lengths_dt <- first_processed_file[, .(Geneid, Length)]
        rv_gene_lengths_dt <- rv_gene_lengths_dt[!duplicated(Geneid)]
        rv$gene_lengths <- as.data.frame(rv_gene_lengths_dt)
        
        count_data_list_processed <- map(processed_raw_data_list, ~ .x[, .(Geneid, Count)])
        
        merged_dt_processed <- data.table::copy(count_data_list_processed[[1]])
        setnames(merged_dt_processed, "Count", initial_sample_names[1])
        
        if (length(count_data_list_processed) > 1) {
          merged_df_processed_temp <- as.data.frame(merged_dt_processed)
          for(i in 2:length(count_data_list_processed)) {
            current_dt_processed_temp <- data.table::copy(count_data_list_processed[[i]])
            current_sample_name_processed <- initial_sample_names[i]
            setnames(current_dt_processed_temp, "Count", current_sample_name_processed)
            merged_df_processed_temp <- dplyr::full_join(merged_df_processed_temp, as.data.frame(current_dt_processed_temp), by = "Geneid")
          }
          setDT(merged_df_processed_temp)
          merged_dt_processed <- merged_df_processed_temp
        }
        
        numeric_cols <- setdiff(colnames(merged_dt_processed), "Geneid")
        for (col in numeric_cols) {
          if(is.numeric(merged_dt_processed[[col]])) {
            set(merged_dt_processed, which(is.na(merged_dt_processed[[col]])), col, 0)
          }
        }
        
        if(any(duplicated(merged_dt_processed$Geneid))) {
          count_cols_for_sum_final <- setdiff(colnames(merged_dt_processed), "Geneid")
          merged_dt_processed_final <- merged_dt_processed[, lapply(.SD, sum, na.rm = TRUE), by = Geneid, .SDcols = count_cols_for_sum_final]
          showNotification("重複するGeneIDがあったため、カウント値を集約しました。", type="warning", duration=7)
        } else {
          merged_dt_processed_final <- merged_dt_processed
        }
        
        rv$merged_data <- as.data.frame(merged_dt_processed_final)
        message("--- rv$merged_data (", rv$current_gene_id_type(), "ベース) 作成完了. Dimensions: ", paste(dim(rv$merged_data), collapse="x"), " ---")
        
        ### ★ 追加 ★ rv$sample_metadataに 'active' 列を追加
        rv$sample_metadata <- data.frame(
          id = paste0("sample_", 1:length(initial_sample_names)),
          current_name = initial_sample_names,
          group = rep("GroupA", length(initial_sample_names)),
          active = rep(TRUE, length(initial_sample_names)), # デフォルトで全てアクティブ
          stringsAsFactors = FALSE
        )
        
      }, error = function(e) {
        error_msg <- paste("ファイル処理中にエラー:", e$message)
        showNotification(error_msg, type = "error", duration = 15)
        rv$merged_data <- NULL; rv$sample_metadata <- NULL;
        rv$gene_lengths <- NULL; rv$filtered_keep <- NULL; rv$deg_results <- NULL
        rv$background_genes_original <- NULL
        validate(error_msg)
      })
      message("--- File Upload/Species Change Event End ---")
    })
    
    output$sampleMetadataInput <- renderUI({
      if (is.null(rv$sample_metadata)) { return(tags$p("ファイルをアップロードしてください。")) }
      input_tags <- map(1:nrow(rv$sample_metadata), ~{
        sample_info <- rv$sample_metadata[.x, ]
        tagList(
          ### ★ 変更 ★ チェックボックスを追加し、レイアウトを調整
          fluidRow(
            column(3, checkboxInput(inputId = ns(paste0("active_sample_", sample_info$id)), 
                                    label = "解析対象", 
                                    value = sample_info$active)),
            column(5, textInput(inputId = ns(paste0("sample_name_", sample_info$id)), 
                                label = paste0("サンプル ", .x, " 名前"), 
                                value = sample_info$current_name)),
            column(4, textInput(inputId = ns(paste0("group_name_", sample_info$id)), 
                                label = "グループ", 
                                value = sample_info$group))
          )
        )
      })
      do.call(tagList, input_tags)
    })
    
    observe({
      # Reactive dependency on all relevant inputs
      # This will trigger the debounced reactive when any of these change
      all_inputs <- reactive({
        req(rv$sample_metadata)
        # Create a list of dependencies on all dynamic inputs
        lapply(1:nrow(rv$sample_metadata), function(i) {
          sample_id <- rv$sample_metadata$id[i]
          list(
            input[[paste0("sample_name_", sample_id)]],
            input[[paste0("group_name_", sample_id)]],
            input[[paste0("active_sample_", sample_id)]]
          )
        })
      })

      # Debounce the reaction to input changes
      debounced_inputs <- debounce(all_inputs, 500) # 500ms delay

      observeEvent(debounced_inputs(), {
        req(rv$sample_metadata, rv$merged_data)

        n_samples <- nrow(rv$sample_metadata)
        new_names <- character(n_samples)
        new_groups <- character(n_samples)
        new_actives <- logical(n_samples)
        valid_inputs <- TRUE

        # Use isolated inputs to prevent re-triggering
        isolate({
          isolate_inputs <- reactiveValuesToList(input)
          for (i in 1:n_samples) {
            sample_id <- rv$sample_metadata$id[i]
            name_val <- isolate_inputs[[paste0("sample_name_", sample_id)]]
            group_val <- isolate_inputs[[paste0("group_name_", sample_id)]]
            active_val <- isolate_inputs[[paste0("active_sample_", sample_id)]]

            if (is.null(name_val) || is.null(group_val) || is.null(active_val)) {
              valid_inputs <- FALSE
              break
            }
            new_names[i] <- name_val
            new_groups[i] <- group_val
            new_actives[i] <- active_val
          }
        })

        if (valid_inputs) {
          # Check for actual changes before updating
          if (!identical(new_names, rv$sample_metadata$current_name) ||
              !identical(new_groups, rv$sample_metadata$group) ||
              !identical(new_actives, rv$sample_metadata$active)) {

            unique_new_names <- make.unique(new_names)
            if(any(unique_new_names != new_names)) {
              showNotification("サンプル名に重複があったため変更されました。", type = "warning", duration = 5)
              new_names <- unique_new_names
              for(i_update in 1:n_samples) {
                if (new_names[i_update] != isolate(input[[paste0("sample_name_", rv$sample_metadata$id[i_update])]])) {
                  updateTextInput(session, paste0("sample_name_", rv$sample_metadata$id[i_update]), value = new_names[i_update])
                }
              }
            }

            if (!is.null(rv$merged_data) && ("Geneid" %in% colnames(rv$merged_data)) && (ncol(rv$merged_data) - 1 == length(new_names))) {
              # Update metadata
              rv$sample_metadata$current_name <- new_names
              rv$sample_metadata$group <- new_groups
              rv$sample_metadata$active <- new_actives

              # Update merged_data column names
              current_merged_data_colnames <- colnames(rv$merged_data)
              new_colnames_for_merged <- c("Geneid", new_names)

              if(!identical(current_merged_data_colnames, new_colnames_for_merged)) {
                colnames(rv$merged_data) <- new_colnames_for_merged
                message("--- サンプル名とグループ、active状態が更新されました ---")
              } else {
                message("--- サンプル情報（グループ or active状態）が更新されました ---")
              }
            }
          }
        }
      })
    })
    
    plot_data <- reactive({
      req(rv$merged_data, rv$sample_metadata)
      
      ### ★ 変更 ★ アクティブなサンプルのみを対象にする
      active_metadata <- rv$sample_metadata[rv$sample_metadata$active, , drop = FALSE]
      validate(need(nrow(active_metadata) > 0, "プロット対象のサンプルがありません。少なくとも1つはチェックを入れてください。"))
      
      expected_sample_colnames <- active_metadata$current_name
      
      current_merged_data <- rv$merged_data
      
      cols_for_summing <- intersect(expected_sample_colnames, colnames(current_merged_data))
      if(length(cols_for_summing) == 0) return(NULL)
      numeric_data_for_summing <- current_merged_data[, cols_for_summing, drop = FALSE]
      
      are_numeric <- sapply(numeric_data_for_summing, is.numeric)
      if(!all(are_numeric)) return(NULL)
      total_reads_vector <- colSums(numeric_data_for_summing, na.rm = TRUE)
      
      if(length(total_reads_vector) == 0) return(NULL)
      
      reads_df <- data.frame(current_name = names(total_reads_vector), TotalReads = total_reads_vector, row.names = NULL, stringsAsFactors = FALSE)
      
      # active_metadataからプロット用のメタデータを取得
      plot_metadata_df <- active_metadata[, c("current_name", "group"), drop = FALSE]
      
      combined_plot_data <- dplyr::left_join(reads_df, plot_metadata_df, by = "current_name")
      
      # 表示順を元のサンプル順に合わせる
      ordered_levels <- intersect(active_metadata$current_name, combined_plot_data$current_name)
      combined_plot_data$current_name <- factor(combined_plot_data$current_name, levels = ordered_levels)
      
      combined_plot_data <- combined_plot_data[!is.na(combined_plot_data$current_name), ]
      if(nrow(combined_plot_data) == 0) return(NULL)
      
      return(combined_plot_data)
    })
    
    output$librarySizePlot <- renderPlotly({
      data_to_plot <- plot_data()
      req(data_to_plot, nrow(data_to_plot) > 0)
      
      p <- plot_ly(data = data_to_plot, x = ~current_name, y = ~TotalReads, color = ~group, type = 'bar',
                   text = ~paste("サンプル:", current_name, "<br>合計リード数:", format(TotalReads, big.mark = ",", scientific = FALSE), "<br>グループ:", group),
                   hoverinfo = 'text') %>%
        layout(
          xaxis = list(title = "サンプル", tickangle = -45, categoryorder = "array", categoryarray = ~current_name),
          yaxis = list(title = "合計リード数", tickformat = ','),
          legend = list(title = list(text = '<b> グループ </b>')),
          margin = list(b = 150)
        )
      return(p)
    })
    
    # --- セッションの保存 (RDS) ---
    output$downloadRDS <- downloadHandler(
      filename = function() {
        paste0("shiny_session_", Sys.Date(), ".rds")
      },
      content = function(file) {
        showNotification("セッションデータを準備しています...時間がかかる場合があります。",
                         duration = NULL, id = "save_session_message", type = "message")
        session_data <- reactiveValuesToList(rv)
        saveRDS(session_data, file)
        removeNotification("save_session_message")
        showNotification("セッションの保存が完了しました。", type = "message")
      }
    )
    
    # --- セッションの復元 (RDS) ---
    observeEvent(input$uploadRDS, {
      rds_file <- input$uploadRDS
      if (is.null(rds_file)) {
        return(NULL)
      }
      showNotification("セッションファイルを読み込んでいます...", duration = 5, type = "message")
      tryCatch({
        loaded_session_data <- readRDS(rds_file$datapath)
        loaded_names <- names(loaded_session_data)
        current_rv_names <- names(reactiveValuesToList(rv))
        for (name in loaded_names) {
          if (name %in% current_rv_names) {
            rv[[name]] <- loaded_session_data[[name]]
          }
        }
        showNotification("セッションの復元が完了しました！各タブの表示が更新されるまでお待ちください。", type = "success", duration = 10)
      }, error = function(e) {
        showNotification(paste("セッションの復元に失敗しました:", e$message), type = "error", duration = 15)
      })
    })
    
    `%||%` <- function(a, b) {
      if (!is.null(a)) a else b
    }
    
  })
}