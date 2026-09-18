# gtf_utils.R
# GTF/GFF アノテーションのパースと、全モジュール共通の遺伝子ID表示変換ヘルパー。
# OrgDb が存在しない生物種(例: Lotus japonicus)でも、アップロードされた GTF から
# gene_id -> gene_name / biotype / length を取得して Symbol 表示・TPM/FPKM を可能にする。

# 属性文字列から値を抽出する (GTF: key "value";  /  GFF3: key=value; の両対応)
.extract_gtf_attr <- function(attr_vec, keys) {
  out <- rep(NA_character_, length(attr_vec))
  for (key in keys) {
    need <- is.na(out)
    if (!any(need)) break
    pat_full <- paste0('.*\\b', key, '[ =]"?([^";]+)"?.*')
    has <- grepl(paste0('\\b', key, '[ =]'), attr_vec[need])
    vals <- rep(NA_character_, sum(need))
    if (any(has)) {
      vals[has] <- sub(pat_full, "\\1", attr_vec[need][has], perl = TRUE)
    }
    out[need] <- vals
  }
  trimws(out)
}

# GTF/GFF ファイルをパースして遺伝子アノテーション data.frame を返す。
# 返り値: data.frame(gene_id, gene_name, biotype, gene_length)  (gene_length は union exon 長, bp)
parse_gtf_annotation <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)

  # 先頭のコメント行 (#) 数を検出 (fread がコメント行で列数を誤認するのを防ぐ)
  is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  peek_con <- if (is_gz) gzfile(path, "r") else file(path, "r")
  hdr <- tryCatch(readLines(peek_con, n = 5000), error = function(e) character(0))
  close(peek_con)
  n_skip <- 0L
  for (ln in hdr) { if (startsWith(ln, "#")) n_skip <- n_skip + 1L else break }

  dt <- tryCatch(
    data.table::fread(
      path, header = FALSE, sep = "\t", quote = "", skip = n_skip,
      showProgress = FALSE, fill = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0 || ncol(dt) < 9) {
    stop("GTF/GFFファイルを読み込めませんでした。タブ区切り9列の標準GTF/GFF形式か確認してください。")
  }
  # 先頭9列を標準名に割り当て
  dt <- dt[, 1:9]
  data.table::setnames(dt, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr"))
  dt[, start := suppressWarnings(as.integer(start))]
  dt[, end := suppressWarnings(as.integer(end))]

  # 残存コメント行・不正行を除去
  dt <- dt[!startsWith(as.character(seqname), "#")]
  dt <- dt[!is.na(feature) & feature != ""]

  # --- gene_id ---
  gene_id <- .extract_gtf_attr(dt$attr, c("gene_id", "geneID"))
  # GFF3 では gene 行の ID=gene-XXX 形式になる場合があるため補完
  na_gid <- is.na(gene_id)
  if (any(na_gid)) {
    gff_id <- .extract_gtf_attr(dt$attr[na_gid], c("ID"))
    gene_id[na_gid] <- sub("^gene[-:]", "", gff_id)
  }
  dt[, gene_id := gene_id]
  dt <- dt[!is.na(gene_id) & gene_id != ""]
  if (nrow(dt) == 0) stop("GTF/GFFから gene_id を抽出できませんでした。")

  # --- メタデータ (gene_name, biotype): gene 行を優先 ---
  meta_src <- dt[feature == "gene"]
  if (nrow(meta_src) == 0) meta_src <- dt  # gene 行が無い GTF は全行から拾う
  meta_src <- meta_src[!duplicated(gene_id)]

  gene_name <- .extract_gtf_attr(meta_src$attr, c("gene_name", "gene", "Name", "gene_symbol"))
  biotype   <- .extract_gtf_attr(meta_src$attr, c("gene_biotype", "gene_type", "biotype"))
  # gene_name が無ければ gene_id をシンボルとして使う
  gene_name[is.na(gene_name) | gene_name == ""] <- meta_src$gene_id[is.na(gene_name) | gene_name == ""]

  meta_df <- data.frame(
    gene_id = meta_src$gene_id,
    gene_name = gene_name,
    biotype = ifelse(is.na(biotype), "unknown", biotype),
    stringsAsFactors = FALSE
  )

  # --- 機能記述 (product / description): transcript / mRNA / CDS 行から取得 ---
  # NCBI GTF では LOC遺伝子の gene_name は "LOCxxxx" だが product に
  # 人間可読な機能記述 (例 "disease resistance protein ...") が入る。
  # product は transcript/mRNA/CDS に加え exon 行にも入る GTF があるため広めに拾う
  desc_df <- NULL
  prod_src <- dt[feature %in% c("transcript", "mRNA", "CDS", "exon")]
  if (nrow(prod_src) > 0) {
    prod_vals <- .extract_gtf_attr(prod_src$attr, c("product", "description", "Note"))
    prod_src[, product := prod_vals]
    prod_src <- prod_src[!is.na(product) & product != ""]
    if (nrow(prod_src) > 0) {
      # アイソフォーム接尾辞 (", transcript variant X1" 等) を除去し遺伝子レベルに正規化
      prod_src[, product := sub(",?\\s*transcript variant [^;]*$", "", product)]
      prod_src <- prod_src[!duplicated(gene_id)]
      desc_df <- as.data.frame(prod_src[, list(gene_id, gene_description = product)])
    }
  }

  # --- length: exon の union 長 ---
  exon_dt <- dt[feature == "exon"]
  len_df <- NULL
  if (nrow(exon_dt) > 0) {
    len_df <- tryCatch({
      if (requireNamespace("IRanges", quietly = TRUE)) {
        # gene 毎に exon 区間を union して幅を合計 (= featureCounts の meta-feature 長相当)
        agg <- exon_dt[, {
          ir <- IRanges::reduce(IRanges::IRanges(start = start, end = end))
          list(gene_length = sum(IRanges::width(ir)))
        }, by = gene_id]
        as.data.frame(agg)
      } else {
        # フォールバック: exon 長の単純合計 (重複区間は過大評価しうる)
        agg <- exon_dt[, list(gene_length = sum(end - start + 1L)), by = gene_id]
        as.data.frame(agg)
      }
    }, error = function(e) NULL)
  }

  ann <- meta_df
  if (!is.null(len_df)) {
    ann <- merge(ann, len_df, by = "gene_id", all.x = TRUE)
  } else {
    ann$gene_length <- NA_real_
  }
  if (!is.null(desc_df)) {
    ann <- merge(ann, desc_df, by = "gene_id", all.x = TRUE)
  } else {
    ann$gene_description <- NA_character_
  }
  # 機能記述が無い遺伝子は gene_name で補完
  ann$gene_description[is.na(ann$gene_description) | ann$gene_description == ""] <-
    ann$gene_name[is.na(ann$gene_description) | ann$gene_description == ""]
  ann[!duplicated(ann$gene_id), c("gene_id", "gene_name", "biotype", "gene_length", "gene_description")]
}

# ─────────────────────────────────────────────────────────────────────────
# 同梱 gene2go レジストリ (非モデル生物の GO エンリッチメント用)
# ─────────────────────────────────────────────────────────────────────────
# OrgDb の無い生物種でも、アプリに同梱した gene2go 注釈 (GeneID, GO, Term, Category)
# を使えば、ユーザーは種を選ぶだけで GO 解析できる (ファイルのアップロード不要)。
#
# ★ 新しい非モデル生物を「種を選ぶだけ」で使えるようにする手順 (管理者向け):
#   1. scripts/build_gene2go.sh <tax_id> <prefix> で {prefix}_gene2go.tsv.gz を生成
#   2. それを after_count/data/ に配置
#   3. 下の bundled_gene2go_registry に  "<種コード>" = "{prefix}_gene2go.tsv.gz"  を1行追加
#   4. module_data_upload_metadata_new.R の species_choices_ui に種を追加
#   5. (KEGGを使うなら) module_go_enrichment_integrated.R の kegg_species_map にコードを追加
#   → これでエンドユーザーはコマンド不要、種選択だけで GO/KEGG が回る。
bundled_gene2go_registry <- c(
  "Lotus_japonicus" = "lja_gene2go.tsv.gz"
)

# 種コードに対応する同梱 gene2go ファイルの実パスを解決する (getwd 非依存)。
# 登録が無い / ファイルが見つからない場合は NULL。
resolve_bundled_gene2go <- function(species_code) {
  if (is.null(species_code) || !nzchar(species_code)) return(NULL)
  fname <- bundled_gene2go_registry[[species_code]]
  if (is.null(fname) || is.na(fname) || !nzchar(fname)) return(NULL)
  cands <- c(
    file.path("after_count", "data", fname),
    file.path(getwd(), "after_count", "data", fname),
    file.path(dirname(getwd()), "after_count", "data", fname)
  )
  hit <- cands[file.exists(cands)]
  if (length(hit) > 0) hit[1] else NULL
}

# ─────────────────────────────────────────────────────────────────────────
# アプリ内 gene2go 自動取得 (Tax ID 指定。コマンド不要)
# ─────────────────────────────────────────────────────────────────────────
# build_gene2go.sh の処理を R 関数化したもの。ユーザーは GUI で Tax ID を入れるだけ。
# NCBI gene2go.gz (全生物 約1.3GB) を「初回のみ」永続キャッシュにダウンロードし、
# 以降は任意の生物を即座に切り出す。Tax ID 別スライスもキャッシュする。

# gene2go.gz の永続キャッシュ先 (Google Drive 同期外のローカル R ユーザーキャッシュ)
ncbi_gene2go_cache_path <- function() {
  base <- tryCatch(tools::R_user_dir("EasyRNASeq", which = "cache"),
                   error = function(e) file.path(tempdir(), "EasyRNASeq_cache"))
  dir.create(base, showWarnings = FALSE, recursive = TRUE)
  file.path(base, "gene2go.gz")
}

# 指定 Tax ID の gene2go 注釈 data.frame(GeneID, GO, Term, Category) を返す。
#   tax_id   : NCBI Taxonomy ID (数値文字列。例 ミヤコグサ=34305, シロイヌナズナ=3702)
#   progress : function(fraction, message) 進捗通知 (Shiny の setProgress 等)。NULL可。
# 例外時は stop()。
fetch_gene2go_by_taxid <- function(tax_id, progress = NULL) {
  tax_id <- trimws(as.character(tax_id))
  if (!grepl("^[0-9]+$", tax_id)) {
    stop("Tax ID は数値で指定してください (例: ミヤコグサ=34305, シロイヌナズナ=3702)。")
  }
  pr <- function(f, m) if (!is.null(progress)) try(progress(f, m), silent = TRUE)

  cache_gz <- ncbi_gene2go_cache_path()
  slice_gz <- file.path(dirname(cache_gz), paste0("gene2go_", tax_id, ".tsv.gz"))

  # 1) Tax ID 別スライスがキャッシュ済みなら即読み込み
  if (file.exists(slice_gz) && file.info(slice_gz)$size > 0) {
    pr(0.9, "キャッシュから注釈を読み込み中...")
    df <- tryCatch(read.table(gzfile(slice_gz), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) return(df)
  }

  # 2) 全生物 gene2go.gz を初回のみダウンロード (約1.3GB, 数分)
  if (!file.exists(cache_gz) || file.info(cache_gz)$size < 1e8) {
    pr(0.05, "NCBI gene2go (初回のみ 約1.3GB) をダウンロード中... 数分かかります")
    old_to <- getOption("timeout"); on.exit(options(timeout = old_to), add = TRUE)
    options(timeout = 7200)
    ok <- tryCatch({
      utils::download.file("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",
                           destfile = cache_gz, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) { stop(paste("gene2go.gz のダウンロードに失敗:", conditionMessage(e))) })
  }

  # 3) ストリーミングで tax_id 行のみ抽出 (メモリ安全 + sorted前提で早期終了)
  #    gene2go 列順: 1 tax_id 2 GeneID 3 GO_ID 4 Evidence 5 Qualifier 6 GO_term 7 PubMed 8 Category
  pr(0.7, "該当生物の GO 注釈を抽出中...")
  con <- gzfile(cache_gz, "r"); on.exit(close(con), add = TRUE)
  invisible(readLines(con, n = 1))  # ヘッダ行を読み飛ばす
  rows <- list(); seen <- FALSE; chunk <- 200000L
  repeat {
    ls <- readLines(con, n = chunk)
    if (length(ls) == 0) break
    f1 <- sub("\t.*$", "", ls)
    hit <- f1 == tax_id
    if (any(hit)) { seen <- TRUE; rows[[length(rows) + 1]] <- ls[hit] }
    if (seen && !any(hit)) {
      nums <- suppressWarnings(as.numeric(f1))
      if (all(!is.na(nums)) && min(nums) > as.numeric(tax_id)) break  # tax_id を通過したら終了
    }
  }
  if (length(rows) == 0) {
    stop(paste0("Tax ID=", tax_id, " の GO 注釈が NCBI gene2go に見つかりませんでした。",
                "NCBIにGO注釈の無い生物の可能性があります (eggNOG-mapper 経路を検討してください)。"))
  }
  all_lines <- unlist(rows)
  parts <- strsplit(all_lines, "\t", fixed = TRUE)
  keep  <- lengths(parts) >= 8                      # 8列揃った行のみ採用 (ragged行を除外)
  parts <- parts[keep]
  if (length(parts) == 0) stop(paste0("Tax ID=", tax_id, " の有効な gene2go 行がありませんでした。"))
  m <- do.call(rbind, lapply(parts, function(p) p[1:8]))
  df <- data.frame(GeneID = m[, 2], GO = m[, 3], Term = m[, 6], Category = m[, 8], stringsAsFactors = FALSE)
  df <- unique(df[nzchar(df$GeneID) & nzchar(df$GO), , drop = FALSE])

  # 4) スライスをキャッシュ保存 (次回以降は即時)
  tryCatch({
    gzf <- gzfile(slice_gz, "w")
    utils::write.table(df, gzf, sep = "\t", quote = FALSE, row.names = FALSE)
    close(gzf)
  }, error = function(e) NULL)

  pr(1, "完了")
  df
}

# アップロードされた GTF アノテーションを調べ、「表示する遺伝子IDタイプ」の
# selectInput 用 choices と推奨デフォルトを返す。GTF の中身に応じて、意味のある
# 選択肢だけを提示し、ラベル(表示名)も内容に合わせて変える。
#   gtf_df : parse_gtf_annotation() の返り値 (gene_id, gene_name, gene_description, ...)
#            NULL / 無効なときは GTF 未適用時のデフォルト3択を返す。
# 返り値: list(choices = 名前付きベクトル(ラベル=値), selected = 既定値)
#   - SYMBOL  : gene_name が gene_id と実質的に異なる遺伝子が一定割合あるときのみ提示
#   - GENENAME: product/description が gene_name/gene_id と異なる遺伝子が一定割合あるとき提示
#   - ENTREZID: 生の内部ID。常に提示。NCBI由来(LOC+数字/数字)ならラベルを "NCBI Gene ID" に変える
gtf_display_id_choices <- function(gtf_df) {
  default <- list(
    choices  = c("Gene Symbol" = "SYMBOL", "Entrez ID (内部ID)" = "ENTREZID", "Gene Name" = "GENENAME"),
    selected = "SYMBOL"
  )
  if (is.null(gtf_df) || !is.data.frame(gtf_df) || nrow(gtf_df) == 0 ||
      !all(c("gene_id", "gene_name") %in% colnames(gtf_df))) {
    return(default)
  }

  gid  <- as.character(gtf_df$gene_id)
  gnm  <- as.character(gtf_df$gene_name)
  gdsc <- if ("gene_description" %in% colnames(gtf_df)) {
    as.character(gtf_df$gene_description)
  } else {
    rep(NA_character_, length(gid))
  }
  n <- length(gid)
  frac <- function(x) if (n == 0) 0 else sum(x, na.rm = TRUE) / n

  # gene_name が gene_id と実質的に異なる遺伝子の割合 (>=5% で「実シンボルあり」と判定)
  symbol_informative <- frac(!is.na(gnm) & nzchar(gnm) & gnm != gid) >= 0.05
  # 機能記述(product)が gene_name/gene_id と異なる遺伝子の割合
  desc_informative <- frac(!is.na(gdsc) & nzchar(gdsc) & gdsc != gnm & gdsc != gid) >= 0.05

  # gene_id が NCBI由来 (LOC+数字 / 純数字) かでラベルを変える
  id_is_ncbi <- frac(grepl("^(LOC)?[0-9]+$", gid)) >= 0.5
  id_label   <- if (id_is_ncbi) "NCBI Gene ID" else "Gene ID (内部ID)"

  choices <- character(0)
  if (symbol_informative) choices <- c(choices, setNames("SYMBOL", "Gene Symbol"))
  if (desc_informative)   choices <- c(choices, setNames("GENENAME", "遺伝子名 (説明 / product)"))
  choices <- c(choices, setNames("ENTREZID", id_label))

  # 既定の表示優先順位: 実シンボル > 機能記述 > 生ID
  selected <- if (symbol_informative) "SYMBOL" else if (desc_informative) "GENENAME" else "ENTREZID"

  list(choices = choices, selected = selected)
}

# 内部 Geneid を表示用ラベルに変換する共通ヘルパー。
# 優先順位: (1) アップロード GTF アノテーション -> (2) OrgDb -> (3) 生ID
#   internal_ids        : rv$merged_data の Geneid (変換後の内部ID)
#   target_display_type : "SYMBOL" / "ENTREZID" / "GENENAME"
#   species_code        : rv$selected_species
#   gene_annotation     : rv$gene_annotation (data.frame: Geneid, gene_name, biotype, length) or NULL
#   orgdb_map           : orgdb_species_map (名前付きベクトル)
annotate_display_ids <- function(internal_ids, target_display_type, species_code, gene_annotation = NULL, orgdb_map = NULL) {
  internal_ids <- as.character(internal_ids)
  if (is.null(target_display_type) || !nzchar(target_display_type) || target_display_type == "ENTREZID") {
    return(internal_ids)
  }

  # (1) GTF アノテーション優先
  #   SYMBOL  -> gene_name (例 LOC130709451 / 命名済みなら実シンボル)
  #   GENENAME-> gene_description (product 機能記述, 例 "disease resistance protein ...")
  if (!is.null(gene_annotation) && is.data.frame(gene_annotation) &&
      all(c("Geneid", "gene_name") %in% colnames(gene_annotation))) {
    src_col <- if (target_display_type == "GENENAME" && "gene_description" %in% colnames(gene_annotation)) {
      "gene_description"
    } else {
      "gene_name"
    }
    lut <- setNames(as.character(gene_annotation[[src_col]]), as.character(gene_annotation$Geneid))
    mapped <- lut[internal_ids]
    if (any(!is.na(mapped))) {
      return(ifelse(is.na(mapped) | mapped == "", internal_ids, as.character(mapped)))
    }
  }

  # (2) Lotus japonicus: 専用 OrgDb 無し。SYMBOL は LOC+番号 で復元
  if (!is.null(species_code) && species_code == "Lotus_japonicus") {
    if (target_display_type == "SYMBOL") return(paste0("LOC", internal_ids))
    return(internal_ids)
  }

  # (3) OrgDb 変換
  if (is.null(orgdb_map) || is.null(species_code) || !(species_code %in% names(orgdb_map))) {
    return(internal_ids)
  }
  orgdb_pkg_name <- orgdb_map[[species_code]]
  if (is.null(orgdb_pkg_name) || !requireNamespace(orgdb_pkg_name, quietly = TRUE)) {
    return(internal_ids)
  }
  org_db <- get(orgdb_pkg_name)
  if (!target_display_type %in% AnnotationDbi::columns(org_db) ||
      !"ENTREZID" %in% AnnotationDbi::keytypes(org_db)) {
    return(internal_ids)
  }
  unique_keys <- unique(internal_ids)
  converted_map <- tryCatch(
    suppressMessages(AnnotationDbi::mapIds(org_db, keys = unique_keys, column = target_display_type, keytype = "ENTREZID", multiVals = "first")),
    error = function(e) NULL
  )
  if (is.null(converted_map)) return(internal_ids)
  final_ids <- converted_map[internal_ids]
  na_idx <- is.na(final_ids)
  if (any(na_idx)) final_ids[na_idx] <- paste0(internal_ids[na_idx], " (変換不可)")
  as.character(final_ids)
}
