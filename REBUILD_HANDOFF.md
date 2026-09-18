# toll_rnaseq 再構築 引き継ぎ書（2026-06-23）

次セッションで toll の「タブ10: 作図 & 非モデル生物エンリッチメント」を作り直すための完全仕様。
これ1枚で、文脈ゼロからでも再実装できる。

---
## 0. このドキュメントの使い方
- 前回実装した `after_count/R/module_figure_enrichment.R` は yakuri 上に残存（`app.R.bak_*` も両app横に）。
  作り直す場合は参考にしつつ置換、または一旦退避してゼロから。
- 「もう一度作る」の意図（タブ10だけ作り直す / アプリ全体を再設計 / agyで作る 等）はセッション冒頭で確認すること。

## 1. toll とは / 場所
- **EasyRNA-Seq**: ウェット実験者がGUI操作のみでFASTQ→論文図表まで完結する R Shiny 統合RNA-seqプラットフォーム（Bioinformatics誌Applications Note投稿予定）。
- repo: yakuri `/home/masakazu/gdrive/toll_rnaseq/`（git）。※Google Drive実体。ローカル `/mnt/g/マイドライブ/`（要マウント）と同一になるはず。
- 起動: **Docker** `rocker/shiny-verse:4.4.1`（`docker-compose up --build`）。**conda toll env は存在しない**。
- app.R が2つある:
  - 統合（上流+下流）: `app.R`（774行, 英語タブ1-8）← **ユーザーが起動する本体**
  - 下流単体: `after_count/app.R`（206行, 日本語タブ1-9）
  - 上流単体: `before_count/EasyRNASeq_Preprocessor/app.R`
- 状態記録: `toll_project.md`(v1.2.0), `toll_log.md`, 図表案 `manuscript_figures.md`

## 2. アーキテクチャ規約（必ず踏襲）
- Shinyモジュール: `xxxUI(id)` + `xxxServer(id, rv)`、`moduleServer(id, function(input,output,session){ ns<-session$ns; ... })`, UIは `ns<-NS(id)`。
- 共有状態 `rv <- reactiveValues(merged_data, sample_metadata, file_info, gene_lengths, filtered_keep, background_genes_original, deg_results, selected_species, current_gene_id_type="ENTREZID")`。
- DEG結果は `rv$deg_results$top_tags$table`（edgeRベース, 列: Geneid=EntrezID, logFC, PValue, FDR, F/LR...）。`rv$deg_results$comparison` は2要素。
- 既存GO/GSEAモジュールは **msigdbr + OrgDb(org.Mm.eg.db等)＝モデル生物専用**。非モデル生物経路が無いのがタブ10の存在理由。
- after_count/R/ 既存モジュール: data_upload_metadata_new, filtering, processing, dimension_reduction, deg_analysis, gsea, go_enrichment_integrated, timeseries_analysis, gene_barplot_swap。
- パッケージは app.R 冒頭で全 `library()` 済 → **新モジュールで追加依存不要**:
  shiny, shinycssloaders(withSpinner), shinyjs, DT, plotly, dplyr, purrr, data.table, tibble, edgeR, ggplot2, matrixStats, pheatmap, RColorBrewer, igraph, Rtsne, umap, fgsea, msigdbr, clusterProfiler, enrichplot, AnnotationDbi, ggrepel, maSigPro, rhandsontable, org.Mm.eg.db。

## 3. タブ10 仕様（作るもの）
モジュール `after_count/R/module_figure_enrichment.R` = `figureEnrichmentUI(id)` / `figureEnrichmentServer(id)`。
**rv非依存・自己完結**（既存モジュールは無改変）。内部 `tabsetPanel` で5セクション、各々ファイルアップロード→作図→PDF/TSV DL:
1. **① DEG → Volcano/MA**: 列自動検出(Geneid/log2FoldChange/padj, MAはbaseMean)。padj・|lfc|閾値, Up/Down/NS色分け, 上位ラベル(ggrepel)。
2. **② GO/KEGG ORA → dotplot/barplot**: 列 ID/Description/p.adjust/Count/GeneRatio。enrichResultを再構築せず**TSVからggplot2で直接**描画（堅牢）。
3. **③ GSEA → NESバー**: 列 ID/Description/NES/p.adjust。|NES|上位, Activated/Suppressed色分け。
4. **④ カウント行列 → Heatmap**: pheatmap。選択遺伝子 or 分散上位N, log2(x+1)・行Zスコア option。
5. **⑤ 非モデル生物 eggNOG**: DEG表 + `*.emapper.annotations`(+任意 query2gene) → GO(BP/MF/CC)/KEGG の **ORA+GSEA**（clusterProfiler `enricher`/`GSEA`, **term2gene方式, OrgDb不要**）→ 結果表 + dotplot/NESバー + TSV/PDF DL。

## 4. 入力ファイル形式（Nextflow rnaseq `--run_func_enrich` 出力をそのまま受理）
- **DEG(DESeq2)**: `Geneid, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj`
- **ORA**: `ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count`
- **GSEA**: `ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, ..., core_enrichment`
- **eggNOG**: `*.emapper.annotations` → `read.delim(comment.char="#", header=FALSE, quote="")` で#行全除去→列名手動付与。GOs=10列目, KEGG_Pathway=13列目(`ko#####`)。
- 共通: CSV/TSV自動判別 + 先頭`#`コメント除去ヘルパ（`readLines`→#除去→`fread`）。列名は大小無視で候補マッチ。

## 5. eggNOG term2gene ロジック（⑤の核, Nextflow clusterprofiler_enrich.nf と同一）
- GO: GOs を comma split → (GO, gene) → `GO.db` で ONTOLOGY/TERM 取得しBP/MF/CC絞込 → `enricher`/`GSEA`(TERM2GENE, TERM2NAME)。
- KEGG: KEGG_Pathway split → `^ko[0-9]{5}$` のみ → 名称は `clusterProfiler:::kegg_list("pathway","ko")`（要ネット, 失敗時はID）。
- universe = 注釈済 `unique(em$gene)`。query2gene で query→gene ID マップ（無ければ query=gene）。
- ORA: DEG(padj<th & |lfc|>th) ∩ universe を gene に。GSEA: stat(or logFC)降順ランキング。

## 6. 配線（**app.R 2つとも**に反映必須）
- `required_after`/source リストに `"after_count/R/module_figure_enrichment.R"` 追加 → `source()`。
- `tabPanel(...figureEnrichmentUI("figEnrichTab"))` 追加。
- server に `figureEnrichmentServer("figEnrichTab")` 追加（統合appは「DOWNSTREAM SERVER LOGIC」節, timeseries呼び出しの後）。
- **更新前に `cp app.R app.R.bak_$(date)` バックアップ**。

## 7. 検証手順
- 構文: `/home/masakazu/miniconda3/bin/Rscript -e 'parse("....R")'`（base R, shiny不要）。
- ロジック（実データ）: env `gsea`（clusterProfiler/enrichplot/GO.db/fgsea/data.table/ggrepel あり, pheatmapは無し）で `.fe_*`ヘルパ + ⑤核を実行。
  期待値（Lotus C1 = `Data/260615_gale_drnomura_rnaseq/downstream/de/C1_*_DESeq2_all.tsv` + annotation/）: **eggNOG GO_BP ORA=7, GSEA=62**（前回一致確認済）。
- 実UI: `cd /home/masakazu/gdrive/toll_rnaseq && docker-compose up --build` → タブ「10. Figures & Non-model Enrichment」で各セクションにファイル投入。

## 8. 落とし穴 / 注意
- **GSEA**: `clusterProfiler::GSEA(seed=TRUE)` は非対話Rscriptで `.Random.seed not found` エラー → スクリプト冒頭 `set.seed(1234)` 必須。
- `%||%` はserver関数**末尾**定義でOK（reactive実行時には定義済, 既存gseaモジュールと同パターン）。
- タンパク質抽出は **AGAT**（gffreadはNCBI GTFを読めない）。
- KEGG KO経由は動物疾患pathway(Toxoplasmosis等)混入あり → 対象生物の関連pathwayを主に解釈。
- 統合appはタブ9(swap)未統合のため番号にgap。識別子は `figEnrichTab` で統一。
- Shinyをagyに委任する場合の**禁止パターン5点**: ExtendedTask / future_promise / status / open_export_modal / set_status配置 をプロンプトに明記。
- yakuri: conda/git はPATH未設定 → フルパス。`/home/masakazu/miniconda3/bin/{Rscript,nextflow}`, env は `envs/<env>/bin/Rscript`。複雑スクリプトはscp、待機はuntilループ。

## 9. 関連資産（この機能の供給元）
- **Nextflowエンリッチメント資産**: `/home/masakazu/src/nextflow` の rnaseq `--run_func_enrich`
  （`modules/{extract_proteins,deseq2_contrasts,clusterprofiler_enrich}.nf`, `rnaseq/FUNCTIONAL_ENRICHMENT.md`）。git未コミット(main)。
- **テストデータ**: `/home/masakazu/Data/260615_gale_drnomura_rnaseq/downstream/`（de/ annotation/ enrichment/, README.md）。
- 全体引き継ぎ: `/home/masakazu/SESSION_HANDOFF_2026-06-23.md`。
