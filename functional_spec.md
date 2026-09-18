# EasyRNA-Seq — Downstream Analysis Functional Specification

## Overview

本アプリは、RNA-Seqのリードカウントデータから、フィルタリング、正規化、差次発現解析（DEG）、パスウェイ解析、および各種可視化までを一気通貫で行うShinyアプリケーションである。Dockerを用いた上流解析（FASTQからのマッピング・カウント）機能と、豊富な統計手法（edgeR, DESeq2, maSigPro等）を備えた下流解析機能の両方を提供する。

---

## Module Structure

| # | モジュール名 | 対応する主な機能 / タブ | 関連ファイル |
|:--|:---|:---|:---|
| 1 | `dataUploadMetadata` | データアップロード、メタデータ設定、ID変換 | `module_data_upload_metadata_new.R` |
| 2 | `filtering` | 低発現遺伝子の除外（フィルタリング）、QC可視化 | `module_filtering.R` |
| 3 | `dimensionReduction` + `processing` | **Visualization**: 次元圧縮・ヒートマップ・個別遺伝子発現 | `module_dimension_reduction.R`, `module_processing.R` |
| 4 | `degAnalysis` | 差次発現解析（2群比較、相互作用、多群比較） | `module_deg_analysis.R` |
| 5 | `gsea` | 遺伝子セットエンリッチメント解析 (GSEA) | `module_gsea.R` |
| 6 | `goEnrichmentIntegrated` | GO/KEGG/Reactome エンリッチメント解析 | `module_go_enrichment_integrated.R` |
| 7 | `timeseriesAnalysis` | 時系列解析 (maSigPro) | `module_timeseries_analysis.R` |
| 8 | `deconvolution` | 細胞画分推定 (immunedeconv) | `module_deconvolution.R` |

**タブ順序 (サイドバーナビ):**
`Data Upload` → `QC & Filtering` → `Visualization` → `DEG Analysis` → `GSEA` → `GO Enrichment` → `Time-series` → `Deconvolution`

---

## Tab 1: Data Upload & Metadata

### 1.1 File Upload

- **入力形式**:
  - `individual`: featureCounts の個別ファイル（複数選択可）。1列目(Geneid)、6列目(Length)、7列目(Count)を抽出。
  - `merged`: 結合済みのカウント行列 (CSV/TSV)。1列目を遺伝子ID、2列目以降をカウント値として自動認識。
- **自動トリミング**: ファイル名から共通の接頭辞・接尾辞を削除してサンプル名を自動生成。

### 1.2 Species Selection & Gene ID Conversion

- **生物種選択**: ヒト、マウス、ラット、ゼブラフィッシュ等、20種類以上の主要生物種に対応。
- **ID自動変換**:
  - 入力ID（Ensembl, Symbol, RefSeq, Entrez）を自動検出し、解析用の **Entrez ID** に変換。
  - `OrgDb` パッケージを使用。Symbol変換時は「Alias」検索も行い、変換率を向上させる。
  - Ensembl IDのパターンから生物種を自動推測し、不一致がある場合は警告を表示して自動補正。

### 1.3 Metadata Settings (Sample Labels, Conditions, Batches)

- **編集**: `rhandsontable` によるスプレッドシート形式の編集。
- **設定項目**: `active`（解析対象外の設定）、`current_name`（表示名）、`group`（主要群）のほか、任意の因子（Batch等）を「＋ Add Group」ボタンで追加可能。
- **バリデーション**: 未入力セルやアクティブサンプルなしの状態を検知し、バナーで警告。

### 1.4 Session Management (Save/Restore RDS)

- 現在の全 `reactiveValues` を `.rds` ファイルとして保存・復元可能。

---

## Tab 2: Quality Control & Filtering

### 2.1 Summary Statistics

- **Value Boxes**: サンプル数、総リード数、総遺伝子数、マッピング率（固定値）を表示。
- **Library Size**: サンプルごとの総カウント数を棒グラフ（Plotly）で表示。

### 2.2 edgeR Low Count Filtering

- **アルゴリズム**: `edgeR::filterByExpr` を使用。
- **パラメータ**: `min.count` (default: 10), `min.total.count` (15), `min.prop` (0.7), `large.n` (10)。
- **グループ考慮**: メタデータの `group` 情報を自動的に考慮してフィルタリングを実行。

### 2.3 QC Plots

- **LogCPM Density**: フィルタリング前後のカウント分布を比較表示。
- **Correlation Heatmap**: サンプル間のピアソン相関係数ヒートマップ。
- **PCA Projection**: 変動上位500遺伝子に基づくサンプル間距離の俯瞰。

---

## Tab 3: Differential Expression Analysis

### 3.1 Model Settings

- **解析タイプ**:
  - `std`: ペアワイズ比較（2群間）または相互作用モデル（2要素）。
  - `lrt`: 多群比較（ANOVA-like）。
- **統計アルゴリズム**: `edgeR` (QLFTest) または `DESeq2` (Wald test / LRT)。
- **バッチ補正**: モデル式にバッチ項（Batch）を含めることが可能。

### 3.2 Results & Plots

- **閾値設定**: FDR または P-value (default: 0.05)、および Log2 Fold Change (default: 1.0)。
- **MAプロット**: 平均発現量とLogFCの関係。
- **Volcano Plot**: LogFCと有意性の散布図（Plotly）。指定遺伝子のハイライト機能。
- **K-means Clustering**: DEGsの発現パターンをクラスタリング（default: k=4）。トレンドを可視化し、クラスターリストをGO解析へ送信可能。

### 3.3 Export

- 全遺伝子の解析結果テーブル（LogFC, PValue, FDR等）を Excel/CSV で出力。

---

## Tab 4: Pathway Analysis (GSEA)

### 4.1 GSEA Settings

- **アルゴリズム**: `fgsea` を使用。
- **データベース**: `msigdbr` から Hallmark, KEGG, Reactome, GO 等を選択。キーワード検索による特定セットの指定も可能。
- **ランキング指標**: `logFC * -log10(PValue)` [推奨], `signed -log10(PValue)`, `log2FC`, `stat` から選択。

### 4.2 Results & Plots

- **Table**: NES（Normalized Enrichment Score）、padj（有意性）を表示。
- **Enrichment Plot**: 選択したパスウェイのES推移グラフを出力。

---

## Tab 5: GO/Pathway Enrichment

### 5.1 Settings

- **アルゴリズム**: `clusterProfiler` による過剰表現解析 (ORA)。
- **対象リスト**: Up-regulated, Down-regulated, All Significant、または DEGタブから送信された k-means クラスター。
- **カテゴリ**: GO (BP, MF, CC), KEGG, Reactome。

### 5.2 Results & Plots

- **可視化**: Bar Plot, Dot Plot, Cnet Plot（遺伝子とタームのネットワーク図）。
- **ID変換**: 結果テーブル内の `geneID` 列（Entrez）を表示設定に合わせて Symbol 等に変換。

---

## Tab 6: Visualization

### 6.1 Dimension Reduction (PCA/tSNE/UMAP)

- **手法**: MDS, PCA, t-SNE, UMAP, Dendrogram, Network Graph 等。
- **バッチ補正**: `limma::removeBatchEffect` を logCPM データに適用可能。
- **カスタマイズ**: PCAの色（Color）や形（Shape）をメタデータの任意の列に紐付け可能。

### 6.2 Heatmaps

- **正規化**: Raw, CPM, TPM, FPKM, logCPM (Z-score)。
- **カスタム指定**: 任意の遺伝子名（カンマ区切り）を入力して描画。
- **オプション**: クラスタリングの有無、スケーリング（行/なし）、フォントサイズ調整、サンプル順序の任意変更。

### 6.3 Individual Gene Expression

- **検索**: サーバーサイド `selectizeInput` による高速な遺伝子検索。
- **スタイル**: サンプル別（Bar）、グループ別（Bar+ErrorBar）、ミニヒートマップ。
- **統計検定**: グループ別プロット時に、2群なら t-検定、3群以上なら ANOVA のP値を自動表示。

---

## Tab 7: Time-series Analysis (maSigPro)

### 7.1 Settings

- **実験デザイン**: `edesign`（Time, Replicate）を `rhandsontable` で定義。
- **モデル**: 回帰モデルの次数（Degree）および有意水準（Alpha）、R-squared カットオフを設定。

### 7.2 Results

- **Summary**: 有意に変容する遺伝子セットを抽出。
- **Profiles**: `see.genes` による発現パターンのクラスタリング可視化。
- **Groups**: `PlotGroups` によるクラスターごとの平均トレンド表示。

---

## Tab 8: Deconvolution

### 8.1 Settings

- **手法**: `immunedeconv` (quanTIseq, EPIC, xCell, MCP-counter)。
- **データ準備**: カウントデータを内部で TPM（遺伝子長情報がある場合）または CPM に正規化し、HGNC Symbol に変換。

### 8.2 Results

- **Stack Barplot**: 各サンプルの推定細胞画分を表示。
- **Heatmap**: 細胞タイプごとのスコア分布（Row Z-score）。

---

## Reactive State (rv$*)

`app.R` の `rv` (reactiveValues) で管理される主要な状態:

| 変数名 | 型 | 用途 |
|:---|:---|:---|
| `merged_data` | data.frame | 全遺伝子×全サンプルのカウント行列。1列目は `Geneid` (Entrez)。 |
| `sample_metadata` | data.frame | サンプル名、グループ、Batch、アクティブ状態等の属性情報。 |
| `filtered_keep` | logical vector | `filtering` モジュールで確定した「保持すべき遺伝子」のフラグ。 |
| `deg_results` | list | DEG解析結果（統計量テーブル、有意遺伝子リスト、k-means結果等）。 |
| `gene_lengths` | data.frame | 遺伝子IDごとの長さ情報（TPM/FPKM計算用）。 |
| `selected_species` | string | `Homo_sapiens`, `Mus_musculus` 等の生物種コード。 |
| `current_gene_id_type` | string | 入力データのID形式（ENSEMBL, SYMBOL 等）。 |
| `background_genes_original` | character vector | フィルタリングを通過した全遺伝子のEntrez IDリスト。 |

---

## UI Input/Output ID Table

各タブの主要 input/output ID（module namespace含む）:

| タブ | ID (Namespace含む) | 種別 | 説明 |
|:---|:---|:---|:---|
| Data Upload | `dataTab-countFiles` | input | カウントデータのアップロード |
| | `dataTab-species` | input | 生物種選択 |
| | `dataTab-sampleMetadataTable` | output | メタデータ編集テーブル |
| QC/Filtering | `filterTab-perform_filtering` | input | フィルタリング実行の有無 |
| | `filterTab-min_count` | input | 最小カウント閾値 |
| | `dataTab-librarySizePlot` | output | ライブラリサイズプロット |
| DEG Analysis | `degTab-runDEG` | input | 解析実行ボタン |
| | `degTab-analysis_type` | input | モデル選択 (std / lrt) |
| | `degTab-degResultTable` | output | 結果テーブル |
| Pathway | `go_module-run_analysis` | input | エンリッチメント解析実行 |
| | `go_module-gene_set` | input | 解析対象リスト (up/down/cluster) |
| | `go_module-goDotPlot` | output | ドットプロット |
| Visualization | `dimRedTab-dim_reduction_method` | input | 可視化手法 (PCA/UMAP等) |
| | `procTab-genes_to_plot` | input | 検索・プロット対象遺伝子 |
| | `procTab-custom_heatmap_plot` | output | カスタムヒートマップ |
| Time-series | `timeseriesTab-run_analysis` | input | maSigPro 実行 |
| Deconv | `deconvTab-run_deconv` | input | Deconvolution 実行 |
