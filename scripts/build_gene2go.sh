#!/usr/bin/env bash
# build_gene2go.sh  【アプリ管理者用ツール — エンドユーザーは実行しない】
#
# 任意の生物の gene2go アノテーション (toll_rnaseq の GO エンリッチメント入力) を
# NCBI の全生物統合ファイルから生成する。生成した {prefix}_gene2go.tsv.gz を
# after_count/data/ に置き、gtf_utils.R の bundled_gene2go_registry に登録すれば、
# エンドユーザーは「種を選ぶだけ」(コマンド不要)で GO 解析できる。
# (登録手順は gtf_utils.R の bundled_gene2go_registry 付近のコメント参照)
#
# 出力形式 (タブ区切り, gzip):  GeneID <TAB> GO <TAB> Term <TAB> Category
#   Category = Process / Function / Component (= GO BP / MF / CC)
#   → アプリの「2d. gene2go アノテーション」にそのままアップロードできる。
#
# 用法:
#   ./build_gene2go.sh <tax_id> <output_prefix>
# 例 (ミヤコグサ Lotus japonicus, tax_id=34305):
#   ./build_gene2go.sh 34305 lja      →  lja_gene2go.tsv.gz
# 例 (シロイヌナズナ Arabidopsis thaliana, tax_id=3702):
#   ./build_gene2go.sh 3702 ath       →  ath_gene2go.tsv.gz
#
# tax_id の調べ方:
#   - NCBI Taxonomy で種名検索: https://www.ncbi.nlm.nih.gov/taxonomy
#   - またはアセンブリレポート(GCF_*_assembly_report.txt)の "Taxid:" 行
#
# KEGG 生物種コードの調べ方 (アプリの「KEGG生物種コード」入力用):
#   curl -s https://rest.kegg.jp/list/organism | grep -i "<species name>"
#   例: Lotus japonicus → lja
#
# 注意:
#   - NCBI gene2go.gz は全生物を含むため大きい (~1GB前後)。一度DLしたら GENE2GO_CACHE で使い回す。
#   - NCBI が Gnomon 等で GO 注釈を付けている生物のみ抽出可能。0件なら NCBI に GO 注釈が無い
#     生物 → eggNOG-mapper 経路 (アプリの「2c. カスタムアノテーション」) を使う。
#   - 遺伝子IDはカウント行列の Geneid と一致させること。NCBI RefSeq(LOC######)を使った
#     カウントなら GeneID は LOC の数値部分 (= gene2go の GeneID 列) と一致する。

set -euo pipefail

TAXID="${1:?Usage: build_gene2go.sh <tax_id> <output_prefix>   (例: build_gene2go.sh 34305 lja)}"
PREFIX="${2:?Usage: build_gene2go.sh <tax_id> <output_prefix>   (例: build_gene2go.sh 34305 lja)}"

SRC_URL="https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
CACHE="${GENE2GO_CACHE:-./gene2go.gz}"   # 既存ファイルがあれば再DLしない。環境変数で場所指定可
OUT="${PREFIX}_gene2go.tsv.gz"

if [ ! -s "$CACHE" ]; then
  echo "[1/2] NCBI gene2go.gz をダウンロード中 (全生物・大容量)..."
  echo "      保存先: $CACHE  (次回以降は GENE2GO_CACHE=$CACHE で再利用)"
  curl -L --fail --retry 3 -o "$CACHE" "$SRC_URL"
else
  echo "[1/2] キャッシュ済み gene2go.gz を使用: $CACHE"
fi

echo "[2/2] tax_id=$TAXID を抽出 → $OUT を生成中..."
# NCBI gene2go 列順: 1=tax_id 2=GeneID 3=GO_ID 4=Evidence 5=Qualifier 6=GO_term 7=PubMed 8=Category
zcat "$CACHE" | awk -F'\t' -v tax="$TAXID" '
  BEGIN { OFS="\t"; print "GeneID","GO","Term","Category" }
  $1 == tax { print $2, $3, $6, $8 }
' | gzip > "$OUT"

N=$(( $(zcat "$OUT" | wc -l) - 1 ))
echo "完了: $OUT  (${N} 注釈)"
if [ "$N" -eq 0 ]; then
  echo "警告: tax_id=$TAXID の GO 注釈が 0 件でした。" >&2
  echo "      NCBI に GO 注釈が無い生物の可能性 → eggNOG-mapper 経路を検討してください。" >&2
  exit 2
fi
echo ""
echo "次の手順: アプリの「2. データアップロード」→「2d. gene2go アノテーション」に $OUT をアップロードし、"
echo "          (KEGG解析する場合は「KEGG生物種コード」も入力して) GOエンリッチメントを実行してください。"
