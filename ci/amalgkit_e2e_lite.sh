#!/usr/bin/env bash
set -euo pipefail

echo "== amalgkit e2e-lite start =="

WORK="$PWD/.amalgkit_e2e_work"
FASTQ="$WORK/fastq"
FASTA="$WORK/fasta"
mkdir -p "$FASTQ" "$FASTA"

# --- 1) ダミー FASTQ を2サンプル作成（ペアエンド）
mkfq () {
  local prefix="$1"
  printf "@%s/1\nACGTACGT\n+\nFFFFFFFF\n" "$prefix" | gzip -c > "$FASTQ/${prefix}_1.fq.gz"
  printf "@%s/2\nACGTACGT\n+\nFFFFFFFF\n" "$prefix" | gzip -c > "$FASTQ/${prefix}_2.fq.gz"
}
mkfq "S1"
mkfq "S2"

# --- 2) integrate：metadata.tsv を生成（seqkit が必要）
amalgkit integrate \
  --fastq_dir "$FASTQ" \
  --out_dir "$WORK"

# === resolve metadata tsv path robustly ===
# ローカルFASTQのみの場合は metadata_private_fastq.tsv、
# 既存メタと統合した場合は metadata_updated_for_private_fastq.tsv、
# それ以外は metadata.tsv が出力される可能性あり。
META="$(
  find "$WORK" -maxdepth 2 -type f \
    \( -name 'metadata_private_fastq.tsv' \
       -o -name 'metadata_updated_for_private_fastq.tsv' \
       -o -name 'metadata.tsv' \) \
  | head -n1 || true
)"

if [ -z "$META" ] || [ ! -s "$META" ]; then
  echo "metadata.tsv not found under $WORK" >&2
  echo "-- debug: list candidates --" >&2
  find "$WORK" -maxdepth 3 -type f -name 'metadata*.tsv' -ls >&2 || true
  exit 1
fi

echo "Found metadata sheet: $META"

# ★ 追加: 以降のコマンドが期待する既定パスに合わせてコピー
mkdir -p "$WORK/metadata"
cp -f "$META" "$WORK/metadata/metadata.tsv"
META="$WORK/metadata/metadata.tsv"

# scientific_name / curate_group を埋める（全行を対象に安全に追記）
python - "$META" <<'PY'
import csv, sys, pathlib
meta_path = pathlib.Path(sys.argv[1])

with meta_path.open("r", newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    fieldnames = list(reader.fieldnames or [])

# 無ければ列を追加、空なら埋める（テスト用の固定値）
for r in rows:
    if 'scientific_name' not in r or not r['scientific_name']:
        r['scientific_name'] = 'Testus testus'
    if 'curate_group' not in r or not r['curate_group']:
        r['curate_group'] = 'group1'

# DictWriter 用に全列のユニオン
all_fields, seen = [], set()
for k in fieldnames + [k for row in rows for k in row.keys()]:
    if k not in seen:
        seen.add(k); all_fields.append(k)

with meta_path.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=all_fields, delimiter="\t")
    w.writeheader(); w.writerows(rows)

print(f"Updated metadata: {meta_path}")
PY

# --- 4) config / select：最低限動作
amalgkit config --config base --out_dir "$WORK"   # $WORK/config_base を作る

# ★ 追加：select が期待する $WORK/config を用意（リンク or コピー）
if [ ! -d "$WORK/config" ]; then
  ln -s "$WORK/config_base" "$WORK/config" 2>/dev/null || {
    mkdir -p "$WORK/config"
    cp -f "$WORK"/config_base/*.config "$WORK/config/" || true
  }
fi
ls -l "$WORK"/config* || true

amalgkit select --metadata "$META" --out_dir "$WORK" --max_sample 99

# --- 5) quant：モックの kallisto を使って軽量実行
# 参照配列（FASTA）を用意。ファイル名は "Genus_species*.fasta" 形式が望ましい
cat > "$FASTA/Testus_testus_dummy.fasta" <<FA
>tx1
ACGTACGT
>tx2
ACGTACGT
FA

# インデックス作成も quant 側に任せる（--build_index yes）
# バッチ指定で 1行ずつ処理
amalgkit quant --metadata "$META" --out_dir "$WORK" --fasta_dir "$FASTA" --build_index yes --threads 1 --batch 1
amalgkit quant --metadata "$META" --out_dir "$WORK" --fasta_dir "$FASTA" --build_index yes --threads 1 --batch 2
echo "[ok] quant (mocked)"

# --- 6) merge：種ごとに統合テーブル生成
amalgkit merge --out_dir "$WORK"

# 生成物の存在確認（どこに出ても拾えるよう find）
TPM_FILE="$(find "$WORK" -type f -name 'Testus_testus_tpm.tsv' | head -n1 || true)"
EFF_FILE="$(find "$WORK" -type f -name 'Testus_testus_eff_length.tsv' | head -n1 || true)"
EST_FILE="$(find "$WORK" -type f -name 'Testus_testus_est_counts.tsv' | head -n1 || true)"
test -s "${TPM_FILE:-}" && echo "[ok] merge -> $TPM_FILE" || (echo "tpm table not found" >&2; exit 1)
test -s "${EFF_FILE:-}" && echo "[ok] merge -> $EFF_FILE" || (echo "eff_length table not found" >&2; exit 1)
test -s "${EST_FILE:-}" && echo "[ok] merge -> $EST_FILE" || (echo "est_counts table not found" >&2; exit 1)

# 列数が target_id + 2サンプル以上になっているか軽く検査
awk -F'\t' 'NR==1{ if (NF<3) { printf("header columns=%d (expected>=3)\n", NF) > "/dev/stderr"; exit 1 } }' "$TPM_FILE"
echo "[ok] merge header columns >= 3"

# --- 7) sanity：一通り出力の存在チェック
amalgkit sanity --metadata "$META" --out_dir "$WORK" --all
echo "[ok] sanity"

echo "== amalgkit e2e-lite done =="
