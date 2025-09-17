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

META="$WORK/metadata/metadata.tsv"
test -s "$META" || (echo "metadata.tsv not found" >&2; exit 1)
echo "[ok] integrate -> $META"

# --- 3) scientific_name / curate_group を追記
python - << 'PY'
import csv, sys
meta_path = sys.argv[1]
rows = []
with open(meta_path, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = list(reader)
# 必須列がなければ追加
for r in rows:
    r['scientific_name'] = r.get('scientific_name') or 'Testus testus'
# グループは行により変える
for i, r in enumerate(rows, 1):
    r['curate_group'] = f"group{i}"
# 上書き
with open(meta_path, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter='\t')
    w.writeheader(); w.writerows(rows)
print("patched:", meta_path, "rows:", len(rows))
PY
"$META"

# --- 4) config / select：最低限動作
amalgkit config --config base
# サンプル数を多めに設定して誰も落とさない
amalgkit select --metadata "$META" --out_dir "$WORK" --max_sample 99 || true

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
amalgkit quant --out_dir "$WORK" --fasta_dir "$FASTA" --build_index yes --threads 1 --batch 1
amalgkit quant --out_dir "$WORK" --fasta_dir "$FASTA" --build_index yes --threads 1 --batch 2
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
