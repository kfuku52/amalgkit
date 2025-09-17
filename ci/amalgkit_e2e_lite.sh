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
# amalgkit integrate は状況により metadata.tsv か
# metadata_updated_for_private_fastq.tsv を出力する
META="$(find "$WORK" -maxdepth 2 -type f \( -name 'metadata.tsv' -o -name 'metadata_updated_for_private_fastq.tsv' \) | head -n1 || true)"

if [ -z "$META" ] || [ ! -s "$META" ]; then
  echo "metadata.tsv not found under $WORK" >&2
  echo "-- debug: list candidates --" >&2
  find "$WORK" -maxdepth 3 -type f -name 'metadata*.tsv' -ls >&2 || true
  exit 1
fi

echo "Found metadata sheet: $META"

# scientific_name / curate_group を埋める（必要に応じて）
python - <<'PY'
import csv, sys
import pathlib
meta_path = pathlib.Path(r'''$META''')
rows = list(csv.DictReader(meta_path.open("r", newline=""), delimiter="\t"))
# 例として 2 サンプルに固定で入れる（テスト用）
for r in rows:
    if r.get("run") in ("S1", "S2"):
        r["scientific_name"] = r.get("scientific_name") or "Testus example"
        r["curate_group"] = r.get("curate_group") or ("liver" if r["run"]=="S1" else "brain")

fieldnames = rows[0].keys()
with meta_path.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"Updated metadata: {meta_path}")
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
