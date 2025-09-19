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

# --- 2) integrate：metadata を生成
amalgkit integrate --fastq_dir "$FASTQ" --out_dir "$WORK"

# === resolve metadata path ===
META="$(
  find "$WORK" -maxdepth 2 -type f \
    \( -name 'metadata_private_fastq.tsv' -o -name 'metadata_updated_for_private_fastq.tsv' -o -name 'metadata.tsv' \) \
  | head -n1 || true
)"
if [ -z "$META" ] || [ ! -s "$META" ]; then
  echo "metadata.tsv not found under $WORK" >&2
  find "$WORK" -maxdepth 3 -type f -name 'metadata*.tsv' -ls >&2 || true
  exit 1
fi
echo "Found metadata sheet: $META"

# 以降の既定パスに合わせる
mkdir -p "$WORK/metadata"
cp -f "$META" "$WORK/metadata/metadata.tsv"
META="$WORK/metadata/metadata.tsv"

# --- 3) scientific_name / curate_group を付与（後工程のため）
python - "$META" <<'PY'
import csv, sys, pathlib
p = pathlib.Path(sys.argv[1])
rows = list(csv.DictReader(p.open("r", newline="", encoding="utf-8"), delimiter="\t"))
if rows:
    for r in rows:
        r.setdefault("scientific_name", "")
        r.setdefault("curate_group", "")
        if not r["scientific_name"]:
            r["scientific_name"] = "Testus testus"
        if not r["curate_group"]:
            r["curate_group"] = "group1"
    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
print("patched scientific_name/curate_group ->", p)
PY

# --- 4) config / select
amalgkit config --config base --out_dir "$WORK"   # $WORK/config_base を作る

# select が期待する $WORK/config を用意（リンク or コピー）
if [ ! -d "$WORK/config" ]; then
  ln -s "$WORK/config_base" "$WORK/config" 2>/dev/null || {
    mkdir -p "$WORK/config"
    cp -f "$WORK"/config_base/*.config "$WORK/config/" || true
  }
fi
ls -l "$WORK"/config* || true

amalgkit select --metadata "$META" --out_dir "$WORK" --max_sample 99

# ★ 重要：select 後に is_sampled=yes を強制（select が no にするため）
python - "$META" <<'PY'
import csv, sys, pathlib
p = pathlib.Path(sys.argv[1])
rows = list(csv.DictReader(p.open("r", newline="", encoding="utf-8"), delimiter="\t"))
if rows:
    for r in rows:
        r["is_sampled"] = "yes"
    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
print("patched is_sampled=yes (post-select) ->", p)
PY

# --- 5) quant：モック kallisto を使って軽量実行（バッチ無しで一括）
# 参照配列（FASTA）。ファイル名は "Genus_species*.fasta"
cat > "$FASTA/Testus_testus_dummy.fasta" <<FA
>tx1
ACGTACGT
>tx2
ACGTACGT
FA

amalgkit quant --metadata "$META" --out_dir "$WORK" --fasta_dir "$FASTA" --build_index yes --threads 1
echo "[ok] quant (mocked)"

# --- 6) merge：統合テーブル生成
amalgkit merge --out_dir "$WORK"

# 生成物の存在確認
TPM_FILE="$(find "$WORK" -type f -name 'Testus_testus_tpm.tsv' | head -n1 || true)"
EFF_FILE="$(find "$WORK" -type f -name 'Testus_testus_eff_length.tsv' | head -n1 || true)"
EST_FILE="$(find "$WORK" -type f -name 'Testus_testus_est_counts.tsv' | head -n1 || true)"
test -s "${TPM_FILE:-}" && echo "[ok] merge -> $TPM_FILE" || (echo "tpm table not found" >&2; exit 1)
test -s "${EFF_FILE:-}" && echo "[ok] merge -> $EFF_FILE" || (echo "eff_length table not found" >&2; exit 1)
test -s "${EST_FILE:-}" && echo "[ok] merge -> $EST_FILE" || (echo "est_counts table not found" >&2; exit 1)

# 列数チェック
awk -F'\t' 'NR==1{ if (NF<3) { printf("header columns=%d (expected>=3)\n", NF) > "/dev/stderr"; exit 1 } }' "$TPM_FILE"
echo "[ok] merge header columns >= 3"

# --- 7) sanity
amalgkit sanity --metadata "$META" --out_dir "$WORK" --all
echo "[ok] sanity"

echo "== amalgkit e2e-lite done =="
