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

# --- 3) scientific_name / curate_group を強制上書き
python - "$META" <<'PY'
import csv, sys, pathlib
p = pathlib.Path(sys.argv[1])
rows = list(csv.DictReader(p.open("r", newline="", encoding="utf-8"), delimiter="\t"))
if rows:
    # ヘッダに列を追加（無ければ）
    fields = list(rows[0].keys())
    for col in ("scientific_name","curate_group"):
        if col not in fields: fields.append(col)

    for r in rows:
        # ★ 常に上書き（placeholderを確実に消す）
        r["scientific_name"] = "Testus_testus"
        if not r.get("curate_group"):
            r["curate_group"] = "group1"

    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
print("force-set scientific_name=Testus_testus, curate_group=group1 ->", p)
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

# post-select patch: is_sampled=yes & scientific_name=Testus testus を強制
python - "$META" <<'PY'
import csv, sys, pathlib
p = pathlib.Path(sys.argv[1])
with p.open("r", newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
if rows:
    fields = list(rows[0].keys())
    for col in ("is_sampled", "scientific_name"):
        if col not in fields: fields.append(col)
    for r in rows:
        r["is_sampled"] = "yes"
        r["scientific_name"] = "Testus_testus"
    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
print("post-select patch ->", p)
PY

python - "$META" <<'PY'
import csv, pathlib, sys
p = pathlib.Path(sys.argv[1])
rows = list(csv.DictReader(p.open("r", encoding="utf-8"), delimiter="\t"))
print("[debug] unique scientific_name:", sorted({r.get("scientific_name","") for r in rows}))
PY

echo "[debug] FASTA exact hits for Testus_testus:"
ls -l "$FASTA"/Testus_testus.* || true

# ---- FASTA を先に作る → getfastq リンク → quant ----

# 参照配列（FASTA）を検出されやすい配置で用意
mkdir -p "$FASTA"
cat > "$FASTA/Testus_testus.fasta" <<'FA'
>tx1
ACGTACGT
>tx2
ACGTACGT
FA
# 互換エイリアス（拡張子/別名）
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/Testus_testus.fa"
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/Testus_testus_dummy.fasta"

# ★サブディレクトリにも同名を用意（実装が <fasta_dir>/<Genus_species>/*.fa を探す場合に対応）
mkdir -p "$FASTA/Testus_testus"
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/Testus_testus/Testus_testus.fasta"
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/Testus_testus/Testus_testus.fa"

# ★ 追加：空白名・小文字・別名も全て用意（実装の探索パターン総当たり）
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/Testus testus.fasta"
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/Testus testus.fa"
mkdir -p "$FASTA/Testus testus"
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/Testus testus/Testus testus.fasta"
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/Testus testus/Testus testus.fa"

ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/testus_testus.fasta"
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/testus_testus.fa"
mkdir -p "$FASTA/testus_testus"
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/testus_testus/testus_testus.fasta"
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/testus_testus/testus_testus.fa"

# ありがちな別名も置いておく
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/Testus_testus_transcripts.fa"
ln -sf "$FASTA/Testus_testus.fasta" "$FASTA/Testus_testus_cdna.fa"

# FASTA を作った直後に追加
gzip -c "$FASTA/Testus_testus.fasta" > "$FASTA/Testus_testus.fa.gz"
ln -sf "$FASTA/Testus_testus.fa.gz" "$FASTA/Testus_testus.fasta.gz" || true

# quant が期待する getfastq 出力の場所へ、ダミー FASTQ をリンク
for id in S1 S2; do
  d="$WORK/getfastq/$id"
  mkdir -p "$d"
  ln -sf "$FASTQ/${id}_1.fq.gz" "$d/${id}_1.fq.gz"
  ln -sf "$FASTQ/${id}_2.fq.gz" "$d/${id}_2.fq.gz"
  ln -sf "$FASTQ/${id}_1.fq.gz" "$d/${id}_1.fastq.gz"
  ln -sf "$FASTQ/${id}_2.fq.gz" "$d/${id}_2.fastq.gz"
done

# デバッグ（シンボリックリンクも見えるように -type f/-type l 両方）
echo "[debug] FASTA dir tree (wide):"
find "$FASTA" -maxdepth 2 \( -type f -o -type l \) -printf "%p -> %l\n" | sort
echo "[debug] getfastq tree:"; find "$WORK/getfastq" -maxdepth 2 \( -type f -o -type l \) -ls | sed -n '1,20p' || true

# ★ quant はここで1回だけ実行（バッチ無し）
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
