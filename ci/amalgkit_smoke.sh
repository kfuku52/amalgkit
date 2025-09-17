#!/usr/bin/env bash
set -Eeuo pipefail
(set -o pipefail 2>/dev/null) || true

echo "[amalgkit] smoke start"
export LC_ALL=C LANG=C

# ---------- helpers ----------
assert_contains () {
  # $1: haystack (string), $2: regex
  grep -Eiq "$2" <<<"$1" || {
    echo "ASSERT FAIL: output does not contain /$2/"
    echo "----- output begin -----"
    echo "$1"
    echo "----- output end -----"
    exit 1
  }
}
require_cmd () {
  command -v "$1" >/dev/null || { echo "MISSING command: $1"; exit 1; }
}

# ---------- versions ----------
python --version
pip --version

# ---------- CLI presence ----------
require_cmd amalgkit

# Prefer --version if実装あり。無ければ -h の冒頭を確認
if amalgkit --version >/tmp/amalgkit.version 2>&1; then
  echo "[amalgkit] --version:"
  cat /tmp/amalgkit.version
else
  echo "[amalgkit] --version not available; showing -h head"
  amalgkit -h | head -n 20 || true
fi

# ---------- top-level help ----------
top_help="$(amalgkit -h || true)"
assert_contains "$top_help" "amalgkit"

# ---------- subcommands help ----------
subs=(
  metadata
  integrate
  config
  select
  getfastq
  quant
  merge
  cstmm
  curate
  csca
  sanity
)

for s in "${subs[@]}"; do
  echo "---- [help] amalgkit $s ----"
  out="$(amalgkit "$s" -h || true)"
  # helpが返る＆そのサブコマンド名が含まれていることを最低限確認
  assert_contains "$out" "$s"
done

# ---------- python import ----------
python - <<'PY'
import importlib, sys
m = importlib.import_module("amalgkit")
print("import OK; version:", getattr(m, "__version__", "unknown"))
PY

echo "[amalgkit] smoke OK"
