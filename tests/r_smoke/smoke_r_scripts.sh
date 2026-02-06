#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
FIXTURE_DIR="$REPO_DIR/tests/r_smoke/fixtures"
TMP_ROOT="${1:-$(mktemp -d /tmp/amalgkit_r_smoke.XXXXXX)}"

echo "[r-smoke] repo: $REPO_DIR"
echo "[r-smoke] fixture: $FIXTURE_DIR"
echo "[r-smoke] tmp:  $TMP_ROOT"

test -d "$FIXTURE_DIR" || { echo "[r-smoke] ERROR: fixture directory not found: $FIXTURE_DIR"; exit 1; }

run_rscript_checked() {
    local name="$1"
    shift
    local log_file="$TMP_ROOT/${name}.log"
    echo "[r-smoke] running: $name"
    if ! Rscript "$@" >"$log_file" 2>&1; then
        echo "[r-smoke] ERROR: $name failed"
        cat "$log_file"
        exit 1
    fi
    if grep -q "Execution halted" "$log_file"; then
        echo "[r-smoke] ERROR: 'Execution halted' detected in $name"
        cat "$log_file"
        exit 1
    fi
    tail -n 10 "$log_file"
}

# 1) merge.r smoke
cp -R "$FIXTURE_DIR/merge" "$TMP_ROOT/merge"
run_rscript_checked \
    merge \
    "$REPO_DIR/amalgkit/merge.r" \
    "$TMP_ROOT/merge" \
    "$TMP_ROOT/merge/metadata.tsv" \
    "$REPO_DIR/amalgkit/util.r"
test -f "$TMP_ROOT/merge/merge_mapping_rate.pdf"

# 2) curate.r smoke
mkdir -p "$TMP_ROOT/curate_out"
run_rscript_checked \
    curate \
    "$REPO_DIR/amalgkit/curate.r" \
    "$FIXTURE_DIR/cstmm/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae_cstmm_counts.tsv" \
    "$FIXTURE_DIR/cstmm/metadata.tsv" \
    "$TMP_ROOT/curate_out" \
    "$FIXTURE_DIR/cstmm/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae_eff_length.tsv" \
    "pearson" \
    "0.2" \
    "0" \
    "0" \
    "h2bk119r|wt|ire1-delta|wild_type_genotype" \
    "DEFAULT" \
    "log2p1-fpkm" \
    "0" \
    "0.3" \
    "no" \
    "1" \
    "1" \
    "$REPO_DIR/amalgkit/util.r" \
    "1"
test -f "$TMP_ROOT/curate_out/curate/Saccharomyces_cerevisiae/tables/Saccharomyces_cerevisiae.no.curation_final_summary.tsv"

# 3) csca.r smoke
cp -R "$FIXTURE_DIR/csca" "$TMP_ROOT/csca"
run_rscript_checked \
    csca \
    "$REPO_DIR/amalgkit/csca.r" \
    "h2bk119r|wt|ire1-delta|wild_type_genotype" \
    "DEFAULT" \
    "$TMP_ROOT/csca" \
    "$TMP_ROOT/csca/csca_input" \
    "$TMP_ROOT/csca/multispecies_busco_table.tsv" \
    "$TMP_ROOT/csca/multispecies_genecount.tsv" \
    "$REPO_DIR/amalgkit/util.r" \
    "$TMP_ROOT/csca" \
    "no"
test -f "$TMP_ROOT/csca/csca_exclusion.pdf"

echo "[r-smoke] all checks passed"
