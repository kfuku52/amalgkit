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

write_dcf() {
    local out_file="$1"
    shift
    : >"$out_file"
    while (($# >= 2)); do
        printf '%s: %s\n' "$1" "$2" >>"$out_file"
        shift 2
    done
}

# 1) merge.r smoke
cp -R "$FIXTURE_DIR/merge" "$TMP_ROOT/merge"
write_dcf \
    "$TMP_ROOT/merge.dcf" \
    "dir_merge" "$TMP_ROOT/merge" \
    "file_metadata" "$TMP_ROOT/merge/metadata.tsv" \
    "r_util_path" "$REPO_DIR/amalgkit/util.r"
run_rscript_checked \
    merge \
    "$REPO_DIR/amalgkit/merge.r" \
    "$TMP_ROOT/merge.dcf"
test -f "$TMP_ROOT/merge/merge_mapping_rate.pdf"
test -f "$TMP_ROOT/merge/merge_mean_expression_boxplot.pdf"

# 2) prepare_tables.r smoke
mkdir -p "$TMP_ROOT/per_species_out"
write_dcf \
    "$TMP_ROOT/per_species.dcf" \
    "est_counts_path" "$FIXTURE_DIR/cstmm/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae_cstmm_counts.tsv" \
    "metadata_path" "$FIXTURE_DIR/cstmm/metadata.tsv" \
    "out_dir" "$TMP_ROOT/per_species_out" \
    "eff_length_path" "$FIXTURE_DIR/cstmm/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae_eff_length.tsv" \
    "dist_method" "pearson" \
    "mapping_rate_cutoff" "0.2" \
    "min_dif" "0" \
    "plot_intermediate" "0" \
    "selected_sample_groups" "h2bk119r|wt|ire1-delta|wild_type_genotype" \
    "sample_group_colors" "DEFAULT" \
    "transform_method" "log2p1-fpkm" \
    "one_outlier_per_iteration" "0" \
    "correlation_threshold" "0.3" \
    "batch_effect_alg" "no" \
    "clip_negative" "1" \
    "maintain_zero" "1" \
    "r_util_path" "$REPO_DIR/amalgkit/util.r" \
    "skip_curation_flag" "1" \
    "outlier_method" "legacy" \
    "robust_margin_threshold" "0" \
    "robust_z_threshold" "-2.5" \
    "disable_auto_outlier_filter_flag" "0"
run_rscript_checked \
    per_species \
    "$REPO_DIR/amalgkit/prepare_tables.r" \
    "$TMP_ROOT/per_species.dcf"
test -f "$TMP_ROOT/per_species_out/per_species/Saccharomyces_cerevisiae/tables/Saccharomyces_cerevisiae.no.curation_final_summary.tsv"

# 3) cross_species_filter.r smoke
cp -R "$FIXTURE_DIR/csca" "$TMP_ROOT/cross_species"
if [ -d "$TMP_ROOT/cross_species/csca_input" ]; then
    mv "$TMP_ROOT/cross_species/csca_input" "$TMP_ROOT/cross_species/cross_species_input_symlinks"
fi
write_dcf \
    "$TMP_ROOT/cross_species.dcf" \
    "selected_sample_groups" "h2bk119r|wt|ire1-delta|wild_type_genotype" \
    "sample_group_colors" "DEFAULT" \
    "dir_work" "$TMP_ROOT/cross_species" \
    "dir_cross_species_input_table" "$TMP_ROOT/cross_species/cross_species_input_symlinks" \
    "file_orthogroup" "$TMP_ROOT/cross_species/multispecies_busco_table.tsv" \
    "file_genecount" "$TMP_ROOT/cross_species/multispecies_genecount.tsv" \
    "r_util_path" "$REPO_DIR/amalgkit/util.r" \
    "dir_cross_species" "$TMP_ROOT/cross_species" \
    "batch_effect_alg" "no" \
    "missing_strategy" "em_pca" \
    "cross_species_outlier_method" "none" \
    "cross_species_margin_threshold" "0" \
    "cross_species_robust_z_threshold" "-2.5" \
    "cross_species_plot_mode" "dual"
run_rscript_checked \
    cross_species \
    "$REPO_DIR/amalgkit/cross_species_filter.r" \
    "$TMP_ROOT/cross_species.dcf"
test -f "$TMP_ROOT/cross_species/cross_species_exclusion.pdf"

echo "[r-smoke] all checks passed"
