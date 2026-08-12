"""Shared, streaming validation for files exchanged between workflows."""

from __future__ import annotations

import json
import math
import os
from collections.abc import Iterable, Sequence

import numpy
import pandas

# Quant abundance table contract. Every producer (kallisto, oarfish) must emit
# these columns. For kallisto, "eff_length" is a true effective length. oarfish
# has no notion of effective length (long-read protocols have no fragment-length
# bias) and emits the annotated transcript length in that column instead; the
# run-info JSON records "eff_length_source" as "annotated_length" (oarfish) or
# "effective_length" (a producer that provided a genuine effective length).
QUANT_ABUNDANCE_REQUIRED_COLUMNS = (
    "target_id",
    "length",
    "eff_length",
    "est_counts",
    "tpm",
)
QUANT_ABUNDANCE_NUMERIC_COLUMNS = ("length", "eff_length", "est_counts", "tpm")
BUSCO_REQUIRED_COLUMNS = (
    "busco_id",
    "status",
    "sequence",
    "score",
    "length",
    "orthodb_url",
    "description",
)


def validate_nonempty_table(
    path: str,
    required_columns: Sequence[str],
    context: str,
    comment: str | None = None,
    require_data_rows: bool = True,
    require_non_target_columns: bool = False,
    numeric_nonnegative_columns: Iterable[str] = (),
    chunk_size: int = 10_000,
) -> str:
    """Return an error string for an invalid TSV, or an empty string if valid.

    Only the header and a small preview are loaded initially. Columns requiring
    row-level checks are then streamed in bounded chunks, so validation memory
    usage does not grow with an abundance matrix's row count.
    """
    try:
        preview = pandas.read_csv(
            path,
            sep="\t",
            header=0,
            comment=comment,
            nrows=5,
            low_memory=False,
        )
    except Exception as exc:
        return f"Failed to read {context}: {exc}"

    missing_columns = [column for column in required_columns if column not in preview.columns]
    if missing_columns:
        return "Missing required column(s) in {}: {}".format(context, ", ".join(missing_columns))
    if require_non_target_columns and len(preview.columns) <= len(required_columns):
        return "{} did not include any data columns beyond {}.".format(
            context,
            ", ".join(required_columns),
        )
    if require_data_rows and preview.shape[0] == 0:
        return f"{context} did not contain any data rows."

    numeric_columns = [column for column in numeric_nonnegative_columns if column in preview.columns]
    scan_columns = list(numeric_columns)
    if "target_id" in preview.columns:
        scan_columns.append("target_id")
    if not scan_columns:
        return ""

    saw_data_row = False
    saw_valid_target = False
    seen_target_ids: set[str] = set()
    try:
        with pandas.read_csv(
            path,
            sep="\t",
            header=0,
            comment=comment,
            usecols=scan_columns,
            chunksize=chunk_size,
            low_memory=False,
        ) as chunks:
            for chunk in chunks:
                saw_data_row = saw_data_row or chunk.shape[0] > 0
                if "target_id" in chunk.columns:
                    target_ids = chunk["target_id"].fillna("").astype(str).str.strip()
                    if target_ids.eq("").any():
                        return f"{context} contains missing target_id values."
                    duplicate_target_ids = set(target_ids.loc[target_ids.duplicated()].tolist())
                    duplicate_target_ids.update(set(target_ids.tolist()).intersection(seen_target_ids))
                    if duplicate_target_ids:
                        return "{} contains duplicate target_id values: {}.".format(
                            context,
                            ", ".join(sorted(duplicate_target_ids)[:5]),
                        )
                    seen_target_ids.update(target_ids.tolist())
                    saw_valid_target = saw_valid_target or bool(target_ids.ne("").any())
                for column in numeric_columns:
                    numeric = pandas.to_numeric(chunk[column], errors="coerce")
                    values = numeric.to_numpy(dtype=float)
                    if numeric.isna().any() or not numpy.isfinite(values).all():
                        return f'{context} contains non-finite values in "{column}".'
                    if (numeric < 0).any():
                        return f'{context} contains negative values in "{column}".'
    except Exception as exc:
        return f"Failed to scan {context}: {exc}"

    if require_data_rows and not saw_data_row:
        return f"{context} did not contain any data rows."
    if require_data_rows and "target_id" in preview.columns and not saw_valid_target:
        return f"{context} did not contain valid target_id values."
    return ""


def validate_quant_run_info_json(path: str) -> str:
    """Return an error string unless a quant run-info JSON is well formed."""
    try:
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception as exc:
        return f"Failed to read quant run info JSON: {exc}"
    if not isinstance(payload, dict):
        return "quant run info JSON must contain an object."
    if "p_pseudoaligned" not in payload:
        return 'quant run info JSON is missing "p_pseudoaligned".'
    try:
        value = float(payload["p_pseudoaligned"])
    except (TypeError, ValueError):
        return 'quant run info JSON has an invalid "p_pseudoaligned" value.'
    if not math.isfinite(value) or value < 0.0 or value > 100.0:
        return f'quant run info JSON has out-of-range "p_pseudoaligned": {value}'
    eff_length_source = payload.get("eff_length_source")
    if eff_length_source is not None and eff_length_source not in {
        "annotated_length",
        "effective_length",
    }:
        return 'quant run info JSON has an invalid "eff_length_source" value.'
    return ""


def validate_quant_output_files(sra_id: str, output_dir: str) -> tuple[bool, str]:
    """Validate the abundance TSV and run-info JSON produced for one run."""
    abundance_path = os.path.join(output_dir, sra_id + "_abundance.tsv")
    run_info_path = os.path.join(output_dir, sra_id + "_run_info.json")
    if not os.path.isfile(abundance_path) or os.path.getsize(abundance_path) <= 0:
        return False, f"Missing or empty quant abundance table: {abundance_path}"
    if not os.path.isfile(run_info_path) or os.path.getsize(run_info_path) <= 0:
        return False, f"Missing or empty quant run-info JSON: {run_info_path}"

    error = validate_nonempty_table(
        abundance_path,
        required_columns=QUANT_ABUNDANCE_REQUIRED_COLUMNS,
        context="quant abundance table",
        numeric_nonnegative_columns=QUANT_ABUNDANCE_NUMERIC_COLUMNS,
    )
    if error:
        return False, error
    error = validate_quant_run_info_json(run_info_path)
    if error:
        return False, error
    return True, ""
