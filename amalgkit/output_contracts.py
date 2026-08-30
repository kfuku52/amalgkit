"""Shared, streaming validation for files exchanged between workflows."""

from __future__ import annotations

import json
import math
import os
from collections.abc import Iterable, Sequence

import numpy
import pandas

from amalgkit.identifier_validation import TargetIdTracker
from amalgkit.table_io import read_identifier_tsv

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


def validate_table_frame(
    frame: pandas.DataFrame,
    required_columns: Sequence[str],
    context: str,
    *,
    require_data_rows: bool = True,
    numeric_nonnegative_columns: Iterable[str] = (),
) -> str:
    """Validate already-loaded rows without reading a file a second time."""
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        return f"Missing required column(s) in {context}: {', '.join(missing)}"
    if require_data_rows and frame.empty:
        return f"{context} did not contain any data rows."
    if "target_id" in frame.columns:
        ids = frame["target_id"].fillna("").astype(str).str.strip()
        if ids.eq("").any():
            return f"{context} contains missing target_id values."
        duplicates = ids.loc[ids.duplicated()]
        if not duplicates.empty:
            return f"{context} contains duplicate target_id values: {', '.join(sorted(set(duplicates))[:5])}."
    for column in numeric_nonnegative_columns:
        values = pandas.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        if not numpy.isfinite(values).all():
            return f'{context} contains non-finite values in "{column}".'
        if (values < 0).any():
            return f'{context} contains negative values in "{column}".'
    return ""


def read_quant_abundance(path: str, value_columns: Sequence[str], context: str) -> pandas.DataFrame:
    """Read only needed quant columns once, preserving IDs and enforcing values."""
    columns = ["target_id", *value_columns]
    try:
        frame = read_identifier_tsv(path, identifier_columns=("target_id",), usecols=columns)
    except Exception as exc:
        raise ValueError(f"Failed to read {context}: {exc}") from exc
    error = validate_table_frame(frame, columns, context, numeric_nonnegative_columns=value_columns)
    if error:
        raise ValueError(error)
    frame["target_id"] = frame["target_id"].str.strip()
    return frame


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
    try:
        with (
            TargetIdTracker() as seen_target_ids,
            pandas.read_csv(
                path,
                sep="\t",
                header=0,
                comment=comment,
                usecols=scan_columns,
                chunksize=chunk_size,
                low_memory=False,
                converters={"target_id": str},
            ) as chunks,
        ):
            for chunk in chunks:
                saw_data_row = saw_data_row or chunk.shape[0] > 0
                error = validate_table_frame(
                    chunk, (), context, require_data_rows=False, numeric_nonnegative_columns=numeric_columns
                )
                if error:
                    return error
                if "target_id" in chunk.columns:
                    target_ids = chunk["target_id"].fillna("").astype(str).str.strip()
                    duplicate_target_ids = seen_target_ids.add(target_ids.tolist())
                    if duplicate_target_ids:
                        return "{} contains duplicate target_id values: {}.".format(
                            context,
                            ", ".join(sorted(duplicate_target_ids)[:5]),
                        )
                    saw_valid_target = saw_valid_target or bool(target_ids.ne("").any())
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
