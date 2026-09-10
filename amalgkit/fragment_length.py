"""Single-end fragment distributions, distinct from read/spot length statistics."""

from __future__ import annotations

import csv
import math
import sys
from collections.abc import Mapping
from typing import Any

import pandas

PROVENANCE_KEY = "amalgkit_fragment_length"
FRAGMENT_METADATA_COLUMNS = (
    "fragment_length_mean",
    "fragment_length_sd",
    "fragment_length_source",
    "fragment_length_source_detail",
)


def _text(value: Any) -> str:
    if not pandas.api.types.is_scalar(value):
        raise ValueError(f"Expected a scalar fragment field, received {value!r}.")
    return "" if value is None or pandas.isna(value) else str(value).strip()


def positive_fragment_value(value: Any, label: str) -> float | None:
    """Only blank/NA cells are missing; invalid supplied numbers are errors."""
    if not _text(value):
        return None
    try:
        number = float(_text(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite number > 0; received {value!r}.") from exc
    if isinstance(value, bool) or not math.isfinite(number) or number <= 0:
        raise ValueError(f"{label} must be a finite number > 0; received {value!r}.")
    return number


def validate_fragment_args(args: Any) -> None:
    policy = getattr(args, "fragment_length_policy", "assume")
    if policy not in ("assume", "error"):
        raise ValueError("--fragment_length_policy must be assume or error.")
    raw_mean = getattr(args, "fragment_length_mean", None)
    raw_sd = getattr(args, "fragment_length_sd", None)
    mean = positive_fragment_value(None if raw_mean is None else str(raw_mean), "--fragment_length_mean")
    sd = positive_fragment_value(None if raw_sd is None else str(raw_sd), "--fragment_length_sd")
    if (mean is None) != (sd is None):
        raise ValueError("Supply both --fragment_length_mean and --fragment_length_sd.")


def load_fragment_length_file(path: str, known_runs: set[str]) -> dict[str, dict[str, str]]:
    """Read once before batch selection; keep run identifiers lexical."""
    required = {"run", "fragment_length_mean", "fragment_length_sd", "source", "source_detail"}
    records = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        if len(columns) != len(set(columns)) or not required.issubset(columns):
            raise ValueError(f"Fragment length TSV requires unique columns: {', '.join(sorted(required))}.")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed fragment length TSV row {reader.line_num}.")
            run = row["run"].strip()
            if not run or run not in known_runs:
                raise ValueError(f"Unknown or empty run in fragment length TSV: {run!r}.")
            if run in records:
                raise ValueError(f"Duplicate run in fragment length TSV: {run}.")
            for field in ("fragment_length_mean", "fragment_length_sd"):
                if positive_fragment_value(row[field], f"Run {run}: {field}") is None:
                    raise ValueError(f"Run {run}: fragment length TSV requires both mean and SD.")
            if not row["source"].strip() or not row["source_detail"].strip():
                raise ValueError(f"Run {run}: fragment length TSV requires source and source_detail.")
            records[run] = {key: value.strip() for key, value in row.items()}
    return records


def fragment_file_records(args: Any, known_runs: set[str]) -> dict[str, dict[str, str]]:
    cached = getattr(args, "_fragment_length_by_run", None)
    if cached is not None:
        return cached
    path = getattr(args, "fragment_length_file", None)
    return load_fragment_length_file(path, known_runs) if path else {}


def has_explicit_fragment_input(args: Any, row: Mapping[str, Any], records: Mapping[str, Any], run: str) -> bool:
    return (
        run in records
        or getattr(args, "fragment_length_mean", None) is not None
        or any(_text(row.get(field)) for field in FRAGMENT_METADATA_COLUMNS[:2])
    )


def _nominal_candidate(row: Mapping[str, Any], run: str) -> tuple[Any, Any, str, str]:
    mean = row.get("nominal_length")
    sd = row.get("nominal_sdev")
    insert = row.get("mean_insert_size")
    nominal = positive_fragment_value(mean, f"Run {run}: nominal_length")
    alternative = positive_fragment_value(insert, f"Run {run}: mean_insert_size")
    if nominal is not None and alternative is not None and nominal != alternative:
        raise ValueError(f"Run {run}: conflicting nominal_length and mean_insert_size; supply a run-specific override.")
    if nominal is None and alternative is not None:
        # NOMINAL_SDEV belongs to NOMINAL_LENGTH, not an unrelated sample attribute.
        if positive_fragment_value(sd, f"Run {run}: nominal_sdev") is not None:
            raise ValueError(f"Run {run}: mean_insert_size cannot be combined with an unassociated nominal_sdev.")
        return insert, None, "metadata_mean_insert_size", "mean_insert_size sample attribute; measurement unverified"
    source = _text(row.get("nominal_length_source")) or "metadata_nominal"
    detail = _text(row.get("nominal_length_source_detail")) or "nominal_length/nominal_sdev; measurement unverified"
    return mean, sd, source, detail


def resolve_fragment_distribution(
    args: Any, row: Mapping[str, Any], run: str, records: Mapping[str, Mapping[str, str]], *, emit_warning: bool = True
) -> dict[str, Any]:
    """Choose one source pair, then fill only missing components if explicitly allowed."""
    validate_fragment_args(args)
    policy = getattr(args, "fragment_length_policy", "assume")
    mean: Any
    sd: Any
    if run in records:
        record = records[run]
        mean, sd = record["fragment_length_mean"], record["fragment_length_sd"]
        source, detail = "run_file:" + record["source"], record["source_detail"]
    elif any(_text(row.get(field)) for field in FRAGMENT_METADATA_COLUMNS[:2]):
        mean, sd = row.get("fragment_length_mean"), row.get("fragment_length_sd")
        declared_source = _text(row.get("fragment_length_source")) or "user"
        detail = _text(row.get("fragment_length_source_detail"))
        if declared_source != "user" and not detail:
            raise ValueError(f"Run {run}: fragment_length_source_detail is required for source {declared_source!r}.")
        source = "metadata:" + declared_source
        detail = detail or "run-specific fragment_length_mean/fragment_length_sd"
    elif getattr(args, "fragment_length_mean", None) is not None:
        mean, sd = args.fragment_length_mean, args.fragment_length_sd
        source, detail = "cli", "--fragment_length_mean/--fragment_length_sd"
    else:
        mean, sd, source, detail = _nominal_candidate(row, run)
    supplied = {"mean": _text(mean), "sd": _text(sd)}
    values = {
        "mean": positive_fragment_value(mean, f"Run {run}: fragment length mean ({source})"),
        "sd": positive_fragment_value(sd, f"Run {run}: fragment length SD ({source})"),
    }
    missing = [field for field, value in values.items() if value is None]
    if missing and policy == "error":
        raise ValueError(
            f"Run {run}: missing fragment length {', '.join(missing)}; supply mean/SD or use --fragment_length_policy assume."
        )
    defaults = {"mean": 200.0, "sd": 20.0}
    parameters = {}
    for field, value in values.items():
        parameters[field] = {
            "value": defaults[field] if value is None else value,
            "source": "assumed" if value is None else source,
            "source_detail": f"assume policy default {defaults[field]:g} bp" if value is None else detail,
            "supplied_value": supplied[field],
        }
    if missing and emit_warning:
        adopted = ", ".join(f"{field}={defaults[field]:g} bp" for field in missing)
        print(
            f"WARNING: Run {run}: assuming fragment length {adopted}; not measured. Supply run-specific values or use --fragment_length_policy error.",
            file=sys.stderr,
        )
    return {
        "schema_version": 1,
        "layout": "single",
        "policy": policy,
        "input_source": source,
        "input_source_detail": detail,
        **parameters,
    }


def validate_fragment_provenance(value: Any) -> str:
    """Shared optional JSON contract; old kallisto outputs have no such field."""
    if not isinstance(value, dict) or type(value.get("schema_version")) is not int or value["schema_version"] != 1:
        return "Invalid AMALGKIT fragment length provenance schema."
    layout = value.get("layout")
    if layout == "paired":
        return (
            "" if value.get("source") == "kallisto_paired_estimation" else "Invalid paired fragment length provenance."
        )
    if layout != "single" or value.get("policy") not in ("assume", "error"):
        return "Invalid fragment length layout or policy."
    for field in ("mean", "sd"):
        parameter = value.get(field)
        if not isinstance(parameter, dict):
            return f"Missing fragment length {field} provenance."
        try:
            if positive_fragment_value(parameter.get("value"), field) is None:
                return f"Missing fragment length {field} value."
        except ValueError as exc:
            return str(exc)
        if not isinstance(parameter.get("source"), str) or not parameter["source"].strip():
            return f"Missing fragment length {field} source."
        if not isinstance(parameter.get("source_detail"), str) or not parameter["source_detail"].strip():
            return f"Missing fragment length {field} source_detail."
        if not isinstance(parameter.get("supplied_value"), str):
            return f"Missing fragment length {field} supplied_value."
        try:
            supplied = positive_fragment_value(parameter["supplied_value"], f"fragment length {field} supplied_value")
        except ValueError as exc:
            return str(exc)
        if parameter["source"] == "assumed":
            expected = 200.0 if field == "mean" else 20.0
            if supplied is not None or parameter["value"] != expected:
                return f"Inconsistent assumed fragment length {field}."
        elif supplied is None or supplied != float(parameter["value"]):
            return f"Fragment length {field} differs from its supplied value."
        if value["policy"] == "error" and parameter["source"] == "assumed":
            return "Strict fragment length provenance cannot contain assumptions."
    return ""
