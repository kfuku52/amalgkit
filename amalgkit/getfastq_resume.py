"""Getfastq resume fingerprints, state serialization, and file identities.

Keep output-affecting options in the fingerprint catalog below. Orchestration,
FASTQ validation, and decisions to restart remain in the command module.
"""

from __future__ import annotations

import hashlib
import json
import os
import stat
from collections.abc import Mapping, Sequence
from typing import Any

import numpy
import pandas

from amalgkit.metadata_utils import get_metadata_row_index_by_run
from amalgkit.output_utils import atomic_output_path

GETFASTQ_RUN_STATE_FILENAME = "getfastq_run_state.json"
GETFASTQ_COMPLETION_FILENAME = "getfastq_completion.json"
# Version 3 covers mate conversion; older fingerprints cannot prove its setting.
GETFASTQ_RESUME_SCHEMA_VERSION = 3
GETFASTQ_PHASE_FIRST_ROUND = "first_round"
GETFASTQ_PHASE_SECOND_ROUND_IN_PROGRESS = "second_round_in_progress"
GETFASTQ_PHASE_COMPLETE = "complete"


def _normalize_getfastq_resume_value(value: Any) -> Any:
    if isinstance(value, numpy.generic):
        value = value.item()
    if pandas.isna(value):
        return None
    if isinstance(value, (bool, int, float, str)) or value is None:
        return value
    return str(value)


def build_getfastq_run_fingerprint(
    args: Any, sra_stat: Mapping[str, Any], g: Mapping[str, Any], run_metadata: Any
) -> str:
    option_names = [
        "max_bp",
        "min_read_length",
        "treat_identical_paired_as_single",
        "fastp",
        "fastp_option",
        "rrna_filter",
        "rrna_filter_sensitivity",
        "rrna_filter_max_seqs",
        "rrna_filter_chunk_spots",
        "rrna_filter_memory_limit",
        "filter_order",
        "contam_filter",
        "contam_filter_rank",
        "contam_filter_db_name",
        "contam_filter_db",
        "contam_filter_sensitivity",
        "contam_filter_max_seqs",
        "read_name",
        "tol",
    ]
    metadata_names = [
        "total_bases",
        "total_spots",
        "spot_length",
        "scientific_name",
        "taxid",
        "private_file",
        "read1_path",
        "read2_path",
    ]
    ind_sra = sra_stat.get("metadata_idx")
    if ind_sra is None:
        ind_sra = get_metadata_row_index_by_run(run_metadata, sra_stat["sra_id"])
    payload = {
        "schema_version": GETFASTQ_RESUME_SCHEMA_VERSION,
        "run": sra_stat["sra_id"],
        "layout": sra_stat["layout"],
        "target_bp_per_run": int(g["num_bp_per_sra"]),
        "options": {
            name: _normalize_getfastq_resume_value(
                getattr(args, name, False if name == "treat_identical_paired_as_single" else None)
            )
            for name in option_names
        },
        "metadata": {
            name: _normalize_getfastq_resume_value(run_metadata.df.at[ind_sra, name])
            for name in metadata_names
            if name in run_metadata.df.columns
        },
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def get_getfastq_run_state_path(run_dir: str) -> str:
    return os.path.join(run_dir, GETFASTQ_RUN_STATE_FILENAME)


def _atomic_write_json(payload: Mapping[str, Any], output_path: str) -> None:
    with atomic_output_path(outpath=output_path, suffix=".json") as tmp_path:
        with open(tmp_path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")


def read_getfastq_run_state(run_dir: str) -> dict[str, Any] | None:
    state_path = get_getfastq_run_state_path(run_dir)
    if not os.path.isfile(state_path):
        return None
    try:
        with open(state_path, encoding="utf-8") as handle:
            state = json.load(handle)
    except Exception as exc:
        raise ValueError("Failed to read getfastq resume state {}: {}".format(state_path, exc)) from exc
    if not isinstance(state, dict):
        raise ValueError("Invalid getfastq resume state (expected an object): {}".format(state_path))
    return state


def _get_getfastq_output_candidates(sra_stat: Mapping[str, Any]) -> list[tuple[str, str]]:
    sra_id = sra_stat["sra_id"]
    if sra_stat["layout"] == "single":
        suffixes = [""]
    elif sra_stat["layout"] == "paired":
        suffixes = ["_1", "_2"]
    else:
        raise ValueError("Unsupported library layout for getfastq resume: {}".format(sra_stat["layout"]))
    return [
        (sra_id + suffix + ".amalgkit.fastq.gz", sra_id + suffix + ".amalgkit.fastq.gz.safely_removed")
        for suffix in suffixes
    ]


def _snapshot_getfastq_outputs(run_dir: str, output_names: Sequence[str]) -> list[dict[str, Any]]:
    snapshots = []
    for output_name in output_names:
        output_path = os.path.join(run_dir, output_name)
        try:
            path_stat = os.stat(output_path, follow_symlinks=False)
            stat_result = os.stat(output_path, follow_symlinks=True)
        except FileNotFoundError as exc:
            raise FileNotFoundError("getfastq output is missing or is not a file: {}".format(output_path)) from exc
        if not stat.S_ISREG(stat_result.st_mode):
            raise OSError("getfastq output does not resolve to a regular file: {}".format(output_path))
        snapshots.append(
            {
                "name": output_name,
                "dev": int(path_stat.st_dev),
                "inode": int(path_stat.st_ino),
                "size": int(stat_result.st_size),
                "mtime_ns": int(stat_result.st_mtime_ns),
                "ctime_ns": int(path_stat.st_ctime_ns),
                "target_dev": int(stat_result.st_dev),
                "target_inode": int(stat_result.st_ino),
                "target_ctime_ns": int(stat_result.st_ctime_ns),
            }
        )
    return snapshots
