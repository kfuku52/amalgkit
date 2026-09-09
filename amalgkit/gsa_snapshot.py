"""Publish measured GSA metadata without changing source tables or batch indices."""

import hashlib
import json
import os
import re

import pandas

from amalgkit.download_utils import acquire_exclusive_lock
from amalgkit.gsa_fastq import is_gsa_row, manifest_fingerprint
from amalgkit.gsa_select import resolve_deferred_selection, validate_selection_ready
from amalgkit.metadata_utils import Metadata
from amalgkit.output_utils import atomic_output_path, atomic_write_dataframe

MEASURED_COLUMNS = (
    "total_spots",
    "total_bases",
    "spot_length",
    "read_count_status",
    "gsa_input_seconds",
    "gsa_input_fingerprint",
)


def capture_gsa_file_source(args, metadata, path, digest, source_bytes=None):
    """Keep the table and digest from the same read, before filtering/batching."""
    args._gsa_snapshot_source = {"table": metadata.df.copy(deep=True), "path": path, "digest": digest}
    _, snapshot_path, state_path = _paths(args)
    if os.path.realpath(path) == os.path.realpath(snapshot_path):
        _capture_snapshot_input(args, path, digest, source_bytes, state_path)


def _capture_snapshot_input(args, path, digest, source_bytes, state_path):
    directory = os.path.dirname(path)
    with acquire_exclusive_lock(os.path.join(directory, "gsa_metadata.lock"), lock_label="GSA measured metadata"):
        if _digest(path) != digest:
            raise ValueError("GSA metadata source changed since loading; rerun getfastq with the current source")
        state = {}
        if os.path.isfile(state_path):
            with open(state_path) as handle:
                state = json.load(handle)
            if not isinstance(state, dict):
                raise ValueError("Invalid GSA metadata snapshot state")
        previous = state.get("snapshot_input_digest")
        reuse = (
            isinstance(previous, str)
            and re.fullmatch(r"[0-9a-f]{64}", previous)
            and state.get("snapshot_digest") == digest
        )
        frozen_digest = digest
        if reuse and isinstance(previous, str):
            frozen_digest = previous
        frozen_dir = os.path.join(directory, "gsa_input_metadata")
        frozen_path = os.path.join(frozen_dir, frozen_digest + ".tsv")
        if os.path.islink(frozen_dir) or os.path.islink(frozen_path):
            raise ValueError("GSA frozen metadata input must not be a symbolic link")
        os.makedirs(frozen_dir, exist_ok=True)
        if reuse and not os.path.isfile(frozen_path):
            raise ValueError("GSA frozen metadata input is missing")
        if not os.path.exists(frozen_path):
            if source_bytes is None:
                with open(path, "rb") as handle:
                    source_bytes = handle.read()
            if hashlib.sha256(source_bytes).hexdigest() != frozen_digest:
                raise ValueError("GSA metadata source changed since loading")
            with atomic_output_path(frozen_path) as temporary:
                with open(temporary, "wb") as handle:
                    handle.write(source_bytes)
        if _digest(frozen_path) != frozen_digest:
            raise ValueError("GSA frozen metadata input is invalid")
        args._gsa_snapshot_source = {
            "table": Metadata.from_DataFrame(_read_table(frozen_path)).df,
            "path": frozen_path,
            "digest": frozen_digest,
            "snapshot_input": path,
            "snapshot_state": state_path,
            "snapshot_input_digest": frozen_digest,
        }


def capture_gsa_accession_source(args, metadata, requested_ids):
    """All resolved rows define one input generation shared by accession jobs."""
    table = metadata.df.copy(deep=True)
    canonical = table.drop(columns=["gsa_retrieved_at"], errors="ignore").fillna("").astype(str)
    payload = {
        "ids": list(requested_ids),
        "layout": getattr(args, "layout", None),
        "scientific_name": getattr(args, "sci_name", None),
        "rows": canonical.reindex(sorted(canonical.columns), axis=1).to_dict("records"),
    }
    digest = "accessions:" + hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    args._gsa_snapshot_source = {
        "table": table,
        "path": None,
        "digest": digest,
        "id_list": getattr(args, "id_list", None),
        "requested_ids": list(requested_ids),
    }


def _validate_source(source):
    if source["path"] is not None and _digest(source["path"]) != source["digest"]:
        raise ValueError("GSA metadata source changed since loading; rerun getfastq with the current source")
    if source.get("snapshot_input"):
        current = _digest(source["snapshot_input"])
        if current != source["digest"]:
            with open(source["snapshot_state"]) as handle:
                state = json.load(handle)
            if (
                not isinstance(state, dict)
                or state.get("snapshot_input_digest") != source["digest"]
                or state.get("snapshot_digest") != current
            ):
                raise ValueError("GSA metadata source changed since loading; rerun getfastq with the current source")
    if source.get("id_list") is not None:
        with open(source["id_list"]) as handle:
            ids = [line.strip() for line in handle if line.strip() and not line.lstrip().startswith("#")]
        if ids != source["requested_ids"]:
            raise ValueError("GSA accession list changed since loading; rerun getfastq with the current source")


def _digest(path):
    with open(path, "rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _paths(args):
    directory = os.path.join(args.out_dir, "getfastq")
    paths = (directory, os.path.join(directory, "metadata.tsv"), os.path.join(directory, "gsa_metadata_state.json"))
    if any(os.path.islink(path) for path in paths):
        raise ValueError("GSA metadata snapshot paths must not be symbolic links")
    return paths


def _read_table(path):
    return pandas.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _valid_state(state_path, source_digest, snapshot_path):
    if not os.path.isfile(state_path):
        return False
    with open(state_path) as handle:
        state = json.load(handle)
    if not isinstance(state, dict):
        raise ValueError("Invalid GSA metadata snapshot state; rerun getfastq to rebuild it")
    if state.get("source_digest") != source_digest:
        return False
    if not os.path.isfile(snapshot_path) or state.get("snapshot_digest") != _digest(snapshot_path):
        raise ValueError("GSA measured metadata snapshot is invalid; rerun getfastq to rebuild it")
    return True


def publish_gsa_snapshot(args, metadata):
    source = getattr(args, "_gsa_snapshot_source", None)
    if source is None:
        raise ValueError("GSA snapshot requires the metadata source captured before processing")
    directory, snapshot_path, state_path = _paths(args)
    os.makedirs(directory, exist_ok=True)
    source_digest = source["digest"]
    with acquire_exclusive_lock(os.path.join(directory, "gsa_metadata.lock"), lock_label="GSA measured metadata"):
        _validate_source(source)
        combined = source["table"].copy(deep=True)
        if combined["run"].duplicated().any():
            raise ValueError("Duplicate Run IDs in GSA metadata snapshot source")
        combined = combined.set_index("run", drop=False)
        try:
            valid_snapshot = _valid_state(state_path, source_digest, snapshot_path)
        except (OSError, ValueError, TypeError):
            valid_snapshot = False
        if valid_snapshot:
            _overlay_measurements(combined, _read_table(snapshot_path))
        _overlay_measurements(combined, metadata.df)
        combined = Metadata.from_DataFrame(combined.reset_index(drop=True))
        combined = resolve_deferred_selection(combined)
        _validate_source(source)
        atomic_write_dataframe(combined.df, snapshot_path, sep="\t", index=False)
        with atomic_output_path(state_path) as temporary:
            with open(temporary, "w") as handle:
                json.dump(
                    {
                        "source_digest": source_digest,
                        "snapshot_digest": _digest(snapshot_path),
                        "snapshot_input_digest": source.get("snapshot_input_digest"),
                    },
                    handle,
                    sort_keys=True,
                )
        by_run = combined.df.set_index("run")
        for index, row in metadata.df.iterrows():
            for column in ("exclusion", "gsa_deferred_select", "gsa_selection_status"):
                if column in by_run:
                    if column not in metadata.df:
                        metadata.df[column] = ""
                    metadata.df.at[index, column] = by_run.at[row["run"], column]
    print("Measured GSA metadata: {}".format(snapshot_path), flush=True)
    return metadata


def _overlay_measurements(combined, measured):
    """Reapply original annotations/rules; only measured inputs survive jobs."""
    for _, row in measured.iterrows():
        if row["run"] not in combined.index:
            raise ValueError("GSA snapshot source does not contain {}".format(row["run"]))
        if not is_gsa_row(row) or row.get("read_count_status") != "measured":
            continue
        original = combined.loc[row["run"]]
        if not is_gsa_row(original) or manifest_fingerprint(original) != manifest_fingerprint(row):
            raise ValueError("GSA measurement does not match the captured input manifest: {}".format(row["run"]))
        for column in MEASURED_COLUMNS:
            if column not in row:
                continue
            if column not in combined:
                combined[column] = ""
            combined[column] = combined[column].astype(object)
            combined.at[row["run"], column] = row[column]


def preferred_gsa_snapshot(args, source_path, read_table=False):
    """Use a verified snapshot only for an unchanged inferred metadata source."""
    directory, snapshot_path, state_path = _paths(args)
    if not os.path.isfile(state_path):
        return None if read_table else source_path
    with acquire_exclusive_lock(os.path.join(directory, "gsa_metadata.lock"), lock_label="GSA measured metadata"):
        if not _valid_state(state_path, _digest(source_path), snapshot_path):
            return None if read_table else source_path
        table = _read_table(snapshot_path)
        validate_selection_ready(table)
        if read_table:
            print("Using measured GSA metadata snapshot: {}".format(snapshot_path), flush=True)
            return table
        return snapshot_path
