"""Transactional retirement of managed FASTQ entries, including private links."""

from __future__ import annotations

import os
import shutil
import tempfile
import warnings
from collections.abc import Callable, Sequence
from contextlib import AbstractContextManager
from dataclasses import dataclass

from amalgkit.output_utils import atomic_output_path


@dataclass
class _Removal:
    input_path: str
    marker_path: str
    had_marker: bool
    quarantine_path: str = ""
    backup_path: str = ""
    marker_backup: str = ""


def _rollback_fastq_removal(records: Sequence[_Removal], quarantine_dir: str) -> None:
    failures = []
    for record in records:
        try:
            if os.path.lexists(record.marker_backup):
                os.replace(record.marker_backup, record.marker_path)
            elif not record.had_marker and os.path.lexists(record.marker_path):
                os.remove(record.marker_path)
        except OSError as exc:
            failures.append(exc)
    for record in reversed(records):
        if os.path.lexists(record.input_path):
            continue
        source = record.quarantine_path if os.path.lexists(record.quarantine_path) else record.backup_path
        try:
            os.replace(source, record.input_path)
        except OSError as exc:
            failures.append(exc)
    if failures:
        # Keep every remaining input/marker backup available for recovery.
        raise ExceptionGroup(f"FASTQ rollback failed; recover files from {quarantine_dir}", failures)
    shutil.rmtree(quarantine_dir)


def safely_remove_quant_fastq_files(
    in_files: Sequence[str],
    *,
    atomic_writer: Callable[..., AbstractContextManager[str]] = atomic_output_path,
) -> None:
    """Replace managed input entries with markers without following symlinks.

    Private gzip inputs are links owned by getfastq, not owned source files.
    Renaming/unlinking the link retires only the managed entry. Backups retain
    the link text (including relative targets) so rollback restores its identity.
    """
    inputs = [os.path.abspath(path) for path in in_files]
    if not inputs:
        return
    if len(inputs) != len(set(inputs)):
        raise ValueError("FASTQ cleanup inputs must be distinct paths.")
    parents = {os.path.dirname(path) for path in inputs}
    if len(parents) != 1:
        raise ValueError("FASTQ cleanup requires all inputs to share one directory.")
    records = []
    for path in inputs:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"FASTQ cleanup input does not resolve to a regular file: {path}")
        marker = path + ".safely_removed"
        if os.path.lexists(marker) and (os.path.islink(marker) or not os.path.isfile(marker)):
            raise ValueError(f"Refusing non-regular safe-removal marker: {marker}")
        records.append(_Removal(path, marker, os.path.lexists(marker)))
    quarantine_dir = tempfile.mkdtemp(prefix=".amalgkit_quant_cleanup_", dir=parents.pop())
    for index, record in enumerate(records):
        record.quarantine_path = os.path.join(quarantine_dir, f"input_{index}")
        record.backup_path = os.path.join(quarantine_dir, f"backup_{index}")
        record.marker_backup = os.path.join(quarantine_dir, f"marker_{index}")
    try:
        for record in records:
            os.replace(record.input_path, record.quarantine_path)
        for record in records:
            if os.path.islink(record.quarantine_path):
                os.symlink(os.readlink(record.quarantine_path), record.backup_path)
            else:
                os.link(record.quarantine_path, record.backup_path, follow_symlinks=False)
        for record in records:
            if record.had_marker:
                os.replace(record.marker_path, record.marker_backup)
        for record in records:
            with atomic_writer(record.marker_path, suffix=".safely_removed") as temporary_marker:
                with open(temporary_marker, "w", encoding="utf-8") as handle:
                    handle.write("This fastq file was safely removed after `amalgkit quant`.")
        for record in records:
            os.remove(record.quarantine_path)
    except BaseException as exc:
        try:
            _rollback_fastq_removal(records, quarantine_dir)
        except BaseException as rollback_exc:
            exc.add_note(f"FASTQ rollback failed: {rollback_exc}")
            exc.add_note(f"FASTQ recovery files preserved at: {quarantine_dir}")
        raise
    try:
        shutil.rmtree(quarantine_dir)
    except OSError:
        warnings.warn(
            f"FASTQ cleanup committed, but quarantine cleanup was incomplete: {quarantine_dir}",
            RuntimeWarning,
            stacklevel=2,
        )
