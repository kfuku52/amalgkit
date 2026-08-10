"""In-memory file inventory for one getfastq run directory."""

from __future__ import annotations

import os
from collections.abc import Callable, Iterable


def list_run_dir_files(work_dir: str) -> set[str]:
    """List regular files in a run directory without failing if it vanished."""
    try:
        with os.scandir(work_dir) as entries:
            return {entry.name for entry in entries if entry.is_file()}
    except FileNotFoundError:
        return set()


class RunFileState:
    """Mutable snapshot that avoids repeated directory scans within a run."""

    def __init__(
        self,
        work_dir: str,
        files: Iterable[str] | None = None,
        file_lister: Callable[[str], set[str]] = list_run_dir_files,
    ) -> None:
        self.work_dir = work_dir
        self.files = file_lister(work_dir) if files is None else set(files)

    def has(self, filename: str) -> bool:
        return filename in self.files

    def path(self, filename: str) -> str:
        return os.path.join(self.work_dir, filename)

    def add(self, filename: str) -> None:
        self.files.add(filename)

    def discard(self, filename: str) -> None:
        self.files.discard(filename)

    def to_set(self) -> set[str]:
        return set(self.files)
