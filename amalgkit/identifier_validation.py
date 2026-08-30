"""Exact duplicate detection with a bounded in-memory identifier inventory."""

from __future__ import annotations

import os
import sqlite3
import tempfile
from collections.abc import Sequence
from types import TracebackType


class TargetIdTracker:
    """Keep small tables in memory and spill large inventories to a private DB.

    The SQLite file is scratch data, never an output. Both its page cache and
    query batches are bounded; no probabilistic duplicate checks are used.
    """

    def __init__(self, memory_limit: int = 50_000) -> None:
        self.memory_limit = memory_limit
        self._seen: set[str] = set()
        self._scratch: tempfile.TemporaryDirectory[str] | None = None
        self._db: sqlite3.Connection | None = None

    def __enter__(self) -> TargetIdTracker:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        if self._db is not None:
            self._db.close()
        if self._scratch is not None:
            self._scratch.cleanup()

    def _spill(self) -> None:
        self._scratch = tempfile.TemporaryDirectory(prefix="amalgkit-ids-")
        self._db = sqlite3.connect(os.path.join(self._scratch.name, "ids.sqlite3"))
        self._db.execute("PRAGMA cache_size=-2048")
        # This disposable inventory needs no crash recovery or rollback.
        self._db.execute("PRAGMA journal_mode=OFF")
        self._db.execute("PRAGMA synchronous=OFF")
        self._db.execute("CREATE TABLE ids (value TEXT PRIMARY KEY) WITHOUT ROWID")
        self._db.executemany("INSERT INTO ids VALUES (?)", ((value,) for value in self._seen))
        self._seen.clear()

    def add(self, values: Sequence[str]) -> set[str]:
        """Return identifiers already seen; callers check within-chunk duplicates."""
        if self._db is None and len(self._seen) + len(values) <= self.memory_limit:
            duplicates = self._seen.intersection(values)
            self._seen.update(values)
            return duplicates
        if self._db is None:
            self._spill()
        assert self._db is not None  # noqa: S101 - established by _spill
        duplicates = set()
        for start in range(0, len(values), 500):
            batch = values[start : start + 500]
            placeholders = ",".join("?" for _ in batch)
            query = f"SELECT value FROM ids WHERE value IN ({placeholders})"  # noqa: S608 - placeholders only
            duplicates.update(row[0] for row in self._db.execute(query, batch))
        if not duplicates:
            self._db.executemany("INSERT INTO ids VALUES (?)", ((value,) for value in values))
        return duplicates
