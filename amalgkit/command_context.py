import os
from dataclasses import dataclass, field


@dataclass
class PerSpeciesTableContext:
    metadata: object = None
    input_dir: str | None = None

    def resolve(self):
        if self.metadata is None:
            raise ValueError('PerSpeciesTableContext.metadata is required.')
        if self.input_dir is None:
            raise ValueError('PerSpeciesTableContext.input_dir is required.')
        return self.metadata, self.input_dir


@dataclass
class CrossSpeciesFilterContext:
    metadata: object = None


@dataclass
class PrefetchedDirEntries:
    entries: object = None
    entries_sorted: list | tuple | None = None
    path_dir: str | None = None

    @classmethod
    def from_entries(cls, entries=None, path_dir=None):
        sorted_entries = None
        if isinstance(entries, (set, list, tuple)):
            sorted_entries = sorted(entries)
        real_path = os.path.realpath(path_dir) if isinstance(path_dir, str) else None
        return cls(entries=entries, entries_sorted=sorted_entries, path_dir=real_path)

    def resolve_entries(self, path_dir):
        if (not isinstance(self.path_dir, str)) or (os.path.realpath(path_dir) != self.path_dir):
            return None
        if isinstance(self.entries_sorted, (list, tuple)):
            return self.entries_sorted
        if isinstance(self.entries, (set, list, tuple)):
            return self.entries
        return None


@dataclass
class QuantRuntimeContext:
    run_files_by_run: dict = field(default_factory=dict)
    quant_backend_by_run: dict = field(default_factory=dict)
    oarfish_seq_tech_by_run: dict = field(default_factory=dict)
    resolved_index_cache: dict = field(default_factory=dict)
    prefetched_fasta: PrefetchedDirEntries = field(default_factory=PrefetchedDirEntries)
    prefetched_index: PrefetchedDirEntries = field(default_factory=PrefetchedDirEntries)


@dataclass
class GetfastqRuntimeContext:
    mmseqs_dbtype_cache: dict = field(default_factory=dict)
