"""Validated, reusable public GSA FASTQ inputs and paired spot-range extraction."""

import bz2
from contextlib import ExitStack, closing, contextmanager
import gzip
import hashlib
import itertools
import json
import os
import re
import stat
import time
import urllib.parse
import zlib

import pandas

from amalgkit.fastq_utils import is_private_file_value
from amalgkit.download_utils import acquire_exclusive_lock, calculate_file_md5, resolve_download_dir
from amalgkit.fastq_download_integrity import validate_integrity_metadata
from amalgkit.output_utils import atomic_output_path
from amalgkit.runtime_utils import safe_join_component, validate_safe_path_component

CACHE_SCHEMA = 1


class GsaCorruptInputError(ValueError):
    def __init__(self, path, reason):
        self.path = path
        super().__init__("{}: {}".format(reason, path))


def is_gsa_row(row):
    return str(row.get("data_source", "")).strip().lower() == "gsa"


def read_manifest(row):
    run = str(row["run"])
    if not re.fullmatch(r"CRR\d+", run):
        raise ValueError("GSA metadata requires a CRR Run accession: {}".format(run))
    if is_private_file_value(row.get("private_file", "")):
        raise ValueError("GSA public metadata cannot be marked private: {}".format(run))
    if str(row.get("lib_layout", "")).lower() not in {"single", "paired"}:
        raise ValueError("Unsupported GSA library layout for {}".format(run))
    instrument = str(row.get("instrument", "")) + " " + str(row.get("platform", ""))
    if re.search(r"pacbio|pacific biosciences|nanopore|\bONT\b", instrument, re.I):
        raise ValueError("GSA native FASTQ currently supports short-read data only: {}".format(run))
    try:
        files = json.loads(row["gsa_fastq_files"])
    except (KeyError, ValueError, TypeError) as exc:
        raise ValueError("Missing or invalid GSA FASTQ manifest for {}".format(run)) from exc
    if not isinstance(files, list) or not files:
        raise ValueError("Empty GSA FASTQ manifest for {}".format(run))
    groups: dict[int, set[int]] = {}
    filenames = set()
    for entry in files:
        if not isinstance(entry, dict):
            raise ValueError("Invalid GSA FASTQ file entry")
        validate_integrity_metadata(entry)
        name = entry.get("filename")
        if not isinstance(name, str):
            raise ValueError("GSA FASTQ filename must be a string")
        validate_safe_path_component(name, label="GSA FASTQ filename")
        if name in filenames or not re.search(r"\.(fastq|fq)(\.(gz|bz2))?$", name, re.I):
            raise ValueError("Duplicate or unsupported GSA FASTQ filename: {}".format(name))
        filenames.add(name)
        group, mate = entry.get("group"), entry.get("mate")
        if type(group) is not int or group < 0 or type(mate) is not int or mate not in {0, 1, 2}:
            raise ValueError("Invalid GSA file group/mate")
        if mate in groups.setdefault(group, set()):
            raise ValueError("Duplicate GSA file group/mate")
        groups[group].add(mate)
        sources = entry.get("sources")
        if not isinstance(sources, list) or not sources:
            raise ValueError("Missing GSA download source")
        for source in sources:
            url = urllib.parse.urlparse(source.get("url", ""))
            if (
                source.get("source_name") != "GSA"
                or url.scheme != "https"
                or url.hostname not in {"download.cncb.ac.cn", "download.big.ac.cn"}
                or url.username
                or url.password
                or url.port not in {None, 443}
                or run not in url.path.split("/")
                or urllib.parse.unquote(url.path.rsplit("/", 1)[-1]) != name
            ):
                raise ValueError("Invalid GSA FASTQ source for {}".format(run))
    expected = {1, 2} if str(row["lib_layout"]).lower() == "paired" else {0}
    if set(groups) != set(range(len(groups))) or any(mates != expected for mates in groups.values()):
        raise ValueError("Incomplete GSA FASTQ file groups for {}".format(run))
    return sorted(files, key=lambda entry: (entry["group"], entry["mate"]))


def manifest_fingerprint(row):
    manifest = {"schema": CACHE_SCHEMA, "layout": str(row["lib_layout"]).lower(), "files": read_manifest(row)}
    return hashlib.sha256(json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def cache_directory(args, row):
    root = resolve_download_dir(args)
    for name in ("gsa", str(row["run"]), manifest_fingerprint(row)):
        root = safe_join_component(root, name, label="GSA cache component")
        if os.path.islink(root):
            raise ValueError("GSA cache directory must not be a symbolic link: {}".format(root))
        os.makedirs(root, exist_ok=True)
    return root


def _file_identity(path):
    item = os.stat(path, follow_symlinks=False)
    if not stat.S_ISREG(item.st_mode):
        raise ValueError("GSA cache payload must be a regular file: {}".format(path))
    return {"size": item.st_size, "mtime_ns": item.st_mtime_ns, "ctime_ns": item.st_ctime_ns}


def _open_fastq(path):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rb")
    if path.lower().endswith(".bz2"):
        return bz2.open(path, "rb")
    return open(path, "rb")


def _records(handle, label):
    try:
        while header := handle.readline():
            sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
            seq, qual = sequence.rstrip(b"\r\n"), quality.rstrip(b"\r\n")
            if (
                not header.startswith(b"@")
                or len(header.split(None, 1)[0]) <= 1
                or not plus.startswith(b"+")
                or not seq
                or len(seq) != len(qual)
                or not quality
            ):
                raise GsaCorruptInputError(label, "Malformed or truncated GSA FASTQ")
            # A missing final newline is harmless within one source file, but
            # would join its quality to the next group's header on extraction.
            yield (header, sequence, plus, quality if quality.endswith(b"\n") else quality + b"\n")
    except (OSError, EOFError, zlib.error) as exc:
        raise GsaCorruptInputError(label, "Invalid GSA FASTQ compression stream") from exc


def _pair_key(header):
    token = header.split(None, 1)[0]
    return re.sub(rb"/[12]$", b"", token)


def iter_spots(cache_dir, files):
    for _, entry_group in itertools.groupby(files, key=lambda entry: entry["group"]):
        entries = list(entry_group)
        with ExitStack() as stack:
            paths = [os.path.join(cache_dir, entry["filename"]) for entry in entries]
            streams = [_records(stack.enter_context(_open_fastq(path)), path) for path in paths]
            for spot in itertools.zip_longest(*streams):
                if any(record is None for record in spot):
                    raise ValueError("GSA paired FASTQ record counts differ: {}".format(paths))
                if len(spot) == 2:
                    if _pair_key(spot[0][0]) != _pair_key(spot[1][0]):
                        raise ValueError("GSA paired FASTQ read IDs/order differ: {}".format(paths))
                    for mate, record in enumerate(spot, 1):
                        fields = record[0].split()
                        marker = re.search(rb"/([12])$", fields[0])
                        if marker and int(marker[1]) != mate:
                            raise ValueError("GSA FASTQ header contradicts its mate assignment: {}".format(paths))
                        if len(fields) > 1:
                            marker = re.match(rb"([12]):", fields[1])
                            if marker and int(marker[1]) != mate:
                                raise ValueError("GSA FASTQ header contradicts its mate assignment: {}".format(paths))
                yield spot


def _validate_file_bytes(path, entry):
    identity = _file_identity(path)
    if identity["size"] <= 0:
        raise GsaCorruptInputError(path, "Empty GSA FASTQ file")
    if entry.get("expected_bytes") is not None and identity["size"] != entry["expected_bytes"]:
        raise GsaCorruptInputError(path, "GSA FASTQ byte count mismatch")
    if entry.get("expected_md5") is not None and calculate_file_md5(path) != entry["expected_md5"]:
        raise GsaCorruptInputError(path, "GSA FASTQ checksum mismatch")


def _scan_inputs(directory, files):
    spots, bases = 0, 0
    for spot in iter_spots(directory, files):
        spots += 1
        bases += sum(len(record[1].rstrip(b"\r\n")) for record in spot)
    if not spots or not bases:
        raise ValueError("GSA Run contains no FASTQ records")
    return {"total_spots": spots, "total_bases": bases, "spot_length": bases / spots, "read_count_status": "measured"}


def _read_cached_stats(directory, files):
    path = os.path.join(directory, "validated.json")
    if not os.path.exists(path):
        return None
    _file_identity(path)
    try:
        with open(path) as handle:
            result = json.load(handle)
        identities = {entry["filename"]: _file_identity(os.path.join(directory, entry["filename"])) for entry in files}
        if result["schema"] != CACHE_SCHEMA or result["files"] != identities:
            return None
        stats = result["stats"]
        if not re.fullmatch(r"[0-9a-f]{64}", str(stats.get("gsa_input_fingerprint", ""))):
            return None
        for key in ("total_spots", "total_bases", "spot_length"):
            if not isinstance(stats[key], (int, float)) or not 0 < stats[key] < float("inf"):
                return None
        return stats
    except (OSError, KeyError, TypeError, ValueError):
        return None


def prepare_run(args, row, download):
    """Download and validate once under a cache lock; retain originals for round 2."""
    files = read_manifest(row)
    directory = cache_directory(args, row)
    with acquire_exclusive_lock(os.path.join(directory, "input.lock"), lock_label="GSA original FASTQ"):
        cached = _read_cached_stats(directory, files)
        if cached is not None:
            print("Validated GSA input cache: {}".format(row["run"]), flush=True)
            return cached
        retried_files = set()
        for _ in range(len(files) + 1):
            try:
                for entry in files:
                    path = os.path.join(directory, entry["filename"])
                    partial = path + ".part"
                    for candidate in (path, partial):
                        if os.path.lexists(candidate):
                            _file_identity(candidate)
                    if not os.path.exists(path):
                        if not download(
                            sra_id=row["run"],
                            source_candidates=entry["sources"],
                            output_path=partial,
                            args=args,
                            artifact_label="GSA original FASTQ",
                            resume_existing=True,
                        ):
                            raise OSError("GSA FASTQ download failed: {}".format(row["run"]))
                        _validate_file_bytes(partial, entry)
                        os.replace(partial, path)
                    _validate_file_bytes(path, entry)
                stats = _scan_inputs(directory, files)
                break
            except GsaCorruptInputError as exc:
                owned_paths = {
                    os.path.join(directory, entry["filename"]) + suffix for entry in files for suffix in ("", ".part")
                }
                if exc.path not in owned_paths:
                    raise
                _file_identity(exc.path)
                os.remove(exc.path)
                filename = exc.path.removesuffix(".part")
                if filename in retried_files:
                    raise
                retried_files.add(filename)
                print("Corrupt GSA FASTQ; retrying the affected file once: {}".format(exc.path), flush=True)
        digests = {}
        for entry in files:
            with open(os.path.join(directory, entry["filename"]), "rb") as handle:
                digests[entry["filename"]] = hashlib.file_digest(handle, "sha256").hexdigest()
        stats["gsa_input_fingerprint"] = hashlib.sha256(json.dumps(digests, sort_keys=True).encode()).hexdigest()
        result = {
            "schema": CACHE_SCHEMA,
            "stats": stats,
            "files": {entry["filename"]: _file_identity(os.path.join(directory, entry["filename"])) for entry in files},
        }
        with atomic_output_path(os.path.join(directory, "validated.json")) as temporary:
            with open(temporary, "w") as handle:
                json.dump(result, handle, sort_keys=True)
        return stats


def prepare_gsa_metadata(args, metadata, download):
    """Resolve true input sizes before getfastq allocates its base budget."""
    if "data_source" not in metadata.df.columns:
        return metadata
    for key in ("total_spots", "total_bases", "spot_length"):
        # Provider metadata leaves counts blank. Coerce before assigning numbers
        # to avoid object/string dtype differences across supported pandas.
        metadata.df[key] = pandas.to_numeric(metadata.df[key], errors="coerce")
    metadata.df["gsa_input_seconds"] = pandas.to_numeric(
        metadata.df.get("gsa_input_seconds", pandas.Series(index=metadata.df.index, dtype=float)), errors="coerce"
    ).astype(float)
    for index, row in metadata.df.iterrows():
        if not is_gsa_row(row):
            continue
        started = time.perf_counter()
        stats = prepare_run(args, row, download)
        for key, value in stats.items():
            if key not in metadata.df:
                metadata.df[key] = ""
            metadata.df.at[index, key] = value
        metadata.df.at[index, "gsa_input_seconds"] = time.perf_counter() - started
        print(
            "Prepared GSA input {}: {:,} spots, {:,} bp; {:.1f} sec (download/cache and validation).".format(
                row["run"], stats["total_spots"], stats["total_bases"], metadata.df.at[index, "gsa_input_seconds"]
            ),
            flush=True,
        )
    return metadata


def _require_validated_cache(row, directory, files):
    cached = _read_cached_stats(directory, files)
    if cached is None:
        raise ValueError("GSA input cache changed or has not been validated: {}".format(row["run"]))
    if row.get("gsa_input_fingerprint") != cached["gsa_input_fingerprint"]:
        raise ValueError("GSA input cache does not match the measured input fingerprint: {}".format(row["run"]))


@contextmanager
def open_validated_spots(args, row):
    """Hold the input lock while streaming every group from the measured cache."""
    directory = cache_directory(args, row)
    with acquire_exclusive_lock(os.path.join(directory, "input.lock"), lock_label="GSA original FASTQ"):
        files = read_manifest(row)
        _require_validated_cache(row, directory, files)
        with closing(iter_spots(directory, files)) as spots:
            yield spots


def extract_run(args, row, work_dir, start, end, return_stats=False):
    directory = cache_directory(args, row)
    with acquire_exclusive_lock(os.path.join(directory, "input.lock"), lock_label="GSA original FASTQ"):
        return _extract_run_locked(args, row, work_dir, start, end, return_stats=return_stats)


def _extract_run_locked(args, row, work_dir, start, end, return_stats=False):
    files = read_manifest(row)
    directory = cache_directory(args, row)
    _require_validated_cache(row, directory, files)
    start, end = int(start), int(end)
    if start < 1 or end < start:
        raise ValueError("Invalid GSA spot range")
    paired = str(row["lib_layout"]).lower() == "paired"
    paths = [
        os.path.join(work_dir, str(row["run"]) + suffix + ".fastq.gz") for suffix in (["_1", "_2"] if paired else [""])
    ]
    spots, bases, dumped_spots, dumped_bases = 0, 0, 0, 0
    min_length = int(getattr(args, "min_read_length", 0))
    if min_length < 0:
        raise ValueError("--min_read_length must be >= 0.")
    with ExitStack() as stack:
        outputs = [
            stack.enter_context(gzip.open(stack.enter_context(atomic_output_path(path)), "wb", compresslevel=1))
            for path in paths
        ]
        for number, spot in enumerate(iter_spots(directory, files), 1):
            if number < start:
                continue
            if number > end:
                break
            lengths = [len(record[1].rstrip(b"\r\n")) for record in spot]
            dumped_spots += 1
            dumped_bases += sum(lengths)
            if any(length < min_length for length in lengths):
                continue
            for output, record in zip(outputs, spot, strict=True):
                output.writelines(record)
                bases += len(record[1].rstrip(b"\r\n"))
            spots += 1
        if dumped_spots != end - start + 1:
            raise ValueError("Requested GSA spot range exceeds validated input")
    if return_stats:
        return {
            "num_dumped": dumped_spots,
            "num_written": spots,
            "num_rejected": dumped_spots - spots,
            "bp_dumped": dumped_bases,
            "bp_written": bases,
            "bp_rejected": dumped_bases - bases,
        }
    return spots, bases
