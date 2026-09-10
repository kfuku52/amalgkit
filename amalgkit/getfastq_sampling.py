"""Reproducible, non-replacement sampling of validated FASTQ spots.

Selection precedes length filtering. A spot is one SE record or an intact PE
pair; it does not imply an independent molecule. The versioned SHA-256 random
stream and partial Fisher-Yates permutation have a prefix independent of the
requested sample size, so subsequent rounds cannot select a spot twice.
"""

from contextlib import ExitStack, closing
import gzip
import hashlib
import itertools
import json
import os
import re

from amalgkit.fastq_utils import open_fastq_binary
from amalgkit.output_utils import atomic_output_path

ALGORITHM = "sha256-fisher-yates-v1"
MANIFEST = "getfastq_sampling.json"
COUNT_COLUMNS = ("num_dumped", "num_written", "num_rejected", "bp_dumped", "bp_written", "bp_rejected")
RANGE_COLUMNS = [
    "getfastq_sampling_start_1st",
    "getfastq_sampling_end_1st",
    "getfastq_sampling_start_2nd",
    "getfastq_sampling_end_2nd",
]
STATS_COLUMNS = [
    "getfastq_sampling_candidate_spots",
    "getfastq_sampling_candidate_bp",
    "getfastq_sampling_method",
    "getfastq_sampling_seed",
    "getfastq_sampling_algorithm",
    "getfastq_sampling_input_sha256",
] + RANGE_COLUMNS


def random_sampling(args, row=None):
    if getattr(args, "sampling_method", "contiguous") != "random":
        return False
    return (
        row is None
        or str(row.get("private_file", "")).lower() != "yes"
        or bool(getattr(args, "sampling_private", False))
    )


def validate_options(args):
    method = getattr(args, "sampling_method", "contiguous")
    if method not in ("contiguous", "random"):
        raise ValueError("--sampling_method must be contiguous or random")
    if getattr(args, "sampling_private", False) and method != "random":
        raise ValueError("--sampling_private yes requires --sampling_method random")
    if method == "random" and getattr(args, "treat_identical_paired_as_single", False):
        raise ValueError("Random sampling requires preserved mates; disable --treat_identical_paired_as_single")
    seed = getattr(args, "sampling_seed", 0)
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("--sampling_seed must be a nonnegative integer")


def validate_private_sources(args, row, run_dir):
    if not random_sampling(args, row) or str(row.get("private_file", "")).lower() != "yes":
        return
    # Run cleanup precedes extraction; never let it remove a declared source.
    for column in ["read1_path", "read2_path"][: 2 if str(row.get("lib_layout", "")).lower() == "paired" else 1]:
        source = str(row.get(column, ""))
        for resolve in (os.path.abspath, os.path.realpath):
            directory = resolve(run_dir)
            if os.path.commonpath([resolve(source), directory]) == directory:
                raise ValueError("Private sampling source must be outside its managed getfastq run directory")


class _RandomStream:
    def __init__(self, seed, run):
        self.key = json.dumps([ALGORITHM, str(seed), str(run)], separators=(",", ":")).encode()
        self.counter = 0

    def below(self, bound):
        # Rejection sampling avoids modulo bias, including non-power-of-two N.
        ceiling = (1 << 256) - ((1 << 256) % bound)
        while True:
            value = int.from_bytes(hashlib.sha256(self.key + b":" + str(self.counter).encode()).digest(), "big")
            self.counter += 1
            if value < ceiling:
                return value % bound


def selected_indices(total, start, end, seed, run):
    """Return 1-based input positions for inclusive ranks in a random permutation."""
    if not 1 <= start <= end <= total or total >= (1 << 256):
        raise ValueError("Invalid sampling rank range or input count")
    if start == 1 and end == total:
        return range(1, total + 1)
    rng = _RandomStream(seed, run)
    swaps = {}
    selected = set()
    for rank in range(end):
        position = rank + rng.below(total - rank)
        value = swaps.get(position, position)
        swaps[position] = swaps.get(rank, rank)
        swaps.pop(rank, None)
        if rank >= start - 1:
            selected.add(value + 1)
    return selected


def _records(handle, path):
    while True:
        header = handle.readline()
        if not header:
            return
        sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
        if (
            not quality
            or not header.startswith(b"@")
            or len(header.split(None, 1)[0]) <= 1
            or not sequence.rstrip(b"\r\n")
            or not plus.startswith(b"+")
            or len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n"))
        ):
            raise ValueError("Malformed FASTQ during sampling: {}".format(path))
        yield header, sequence, plus, quality


def iter_spots(paths):
    """Validate the complete input, including records that are not selected."""
    with ExitStack() as stack:
        streams = [_records(stack.enter_context(open_fastq_binary(path)), path) for path in paths]
        for spot in itertools.zip_longest(*streams):
            if any(record is None for record in spot):
                raise ValueError("Paired FASTQ record counts differ during sampling")
            if len(spot) == 2:
                keys = [re.sub(rb"/[12]$", b"", record[0].split()[0]) for record in spot]
                if keys[0] != keys[1]:
                    raise ValueError("Paired FASTQ read IDs/order differ during sampling")
                for mate, record in enumerate(spot, 1):
                    fields = record[0].split()
                    suffix = re.search(rb"/([12])$", fields[0])
                    marker = re.match(rb"([12]):", fields[1]) if len(fields) > 1 else None
                    if any(match and int(match[1]) != mate for match in (suffix, marker)):
                        raise ValueError("FASTQ header contradicts mate assignment during sampling")
            yield spot


def manifest_digest(run_dir):
    with open(os.path.join(run_dir, MANIFEST), "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def sample_fastqs(paths, output_paths, *, total, start, end, seed, run, min_length, run_dir, provenance=None):
    """Scan all candidates and atomically publish selected FASTQs plus provenance.

    Outputs must be regular paths, never private input symlinks. The outer
    getfastq state machine restarts an incomplete round; the manifest itself is
    not a completion marker.
    """
    if len(paths) not in (1, 2) or len(paths) != len(output_paths):
        raise ValueError("Sampling requires one FASTQ or two aligned mates")
    with closing(iter_spots(paths)) as spots:
        return sample_spots(
            spots, output_paths, total=total, start=start, end=end, seed=seed,
            run=run, min_length=min_length, run_dir=run_dir, provenance=provenance,
        )


def sample_spots(spots, output_paths, *, total, start, end, seed, run, min_length, run_dir, provenance=None):
    """Sample a validated spot stream; its caller owns and closes the input."""
    if len(output_paths) not in (1, 2):
        raise ValueError("Sampling requires one FASTQ or two aligned mates")
    if min_length < 0:
        raise ValueError("--min_read_length must be >= 0")
    selected = selected_indices(total, start, end, seed, run)
    counts = dict.fromkeys(COUNT_COLUMNS, 0)
    input_digest, selection_digest = hashlib.sha256(), hashlib.sha256()
    candidate_count = candidate_bp = 0
    manifest_path = os.path.join(run_dir, MANIFEST)
    previous = None
    if start > 1:
        with open(manifest_path, encoding="utf-8") as handle:
            previous = json.load(handle)
        if (
            previous["algorithm"] != ALGORITHM
            or previous["seed"] != seed
            or previous["run"] != run
            or previous["candidate_spots"] != total
            or previous["rounds"][-1]["end"] != start - 1
        ):
            raise ValueError("Sampling manifest does not match the next round")
    with ExitStack() as stack:
        temporary_paths = [stack.enter_context(atomic_output_path(path, suffix=".fastq.gz")) for path in output_paths]
        # Close the streams before atomic_output_path publishes either mate.
        with ExitStack() as writers:
            outputs = [writers.enter_context(gzip.open(path, "wb", compresslevel=1)) for path in temporary_paths]
            for number, spot in enumerate(spots, 1):
                if len(spot) != len(outputs):
                    raise ValueError("Sampling input layout differs from output layout")
                spot = [tuple(line.rstrip(b"\r\n") + b"\n" for line in record) for record in spot]
                candidate_count += 1
                lengths = [len(record[1].rstrip(b"\r\n")) for record in spot]
                bases = sum(lengths)
                candidate_bp += bases
                for record in spot:
                    for line in record:
                        input_digest.update(line)
                if number not in selected:
                    continue
                counts["num_dumped"] += 1
                counts["bp_dumped"] += bases
                selection_digest.update(str(number).encode() + b"\n")
                if any(length < min_length for length in lengths):
                    counts["num_rejected"] += 1
                    counts["bp_rejected"] += bases
                    continue
                for output, record in zip(outputs, spot, strict=True):
                    output.writelines(record)
                counts["num_written"] += 1
                counts["bp_written"] += bases
        if candidate_count != total:
            raise ValueError(
                "Sampling input count differs from metadata: {} != {}. "
                "Sampling requires one record/pair per declared spot.".format(candidate_count, total)
            )
        if previous is not None and previous["input_sha256"] != input_digest.hexdigest():
            raise ValueError("Sampling input changed between rounds")
    rounds = previous["rounds"] if previous else []
    rounds.append({"start": start, "end": end, "selected_sha256": selection_digest.hexdigest(), **counts})
    manifest = {
        "schema_version": 1,
        "algorithm": ALGORITHM,
        "seed": seed,
        "run": run,
        "unit": "pair" if len(output_paths) == 2 else "record",
        "candidate_spots": candidate_count,
        "candidate_bp": candidate_bp,
        "input_sha256": input_digest.hexdigest(),
        "selection": "permutation ranks, output in input order within each round",
        "rounds": rounds,
        "provenance": provenance or {},
    }
    with atomic_output_path(manifest_path) as tmp:
        with open(tmp, "w", encoding="utf-8") as handle:
            json.dump(manifest, handle, sort_keys=True, indent=2)
            handle.write("\n")
    return counts, manifest
