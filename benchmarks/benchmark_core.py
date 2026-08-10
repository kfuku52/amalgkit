"""Reproducible microbenchmarks for performance-sensitive AMALGKIT kernels."""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

import numpy
import pandas

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from amalgkit.cross_species_computation import resolve_correlation_matrix  # noqa: E402
from amalgkit.getfastq_file_state import RunFileState  # noqa: E402
from amalgkit.output_contracts import validate_quant_output_files  # noqa: E402


def _measure(function, repeats: int = 3) -> dict[str, float]:
    durations = []
    peaks = []
    for _ in range(repeats):
        tracemalloc.start()
        started_at = time.perf_counter()
        function()
        durations.append(time.perf_counter() - started_at)
        _current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        peaks.append(peak / (1024 * 1024))
    return {
        "median_seconds": round(statistics.median(durations), 6),
        "max_peak_mib": round(max(peaks), 3),
    }


def _benchmark_quant_validation(root: Path, row_count: int) -> dict[str, float]:
    run_id = "BENCH001"
    abundance_path = root / f"{run_id}_abundance.tsv"
    frame = pandas.DataFrame(
        {
            "target_id": [f"tx{index}" for index in range(row_count)],
            "length": numpy.full(row_count, 1000, dtype=int),
            "eff_length": numpy.full(row_count, 900, dtype=int),
            "est_counts": numpy.arange(row_count, dtype=float),
            "tpm": numpy.full(row_count, 1.0, dtype=float),
        }
    )
    frame.to_csv(abundance_path, sep="\t", index=False)
    (root / f"{run_id}_run_info.json").write_text(
        json.dumps({"p_pseudoaligned": 75.0}),
        encoding="utf-8",
    )

    def validate() -> None:
        valid, error = validate_quant_output_files(run_id, str(root))
        if not valid:
            raise RuntimeError(error)

    return _measure(validate)


def _benchmark_correlation(row_count: int, sample_count: int) -> dict[str, float]:
    rng = numpy.random.default_rng(42)
    matrix = pandas.DataFrame(
        rng.normal(size=(row_count, sample_count)),
        columns=[f"run{index}" for index in range(sample_count)],
    )
    matrix.iloc[::20, ::7] = numpy.nan
    return _measure(lambda: resolve_correlation_matrix(matrix, missing_strategy="row_mean"))


def _benchmark_file_inventory(root: Path, file_count: int) -> dict[str, float]:
    run_dir = root / "inventory"
    run_dir.mkdir()
    for index in range(file_count):
        (run_dir / f"chunk-{index}.fastq").touch()
    return _measure(lambda: RunFileState(str(run_dir)), repeats=5)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    quant_rows = 10_000 if args.quick else 100_000
    correlation_rows = 500 if args.quick else 2_000
    samples = 30 if args.quick else 100
    file_count = 200 if args.quick else 2_000

    with tempfile.TemporaryDirectory(prefix="amalgkit-benchmark-") as temp_dir:
        root = Path(temp_dir)
        results = {
            "parameters": {
                "quant_rows": quant_rows,
                "correlation_rows": correlation_rows,
                "correlation_samples": samples,
                "file_count": file_count,
            },
            "quant_output_validation": _benchmark_quant_validation(root, quant_rows),
            "cross_species_correlation": _benchmark_correlation(correlation_rows, samples),
            "getfastq_file_inventory": _benchmark_file_inventory(root, file_count),
        }
    rendered = json.dumps(results, indent=2, sort_keys=True)
    print(rendered)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
