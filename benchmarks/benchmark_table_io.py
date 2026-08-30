"""Isolated-process table I/O benchmark; compare JSON hashes before/after edits."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from amalgkit.merge import load_quant_tables_once  # noqa: E402
from amalgkit.output_contracts import validate_quant_output_files  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=int, default=100_000)
    parser.add_argument("--runs", type=int, default=8)
    parser.add_argument("--operation", choices=("validate", "merge"), default="validate")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="amalgkit-table-benchmark-") as scratch:
        root = Path(scratch)
        runs = [f"R{i}" for i in range(args.runs if args.operation == "merge" else 1)]
        paths = [str(root / f"{run}_abundance.tsv") for run in runs]
        for run, path in zip(runs, paths):
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("target_id\tlength\teff_length\test_counts\ttpm\n")
                for i in range(args.rows):
                    handle.write(f"gene_{i:020d}\t1000\t900\t{i % 100}\t1\n")
            (root / f"{run}_run_info.json").write_text('{"p_pseudoaligned":75}', encoding="utf-8")
        tracemalloc.start()
        started = time.perf_counter()
        if args.operation == "validate":
            valid, message = validate_quant_output_files(runs[0], str(root))
            if not valid:
                raise RuntimeError(message)
            result_hash = hashlib.sha256(b"valid").hexdigest()
        else:
            ids, values = load_quant_tables_once(runs, paths, ["eff_length", "est_counts", "tpm"])
            digest = hashlib.sha256("\n".join(ids).encode())
            for column in values.values():
                for array in column:
                    digest.update(array.astype("<f8").tobytes())
            result_hash = digest.hexdigest()
        seconds = time.perf_counter() - started
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform != "darwin":
            rss *= 1024
        result = {
            "operation": args.operation,
            "rows": args.rows,
            "runs": len(runs),
            "seconds_with_tracemalloc": round(seconds, 4),
            "peak_python_mib": round(peak / 2**20, 3),
            "peak_process_rss_mib": round(rss / 2**20, 3),
            "result_sha256": result_hash,
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
