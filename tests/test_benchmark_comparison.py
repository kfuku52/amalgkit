"""Revision comparisons must reject incorrect outputs, including warmup runs."""

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy
import pytest


@pytest.mark.parametrize('mismatch', ['values', 'counts', 'warmup'])
def test_benchmark_comparison_records_and_rejects_incorrect_output(tmp_path, monkeypatch, mismatch):
    path = Path(__file__).resolve().parents[1] / 'benchmarks' / 'benchmark_recent_changes.py'
    spec = importlib.util.spec_from_file_location('benchmark_comparison_under_test', path)
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)
    baseline = tmp_path / 'baseline'
    current = tmp_path / 'current'

    def run_worker(command, **kwargs):
        output = Path(command[command.index('--output') + 1])
        repo = Path(command[command.index('--repo') + 1])
        output.parent.mkdir(parents=True, exist_ok=True)
        wrong = repo == current and (mismatch == 'values' or (mismatch == 'warmup' and output.stem.endswith('-0')))
        numpy.save(output.with_suffix('.npy'), numpy.array([2. if wrong else 1.]))
        result = dict(seconds=1., cpu_seconds=1., peak_rss_mib=1., workload={},
                      corrected_count_sha256={'counts.tsv': 'wrong' if repo == current and mismatch == 'counts' else 'same'})
        output.write_text(json.dumps(result))

    monkeypatch.setattr(benchmark.subprocess, 'run', run_worker)
    args = SimpleNamespace(output=tmp_path / 'results', baseline=baseline, current=current,
                           cases=['cstmm'], repeats=1, genes=20, samples=20, orthologs=20, reads=20,
                           fraction=.1, without_pairs=False, allocation_scan=False)
    with pytest.raises(RuntimeError, match='output equivalence failed'):
        benchmark.compare(args)
    summary = json.loads((args.output / 'summary.json').read_text())
    assert len(summary['cstmm']['equivalence']) == 2
