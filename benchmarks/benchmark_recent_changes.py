"""Compare revisions with isolated processes, wall time, peak RSS and output checks.

Archive the baseline (do not change branches), then run this same harness against
both trees. Input construction/imports and output comparison are outside the
timed interval. RSS is whole-process high-water memory, including input setup;
unlike tracemalloc it includes NumPy/native allocations. No network is used.
"""

from __future__ import annotations

import argparse
import contextlib
import gc
import gzip
import hashlib
import inspect
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

CASES = ('quant_validation', 'fastq_scan', 'wsfilter', 'csfilter', 'tau', 'sva', 'cstmm')
EXTRA_CASES = ('random_sampling', 'imputation')


def digest_file(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def prepare_case(args, root):
    import numpy as np
    import pandas as pd

    if args.case == 'quant_validation':
        from amalgkit.output_contracts import validate_quant_output_files
        rows = args.genes * 5
        pd.DataFrame(dict(target_id=[f'tx{i}' for i in range(rows)], length=1000,
                          eff_length=900, est_counts=np.arange(rows, dtype=float), tpm=1.)).to_csv(
            root / 'RUN_abundance.tsv', sep='\t', index=False)
        (root / 'RUN_run_info.json').write_text('{"p_pseudoaligned":75.0}', encoding='utf-8')

        def run():
            valid, error = validate_quant_output_files('RUN', str(root))
            if not valid:
                raise RuntimeError(error)
            return np.array([rows])

        return run, {'transcripts': rows}
    if args.case in {'fastq_scan', 'random_sampling'}:
        from amalgkit.fastq_utils import count_fastq_records_and_bases
        path = root / 'reads.fastq.gz'
        record = b'@read\n' + b'ACGT' * 25 + b'\n+\n' + b'I' * 100 + b'\n'
        with gzip.open(path, 'wb', compresslevel=1) as handle:
            for first in range(0, args.reads, 10000):
                handle.write(record * min(10000, args.reads - first))
        if args.case == 'fastq_scan':
            return lambda: np.array(count_fastq_records_and_bases(str(path))), {'reads': args.reads, 'read_length': 100}
        from amalgkit.getfastq_sampling import sample_fastqs

        def run():
            selected = max(1, int(args.reads * args.fraction))
            counts, _manifest = sample_fastqs([str(path)], [str(root / 'sample.fastq.gz')], total=args.reads,
                                             start=1, end=selected, seed=719, run='RUN', min_length=0, run_dir=str(root))
            if counts['num_written'] != selected or counts['bp_written'] != selected * 100:
                raise RuntimeError('Unexpected sampling counts')
            return np.array([counts['num_written'], counts['bp_written']])

        return run, {'reads': args.reads, 'read_length': 100, 'fraction': args.fraction, 'seed': 719}

    genes = min(args.genes, 5000) if args.case == 'sva' else args.genes
    samples = min(args.samples, 80) if args.case == 'sva' else args.samples
    if args.case in {'csfilter', 'imputation'}:
        genes = args.orthologs
    rng = np.random.default_rng(719)
    groups = [f'tissue{i}' for i in range(10)]
    runs = [f'R{i:04d}' for i in range(samples)]
    metadata = pd.DataFrame(dict(
        run=runs, sample_group=[groups[i % 10] for i in range(samples)],
        scientific_name=[f'Species {i // (samples // 10)}' for i in range(samples)],
        species_tag=[f'Species_{i // (samples // 10)}' for i in range(samples)],
        bioproject=[f'P{i % 4}' for i in range(samples)], exclusion='no',
    ))
    if args.case in {'wsfilter', 'tau', 'sva'}:
        metadata['scientific_name'] = 'Species 0'
    if args.case == 'imputation':
        values = np.repeat(rng.gamma(2, 100, size=(genes, 1)), samples, axis=1)
    elif args.case == 'tau':
        # Identical within-group replicates give equivalent tau across the old
        # geometric and new arithmetic definitions, isolating implementation cost.
        group_values = rng.gamma(2, 100, size=(genes, 10))
        values = group_values[:, np.arange(samples) % 10]
    else:
        values = rng.poisson(rng.gamma(2, 100, size=(genes, 1)), size=(genes, samples)).astype(float)
        values *= np.geomspace(.2, 5, samples)[None, :]
    counts = pd.DataFrame(values, index=[f'g{i}' for i in range(genes)], columns=runs)
    workload = dict(genes=genes, samples=samples, groups=10, seed=719)
    if args.case == 'imputation':
        from amalgkit.cstmm_python import _get_df_nonzero
        libraries = counts.sum()
        mask = (np.arange(genes)[:, None] // 20 + np.arange(samples)[None, :] // 10) % 5 == 0
        counts[mask] = np.nan

        def run():
            options = {'library_sizes': libraries, 'scale': 'library_size'} if 'scale' in inspect.signature(_get_df_nonzero).parameters else {}
            return _get_df_nonzero(counts, **options).to_numpy()

        return run, dict(workload, missing_fraction=.2, complete_profiles='identical; equal depth')
    if args.case in {'wsfilter', 'csfilter', 'tau', 'sva'}:
        counts = np.log2(counts + 1)
    if args.case == 'wsfilter':
        from amalgkit.per_species_python import _compute_sample_group_correlation_metrics

        def run():
            result = _compute_sample_group_correlation_metrics(counts, metadata, groups, 'pearson')
            return result[['ws_within_group_cor', 'ws_max_nongroup_cor', 'ws_margin']].to_numpy()

        return run, workload
    if args.case == 'csfilter':
        from amalgkit.cross_species_computation import calculate_correlation_within_group
        counts.columns = metadata['species_tag'] + '_' + metadata['run']
        counts.iloc[::11, ::7] = np.nan

        def run():
            result = calculate_correlation_within_group(metadata, counts, 'corrected')
            return result[['within_group_cor_corrected', 'max_nongroup_cor_corrected']].to_numpy()

        return run, dict(workload, species=10, missing='every 11th row / 7th column')
    if args.case == 'tau':
        from amalgkit import per_species_common as common

        def run():
            if hasattr(common, 'linear_sample_group_summary'):
                means = common.linear_sample_group_summary(counts, metadata, groups)['linear_mean']
                tau = common.sample_group_to_tau(means)
            else:
                means = common.sample_group_mean(counts, metadata, groups)['tc_ave']
                tau = common.sample_group_to_tau(means, transform_method='log2p1-fpkm')
            return tau[['tau']].to_numpy()

        return run, dict(workload, within_group_replicates='identical')
    if args.case == 'sva':
        from amalgkit.batch_effect_sva import run_sva_backend

        def run():
            result, _sv, summary = run_sva_backend(counts, metadata, nsv_setting='3', B_setting='5',
                                                   input_scale='transformed', random_seed=719)
            if not summary['corrected_run_ids']:
                raise RuntimeError('SVA benchmark did not correct the input: ' + str(summary))
            return result.to_numpy()

        return run, dict(workload, nsv=3, irw_iterations=5)
    from amalgkit.cstmm_python import _run_cstmm_python
    counts.columns = metadata['species_tag'] + '_' + metadata['run']
    uncorrected = {}
    for species in metadata['species_tag'].unique():
        uncorrected[species] = counts.loc[:, counts.columns.str.startswith(species + '_')].copy()
        (root / 'input' / species).mkdir(parents=True)
    sog = counts.iloc[:args.orthologs].copy()
    kwargs = {}
    if args.without_pairs and 'reference_diagnostics' in inspect.signature(_run_cstmm_python).parameters:
        kwargs['reference_diagnostics'] = False

    def run():
        result = _run_cstmm_python(uncorrected, sog, str(root / 'input'), str(root / 'output'), metadata,
                                   single_copy_threshold=50., **kwargs)
        return result.round2_factors.to_numpy()

    return run, dict(workload, orthologs=args.orthologs, species=10, missing='none',
                    observed_pair_diagnostics=not args.without_pairs)


def worker(args):
    sys.path.insert(0, str(args.repo.resolve()))
    import numpy as np
    import pandas as pd
    import scipy

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='amalgkit-perf-') as directory:
        root = Path(directory)
        function, workload = prepare_case(args, root)
        gc.collect()
        with (args.output.with_suffix('.log')).open('w') as log, contextlib.redirect_stdout(log):
            start = time.perf_counter()
            result = function()
            elapsed = time.perf_counter() - start
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak_mib = peak / (1024 ** 2 if sys.platform == 'darwin' else 1024)
        np.save(args.output.with_suffix('.npy'), result, allow_pickle=False)
        hashes = {path.name: digest_file(path) for path in sorted(root.glob('output/*/*_cstmm_counts.tsv'))}
    args.output.write_text(json.dumps(dict(
        seconds=elapsed, peak_rss_mib=peak_mib, workload=workload, corrected_count_sha256=hashes,
        environment=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
                         machine=platform.machine(), platform=platform.platform(),
                         threads={name: os.environ.get(name) for name in THREAD_VARS}),
    ), indent=2) + '\n', encoding='utf-8')


THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')


def compare(args):
    import numpy as np

    args.output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, **dict.fromkeys(THREAD_VARS, '1'))
    summary = {}
    for case in args.cases:
        measurements = {'baseline': [], 'current': []}
        equivalence = []
        for repeat in range(args.repeats + 1):
            results = {}
            for label in (('baseline', 'current') if repeat % 2 == 0 else ('current', 'baseline')):
                repo = args.baseline if label == 'baseline' else args.current
                path = args.output / case / f'{label}-{repeat}.json'
                command = [sys.executable, str(Path(__file__).resolve()), '--repo', str(repo), '--case', case,
                           '--output', str(path), '--genes', str(args.genes), '--samples', str(args.samples),
                           '--orthologs', str(args.orthologs), '--reads', str(args.reads), '--fraction', str(args.fraction)]
                if args.without_pairs:
                    command.append('--without-pairs')
                subprocess.run(command, env=env, check=True)  # noqa: S603 - explicit local benchmark revisions
                result = json.loads(path.read_text())
                results[label] = (path, result)
                if repeat:
                    measurements[label].append(result)
                print(f'{case} {label} {repeat}: {result["seconds"]:.3f}s, {result["peak_rss_mib"]:.1f} MiB', flush=True)
            left, right = (np.load(results[label][0].with_suffix('.npy'), allow_pickle=False)
                           for label in ('baseline', 'current'))
            same = bool(np.allclose(left, right, rtol=1e-8, atol=1e-10, equal_nan=True))
            finite = np.isfinite(left) & np.isfinite(right)
            difference = float(np.max(np.abs(left[finite] - right[finite]))) if finite.any() else 0.
            hashes_match = results['baseline'][1]['corrected_count_sha256'] == results['current'][1]['corrected_count_sha256']
            equivalence.append(dict(numeric_allclose=same, max_absolute_difference=difference, count_files_identical=hashes_match))
        summary[case] = {
            label: dict(median_seconds=statistics.median(row['seconds'] for row in rows),
                        min_seconds=min(row['seconds'] for row in rows), max_seconds=max(row['seconds'] for row in rows),
                        median_peak_rss_mib=statistics.median(row['peak_rss_mib'] for row in rows),
                        max_peak_rss_mib=max(row['peak_rss_mib'] for row in rows), workload=rows[0]['workload'])
            for label, rows in measurements.items()
        }
        summary[case]['equivalence'] = equivalence
        summary[case]['time_ratio'] = summary[case]['current']['median_seconds'] / summary[case]['baseline']['median_seconds']
        summary[case]['rss_ratio'] = summary[case]['current']['median_peak_rss_mib'] / summary[case]['baseline']['median_peak_rss_mib']
        (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--current', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--case', choices=CASES + EXTRA_CASES)
    parser.add_argument('--cases', nargs='+', choices=CASES + EXTRA_CASES, default=list(CASES))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--genes', type=int, default=20000)
    parser.add_argument('--samples', type=int, default=200)
    parser.add_argument('--orthologs', type=int, default=2000)
    parser.add_argument('--reads', type=int, default=1000000)
    parser.add_argument('--fraction', type=float, default=.1)
    parser.add_argument('--without-pairs', action='store_true')
    args = parser.parse_args()
    if (args.repeats < 1 or args.samples < 20 or args.samples % 10
            or min(args.genes, args.orthologs, args.reads) < 1 or not 0 < args.fraction <= 1):
        parser.error('Positive sizes/repeats and a sample count >=20 divisible by 10 are required.')
    if args.repo:
        if args.case is None:
            parser.error('--repo requires --case')
        worker(args)
    else:
        if args.baseline is None:
            parser.error('--baseline is required for a comparison')
        compare(args)


if __name__ == '__main__':
    main()
