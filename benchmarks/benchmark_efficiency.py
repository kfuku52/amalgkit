"""Compare the second efficiency audit against an archived source tree.

Each revision runs in a separate process with one warmup and three timed runs.
Inputs/imports and result serialization are outside the timer. The allocation
scan is separate; process RSS includes setup, repeats and that scan.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import inspect
import json
import os
from pathlib import Path
import resource
import statistics
import subprocess
import sys
import tempfile
import time
import tracemalloc

CASES = ('summary', 'tau', 'overlay', 'correlation', 'sampling', 'sva', 'mds', 'diagnostics', 'sva_pipeline')


def prepare(case, root):
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(719)
    if case == 'summary':
        from amalgkit.merge_plots import _read_est_count_summary
        path = root / 'Species_A' / 'Species_A_est_counts.tsv'
        path.parent.mkdir()
        with path.open('w') as handle:
            handle.write('target_id\t' + '\t'.join(f'R{i}' for i in range(64)) + '\n')
            for start in range(0, 120000, 5000):
                frame = pd.DataFrame((np.arange(start, start+5000)[:, None] + np.arange(64)[None, :]) % 1001,
                                     columns=[f'R{i}' for i in range(64)])
                frame.insert(0, 'target_id', [f'g{i}' for i in range(start, start+5000)])
                frame.to_csv(handle, sep='\t', index=False, header=False)
        return lambda: _read_est_count_summary(path)[['mean_before', 'lib_size']].to_numpy(), {'genes': 120000, 'samples': 64}
    if case == 'overlay':
        from amalgkit.gsa_snapshot import _overlay_measurements
        rows = []
        for i in range(1500):
            run = f'CRR{i:08d}'
            name = run + '.fastq.gz'
            manifest = [{'filename': name, 'mate': 0, 'group': 0,
                         'sources': [{'source_name': 'GSA', 'url': f'https://download.cncb.ac.cn/gsa2/CRA1/{run}/{name}'}]}]
            rows.append(dict(run=run, data_source='gsa', private_file='no', lib_layout='single',
                             gsa_fastq_files=json.dumps(manifest), total_spots='', total_bases='',
                             spot_length='', read_count_status='unknown'))
        source = pd.DataFrame(rows).set_index('run', drop=False)
        measured = pd.DataFrame([dict(row, total_spots='8', total_bases='160', spot_length='20',
                                      read_count_status='measured', gsa_input_seconds='1',
                                      gsa_input_fingerprint='a'*64) for row in rows])
        def run():
            out = source.copy(deep=True)
            _overlay_measurements(out, measured)
            return out
        return run, {'rows': 1500, 'measured_fraction': 1}
    if case == 'sampling':
        from amalgkit import getfastq_sampling as module
        function = getattr(module, '_sorted_selected_indices', module.selected_indices)
        return lambda: function(2000000, 1, 300000, 719, 'AUDIT'), {'spots': 2000000, 'chosen': 300000}
    if case == 'sva':
        from amalgkit.batch_effect_sva import _supported_right_vectors
        values = rng.normal(size=(20000, 120))
        return lambda: _supported_right_vectors(values, 3), {'genes': 20000, 'samples': 120, 'rank': 'full'}
    if case == 'sva_pipeline':
        from amalgkit.batch_effect_sva import run_sva_backend
        genes, samples = 1500, 60
        groups = np.repeat(['A', 'B'], samples // 2)
        latent = rng.normal(size=samples)
        loading = rng.normal(size=genes)
        loading[:genes//3] = 0
        values = (20 + rng.normal(size=(genes, samples)) * .3
                  + loading[:, None] * latent[None, :] * 2
                  + rng.normal(size=(genes, 1)) * (groups == 'B')[None, :])
        runs = [f'R{i}' for i in range(samples)]
        frame = pd.DataFrame(values, columns=runs)
        metadata = pd.DataFrame({'run': runs, 'sample_group': groups, 'bioproject': [f'P{i%2}' for i in range(samples)]})
        def run():
            corrected, _sv, summary = run_sva_backend(frame, metadata, nsv_setting='1', B_setting='5',
                                                      input_scale='transformed', random_seed=719)
            if not summary['corrected_run_ids']:
                raise RuntimeError('Synthetic latent-factor correction failed: ' + str(summary['skip_reason']))
            return corrected.to_numpy()
        return run, {'genes': genes, 'samples': samples, 'planted_factors': 1, 'irw_iterations': 5}
    if case == 'mds':
        from amalgkit.per_species_finalize_python import _classic_mds_coordinates
        from scipy.spatial.distance import pdist, squareform
        dist = squareform(pdist(rng.normal(size=(800, 8))))
        return lambda: _classic_mds_coordinates(dist), {'samples': 800, 'dimensions': 8}
    if case == 'diagnostics':
        from amalgkit.cstmm_diagnostics import write_observed_pair_diagnostics
        frame = pd.DataFrame(rng.poisson(200, size=(1000, 80)).astype(float), columns=[f'R{i}' for i in range(80)])
        frame.iloc[::13, ::11] = np.nan
        libraries, factors = frame.sum(), pd.Series(np.linspace(.5, 1.5, 80), index=frame.columns)
        path = root / 'pairs.tsv'
        def run():
            write_observed_pair_diagnostics(frame, libraries, factors, path, global_reference='R17')
            return path
        return run, {'genes': 1000, 'samples': 80, 'pairs': 3160}
    samples = 400 if case == 'tau' else 200
    genes = 5000 if case == 'tau' else 10000
    runs = [f'R{i}' for i in range(samples)]
    frame = pd.DataFrame(rng.gamma(2, 100, size=(genes, samples)), columns=runs)
    metadata = pd.DataFrame({'run': runs, 'sample_group': [f'T{i%10}' for i in range(samples)],
                             'bioproject': [f'P{i%7}' for i in range(samples)], 'exclusion': 'no'})
    if case == 'tau':
        from amalgkit.per_species_common import linear_sample_group_summary
        frame.iloc[::101, ::13] = np.nan
        return lambda: linear_sample_group_summary(frame, metadata, transform_method='none-none'), {'genes': genes, 'samples': samples}
    from amalgkit.per_species_outputs import _compute_corr_matrix, save_correlation_statistics, save_state_overview_pdf
    can_reuse = 'precomputed_tc_dist_matrix' in inspect.signature(save_state_overview_pdf).parameters
    def run():
        corr = _compute_corr_matrix(frame, 'pearson') if can_reuse else None
        stats = save_correlation_statistics(frame, metadata, 'pearson', 0, precomputed_tc_dist_matrix=corr)
        if corr is None:
            corr = _compute_corr_matrix(frame, 'pearson')
        return {'correlation': corr, 'statistics': stats}
    return run, {'genes': genes, 'samples': samples, 'scope': 'statistics and plot correlation preparation'}


def worker(args):
    sys.path.insert(0, str(args.repo.resolve()))
    import numpy as np
    import pandas as pd
    with tempfile.TemporaryDirectory(prefix='amalgkit-efficiency-') as directory:
        function, workload = prepare(args.case, Path(directory))
        samples = []
        for i in range(4):
            gc.collect()
            cpu, wall = time.process_time(), time.perf_counter()
            result = function()
            elapsed, used = time.perf_counter() - wall, time.process_time() - cpu
            if i:
                samples.append({'wall_seconds': elapsed, 'cpu_seconds': used})
            del result
        gc.collect()
        tracemalloc.start()
        result = function()
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        if args.case == 'sampling':
            result = np.asarray(sorted(result), dtype=np.uint64)
        elif args.case in {'sva', 'mds'}:
            result = result @ result.T  # Compare subspaces/geometry, independent of eigenvector signs.
        if isinstance(result, Path):
            digest = hashlib.sha256(result.read_bytes()).hexdigest()
            result = np.frombuffer(bytes.fromhex(digest), dtype=np.uint8)
        elif isinstance(result, (dict, pd.DataFrame)):
            frames = result if isinstance(result, dict) else {'frame': result}
            digest = hashlib.sha256()
            for key, frame in frames.items():
                digest.update(key.encode())
                digest.update(frame.to_csv(sep='\t', index=True).encode())
            result = np.frombuffer(digest.digest(), dtype=np.uint8)
        np.save(args.output.with_suffix('.npy'), result, allow_pickle=False)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (2**20 if sys.platform == 'darwin' else 1024)
    args.output.write_text(json.dumps(dict(
        workload=workload, measurements=samples,
        median_cpu_seconds=statistics.median(row['cpu_seconds'] for row in samples),
        median_wall_seconds=statistics.median(row['wall_seconds'] for row in samples),
        peak_python_mib=peak/2**20, process_peak_rss_mib=rss,
        environment=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__,
                         load_average=os.getloadavg(), native_threads=1),
    ), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--current', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--case', choices=CASES)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.repo:
        worker(args)
        return
    if args.baseline is None:
        parser.error('--baseline is required')
    import numpy as np
    args.output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, **dict.fromkeys(('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'), '1'))
    summary = {}
    for index, case in enumerate(args.cases):
        rows = {}
        for label in (('baseline', 'current') if index % 2 == 0 else ('current', 'baseline')):
            path = args.output / f'{case}-{label}.json'
            command = [sys.executable, str(Path(__file__).resolve()), '--repo',
                       str(args.baseline if label == 'baseline' else args.current),
                       '--case', case, '--output', str(path)]
            subprocess.run(command, env=env, check=True)  # noqa: S603 - explicitly selected local source revisions
            rows[label] = json.loads(path.read_text())
        left, right = (np.load((args.output / f'{case}-{label}.npy'), allow_pickle=False) for label in ('baseline', 'current'))
        exact = bool(np.array_equal(left, right, equal_nan=True))
        equivalent = bool(np.allclose(left, right, rtol=1e-8, atol=1e-10, equal_nan=True))
        summary[case] = dict(results=rows, exact=exact, equivalent=equivalent,
                             cpu_ratio=rows['current']['median_cpu_seconds']/rows['baseline']['median_cpu_seconds'])
        (args.output / 'summary.json').write_text(json.dumps(summary, indent=2))
        print(case, 'CPU ratio', round(summary[case]['cpu_ratio'], 3), 'exact', exact, 'equivalent', equivalent, flush=True)
        if not equivalent:
            raise RuntimeError(f'{case} output equivalence failed')


if __name__ == '__main__':
    main()
