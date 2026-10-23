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
import tracemalloc
from pathlib import Path

CASES = ('quant_validation', 'fastq_scan', 'wsfilter', 'csfilter', 'tau', 'sva', 'cstmm')
EXTRA_CASES = ('random_sampling', 'imputation', 'select_normalization', 'gsa_selection',
               'tmm_rounds', 'averaged_inputs', 'select_exclusion', 'select_exclusion_unique', 'provenance_hash',
               'outlier_groups', 'orthogroup_counts', 'orthogroup_counts_long',
               'mapping_rate', 'metadata_merge', 'tau_annotations', 'correlation_stats', 'tpm_transform')


def digest_file(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def prepare_table_case(args, root):
    from types import SimpleNamespace

    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(719)
    if args.case in {'mapping_rate', 'metadata_merge'}:
        runs = [f'R{i:07d}' for i in range(args.genes)]
    if args.case == 'mapping_rate':
        from amalgkit.metadata_utils import get_mapping_rate
        for index, run in enumerate(runs):
            directory = root / run
            directory.mkdir()
            (directory / f'{run}_run_info.json').write_text(
                json.dumps({'p_pseudoaligned': float(index % 100)}), encoding='utf-8')
        frame = pd.DataFrame({'run': runs + [' ' + run + ' ' for run in runs[::10]] + ['0001', 'NA', ''],
                              'mapping_rate': 88.})
        frame.index = [index // 2 for index in range(len(frame))]

        def run():
            return get_mapping_rate(SimpleNamespace(df=frame.copy()), str(root), max_workers=1).df

        return run, dict(runs=args.genes, rows=len(frame), workers=1, duplicate_normalized_runs=True)
    if args.case == 'metadata_merge':
        from amalgkit.filter_utils import merge_metadata_by_run
        frame = pd.DataFrame({'run': runs, **{f'score{i}': rng.normal(size=args.genes) for i in range(12)}})
        for column in range(4):
            frame[f'annotation{column}'] = pd.Series(['original'] * args.genes, dtype='string')
        update = frame.iloc[rng.permutation(args.genes)[:max(1, args.genes * 3 // 4)]].copy()
        for column in range(12):
            update[f'score{column}'] *= 2
            update.iloc[::13, update.columns.get_loc(f'score{column}')] = np.nan
        for column in range(4):
            update[f'annotation{column}'] = pd.Series('updated', index=update.index, dtype='string')
        return (lambda: merge_metadata_by_run(frame, update, overwrite_columns=('score0',)),
                dict(rows=args.genes, updates=len(update), numeric_columns=12, string_columns=4, seed=719))
    if args.case == 'tau_annotations':
        from amalgkit.per_species_common import sample_group_to_tau
        values = rng.gamma(2, 100, size=(args.genes, args.samples))
        values[::19] = 0
        values[1::23, :2] = 500
        values[2::29, 0] = np.nan
        groups = [f'tissue_{i:02d}' for i in range(args.samples)]
        groups[:3] = ['0001', 'NA', r'tissue|escaped\name']
        frame = pd.DataFrame(values, columns=groups, index=[f'gene{i}' for i in range(args.genes)])
        return (lambda: sample_group_to_tau(frame, rich_annotation=True),
                dict(genes=args.genes, groups=args.samples, rich_annotation=True, seed=719))
    if args.case == 'tpm_transform':
        from amalgkit.per_species_finalize_python import _transform_raw_to_tpm
        columns = [f'R{i:07d}' for i in range(args.samples)]
        index = pd.Index(['0001', 'NA'][:args.genes] + [f'gene{i}' for i in range(2, args.genes)])
        counts = pd.DataFrame(rng.gamma(2, 100, size=(args.genes, args.samples)), index=index, columns=columns)
        lengths = pd.DataFrame(rng.uniform(100, 3000, size=counts.shape), index=index, columns=columns)
        return (lambda: _transform_raw_to_tpm(counts, lengths),
                dict(genes=args.genes, samples=args.samples, transform='tpm', seed=719))
    from amalgkit.per_species_outputs import save_correlation_statistics
    runs = [f'R{i:07d}' for i in range(args.samples)]
    frame = pd.DataFrame(rng.normal(size=(16, args.samples)), columns=runs)
    metadata = pd.DataFrame({'run': runs, 'bioproject': [f'project_{i % 130}' for i in range(args.samples)],
                             'sample_group': [f'tissue_{i % 19}' for i in range(args.samples)]})
    correlation = rng.uniform(-1, 1, size=(args.samples, args.samples))
    correlation[::137, ::151] = np.nan
    return (lambda: save_correlation_statistics(frame, metadata, 'pearson', 0,
                                                precomputed_tc_dist_matrix=correlation),
            dict(samples=args.samples, matrix='precomputed', projects=130, groups=19, seed=719))


def prepare_adoption_case(args, root):
    import numpy as np
    import pandas as pd
    from types import SimpleNamespace

    rng = np.random.default_rng(719)
    if args.case.startswith('orthogroup_counts'):
        from amalgkit.orthology_utils import orthogroup2genecount
        species = [f'Species_{i:02d}' for i in range(min(args.samples, 50))]
        data = {'busco_id': ['0001', 'NA'][:args.genes] + [f'OG{i:07d}' for i in range(2, args.genes)]}
        for column, name in enumerate(species):
            data[name] = [(','.join(f'g{column}_{row}_{copy}' for copy in range(1 + row % 3))
                           if row % 11 else '-') for row in range(args.genes)]
        if args.case.endswith('_long'):
            data[species[0]][0] = 'long_gene_' + 'x' * 246
        path, output = root / 'orthogroups.tsv', root / 'gene_counts.tsv'
        pd.DataFrame(data).to_csv(path, sep='\t', index=False)

        def run():
            orthogroup2genecount(str(path), str(output), list(reversed(species)))
            return np.frombuffer(hashlib.sha256(output.read_bytes()).digest(), dtype=np.uint8)

        return run, dict(rows=args.genes, species=len(species), long_id=args.case.endswith('_long'))
    if args.case == 'outlier_groups':
        from amalgkit.outlier_utils import compute_groupwise_robust_z
        values = rng.normal(size=args.genes)
        groups = np.array([f'group{i % args.samples}' for i in range(args.genes)], dtype=object)
        order = rng.permutation(args.genes)
        values, groups = values[order], groups[order]
        values[::137] = np.nan
        values[::251] = np.inf
        groups[::197] = ''
        return (lambda: compute_groupwise_robust_z(values, groups),
                dict(rows=args.genes, groups=args.samples, seed=719, min_group_size=3))
    if args.case == 'provenance_hash':
        from amalgkit.quant_provenance import file_identity
        path = root / 'input.bin'
        block = np.arange(256 * 1024, dtype=np.uint8).tobytes()
        with path.open('wb') as handle:
            for _ in range(1024):
                handle.write(block)

        def run():
            identity = file_identity(str(path), {})
            return np.array([str(identity['size']), identity['sha256']])

        return run, dict(bytes=256 * 2**20, identity_cache='miss', disk_cache='warm after warmup')
    if args.case.startswith('select_exclusion'):
        from amalgkit.select import apply_select_exclude_rules, read_select_rules
        rules = read_select_rules(args.repo / 'amalgkit/select_rule_sets/base/select_rules.tsv')
        rules = [rule for rule in rules if rule['stage'] == 'exclude']
        columns = list(dict.fromkeys(column for rule in rules for column in rule['columns']))
        vocabulary = ['bulk leaf healthy', 'single cell RNA', 'miRNA root', 'CONTROL', '', 'NA', '0001']
        unique = args.case.endswith('_unique')
        data = {}
        for offset, column in enumerate(columns):
            data[column] = [vocabulary[(i + offset) % len(vocabulary)] + (f' accession_{i}' if unique else '')
                            for i in range(args.genes)]
        data.update(run=[f'R{i}' for i in range(args.genes)], exclusion='no')
        frame = pd.DataFrame(data)
        return (lambda: apply_select_exclude_rules(SimpleNamespace(df=frame.copy()), rules).df,
                dict(rows=args.genes, columns=len(columns), rules=len(rules), cardinality='unique' if unique else 'low'))
    values = rng.gamma(2, 100, size=(args.genes, args.samples))
    if args.case == 'tmm_rounds':
        from amalgkit.normalization_tmm import run_tmm_rounds_for_cstmm
        values *= np.geomspace(.2, 5, args.samples)[None, :]
        counts = pd.DataFrame(values, columns=[f'R{i}' for i in range(args.samples)])

        def run():
            result = run_tmm_rounds_for_cstmm(counts)
            return np.concatenate([result.round1_factors.to_numpy(), result.round2_factors.to_numpy(),
                                   result.library_sizes.to_numpy(), [result.round1_reference_column],
                                   result.median_reference_columns])

        return run, dict(genes=args.genes, samples=args.samples, rounds=2, counts='fractional')
    from amalgkit.cross_species_filter import _build_averaged_cross_species_inputs
    # Five groups and three replicates per species; omit the final incomplete
    # block rather than silently changing replicate weights between revisions.
    sample_count = args.samples // 15 * 15
    rows = [dict(species_tag=f'S{i//15}', scientific_name=f'Species {i//15}', run=f'R{i}',
                 sample_group=f'T{i%5}', bioproject=f'P{i%3}', exclusion='no') for i in range(sample_count)]
    metadata = pd.DataFrame(rows)
    values = np.log2(values[:, :sample_count] + 1)
    values[::31, ::7] = np.nan
    matrix = pd.DataFrame(values, columns=[f'S{i//15}_R{i}' for i in range(sample_count)])
    orthologs = dict(uncorrected=matrix, corrected=matrix * 1.1)
    return (lambda: _build_averaged_cross_species_inputs(metadata, orthologs),
            dict(genes=args.genes, samples=sample_count, species=sample_count//15, groups=5, replicates=3))


def prepare_case(args, root):
    import numpy as np
    import pandas as pd

    if args.case in {'mapping_rate', 'metadata_merge', 'tau_annotations', 'correlation_stats', 'tpm_transform'}:
        return prepare_table_case(args, root)
    if args.case in {'tmm_rounds', 'averaged_inputs', 'select_exclusion', 'select_exclusion_unique', 'provenance_hash',
                     'outlier_groups', 'orthogroup_counts', 'orthogroup_counts_long'}:
        return prepare_adoption_case(args, root)
    if args.case == 'select_normalization':
        from amalgkit.select import normalize_select_metadata_frame, read_select_rules
        rules = read_select_rules(args.repo / 'amalgkit/select_rule_sets/plantae/select_rules.tsv')
        tissues = ['leaf', 'root', 'flower', 'leaf and root', 'leaf; treated',
                   'mature leaf', 'root tip', 'young flower', 'unknown', '0001', 'NA']
        frame = pd.DataFrame({
            'sample_group': '',
            'sample_attribute_tissue': [tissues[i % len(tissues)] for i in range(args.genes)],
            'source_name': ['leaf covering inflorescence' if i % 5 == 0 else '' for i in range(args.genes)],
            'sample_title': [f'treatment {i % 20}' for i in range(args.genes)],
        })

        def run():
            result = normalize_select_metadata_frame(frame, rules)
            return result.fillna('').to_numpy(dtype=str)

        return run, dict(rows=args.genes, distinct_tissues=len(tissues), rule_set='plantae')
    if args.case == 'gsa_selection':
        from types import SimpleNamespace
        from amalgkit.gsa_select import DEFERRED_COLUMN, resolve_deferred_selection
        recipe = json.dumps([dict(kind='min_spots', threshold=20, target='exclusion', outcome='low_spots'),
                             dict(kind='dedup', columns=['biosample'], outcome='redundant')], sort_keys=True)
        frame = pd.DataFrame({
            'run': [f'CRR{i:08d}' for i in range(args.genes)],
            'data_source': ['gsa' if i % 2 == 0 else 'sra' for i in range(args.genes)],
            'biosample': [f'B{i // 20}' for i in range(args.genes)],
            'total_spots': [np.nan if i % 40 == 0 else i % 101 for i in range(args.genes)],
            'read_count_status': ['unknown' if i % 40 == 0 else 'measured' for i in range(args.genes)],
            'exclusion': 'no', 'is_sampled': 'yes',
            DEFERRED_COLUMN: [recipe if i % 10 == 0 else '' for i in range(args.genes)],
            'gsa_selection_status': ['pending_input_counts' if i % 10 == 0 else '' for i in range(args.genes)],
            **{f'annotation_{i}': '' for i in range(40)},
        })

        def run():
            result = resolve_deferred_selection(SimpleNamespace(df=frame.copy())).df
            return result[['exclusion', DEFERRED_COLUMN, 'gsa_selection_status']].to_numpy(dtype=str)

        return run, dict(rows=args.genes, deferred_fraction=.1, missing_count_fraction=.025)
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
            return result[['ws_within_group_cor', 'ws_max_nongroup_cor', 'ws_margin',
                           'ws_within_common_genes', 'ws_min_nongroup_common_genes']].to_numpy()

        return run, workload
    if args.case == 'csfilter':
        from amalgkit.cross_species_computation import calculate_correlation_within_group
        counts.columns = metadata['species_tag'] + '_' + metadata['run']
        counts.iloc[::11, ::7] = np.nan

        def run():
            result = calculate_correlation_within_group(metadata, counts, 'corrected')
            return result[['within_group_cor_corrected', 'max_nongroup_cor_corrected',
                           'within_common_genes_corrected', 'min_nongroup_common_genes_corrected']].to_numpy()

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
        load_average = os.getloadavg() if hasattr(os, 'getloadavg') else None
        with (args.output.with_suffix('.log')).open('w') as log, contextlib.redirect_stdout(log):
            cpu_start = time.process_time()
            start = time.perf_counter()
            result = function()
            elapsed = time.perf_counter() - start
            cpu_elapsed = time.process_time() - cpu_start
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak_mib = peak / (1024 ** 2 if sys.platform == 'darwin' else 1024)
        if isinstance(result, (dict, pd.DataFrame)):
            frames = result if isinstance(result, dict) else {'metadata': result}
            digest = hashlib.sha256()
            for key, frame in frames.items():
                digest.update(key.encode())
                digest.update(frame.to_csv(sep='\t', index=True).encode())
            result = np.frombuffer(digest.digest(), dtype=np.uint8)
        np.save(args.output.with_suffix('.npy'), result, allow_pickle=False)
        hashes = {path.name: digest_file(path) for path in sorted(root.glob('output/*/*_cstmm_counts.tsv'))}
        # Preserve the timed call's output before a separate allocation scan can
        # mutate shared return values or replace files in the output directory.
        del result
        allocation_peak = None
        if args.allocation_scan:
            gc.collect()
            tracemalloc.start()
            with args.output.with_suffix('.allocation.log').open('w') as log, contextlib.redirect_stdout(log):
                function()
            allocation_peak = tracemalloc.get_traced_memory()[1] / 2**20
            tracemalloc.stop()
    args.output.write_text(json.dumps(dict(
        seconds=elapsed, cpu_seconds=cpu_elapsed, peak_rss_mib=peak_mib, peak_python_mib=allocation_peak,
        workload=workload, corrected_count_sha256=hashes,
        environment=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
                         machine=platform.machine(), platform=platform.platform(),
                         threads={name: os.environ.get(name) for name in THREAD_VARS}, load_average=load_average),
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
                if args.allocation_scan:
                    command.append('--allocation-scan')
                subprocess.run(command, env=env, check=True)  # noqa: S603 - explicit local benchmark revisions
                result = json.loads(path.read_text())
                results[label] = (path, result)
                if repeat:
                    measurements[label].append(result)
                print(f'{case} {label} {repeat}: {result["seconds"]:.3f}s wall, '
                      f'{result["cpu_seconds"]:.3f}s CPU, {result["peak_rss_mib"]:.1f} MiB', flush=True)
            left, right = (np.load(results[label][0].with_suffix('.npy'), allow_pickle=False)
                           for label in ('baseline', 'current'))
            numeric = left.dtype.kind in 'biufc' and right.dtype.kind in 'biufc'
            if left.shape != right.shape:
                same, difference, identical = False, None, False
            elif numeric:
                same = bool(np.allclose(left, right, rtol=1e-8, atol=1e-10, equal_nan=True))
                finite = np.isfinite(left) & np.isfinite(right)
                difference = float(np.max(np.abs(left[finite] - right[finite]))) if finite.any() else 0.
                identical = bool(np.array_equal(left, right, equal_nan=True))
            else:
                same, difference = None, None
                identical = bool(np.array_equal(left, right))
            hashes_match = results['baseline'][1]['corrected_count_sha256'] == results['current'][1]['corrected_count_sha256']
            equivalence.append(dict(outputs_identical=identical, numeric_allclose=same,
                                    max_absolute_difference=difference, count_files_identical=hashes_match))
        summary[case] = {
            label: dict(median_seconds=statistics.median(row['seconds'] for row in rows),
                        min_seconds=min(row['seconds'] for row in rows), max_seconds=max(row['seconds'] for row in rows),
                        median_cpu_seconds=statistics.median(row['cpu_seconds'] for row in rows),
                        min_cpu_seconds=min(row['cpu_seconds'] for row in rows),
                        max_cpu_seconds=max(row['cpu_seconds'] for row in rows),
                        median_peak_rss_mib=statistics.median(row['peak_rss_mib'] for row in rows),
                        max_peak_rss_mib=max(row['peak_rss_mib'] for row in rows),
                        median_peak_python_mib=(statistics.median(row['peak_python_mib'] for row in rows)
                                                if args.allocation_scan else None), workload=rows[0]['workload'])
            for label, rows in measurements.items()
        }
        summary[case]['equivalence'] = equivalence
        summary[case]['time_ratio'] = summary[case]['current']['median_seconds'] / summary[case]['baseline']['median_seconds']
        summary[case]['cpu_ratio'] = summary[case]['current']['median_cpu_seconds'] / summary[case]['baseline']['median_cpu_seconds']
        summary[case]['rss_ratio'] = summary[case]['current']['median_peak_rss_mib'] / summary[case]['baseline']['median_peak_rss_mib']
        (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
        if not all((row['outputs_identical'] or row['numeric_allclose']) and row['count_files_identical']
                   for row in equivalence):
            raise RuntimeError(f'{case} output equivalence failed; see {args.output / "summary.json"}')


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
    parser.add_argument('--allocation-scan', action='store_true',
                        help='Measure Python allocations separately, after timing and process peak RSS.')
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
