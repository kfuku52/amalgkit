import pandas
import numpy

import os
import warnings
from amalgkit.util import *

FASTP_STATS_COLUMNS = ['fastp_duplication_rate', 'fastp_insert_size_peak']
MERGE_QUANT_READ_MAX_WORKERS = 4


def validate_metadata_columns(metadata, required_columns, context):
    missing = [col for col in required_columns if col not in metadata.df.columns]
    if len(missing) > 0:
        raise ValueError(
            'Missing required metadata column(s) for {}: {}'.format(
                context,
                ', '.join(missing),
            )
        )


def collect_valid_run_ids(run_values):
    run_ids = []
    seen = set()
    for run_id in run_values:
        if pandas.isna(run_id):
            continue
        run_id = str(run_id).strip()
        if run_id == '':
            continue
        if run_id in seen:
            continue
        seen.add(run_id)
        run_ids.append(run_id)
    return run_ids


def scan_quant_abundance_paths(quant_dir, target_runs=None):
    detected_paths = {}
    target_runs = None if target_runs is None else set(target_runs)
    if os.path.exists(quant_dir) and (not os.path.isdir(quant_dir)):
        raise NotADirectoryError('Quant path exists but is not a directory: {}'.format(quant_dir))
    try:
        with os.scandir(quant_dir) as quant_entries:
            for entry in quant_entries:
                if not entry.is_dir():
                    continue
                run_id = entry.name
                if (target_runs is not None) and (run_id not in target_runs):
                    continue
                abundance_path = os.path.join(entry.path, run_id + '_abundance.tsv')
                if os.path.isfile(abundance_path):
                    detected_paths[run_id] = abundance_path
    except FileNotFoundError:
        return {}
    return detected_paths

def _find_fastp_stats_candidates(getfastq_dir, run_ids):
    candidates = []
    for run_id in run_ids:
        fastp_stats_path = os.path.join(getfastq_dir, run_id, 'fastp_stats.tsv')
        if os.path.isfile(fastp_stats_path):
            candidates.append((run_id, fastp_stats_path))
    return candidates

def _read_fastp_stats_file(sra_id, fastp_stats_path):
    try:
        df_fastp = pandas.read_csv(fastp_stats_path, sep='\t', header=0)
    except Exception as e:
        return sra_id, None, False, 'Failed to read fastp stats. Skipping {}: {}'.format(fastp_stats_path, e)
    if df_fastp.shape[0] == 0:
        return sra_id, dict(), False, None
    values = dict()
    for col in FASTP_STATS_COLUMNS:
        if col in df_fastp.columns:
            values[col] = df_fastp.loc[0, col]
    return sra_id, values, True, None

def merge_fastp_stats_into_metadata(metadata, out_dir, max_workers='auto'):
    validate_metadata_columns(
        metadata=metadata,
        required_columns=['run'],
        context='merge fastp stats',
    )
    for col in FASTP_STATS_COLUMNS:
        if col not in metadata.df.columns:
            metadata.df.loc[:, col] = numpy.nan
    getfastq_dir = os.path.realpath(os.path.join(out_dir, 'getfastq'))
    if os.path.exists(getfastq_dir) and (not os.path.isdir(getfastq_dir)):
        raise NotADirectoryError('getfastq path exists but is not a directory: {}'.format(getfastq_dir))
    if not os.path.exists(getfastq_dir):
        print('getfastq directory not found. Skipping fastp stats import: {}'.format(getfastq_dir), flush=True)
        return metadata
    run_ids = collect_valid_run_ids(metadata.df.loc[:, 'run'].values)
    normalized_runs = metadata.df.loc[:, 'run'].fillna('').astype(str).str.strip()
    if len(run_ids) == 0:
        print('No valid run IDs were found in metadata. Skipping fastp stats import.', flush=True)
        return metadata
    candidates = _find_fastp_stats_candidates(getfastq_dir=getfastq_dir, run_ids=run_ids)
    num_detected = 0
    if is_auto_parallel_option(max_workers):
        worker_cap = 8
    else:
        worker_cap = validate_positive_int_option(max_workers, 'threads')
    max_workers = min(worker_cap, len(candidates))
    results_by_candidate, failures = run_tasks_with_optional_threads(
        task_items=candidates,
        task_fn=lambda task: _read_fastp_stats_file(task[0], task[1]),
        max_workers=max_workers,
    )
    for task, exc in failures:
        results_by_candidate[task] = (
            task[0],
            None,
            False,
            'Failed to read fastp stats. Skipping {}: {}'.format(task[1], exc),
        )
    results = [results_by_candidate[task] for task in candidates if task in results_by_candidate]

    for sra_id, values, detected, error in results:
        if error is not None:
            warnings.warn(error)
            continue
        if not detected:
            continue
        is_run = (normalized_runs == sra_id)
        if not bool(is_run.any()):
            warnings.warn('Run ID from fastp stats was not found in metadata. Skipping {}.'.format(sra_id))
            continue
        for col, value in values.items():
            metadata.df.loc[is_run, col] = value
        num_detected += 1
    print('{:,} fastp stats files were detected and merged into metadata.'.format(num_detected), flush=True)
    return metadata


def collect_species_runs(metadata, sp):
    species_series = metadata.df.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
    is_sp = (species_series == sp)
    exclusion_series = metadata.df.loc[:, 'exclusion'].fillna('').astype(str).str.strip().str.lower()
    is_sampled = (exclusion_series == 'no')
    is_target = (is_sp & is_sampled)
    sra_ids = collect_valid_run_ids(metadata.df.loc[is_target, 'run'].values)
    sampled_sra_ids = set(sra_ids)
    return sra_ids, sampled_sra_ids


def collect_species_quant_outputs(sra_ids, sampled_sra_ids, detected_paths, quant_dir):
    quant_out_paths = []
    detected_sra_ids = []
    for sra_id in sra_ids:
        quant_out_path = detected_paths.get(sra_id)
        if quant_out_path is None:
            if sra_id in sampled_sra_ids:
                print('quant outfile not found: {}'.format(os.path.join(quant_dir, sra_id, sra_id + '_abundance.tsv')))
            continue
        quant_out_paths.append(quant_out_path)
        detected_sra_ids.append(sra_id)
    return detected_sra_ids, quant_out_paths


def load_quant_tables_once(detected_sra_ids, quant_out_paths, value_columns):
    table_values = {col: [None] * len(quant_out_paths) for col in value_columns}
    target_ids = None

    def read_one_quant_table(file_idx):
        sra_id = detected_sra_ids[file_idx]
        quant_out_path = quant_out_paths[file_idx]
        usecols = ['target_id'] + value_columns
        try:
            quant_df = pandas.read_csv(
                quant_out_path,
                header=0,
                sep='\t',
                usecols=usecols,
            )
        except Exception as e:
            raise ValueError(
                'Failed to read quant output table for run {} ({}): {}'.format(
                    sra_id,
                    quant_out_path,
                    e,
                )
            ) from e
        row_values = {col: quant_df[col].to_numpy() for col in value_columns}
        current_target_ids = quant_df['target_id'].to_numpy()
        return sra_id, current_target_ids, row_values

    if len(quant_out_paths) <= 1:
        results_by_idx = {}
        for file_idx in range(len(quant_out_paths)):
            results_by_idx[file_idx] = read_one_quant_table(file_idx)
    else:
        max_workers = min(MERGE_QUANT_READ_MAX_WORKERS, len(quant_out_paths))
        print(
            'Reading {:,} quant abundance files with {:,} parallel worker(s).'.format(
                len(quant_out_paths),
                max_workers,
            ),
            flush=True,
        )
        task_indices = list(range(len(quant_out_paths)))
        results_by_idx, failures = run_tasks_with_optional_threads(
            task_items=task_indices,
            task_fn=read_one_quant_table,
            max_workers=max_workers,
        )
        if failures:
            details = '; '.join([
                '{} ({}): {}'.format(detected_sra_ids[file_idx], quant_out_paths[file_idx], exc)
                for file_idx, exc in failures
            ])
            raise ValueError('Failed to read one or more quant output tables. {}'.format(details))

    for file_idx in range(len(quant_out_paths)):
        sra_id, current_target_ids, row_values = results_by_idx[file_idx]
        if file_idx == 0:
            target_ids = current_target_ids
        elif (
            (len(current_target_ids) != len(target_ids))
            or (not numpy.array_equal(current_target_ids, target_ids))
        ):
            first_run = detected_sra_ids[0] if len(detected_sra_ids) > 0 else 'unknown'
            txt = (
                'Mismatched target_id rows across quant files for one species. '
                'First run: {}, current run: {} ({})'
            )
            raise ValueError(txt.format(first_run, sra_id, quant_out_paths[file_idx]))
        for col in value_columns:
            table_values[col][file_idx] = row_values[col]
    return target_ids, table_values


def write_species_merged_quant_tables(merge_species_dir, sp_filled, detected_sra_ids, target_ids, table_values, value_columns):
    for col in value_columns:
        merged = pandas.DataFrame({'target_id': target_ids})
        for file_idx, sra_id in enumerate(detected_sra_ids):
            merged[sra_id] = table_values[col][file_idx]
        outfile_name = sp_filled + '_' + col + '.tsv'
        outfile = os.path.join(merge_species_dir, outfile_name)
        print('Writing output file:', outfile)
        merged.to_csv(outfile, sep='\t', index=False)


def merge_species_quant_tables(sp, metadata, quant_dir, merge_dir, run_abundance_paths=None):
    print('processing: {}'.format(sp), flush=True)
    sp_filled = sp.replace(' ', '_')
    merge_species_dir = os.path.join(merge_dir, sp_filled)
    sra_ids, sampled_sra_ids = collect_species_runs(metadata, sp)
    sra_id_set = set(sra_ids)
    if len(sra_ids) == 0:
        warnings.warn('No SRA Run ID found. Skipping: {}'.format(sp))
        return 0
    if isinstance(run_abundance_paths, dict):
        detected_paths = run_abundance_paths
    else:
        detected_paths = scan_quant_abundance_paths(quant_dir=quant_dir, target_runs=sra_id_set)

    detected_sra_ids, quant_out_paths = collect_species_quant_outputs(
        sra_ids=sra_ids,
        sampled_sra_ids=sampled_sra_ids,
        detected_paths=detected_paths,
        quant_dir=quant_dir,
    )
    print('{:,} quant outfiles were detected.'.format(len(quant_out_paths)))
    if len(quant_out_paths) == 0:
        return 0
    os.makedirs(merge_species_dir, exist_ok=True)

    value_columns = ['eff_length', 'est_counts', 'tpm']
    target_ids, table_values = load_quant_tables_once(
        detected_sra_ids=detected_sra_ids,
        quant_out_paths=quant_out_paths,
        value_columns=value_columns,
    )
    write_species_merged_quant_tables(
        merge_species_dir=merge_species_dir,
        sp_filled=sp_filled,
        detected_sra_ids=detected_sra_ids,
        target_ids=target_ids,
        table_values=table_values,
        value_columns=value_columns,
    )
    return len(quant_out_paths)


def run_merge_species_jobs(metadata, quant_dir, merge_dir, run_abundance_paths, species_jobs):
    spp = (
        metadata.df.loc[:, 'scientific_name']
        .fillna('')
        .astype(str)
        .str.strip()
    )
    spp = spp.loc[spp != ''].drop_duplicates().to_numpy()
    if len(spp) == 0:
        raise ValueError('No valid scientific_name entries were found in metadata for merge.')
    if (species_jobs == 1) or (len(spp) <= 1):
        total_detected = 0
        for sp in spp:
            num_detected = merge_species_quant_tables(
                sp=sp,
                metadata=metadata,
                quant_dir=quant_dir,
                merge_dir=merge_dir,
                run_abundance_paths=run_abundance_paths,
            )
            total_detected += int(num_detected)
        if total_detected == 0:
            raise FileNotFoundError('No quant abundance file was detected for any species.')
        return

    max_workers = min(species_jobs, len(spp))
    print('Running merge for {:,} species with {:,} parallel jobs.'.format(len(spp), max_workers), flush=True)
    counts_by_species, failures = run_tasks_with_optional_threads(
        task_items=spp,
        task_fn=lambda sp: merge_species_quant_tables(
            sp,
            metadata,
            quant_dir,
            merge_dir,
            run_abundance_paths,
        ),
        max_workers=max_workers,
    )
    if failures:
        details = '; '.join(['{}: {}'.format(sp, err) for sp, err in failures])
        raise RuntimeError('merge failed for {}/{} species. {}'.format(len(failures), len(spp), details))
    total_detected = 0
    for sp in spp:
        detected = counts_by_species.get(sp, 0)
        if detected is None:
            detected = 0
        total_detected += int(detected)
    if total_detected == 0:
        raise FileNotFoundError('No quant abundance file was detected for any species.')


def run_merge_plot_rscript(merge_dir, path_metadata_merge):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    r_merge_path = os.path.join(script_dir, 'merge.r')
    r_util_path = os.path.join(script_dir, 'util.r')
    r_command = ['Rscript', r_merge_path, merge_dir, path_metadata_merge, r_util_path]
    print('Starting R script for plot generation: {}'.format(' '.join(r_command)), flush=True)
    subprocess.check_call(r_command)


def merge_main(args):
    check_rscript()
    species_jobs, _ = resolve_worker_allocation(
        requested_workers=getattr(args, 'internal_jobs', 'auto'),
        requested_threads=getattr(args, 'threads', 'auto'),
        internal_cpu_budget=getattr(args, 'internal_cpu_budget', 'auto'),
        worker_option_name='internal_jobs',
        context='merge:',
    )
    postprocess_workers, _ = resolve_worker_allocation(
        requested_workers='auto',
        requested_threads=getattr(args, 'threads', 'auto'),
        internal_cpu_budget=getattr(args, 'internal_cpu_budget', 'auto'),
        worker_option_name='threads',
        context='merge postprocess:',
    )
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    quant_dir = os.path.realpath(os.path.join(out_dir, 'quant'))
    merge_dir = os.path.realpath(os.path.join(out_dir, 'merge'))
    if os.path.exists(merge_dir) and (not os.path.isdir(merge_dir)):
        raise NotADirectoryError('Merge path exists but is not a directory: {}'.format(merge_dir))
    os.makedirs(merge_dir, exist_ok=True)
    metadata = load_metadata(args)
    validate_metadata_columns(
        metadata=metadata,
        required_columns=['run', 'scientific_name', 'exclusion'],
        context='merge',
    )
    run_abundance_paths = scan_quant_abundance_paths(
        quant_dir=quant_dir,
        target_runs=set(collect_valid_run_ids(metadata.df.loc[:, 'run'].values)),
    )
    print('Detected {:,} quant abundance files across all runs.'.format(len(run_abundance_paths)), flush=True)
    run_merge_species_jobs(
        metadata=metadata,
        quant_dir=quant_dir,
        merge_dir=merge_dir,
        run_abundance_paths=run_abundance_paths,
        species_jobs=species_jobs,
    )
    print('Getting mapping rate from quant output and write new metadata file into merge directory.', flush=True)
    metadata = merge_fastp_stats_into_metadata(metadata, out_dir, max_workers=postprocess_workers)
    path_metadata_merge = os.path.realpath(os.path.join(out_dir, 'merge', 'metadata.tsv'))
    write_updated_metadata(metadata, path_metadata_merge, args, max_workers=postprocess_workers)
    run_merge_plot_rscript(merge_dir=merge_dir, path_metadata_merge=path_metadata_merge)
