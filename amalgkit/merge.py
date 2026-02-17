import pandas
import numpy

import os
import warnings
from amalgkit.util import *

FASTP_STATS_COLUMNS = ['fastp_duplication_rate', 'fastp_insert_size_peak']


def scan_quant_abundance_paths(quant_dir, target_runs=None):
    detected_paths = {}
    target_runs = None if target_runs is None else set(target_runs)
    try:
        with os.scandir(quant_dir) as quant_entries:
            for entry in quant_entries:
                if not entry.is_dir():
                    continue
                run_id = entry.name
                if (target_runs is not None) and (run_id not in target_runs):
                    continue
                abundance_path = os.path.join(entry.path, run_id + '_abundance.tsv')
                if os.path.exists(abundance_path):
                    detected_paths[run_id] = abundance_path
    except FileNotFoundError:
        return {}
    return detected_paths

def _find_fastp_stats_candidates(getfastq_dir, run_ids):
    candidates = []
    for run_id in run_ids:
        fastp_stats_path = os.path.join(getfastq_dir, run_id, 'fastp_stats.tsv')
        if os.path.exists(fastp_stats_path):
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

def merge_fastp_stats_into_metadata(metadata, out_dir):
    for col in FASTP_STATS_COLUMNS:
        if col not in metadata.df.columns:
            metadata.df.loc[:, col] = numpy.nan
    getfastq_dir = os.path.realpath(os.path.join(out_dir, 'getfastq'))
    if not os.path.exists(getfastq_dir):
        print('getfastq directory not found. Skipping fastp stats import: {}'.format(getfastq_dir), flush=True)
        return metadata
    run_ids = set(metadata.df.loc[:, 'run'].values)
    candidates = _find_fastp_stats_candidates(getfastq_dir=getfastq_dir, run_ids=run_ids)
    num_detected = 0
    max_workers = min(8, len(candidates))
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
        idx = get_metadata_row_index_by_run(metadata, sra_id)
        for col, value in values.items():
            metadata.df.at[idx, col] = value
        num_detected += 1
    print('{:,} fastp stats files were detected and merged into metadata.'.format(num_detected), flush=True)
    return metadata


def collect_species_runs(metadata, sp):
    is_sp = (metadata.df.loc[:, 'scientific_name'] == sp)
    sra_ids = metadata.df.loc[is_sp, 'run'].values
    is_sampled = (metadata.df.loc[:, 'exclusion'] == 'no')
    sampled_sra_ids = set(metadata.df.loc[is_sampled, 'run'].values)
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
    # Single-pass read: each quant file is read once and split into output arrays.
    for file_idx, (sra_id, quant_out_path) in enumerate(zip(detected_sra_ids, quant_out_paths)):
        usecols = ['target_id'] + value_columns if file_idx == 0 else value_columns
        quant_df = pandas.read_csv(
            quant_out_path,
            header=0,
            sep='\t',
            usecols=usecols,
        )
        if file_idx == 0:
            target_ids = quant_df['target_id'].to_numpy()
        for col in value_columns:
            table_values[col][file_idx] = quant_df[col].to_numpy()
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
    spp = metadata.df.loc[:, 'scientific_name'].dropna().unique()
    if (species_jobs == 1) or (len(spp) <= 1):
        for sp in spp:
            merge_species_quant_tables(
                sp=sp,
                metadata=metadata,
                quant_dir=quant_dir,
                merge_dir=merge_dir,
                run_abundance_paths=run_abundance_paths,
            )
        return

    max_workers = min(species_jobs, len(spp))
    print('Running merge for {:,} species with {:,} parallel jobs.'.format(len(spp), max_workers), flush=True)
    _, failures = run_tasks_with_optional_threads(
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


def run_merge_plot_rscript(merge_dir, path_metadata_merge):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    r_merge_path = os.path.join(script_dir, 'merge.r')
    r_util_path = os.path.join(script_dir, 'util.r')
    r_command = ['Rscript', r_merge_path, merge_dir, path_metadata_merge, r_util_path]
    print('Starting R script for plot generation: {}'.format(' '.join(r_command)), flush=True)
    subprocess.check_call(r_command)


def merge_main(args):
    species_jobs, _ = resolve_worker_allocation(
        requested_workers=getattr(args, 'species_jobs', 1),
        cpu_budget=getattr(args, 'cpu_budget', 0),
        worker_option_name='species_jobs',
        context='merge:',
    )
    quant_dir = os.path.realpath(os.path.join(args.out_dir, 'quant'))
    merge_dir = os.path.realpath(os.path.join(args.out_dir, 'merge'))
    os.makedirs(merge_dir, exist_ok=True)
    metadata = load_metadata(args)
    run_abundance_paths = scan_quant_abundance_paths(
        quant_dir=quant_dir,
        target_runs=set(metadata.df.loc[:, 'run'].values),
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
    metadata = merge_fastp_stats_into_metadata(metadata, args.out_dir)
    path_metadata_merge = os.path.realpath(os.path.join(args.out_dir, 'merge', 'metadata.tsv'))
    write_updated_metadata(metadata, path_metadata_merge, args)
    run_merge_plot_rscript(merge_dir=merge_dir, path_metadata_merge=path_metadata_merge)
