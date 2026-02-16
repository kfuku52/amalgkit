import pandas
import numpy

import os
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from amalgkit.util import *

FASTP_STATS_COLUMNS = ['fastp_duplication_rate', 'fastp_insert_size_peak']


def scan_quant_abundance_paths(quant_dir, target_runs=None):
    detected_paths = {}
    if target_runs is None:
        target_runs = None
    else:
        target_runs = set(target_runs)
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
    candidates = []
    run_ids = set(metadata.df.loc[:, 'run'].values)
    with os.scandir(getfastq_dir) as entries:
        for entry in entries:
            if (not entry.is_dir()) or (entry.name not in run_ids):
                continue
            fastp_stats_path = os.path.join(entry.path, 'fastp_stats.tsv')
            if os.path.exists(fastp_stats_path):
                candidates.append((entry.name, fastp_stats_path))
    num_detected = 0
    results = []
    if len(candidates) <= 1:
        for sra_id, fastp_stats_path in candidates:
            results.append(_read_fastp_stats_file(sra_id, fastp_stats_path))
    else:
        max_workers = min(8, len(candidates))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(_read_fastp_stats_file, sra_id, fastp_stats_path): sra_id
                for sra_id, fastp_stats_path in candidates
            }
            for future in as_completed(futures):
                results.append(future.result())

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


def merge_species_quant_tables(sp, metadata, quant_dir, merge_dir, run_abundance_paths=None):
    print('processing: {}'.format(sp), flush=True)
    sp_filled = sp.replace(' ', '_')
    merge_species_dir = os.path.join(os.path.join(merge_dir, sp_filled))
    is_sp = (metadata.df.loc[:,'scientific_name']==sp)
    sra_ids = metadata.df.loc[is_sp,'run'].values
    sra_id_set = set(sra_ids)
    is_sampled = (metadata.df.loc[:,'exclusion']=='no')
    sampled_sra_ids = set(metadata.df.loc[is_sampled,'run'].values)
    if len(sra_ids)==0:
        warnings.warn('No SRA Run ID found. Skipping: {}'.format(sp))
        return 0
    if isinstance(run_abundance_paths, dict):
        detected_paths = run_abundance_paths
    else:
        detected_paths = scan_quant_abundance_paths(quant_dir=quant_dir, target_runs=sra_id_set)

    quant_out_paths = []
    detected_sra_ids = []
    for sra_id in sra_ids:
        quant_out_path = detected_paths.get(sra_id)
        if quant_out_path is not None:
            os.makedirs(merge_species_dir, exist_ok=True)
            quant_out_paths.append(quant_out_path)
            detected_sra_ids.append(sra_id)
            continue
        if sra_id in sampled_sra_ids:
            print('quant outfile not found: {}'.format(os.path.join(quant_dir, sra_id, sra_id + '_abundance.tsv')))
    print('{:,} quant outfiles were detected.'.format(len(quant_out_paths)))
    if len(quant_out_paths) == 0:
        return 0

    value_columns = ['eff_length', 'est_counts', 'tpm']
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

    for col in value_columns:
        merged = pandas.DataFrame({'target_id': target_ids})
        for file_idx, sra_id in enumerate(detected_sra_ids):
            merged[sra_id] = table_values[col][file_idx]
        outfile_name = sp_filled+'_'+col+'.tsv'
        outfile = os.path.join(merge_species_dir, outfile_name)
        print('Writing output file:', outfile)
        merged.to_csv(outfile, sep='\t', index=False)
    return len(quant_out_paths)


def merge_main(args):
    species_jobs = int(getattr(args, 'species_jobs', 1))
    if species_jobs <= 0:
        raise ValueError('--species_jobs must be > 0.')
    quant_dir = os.path.realpath(os.path.join(args.out_dir, 'quant'))
    merge_dir = os.path.realpath(os.path.join(args.out_dir, 'merge'))
    if not os.path.exists(merge_dir):
        os.makedirs(os.path.join(merge_dir))
    metadata = load_metadata(args)
    run_abundance_paths = scan_quant_abundance_paths(
        quant_dir=quant_dir,
        target_runs=set(metadata.df.loc[:, 'run'].values),
    )
    print('Detected {:,} quant abundance files across all runs.'.format(len(run_abundance_paths)), flush=True)
    spp = metadata.df.loc[:,'scientific_name'].dropna().unique()
    if (species_jobs == 1) or (len(spp) <= 1):
        for sp in spp:
            merge_species_quant_tables(
                sp=sp,
                metadata=metadata,
                quant_dir=quant_dir,
                merge_dir=merge_dir,
                run_abundance_paths=run_abundance_paths,
            )
    else:
        max_workers = min(species_jobs, len(spp))
        print('Running merge for {:,} species with {:,} parallel jobs.'.format(len(spp), max_workers), flush=True)
        failures = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    merge_species_quant_tables,
                    sp,
                    metadata,
                    quant_dir,
                    merge_dir,
                    run_abundance_paths,
                ): sp
                for sp in spp
            }
            for future in as_completed(futures):
                sp = futures[future]
                try:
                    future.result()
                except Exception as exc:
                    failures.append((sp, exc))
        if failures:
            details = '; '.join(['{}: {}'.format(sp, err) for sp, err in failures])
            raise RuntimeError('merge failed for {}/{} species. {}'.format(len(failures), len(spp), details))
    print('Getting mapping rate from quant output and write new metadata file into merge directory.', flush=True)
    metadata = merge_fastp_stats_into_metadata(metadata, args.out_dir)
    path_metadata_merge = os.path.realpath(os.path.join(args.out_dir, 'merge', 'metadata.tsv'))
    write_updated_metadata(metadata, path_metadata_merge, args)
    r_merge_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'merge.r')
    r_util_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'util.r')
    r_command = ['Rscript', r_merge_path, merge_dir, path_metadata_merge, r_util_path]
    print('Starting R script for plot generation: {}'.format(' '.join(r_command)), flush=True)
    subprocess.check_call(r_command)
