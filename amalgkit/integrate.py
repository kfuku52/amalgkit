import pandas as pd

import gzip
import os
import re
import warnings

from amalgkit.sanity import check_getfastq_outputs
from amalgkit.util import *

FASTQ_EXTENSIONS = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
FASTQ_MATE_SUFFIX_PATTERN = re.compile(r'^(.*)_([12])$')
QUICK_MODE_SAMPLE_READS = 1000
QUICK_MODE_PAIRED_COUNT_TOLERANCE_FRACTION = 0.05
QUICK_MODE_PAIRED_COUNT_TOLERANCE_MIN_READS = 10
FASTQ_LINECOUNT_CHUNK_SIZE = 16 * 1024 * 1024
PRIVATE_FASTQ_METADATA_COLUMNS = [
    'scientific_name', 'sample_group', 'run', 'read1_path', 'read2_path', 'is_sampled',
    'is_qualified', 'exclusion', 'lib_layout', 'spot_length', 'total_spots', 'total_bases', 'size', 'private_file',
]

def open_fastq_text(path_fastq):
    path_fastq_lower = path_fastq.lower()
    if path_fastq_lower.endswith(('.fq.gz', '.fastq.gz')):
        return gzip.open(path_fastq, 'rt')
    return open(path_fastq, 'rt')

def open_fastq_binary(path_fastq):
    path_fastq_lower = path_fastq.lower()
    if path_fastq_lower.endswith(('.fq.gz', '.fastq.gz')):
        return gzip.open(path_fastq, 'rb')
    return open(path_fastq, 'rb')

def parse_fastq_basename(basename):
    basename_lower = basename.lower()
    for ext in FASTQ_EXTENSIONS:
        if basename_lower.endswith(ext):
            stem = basename[:-len(ext)]
            matched = FASTQ_MATE_SUFFIX_PATTERN.match(stem)
            if matched is not None:
                return matched.group(1), matched.group(2)
            return stem, None
    return None, None


def scan_fastq_directory(fastq_dir):
    grouped_files = {}
    with os.scandir(fastq_dir) as entries:
        for entry in entries:
            if not entry.is_file():
                continue
            run_id, mate = parse_fastq_basename(entry.name)
            if run_id is None:
                continue
            run_id = str(run_id).strip()
            if run_id == '':
                continue
            grouped_files.setdefault(run_id, []).append({
                'basename': entry.name,
                'run': run_id,
                'mate': mate,
                'path': entry.path,
            })
    return grouped_files

def write_fastq_head(path_fastq, path_out, max_lines=4000):
    with open_fastq_text(path_fastq) as fin, open(path_out, 'wt') as fout:
        for i, line in enumerate(fin):
            if i >= max_lines:
                break
            fout.write(line)

def count_fastq_lines(path_fastq):
    num_lines = 0
    with open_fastq_binary(path_fastq) as f:
        while True:
            chunk = f.read(FASTQ_LINECOUNT_CHUNK_SIZE)
            if not chunk:
                break
            num_lines += chunk.count(b'\n')
    return num_lines

def estimate_total_spots_from_line_count(path_fastq):
    total_lines = count_fastq_lines(path_fastq)
    if (total_lines % 4) != 0:
        raise ValueError('Malformed FASTQ (line count not divisible by 4): {}'.format(path_fastq))
    return int(total_lines / 4)

def _sequence_line_length(line_bytes):
    seq_len = len(line_bytes)
    if seq_len == 0:
        return 0
    if line_bytes.endswith(b'\n'):
        seq_len -= 1
        if (seq_len > 0) and line_bytes.endswith(b'\r\n'):
            seq_len -= 1
    return seq_len

def sample_fastq_reads(path_fastq, max_reads=None):
    num_reads = 0
    total_bases = 0
    total_record_chars = 0
    reached_eof = True
    with open_fastq_binary(path_fastq) as f:
        while True:
            line1 = f.readline()
            if line1 == b'':
                break
            line2 = f.readline()
            line3 = f.readline()
            line4 = f.readline()
            if (line2 == b'') or (line3 == b'') or (line4 == b''):
                raise ValueError('Malformed FASTQ (record truncated): {}'.format(path_fastq))
            num_reads += 1
            total_bases += _sequence_line_length(line2)
            total_record_chars += (len(line1) + len(line2) + len(line3) + len(line4))
            if (max_reads is not None) and (num_reads >= max_reads):
                reached_eof = (f.readline() == b'')
                break
    avg_len = int(total_bases / num_reads) if num_reads > 0 else 0
    return num_reads, avg_len, total_record_chars, reached_eof

def scan_fastq_reads(path_fastq, max_reads=None):
    num_reads, avg_len, _, _ = sample_fastq_reads(path_fastq, max_reads=max_reads)
    return num_reads, avg_len

def get_gzip_isize(path_fastq):
    with open(path_fastq, 'rb') as f:
        try:
            f.seek(-4, os.SEEK_END)
        except OSError as exc:
            raise ValueError('Malformed gzip file (too small for ISIZE footer): {}'.format(path_fastq)) from exc
        isize_bytes = f.read(4)
    if len(isize_bytes) != 4:
        raise ValueError('Malformed gzip file (missing ISIZE footer): {}'.format(path_fastq))
    return int.from_bytes(isize_bytes, byteorder='little', signed=False)

def estimate_total_spots_from_gzip_sample(path_fastq, sampled_reads, sampled_record_chars):
    if sampled_reads <= 0:
        raise ValueError('No sampled reads available to estimate total spots: {}'.format(path_fastq))
    if sampled_record_chars <= 0:
        raise ValueError('No sampled record size available to estimate total spots: {}'.format(path_fastq))
    compressed_size = os.path.getsize(path_fastq)
    if compressed_size >= (1 << 32):
        raise ValueError(
            'gzip file is >= 4 GiB and ISIZE footer may overflow. Re-run with --accurate_size yes: {}'.format(
                path_fastq
            )
        )
    avg_record_chars = sampled_record_chars / sampled_reads
    uncompressed_size = get_gzip_isize(path_fastq)
    if uncompressed_size <= 0:
        raise ValueError('gzip ISIZE footer is not positive. Re-run with --accurate_size yes: {}'.format(path_fastq))
    estimated_reads = int(round(uncompressed_size / avg_record_chars))
    return max(1, estimated_reads)

def is_approximately_equal_count(count_a, count_b, tolerance_fraction, tolerance_min):
    diff = abs(int(count_a) - int(count_b))
    max_count = max(int(count_a), int(count_b))
    tolerance = max(int(tolerance_min), int(round(max_count * float(tolerance_fraction))))
    return diff <= tolerance

def scan_fastq_stats(path_fastq, max_reads_for_average=None):
    total_reads = 0
    sampled_reads = 0
    sampled_bases = 0
    with open_fastq_binary(path_fastq) as f:
        while True:
            line1 = f.readline()
            if line1 == b'':
                break
            line2 = f.readline()
            line3 = f.readline()
            line4 = f.readline()
            if (line2 == b'') or (line3 == b'') or (line4 == b''):
                raise ValueError('Malformed FASTQ (record truncated): {}'.format(path_fastq))
            total_reads += 1
            if (max_reads_for_average is None) or (sampled_reads < max_reads_for_average):
                sampled_reads += 1
                sampled_bases += _sequence_line_length(line2)
    average_length = int(sampled_bases / sampled_reads) if sampled_reads > 0 else 0
    return total_reads, average_length, sampled_reads


def resolve_run_fastq_layout(run_id, fastq_records):
    mate_to_paths = {None: [], '1': [], '2': []}
    for record in fastq_records:
        mate_to_paths[record['mate']].append(record['path'])
    for mate_label in ['1', '2']:
        if len(mate_to_paths[mate_label]) > 1:
            raise ValueError('Multiple FASTQ files detected for {} read {}: {}'.format(
                run_id, mate_label, mate_to_paths[mate_label]
            ))

    if (len(mate_to_paths['1']) == 1) and (len(mate_to_paths['2']) == 1) and (len(mate_to_paths[None]) == 0):
        return 'paired', mate_to_paths['1'][0], mate_to_paths['2'][0]
    if (len(fastq_records) == 1) and ((len(mate_to_paths[None]) == 1) or (len(mate_to_paths['1']) == 1)):
        return 'single', fastq_records[0]['path'], 'unavailable'
    if (len(fastq_records) == 1) and (len(mate_to_paths['2']) == 1):
        raise ValueError(
            'Only one paired-end mate file was detected for {}. '
            'Expected both *_1 and *_2 files: {}'.format(
                run_id,
                [record['basename'] for record in fastq_records],
            )
        )
    raise ValueError('Ambiguous FASTQ file set for {}: {}'.format(
        run_id, [record['basename'] for record in fastq_records]
    ))


def scan_run_fastq_stats(run_spec, accurate_size):
    run_id = run_spec['run']
    lib_layout = run_spec['lib_layout']
    read1_path = run_spec['read1_path']
    read2_path = run_spec['read2_path']
    print("Found {} file(s) for ID {}. Lib-layout: {}".format(run_spec['num_files'], run_id, lib_layout), flush=True)
    print("Getting sequence statistics.", flush=True)
    read1_path_lower = read1_path.lower()
    if read1_path_lower.endswith(('.fq', '.fastq')):
        is_decompressed = True
    elif read1_path_lower.endswith(('.fq.gz', '.fastq.gz')):
        is_decompressed = False
    else:
        warnings.warn("{} is not a fastq file. Skipping.".format(read1_path))
        return None

    quick_gzip_mode = (not accurate_size) and (not is_decompressed)
    if accurate_size or is_decompressed:
        print('--accurate_size set to yes. Running accurate sequence scan with Python FASTQ parser.')
        total_spots, avg_len_read1, _ = scan_fastq_stats(path_fastq=read1_path, max_reads_for_average=None)
    else:
        print('--accurate_size set to no. Running quick sequence scan (first 1,000 reads + gzip size estimate).')
        sampled_reads, avg_len_read1, sampled_record_chars, reached_eof_read1 = sample_fastq_reads(
            path_fastq=read1_path,
            max_reads=QUICK_MODE_SAMPLE_READS,
        )
        if sampled_reads == 0:
            raise ValueError('No reads detected in FASTQ file: {}'.format(read1_path))
        if reached_eof_read1:
            total_spots = sampled_reads
        else:
            total_spots = estimate_total_spots_from_gzip_sample(
                path_fastq=read1_path,
                sampled_reads=sampled_reads,
                sampled_record_chars=sampled_record_chars,
            )
    if total_spots <= 0:
        raise ValueError('No reads detected in FASTQ file: {}'.format(read1_path))

    avg_len_read2 = avg_len_read1
    if lib_layout == 'paired':
        if quick_gzip_mode:
            sampled_reads_read2, avg_len_read2, sampled_record_chars_read2, reached_eof_read2 = sample_fastq_reads(
                path_fastq=read2_path,
                max_reads=QUICK_MODE_SAMPLE_READS,
            )
            if sampled_reads_read2 == 0:
                raise ValueError('No reads detected in FASTQ file: {}'.format(read2_path))
            if reached_eof_read2:
                read2_spots = sampled_reads_read2
            else:
                read2_spots = estimate_total_spots_from_gzip_sample(
                    path_fastq=read2_path,
                    sampled_reads=sampled_reads_read2,
                    sampled_record_chars=sampled_record_chars_read2,
                )
            if read2_spots <= 0:
                raise ValueError('No reads detected in FASTQ file: {}'.format(read2_path))
            if not is_approximately_equal_count(
                total_spots,
                read2_spots,
                tolerance_fraction=QUICK_MODE_PAIRED_COUNT_TOLERANCE_FRACTION,
                tolerance_min=QUICK_MODE_PAIRED_COUNT_TOLERANCE_MIN_READS,
            ):
                raise ValueError(
                    'Estimated paired-end read counts for {} differ (read1={}, read2={}). '
                    'Re-run with --accurate_size yes for exact validation.'.format(
                        run_id, total_spots, read2_spots
                    )
                )
            total_spots = min(total_spots, read2_spots)
        else:
            read2_spots, avg_len_read2, _ = scan_fastq_stats(path_fastq=read2_path, max_reads_for_average=None)
            if read2_spots <= 0:
                raise ValueError('No reads detected in FASTQ file: {}'.format(read2_path))
            if read2_spots != total_spots:
                raise ValueError(
                    'Mismatched paired-end read counts for {}: read1 has {} reads, read2 has {} reads.'.format(
                        run_id, total_spots, read2_spots
                    )
                )

    total_bases = total_spots * int(avg_len_read1)
    total_size = os.path.getsize(read1_path)
    if lib_layout == 'paired':
        total_bases = total_spots * (int(avg_len_read1) + int(avg_len_read2))
        total_size += os.path.getsize(read2_path)
    return {
        'run': run_id,
        'read1_path': os.path.abspath(read1_path),
        'read2_path': os.path.abspath(read2_path) if read2_path != 'unavailable' else 'unavailable',
        'lib_layout': lib_layout,
        'total_spots': total_spots,
        'spot_length': int(avg_len_read1),
        'size': total_size,
        'total_bases': total_bases,
    }


def build_run_specs(grouped_files):
    run_specs = []
    for run_id in sorted(grouped_files.keys()):
        fastq_records = grouped_files[run_id]
        lib_layout, read1_path, read2_path = resolve_run_fastq_layout(run_id, fastq_records)
        run_specs.append({
            'run': run_id,
            'num_files': len(fastq_records),
            'lib_layout': lib_layout,
            'read1_path': read1_path,
            'read2_path': read2_path,
        })
    return run_specs


def scan_all_run_fastq_stats(run_specs, accurate_size, threads=1):
    if is_auto_parallel_option(threads):
        jobs = resolve_detected_cpu_count()
    else:
        jobs = validate_positive_int_option(threads, 'threads')
    if (jobs == 1) or (len(run_specs) <= 1):
        stats_by_run = {}
        for run_spec in run_specs:
            stats = scan_run_fastq_stats(run_spec, accurate_size)
            if stats is not None:
                stats_by_run[run_spec['run']] = stats
        return stats_by_run

    max_workers = min(jobs, len(run_specs))
    print('Scanning FASTQ stats for {:,} runs with {:,} parallel worker(s).'.format(len(run_specs), max_workers), flush=True)
    run_spec_indices = list(range(len(run_specs)))
    stats_by_spec_idx, failures = run_tasks_with_optional_threads(
        task_items=run_spec_indices,
        task_fn=lambda idx: scan_run_fastq_stats(run_specs[idx], accurate_size),
        max_workers=max_workers,
    )
    if failures:
        details = '; '.join(['{}: {}'.format(run_specs[idx]['run'], exc) for idx, exc in failures])
        raise RuntimeError('FASTQ stats scan failed for {}/{} runs. {}'.format(len(failures), len(run_specs), details))

    stats_by_run = {}
    for idx, stats in stats_by_spec_idx.items():
        if stats is not None:
            stats_by_run[run_specs[idx]['run']] = stats
    return stats_by_run


def build_private_fastq_metadata_rows(run_specs, stats_by_run):
    rows = []
    for run_spec in run_specs:
        run_id = run_spec['run']
        if run_id not in stats_by_run:
            continue
        stats = stats_by_run[run_id]
        rows.append({
            'scientific_name': 'Please add in format: Genus species',
            'sample_group': 'Please add',
            'run': run_id,
            'read1_path': stats['read1_path'],
            'read2_path': stats['read2_path'],
            'is_sampled': 'yes',
            'is_qualified': 'yes',
            'exclusion': 'no',
            'lib_layout': stats['lib_layout'],
            'spot_length': stats['spot_length'],
            'total_spots': stats['total_spots'],
            'total_bases': stats['total_bases'],
            'size': stats['size'],
            'private_file': 'yes',
        })
    return rows


def get_fastq_stats(args):
    print("Starting integration of fastq-file metadata...")
    if (getattr(args, 'fastq_dir', None) is None) or (str(args.fastq_dir).strip() == ''):
        raise ValueError('--fastq_dir is required.')
    if not os.path.exists(args.fastq_dir):
        raise ValueError("PATH to fastq directory does not exist: {}".format(args.fastq_dir))
    if not os.path.isdir(args.fastq_dir):
        raise NotADirectoryError('PATH to fastq directory is not a directory: {}'.format(args.fastq_dir))
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    metadata_dir = os.path.join(out_dir, 'metadata')
    if os.path.exists(metadata_dir) and (not os.path.isdir(metadata_dir)):
        raise NotADirectoryError('Metadata output path exists but is not a directory: {}'.format(metadata_dir))
    grouped_files = scan_fastq_directory(args.fastq_dir)
    if len(grouped_files) == 0:
        txt = 'No detected fastq files (extensions: {}) in: {}'.format(', '.join(FASTQ_EXTENSIONS), args.fastq_dir)
        raise ValueError(txt)
    run_specs = build_run_specs(grouped_files)
    stats_by_run = scan_all_run_fastq_stats(
        run_specs=run_specs,
        accurate_size=args.accurate_size,
        threads=getattr(args, 'threads', 1),
    )
    rows = build_private_fastq_metadata_rows(run_specs, stats_by_run)
    tmp_metadata = pd.DataFrame.from_records(rows, columns=PRIVATE_FASTQ_METADATA_COLUMNS)
    os.makedirs(metadata_dir, exist_ok=True)
    tmp_metadata = tmp_metadata.sort_values(by='run', axis=0, ascending=True).reset_index(drop=True)
    tmp_metadata.to_csv(os.path.join(out_dir, 'metadata_private_fastq.tsv'), sep='\t', index=False)
    return tmp_metadata

def integrate_main(args):
    if args.metadata=='inferred':
        relative_path = os.path.join(args.out_dir, 'metadata', 'metadata.tsv')
        metadata_path = os.path.realpath(relative_path)
    else:
        metadata_path = os.path.realpath(args.metadata)
        if not os.path.exists(metadata_path):
            raise FileNotFoundError('Metadata file not found: {}'.format(metadata_path))
    if os.path.exists(metadata_path):
        if not os.path.isfile(metadata_path):
            raise IsADirectoryError('Metadata path exists but is not a file: {}'.format(metadata_path))
        print('Merging existing metadata and private fastq info: {}'.format(metadata_path))
        metadata = load_metadata(args)
        if 'run' not in metadata.df.columns:
            raise ValueError('Column "run" is required in metadata.')
        metadata.df.loc[:,'private_file'] = 'no'
        normalized_runs = metadata.df.loc[:, 'run'].fillna('').astype(str).str.strip()
        metadata.df['run'] = normalized_runs
        if (normalized_runs == '').any():
            raise ValueError('Missing Run ID(s) were detected in metadata.')
        duplicate_mask = normalized_runs.duplicated(keep=False)
        if duplicate_mask.any():
            duplicated_runs = normalized_runs.loc[duplicate_mask].drop_duplicates().tolist()
            raise ValueError('Duplicate Run ID(s) were detected in metadata: {}'.format(', '.join(duplicated_runs)))
        print("scanning for getfastq output")
        sra_ids = normalized_runs
        data_available, data_unavailable = check_getfastq_outputs(args, sra_ids, metadata, args.out_dir)
        metadata.df.loc[metadata.df['run'].isin(data_available), 'data_available'] = 'yes'
        metadata.df.loc[metadata.df['run'].isin(data_unavailable), 'data_available'] = 'no'
        tmp_metadata = get_fastq_stats(args)
        df = pd.concat([metadata.df, tmp_metadata])
        merged_runs = df.loc[:, 'run'].fillna('').astype(str).str.strip()
        duplicate_mask = merged_runs.duplicated(keep=False)
        if duplicate_mask.any():
            duplicated_runs = merged_runs.loc[duplicate_mask].drop_duplicates().tolist()
            raise ValueError(
                'Duplicate Run ID(s) were detected after merging existing metadata and private fastq info: {}'.format(
                    ', '.join(duplicated_runs)
                )
            )
        df['run'] = merged_runs
        df.to_csv(os.path.join(args.out_dir, 'metadata', 'metadata_updated_for_private_fastq.tsv'), sep='\t', index=False)
    else:
        print('Generating a new metadata table.')
        get_fastq_stats(args)
