import pandas as pd

import os
import re
import subprocess
import warnings

from amalgkit.download_utils import get_ete_ncbitaxa
from amalgkit.fastq_utils import (
    count_fastq_lines,
    map_seqkit_stats_rows,
    open_fastq_binary,
    open_fastq_text,
    parse_seqkit_stats_row_num_reads_avg_len as parse_seqkit_stats_row,
    parse_seqkit_stats_rows,
    sequence_line_length as _sequence_line_length,
    sample_fastq_reads,
)
from amalgkit.metadata_utils import Metadata, load_metadata
from amalgkit.output_utils import atomic_write_dataframe
from amalgkit.parallel_utils import (
    is_auto_parallel_option,
    resolve_detected_cpu_count,
    run_tasks_with_optional_threads,
    validate_positive_int_option,
)
from amalgkit.sanity import check_getfastq_outputs
from amalgkit.subprocess_utils import run_logged_command

FASTQ_EXTENSIONS = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
FASTQ_MATE_SUFFIX_PATTERN = re.compile(r'^(.*)_([12])$')
QUICK_MODE_SAMPLE_READS = 1000
QUICK_MODE_PAIRED_COUNT_TOLERANCE_FRACTION = 0.05
QUICK_MODE_PAIRED_COUNT_TOLERANCE_MIN_READS = 10
PRIVATE_FASTQ_METADATA_COLUMNS = [
    'scientific_name', 'sample_group', 'run', 'read1_path', 'read2_path', 'is_sampled',
    'is_qualified', 'exclusion', 'lib_layout', 'spot_length', 'total_spots', 'total_bases', 'size', 'private_file',
]
PRIVATE_FASTQ_SCIENTIFIC_NAME_PLACEHOLDER = 'Please add in format: Genus species'
PRIVATE_FASTQ_SAMPLE_GROUP_PLACEHOLDER = 'Please add'
STANDARD_TAXONOMIC_RANKS = ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

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

def get_relative_path_parts(root_dir, path_dir):
    relative_dir = os.path.relpath(path_dir, root_dir)
    if relative_dir in ('', '.'):
        return []
    return [part for part in relative_dir.split(os.sep) if part not in ('', '.')]

def infer_scientific_name_from_relative_parts(relative_parts):
    if len(relative_parts) == 0:
        return ''
    return normalize_scientific_name_for_lookup(relative_parts[0])

def build_run_id_candidates(run_basename, relative_parts, duplicate_basename):
    candidates = []
    if not duplicate_basename:
        return [run_basename]
    if len(relative_parts) > 0:
        candidates.append(relative_parts[0] + '_' + run_basename)
        candidates.append('_'.join(relative_parts) + '_' + run_basename)
    candidates.append(run_basename)
    seen = set()
    deduplicated = []
    for candidate in candidates:
        normalized = str(candidate).strip()
        if normalized == '' or normalized in seen:
            continue
        seen.add(normalized)
        deduplicated.append(normalized)
    return deduplicated

def assign_unique_run_ids(logical_run_records):
    basename_counts = {}
    for logical_key in logical_run_records.keys():
        basename_counts[logical_key[1]] = basename_counts.get(logical_key[1], 0) + 1
    assigned = {}
    used = set()
    sorted_keys = sorted(logical_run_records.keys(), key=lambda item: (item[1], item[0]))
    for logical_key in sorted_keys:
        relative_parts, run_basename = logical_key
        candidates = build_run_id_candidates(
            run_basename=run_basename,
            relative_parts=list(relative_parts),
            duplicate_basename=(basename_counts.get(run_basename, 0) > 1),
        )
        chosen = None
        for candidate in candidates:
            if candidate not in used:
                chosen = candidate
                break
        if chosen is None:
            base_candidate = candidates[0] if len(candidates) > 0 else run_basename
            suffix = 2
            chosen = base_candidate
            while chosen in used:
                chosen = '{}_{}'.format(base_candidate, suffix)
                suffix += 1
        assigned[logical_key] = chosen
        used.add(chosen)
    return assigned

def scan_fastq_directory(fastq_dir):
    logical_run_records = {}
    for root, dirnames, filenames in os.walk(fastq_dir):
        dirnames.sort()
        relative_parts = get_relative_path_parts(fastq_dir, root)
        scientific_name = infer_scientific_name_from_relative_parts(relative_parts)
        for filename in sorted(filenames):
            path_fastq = os.path.join(root, filename)
            if not os.path.isfile(path_fastq):
                continue
            run_basename, mate = parse_fastq_basename(filename)
            if run_basename is None:
                continue
            run_basename = str(run_basename).strip()
            if run_basename == '':
                continue
            logical_key = (tuple(relative_parts), run_basename)
            logical_run_records.setdefault(logical_key, []).append({
                'basename': filename,
                'run_basename': run_basename,
                'mate': mate,
                'path': path_fastq,
                'relative_parts': list(relative_parts),
                'scientific_name': scientific_name,
            })
    grouped_files = {}
    assigned_run_ids = assign_unique_run_ids(logical_run_records)
    for logical_key, records in logical_run_records.items():
        run_id = assigned_run_ids[logical_key]
        for record in records:
            grouped_files.setdefault(run_id, []).append({
                'basename': record['basename'],
                'run': run_id,
                'mate': record['mate'],
                'path': record['path'],
                'relative_parts': record['relative_parts'],
                'scientific_name': record['scientific_name'],
            })
    return grouped_files

def write_fastq_head(path_fastq, path_out, max_lines=4000):
    with open_fastq_text(path_fastq) as fin, open(path_out, 'wt') as fout:
        for i, line in enumerate(fin):
            if i >= max_lines:
                break
            fout.write(line)

def resolve_seqkit_exe_for_integrate(seqkit_exe):
    if seqkit_exe is None:
        return 'seqkit'
    seqkit_exe = str(seqkit_exe).strip()
    if seqkit_exe == '':
        return 'seqkit'
    return seqkit_exe

def resolve_seqkit_threads_for_integrate(seqkit_threads):
    try:
        seqkit_threads = int(seqkit_threads)
    except (TypeError, ValueError):
        seqkit_threads = 1
    return max(1, seqkit_threads)

def is_decompressed_fastq_path(path_fastq):
    path_fastq_lower = str(path_fastq).lower()
    if path_fastq_lower.endswith(('.fq', '.fastq')):
        return True
    if path_fastq_lower.endswith(('.fq.gz', '.fastq.gz')):
        return False
    return None

def parse_seqkit_stats_tsv_output(stdout_txt, path_fastq):
    rows = parse_seqkit_stats_rows(stdout_txt=stdout_txt)
    if len(rows) != 1:
        raise RuntimeError(
            'seqkit stats returned {} rows for single FASTQ input {}.'.format(len(rows), path_fastq)
        )
    return parse_seqkit_stats_row(row=rows[0], path_fastq=path_fastq)

def scan_fastq_stats_with_seqkit_batch(path_fastq_paths, seqkit_exe='seqkit', seqkit_threads=1):
    if path_fastq_paths is None:
        return {}
    deduplicated_paths = []
    seen = set()
    for path_fastq in path_fastq_paths:
        normalized_path = os.path.abspath(path_fastq)
        if normalized_path in seen:
            continue
        seen.add(normalized_path)
        deduplicated_paths.append(path_fastq)
    if len(deduplicated_paths) == 0:
        return {}
    seqkit_exe = resolve_seqkit_exe_for_integrate(seqkit_exe)
    seqkit_threads = resolve_seqkit_threads_for_integrate(seqkit_threads)
    command = [
        seqkit_exe, 'stats',
        '-T',
        '-a',
        '-j', str(seqkit_threads),
    ] + list(deduplicated_paths)
    out, stdout_txt, stderr_txt = run_logged_command(
        command=command,
        runner=subprocess.run,
        print_command=False,
        print_output=False,
        not_found_label='seqkit',
    )
    if out.returncode != 0:
        raise RuntimeError(
            'seqkit stats failed (exit code {}) for {} input file(s): {}'.format(
                out.returncode,
                len(deduplicated_paths),
                stderr_txt.strip(),
            )
        )
    return map_seqkit_stats_rows(
        stdout_txt=stdout_txt,
        requested_paths=deduplicated_paths,
        row_parser=parse_seqkit_stats_row,
        missing_message='seqkit stats output did not include all requested FASTQ paths: {}',
    )

def scan_fastq_stats_with_seqkit(path_fastq, seqkit_exe='seqkit', seqkit_threads=1):
    stats_by_path = scan_fastq_stats_with_seqkit_batch(
        path_fastq_paths=[path_fastq],
        seqkit_exe=seqkit_exe,
        seqkit_threads=seqkit_threads,
    )
    return stats_by_path[os.path.abspath(path_fastq)]

def estimate_total_spots_from_line_count(path_fastq):
    total_lines = count_fastq_lines(path_fastq)
    if (total_lines % 4) != 0:
        raise ValueError('Malformed FASTQ (line count not divisible by 4): {}'.format(path_fastq))
    return int(total_lines / 4)

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


def scan_run_fastq_stats(run_spec, accurate_size, seqkit_exe='seqkit', seqkit_threads=1, seqkit_stats_cache=None):
    run_id = run_spec['run']
    lib_layout = run_spec['lib_layout']
    read1_path = run_spec['read1_path']
    read2_path = run_spec['read2_path']
    print("Found {} file(s) for ID {}. Lib-layout: {}".format(run_spec['num_files'], run_id, lib_layout), flush=True)
    print("Getting sequence statistics.", flush=True)
    is_decompressed = is_decompressed_fastq_path(read1_path)
    if is_decompressed is None:
        warnings.warn("{} is not a fastq file. Skipping.".format(read1_path))
        return None

    quick_gzip_mode = (not accurate_size) and (not is_decompressed)
    if accurate_size or is_decompressed:
        print('--accurate_size set to yes. Running accurate sequence scan with seqkit stats.')
        read1_stats = None
        if isinstance(seqkit_stats_cache, dict):
            read1_stats = seqkit_stats_cache.get(os.path.abspath(read1_path), None)
        if read1_stats is not None:
            total_spots, avg_len_read1 = read1_stats
        else:
            try:
                total_spots, avg_len_read1 = scan_fastq_stats_with_seqkit(
                    path_fastq=read1_path,
                    seqkit_exe=seqkit_exe,
                    seqkit_threads=seqkit_threads,
                )
            except Exception as exc:
                warnings.warn('seqkit stats failed for {}. Falling back to Python FASTQ parser. {}'.format(read1_path, exc))
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
            read2_stats = None
            if isinstance(seqkit_stats_cache, dict):
                read2_stats = seqkit_stats_cache.get(os.path.abspath(read2_path), None)
            if read2_stats is not None:
                read2_spots, avg_len_read2 = read2_stats
            else:
                try:
                    read2_spots, avg_len_read2 = scan_fastq_stats_with_seqkit(
                        path_fastq=read2_path,
                        seqkit_exe=seqkit_exe,
                        seqkit_threads=seqkit_threads,
                    )
                except Exception as exc:
                    warnings.warn('seqkit stats failed for {}. Falling back to Python FASTQ parser. {}'.format(read2_path, exc))
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
        scientific_names = sorted({
            normalize_scientific_name_for_lookup(record.get('scientific_name', ''))
            for record in fastq_records
            if normalize_scientific_name_for_lookup(record.get('scientific_name', '')) != ''
        })
        scientific_name = scientific_names[0] if len(scientific_names) > 0 else ''
        run_specs.append({
            'run': run_id,
            'num_files': len(fastq_records),
            'lib_layout': lib_layout,
            'read1_path': read1_path,
            'read2_path': read2_path,
            'scientific_name': scientific_name,
        })
    return run_specs

def normalize_text_value(value):
    if pd.isna(value):
        return ''
    text = str(value).strip()
    if text.lower() in ['', 'nan', 'none']:
        return ''
    return re.sub(r'\s+', ' ', text)

def normalize_scientific_name_for_lookup(value):
    normalized = normalize_text_value(value)
    if normalized == '':
        return ''
    return normalize_text_value(normalized.replace('_', ' '))

def is_missing_private_scientific_name(value):
    normalized = normalize_scientific_name_for_lookup(value)
    if normalized == '':
        return True
    return normalized.lower() == PRIVATE_FASTQ_SCIENTIFIC_NAME_PLACEHOLDER.lower()

def is_resolvable_scientific_name(value):
    normalized = normalize_scientific_name_for_lookup(value)
    if normalized == '' or is_missing_private_scientific_name(normalized):
        return False
    tokens = [token for token in normalized.split(' ') if token != '']
    if len(tokens) >= 2:
        return True
    return re.fullmatch(r'[A-Za-z][A-Za-z.-]*', tokens[0]) is not None

def infer_default_private_scientific_name(existing_df):
    if existing_df is None or 'scientific_name' not in existing_df.columns:
        return ''
    unique_names = []
    seen = set()
    for value in existing_df['scientific_name'].tolist():
        normalized = normalize_scientific_name_for_lookup(value)
        if not is_resolvable_scientific_name(normalized):
            continue
        key = normalize_scientific_name_for_lookup(normalized).lower()
        if key in seen:
            continue
        seen.add(key)
        unique_names.append(normalized)
    if len(unique_names) == 1:
        return unique_names[0]
    return ''

def build_lineage_taxid_df(taxid_series, ncbi):
    rank_cols = ['taxid_' + rank for rank in STANDARD_TAXONOMIC_RANKS]
    unique_taxids = [int(taxid) for taxid in taxid_series.dropna().unique().tolist()]
    if len(unique_taxids) == 0:
        return pd.DataFrame({'taxid': pd.Series(dtype='Int64'), **{col: pd.Series(dtype='Int64') for col in rank_cols}})

    row_map = {
        taxid: dict({'taxid': taxid}, **{col: pd.NA for col in rank_cols})
        for taxid in unique_taxids
    }
    lineage_map = {}
    lineage_failures = {}
    missing_taxids = list(unique_taxids)
    if hasattr(ncbi, 'get_lineage_translator'):
        try:
            lineage_map = {
                int(taxid): [int(lineage_taxid) for lineage_taxid in lineage]
                for taxid, lineage in ncbi.get_lineage_translator(unique_taxids).items()
            }
            missing_taxids = [taxid for taxid in unique_taxids if taxid not in lineage_map]
        except KeyboardInterrupt:
            raise
        except Exception:
            lineage_map = {}
            missing_taxids = list(unique_taxids)
    for taxid in missing_taxids:
        try:
            lineage_map[taxid] = [int(lineage_taxid) for lineage_taxid in ncbi.get_lineage(taxid)]
        except KeyboardInterrupt:
            raise
        except Exception as exc:
            lineage_failures.setdefault(taxid, exc)

    resolved_lineage_taxids = sorted({
        lineage_taxid
        for lineage in lineage_map.values()
        for lineage_taxid in lineage
    })
    rank_dict = None
    if len(resolved_lineage_taxids) > 0:
        try:
            rank_dict = ncbi.get_rank(resolved_lineage_taxids)
        except KeyboardInterrupt:
            raise
        except Exception:
            rank_dict = None
    if rank_dict is not None:
        for taxid, lineage in lineage_map.items():
            lineage_row = row_map[taxid]
            for lineage_taxid in lineage:
                rank = rank_dict.get(lineage_taxid)
                if rank in STANDARD_TAXONOMIC_RANKS:
                    lineage_row['taxid_' + rank] = int(lineage_taxid)
    else:
        for taxid, lineage in lineage_map.items():
            lineage_row = row_map[taxid]
            try:
                per_taxid_rank_dict = ncbi.get_rank(lineage)
            except KeyboardInterrupt:
                raise
            except Exception as exc:
                lineage_failures.setdefault(taxid, exc)
                continue
            for lineage_taxid, rank in per_taxid_rank_dict.items():
                if rank in STANDARD_TAXONOMIC_RANKS:
                    lineage_row['taxid_' + rank] = int(lineage_taxid)

    if len(lineage_failures) > 0:
        preview = ', '.join(
            ['{} ({})'.format(taxid, exc.__class__.__name__) for taxid, exc in list(lineage_failures.items())[:5]]
        )
        if len(lineage_failures) > 5:
            preview += ', ...'
        warnings.warn(
            'Failed to resolve NCBI lineage for {} taxid(s): {}'.format(
                len(lineage_failures),
                preview,
            )
        )
    lineage_taxid_df = pd.DataFrame([row_map[taxid] for taxid in unique_taxids])
    lineage_taxid_df = lineage_taxid_df.reindex(columns=['taxid'] + rank_cols, fill_value=pd.NA)
    return lineage_taxid_df.astype('Int64')

def ensure_taxonomy_columns(df, args):
    metadata = Metadata.from_DataFrame(df)
    taxid_values = pd.to_numeric(metadata.df['taxid'], errors='coerce') if 'taxid' in metadata.df.columns else pd.Series(dtype='float64')
    if len(taxid_values) == 0:
        metadata.df['taxid'] = pd.Series([pd.NA] * metadata.df.shape[0], dtype='Int64')
    else:
        metadata.df['taxid'] = taxid_values.astype('Int64')

    missing_taxid_mask = metadata.df['taxid'].isna()
    name_to_rows = {}
    for row_idx in metadata.df.index[missing_taxid_mask].tolist():
        normalized_name = normalize_scientific_name_for_lookup(metadata.df.at[row_idx, 'scientific_name'])
        if not is_resolvable_scientific_name(normalized_name):
            continue
        name_to_rows.setdefault(normalized_name, []).append(row_idx)

    ncbi = None
    if len(name_to_rows) > 0:
        ncbi = get_ete_ncbitaxa(args=args)
        unresolved_names = []
        ambiguous_names = []
        try:
            translated = ncbi.get_name_translator(list(name_to_rows.keys()))
        except KeyboardInterrupt:
            raise
        except Exception as exc:
            warnings.warn('Failed to resolve taxid from scientific_name values. {}'.format(exc))
            translated = {}
        for normalized_name, row_indices in name_to_rows.items():
            taxid_candidates = translated.get(normalized_name, [])
            if len(taxid_candidates) == 0:
                unresolved_names.append(normalized_name)
                continue
            chosen_taxid = int(taxid_candidates[0])
            if len(taxid_candidates) > 1:
                ambiguous_names.append('{} -> {}'.format(normalized_name, taxid_candidates))
            metadata.df.loc[row_indices, 'taxid'] = chosen_taxid
        metadata.df['taxid'] = pd.to_numeric(metadata.df['taxid'], errors='coerce').astype('Int64')
        if len(unresolved_names) > 0:
            preview = ', '.join(unresolved_names[:5])
            if len(unresolved_names) > 5:
                preview += ', ...'
            warnings.warn('Could not resolve taxid for {} scientific_name value(s): {}'.format(len(unresolved_names), preview))
        if len(ambiguous_names) > 0:
            preview = ', '.join(ambiguous_names[:5])
            if len(ambiguous_names) > 5:
                preview += ', ...'
            warnings.warn(
                'Multiple taxids matched {} scientific_name value(s); using the first match: {}'.format(
                    len(ambiguous_names),
                    preview,
                )
            )

    if metadata.df['taxid'].notna().any():
        rank_cols = ['taxid_' + rank for rank in STANDARD_TAXONOMIC_RANKS]
        if ncbi is None:
            ncbi = get_ete_ncbitaxa(args=args)
        lineage_taxid_df = build_lineage_taxid_df(metadata.df['taxid'], ncbi)
        metadata.df = metadata.df.drop(columns=rank_cols, errors='ignore')
        metadata.df = metadata.df.merge(lineage_taxid_df, on='taxid', how='left')
    return metadata.df

def finalize_private_fastq_metadata(tmp_metadata, args, existing_df=None):
    metadata = Metadata.from_DataFrame(tmp_metadata)
    if metadata.df.shape[0] == 0:
        return metadata.df
    metadata.df['data_available'] = 'yes'
    default_scientific_name = infer_default_private_scientific_name(existing_df)
    if default_scientific_name != '':
        missing_name_mask = metadata.df['scientific_name'].apply(is_missing_private_scientific_name)
        metadata.df.loc[missing_name_mask, 'scientific_name'] = default_scientific_name
    return ensure_taxonomy_columns(metadata.df, args=args)

def collect_seqkit_target_paths(run_specs, accurate_size):
    target_paths = []
    seen = set()
    for run_spec in run_specs:
        read1_path = run_spec['read1_path']
        is_decompressed = is_decompressed_fastq_path(read1_path)
        if is_decompressed is None:
            continue
        quick_gzip_mode = (not accurate_size) and (not is_decompressed)
        if quick_gzip_mode:
            continue
        for path_fastq in [read1_path, run_spec['read2_path']]:
            if path_fastq == 'unavailable':
                continue
            abs_path = os.path.abspath(path_fastq)
            if abs_path in seen:
                continue
            seen.add(abs_path)
            target_paths.append(path_fastq)
    return target_paths

def build_seqkit_stats_cache_for_runs(run_specs, accurate_size, seqkit_exe='seqkit', seqkit_threads=1):
    target_paths = collect_seqkit_target_paths(run_specs=run_specs, accurate_size=accurate_size)
    if len(target_paths) == 0:
        return {}
    return scan_fastq_stats_with_seqkit_batch(
        path_fastq_paths=target_paths,
        seqkit_exe=seqkit_exe,
        seqkit_threads=seqkit_threads,
    )


def scan_all_run_fastq_stats(run_specs, accurate_size, threads=1, seqkit_exe='seqkit'):
    if is_auto_parallel_option(threads):
        jobs = resolve_detected_cpu_count()
    else:
        jobs = validate_positive_int_option(threads, 'threads')
    seqkit_stats_cache = {}
    try:
        seqkit_stats_cache = build_seqkit_stats_cache_for_runs(
            run_specs=run_specs,
            accurate_size=accurate_size,
            seqkit_exe=seqkit_exe,
            seqkit_threads=max(1, jobs),
        )
        if len(seqkit_stats_cache) > 0:
            print(
                'Collected seqkit stats cache for {:,} FASTQ file(s) with one command.'.format(
                    len(seqkit_stats_cache)
                ),
                flush=True,
            )
    except Exception as exc:
        warnings.warn(
            'Batch seqkit stats scan failed. Falling back to per-run scans. {}'.format(exc)
        )
        seqkit_stats_cache = {}
    if (jobs == 1) or (len(run_specs) <= 1):
        seqkit_threads = max(1, jobs)
        stats_by_run = {}
        for run_spec in run_specs:
            stats = scan_run_fastq_stats(
                run_spec,
                accurate_size,
                seqkit_exe=seqkit_exe,
                seqkit_threads=seqkit_threads,
                seqkit_stats_cache=seqkit_stats_cache,
            )
            if stats is not None:
                stats_by_run[run_spec['run']] = stats
        return stats_by_run

    max_workers = min(jobs, len(run_specs))
    print('Scanning FASTQ stats for {:,} runs with {:,} parallel worker(s).'.format(len(run_specs), max_workers), flush=True)
    seqkit_threads = 1
    run_spec_indices = list(range(len(run_specs)))
    stats_by_spec_idx, failures = run_tasks_with_optional_threads(
        task_items=run_spec_indices,
        task_fn=lambda idx: scan_run_fastq_stats(
            run_specs[idx],
            accurate_size,
            seqkit_exe=seqkit_exe,
            seqkit_threads=seqkit_threads,
            seqkit_stats_cache=seqkit_stats_cache,
        ),
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
        scientific_name = normalize_scientific_name_for_lookup(run_spec.get('scientific_name', ''))
        if scientific_name == '':
            scientific_name = PRIVATE_FASTQ_SCIENTIFIC_NAME_PLACEHOLDER
        rows.append({
            'scientific_name': scientific_name,
            'sample_group': PRIVATE_FASTQ_SAMPLE_GROUP_PLACEHOLDER,
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


def get_fastq_stats(args, existing_df=None):
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
        seqkit_exe=getattr(args, 'seqkit_exe', 'seqkit'),
    )
    rows = build_private_fastq_metadata_rows(run_specs, stats_by_run)
    tmp_metadata = pd.DataFrame.from_records(rows, columns=PRIVATE_FASTQ_METADATA_COLUMNS)
    os.makedirs(metadata_dir, exist_ok=True)
    tmp_metadata = tmp_metadata.sort_values(by='run', axis=0, ascending=True).reset_index(drop=True)
    tmp_metadata = finalize_private_fastq_metadata(tmp_metadata, args=args, existing_df=existing_df)
    atomic_write_dataframe(
        tmp_metadata,
        os.path.join(out_dir, 'metadata_private_fastq.tsv'),
        sep='\t',
        index=False,
    )
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
        tmp_metadata = get_fastq_stats(args, existing_df=metadata.df)
        df = pd.concat([metadata.df, tmp_metadata], ignore_index=True, sort=False)
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
        df = ensure_taxonomy_columns(df, args=args)
        atomic_write_dataframe(
            df,
            os.path.join(args.out_dir, 'metadata', 'metadata_updated_for_private_fastq.tsv'),
            sep='\t',
            index=False,
        )
    else:
        print('Generating a new metadata table.')
        get_fastq_stats(args)
