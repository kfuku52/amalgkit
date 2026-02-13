import pandas as pd

import gzip
import os
import re
import warnings

from amalgkit.sanity import check_getfastq_outputs
from amalgkit.util import *

FASTQ_EXTENSIONS = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')

def open_fastq_text(path_fastq):
    if path_fastq.endswith(('.fq.gz', '.fastq.gz')):
        return gzip.open(path_fastq, 'rt')
    return open(path_fastq, 'rt')

def parse_fastq_basename(basename):
    for ext in FASTQ_EXTENSIONS:
        if basename.endswith(ext):
            stem = basename[:-len(ext)]
            matched = re.match(r'^(.*)_([12])$', stem)
            if matched is not None:
                return matched.group(1), matched.group(2)
            return stem, None
    return None, None

def write_fastq_head(path_fastq, path_out, max_lines=4000):
    with open_fastq_text(path_fastq) as fin, open(path_out, 'wt') as fout:
        for i, line in enumerate(fin):
            if i >= max_lines:
                break
            fout.write(line)

def count_fastq_lines(path_fastq):
    num_lines = 0
    with open_fastq_text(path_fastq) as f:
        for _ in f:
            num_lines += 1
    return num_lines

def estimate_total_spots_from_line_count(path_fastq):
    total_lines = count_fastq_lines(path_fastq)
    if (total_lines % 4) != 0:
        raise ValueError('Malformed FASTQ (line count not divisible by 4): {}'.format(path_fastq))
    return int(total_lines / 4)

def scan_fastq_reads(path_fastq, max_reads=None):
    num_reads = 0
    total_bases = 0
    with open_fastq_text(path_fastq) as f:
        while True:
            line1 = f.readline()
            if line1 == '':
                break
            line2 = f.readline()
            line3 = f.readline()
            line4 = f.readline()
            if any([line == '' for line in [line2, line3, line4]]):
                raise ValueError('Malformed FASTQ (record truncated): {}'.format(path_fastq))
            num_reads += 1
            total_bases += len(line2.strip())
            if (max_reads is not None) and (num_reads >= max_reads):
                break
    avg_len = int(total_bases / num_reads) if num_reads > 0 else 0
    return num_reads, avg_len

def get_fastq_stats(args):
    print("Starting integration of fastq-file metadata...")
    if not os.path.exists(args.fastq_dir):
        raise ValueError("PATH to fastq directory does not exist: {}".format(args.fastq_dir))
    all_files = os.listdir(args.fastq_dir)
    parsed_files = []
    for basename in all_files:
        run_id, mate = parse_fastq_basename(basename)
        if run_id is None:
            continue
        parsed_files.append({
            'basename': basename,
            'run': run_id,
            'mate': mate,
            'path': os.path.join(args.fastq_dir, basename),
        })
    if len(parsed_files) == 0:
        txt = 'No detected fastq files (extensions: {}) in: {}'.format(', '.join(FASTQ_EXTENSIONS), args.fastq_dir)
        raise ValueError(txt)
    grouped_files = {}
    for record in parsed_files:
        grouped_files.setdefault(record['run'], []).append(record)
    column_names = ['scientific_name', 'sample_group', 'run', 'read1_path','read2_path', 'is_sampled',
                    'is_qualified','exclusion', 'lib_layout', 'spot_length', 'total_spots', 'total_bases', 'size', 'private_file']
    tmp_metadata = pd.DataFrame(columns = column_names)
    row = 0
    for run_id in sorted(grouped_files.keys()):
        fastq_records = grouped_files[run_id]
        mate_to_paths = {None: [], '1': [], '2': []}
        for record in fastq_records:
            mate_to_paths[record['mate']].append(record['path'])
        for mate_label in ['1', '2']:
            if len(mate_to_paths[mate_label]) > 1:
                raise ValueError('Multiple FASTQ files detected for {} read {}: {}'.format(
                    run_id, mate_label, mate_to_paths[mate_label]
                ))

        if (len(mate_to_paths['1']) == 1) and (len(mate_to_paths['2']) == 1) and (len(mate_to_paths[None]) == 0):
            lib_layout = 'paired'
            read1_path = mate_to_paths['1'][0]
            read2_path = mate_to_paths['2'][0]
        elif len(fastq_records) == 1:
            lib_layout = 'single'
            read1_path = fastq_records[0]['path']
            read2_path = 'unavailable'
        else:
            raise ValueError('Ambiguous FASTQ file set for {}: {}'.format(
                run_id, [record['basename'] for record in fastq_records]
            ))
        print("Found {} file(s) for ID {}. Lib-layout: {}".format(len(fastq_records), run_id, lib_layout), flush=True)
        print("Getting sequence statistics.", flush=True)
        if read1_path.endswith(('.fq', '.fastq')):
            is_decompressed = True
        elif read1_path.endswith(('.fq.gz', '.fastq.gz')):
            is_decompressed = False
        else:
            warnings.warn("{} is not a fastq file. Skipping.".format(read1_path))
            continue

        if args.accurate_size or is_decompressed:
            print('--accurate_size set to yes. Running accurate sequence scan with Python FASTQ parser.')
            total_spots, avg_len = scan_fastq_reads(path_fastq=read1_path, max_reads=None)
        else:
            print('--accurate_size set to no. Running quick sequence scan (first 1,000 reads) with Python FASTQ parser.')
            sampled_reads, avg_len = scan_fastq_reads(path_fastq=read1_path, max_reads=1000)
            if sampled_reads == 0:
                raise Exception('No reads detected in FASTQ file: {}'.format(read1_path))
            total_spots = estimate_total_spots_from_line_count(read1_path)
        if lib_layout == 'paired':
            read2_spots = estimate_total_spots_from_line_count(read2_path)
            if read2_spots != total_spots:
                raise ValueError(
                    'Mismatched paired-end read counts for {}: read1 has {} reads, read2 has {} reads.'.format(
                        run_id, total_spots, read2_spots
                    )
                )
        tmp_metadata.at[row, 'scientific_name'] = 'Please add in format: Genus species'
        tmp_metadata.at[row,'sample_group'] = 'Please add'
        tmp_metadata.at[row,'run'] = run_id
        tmp_metadata.at[row,'read1_path'] = os.path.abspath(read1_path)
        tmp_metadata.at[row,'read2_path'] = os.path.abspath(read2_path) if read2_path != 'unavailable' else 'unavailable'
        tmp_metadata.at[row,'is_sampled'] = 'yes'
        tmp_metadata.at[row,'is_qualified'] = 'yes'
        tmp_metadata.at[row, 'exclusion'] = 'no'
        tmp_metadata.at[row,'lib_layout'] = lib_layout
        tmp_metadata.at[row,'total_spots'] = total_spots
        tmp_metadata.at[row,'size'] = os.path.getsize(read1_path)
        tmp_metadata.at[row, 'private_file'] = 'yes'
        tmp_metadata.at[row,'spot_length'] = int(avg_len)
        total_bases = total_spots * int(avg_len)
        if lib_layout == 'paired':
            total_bases = total_bases * 2
        tmp_metadata.at[row,'total_bases'] = total_bases
        row += 1
    if not os.path.exists(os.path.join(args.out_dir, 'metadata')):
        os.makedirs(os.path.join(args.out_dir, 'metadata'))
    tmp_metadata = tmp_metadata.sort_values(by='run', axis=0, ascending=True).reset_index(drop=True)
    tmp_metadata.to_csv(os.path.join(args.out_dir, 'metadata_private_fastq.tsv'), sep='\t', index=False)
    return tmp_metadata

def integrate_main(args):
    if args.metadata=='inferred':
        relative_path = os.path.join(args.out_dir, 'metadata', 'metadata.tsv')
        metadata_path = os.path.realpath(relative_path)
    else:
        metadata_path = os.path.realpath(args.metadata)
    if os.path.exists(metadata_path):
        print('Merging existing metadata and private fastq info: {}'.format(metadata_path))
        metadata = load_metadata(args)
        metadata.df.loc[:,'private_file'] = 'no'
        print("scanning for getfastq output")
        sra_ids = metadata.df.loc[:,'run']
        data_available, data_unavailable = check_getfastq_outputs(args, sra_ids, metadata, args.out_dir)
        metadata.df.loc[metadata.df['run'].isin(data_available), 'data_available'] = 'yes'
        metadata.df.loc[metadata.df['run'].isin(data_unavailable), 'data_available'] = 'no'
        tmp_metadata = get_fastq_stats(args)
        df = pd.concat([metadata.df, tmp_metadata])
        df.to_csv(os.path.join(args.out_dir, 'metadata', 'metadata_updated_for_private_fastq.tsv'), sep='\t', index=False)
    else:
        print('Generating a new metadata table.')
        get_fastq_stats(args)
