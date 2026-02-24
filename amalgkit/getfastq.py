from Bio import Entrez
import itertools
import numpy
import gzip
import pandas
import shlex

from amalgkit.util import *

import xml.etree.ElementTree as ET
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from urllib.error import HTTPError, URLError

IDENTICAL_PAIRED_RATIO_THRESHOLD = 0.99
IDENTICAL_PAIRED_CHECKED_READS = 2000
ID_LIST_METADATA_MAX_WORKERS = 4
ID_LIST_VERBOSE_LIMIT = 20
FASTQ_RECORD_COUNT_CHUNK_BYTES = 16 * 1024 * 1024
SORTMERNA_REFERENCE_SPECS = [
    {
        'label': 'SILVA SSU',
        'url': 'https://ftp.arb-silva.de/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz',
        'gz_filename': 'silva_ssu.fasta.gz',
        'fasta_filename': 'silva_ssu.fasta',
    },
    {
        'label': 'SILVA LSU',
        'url': 'https://ftp.arb-silva.de/release_138.2/Exports/SILVA_138.2_LSURef_tax_silva.fasta.gz',
        'gz_filename': 'silva_lsu.fasta.gz',
        'fasta_filename': 'silva_lsu.fasta',
    },
]
GETFASTQ_STATS_COLUMNS = [
    'num_dumped',
    'num_rejected',
    'num_written',
    'num_fastp_in',
    'num_fastp_out',
    'num_rrna_in',
    'num_rrna_out',
    'bp_dumped',
    'bp_rejected',
    'bp_written',
    'bp_fastp_in',
    'bp_fastp_out',
    'bp_rrna_in',
    'bp_rrna_out',
    'bp_discarded',
    'fastp_duplication_rate',
    'fastp_insert_size_peak',
]

def get_or_detect_intermediate_extension(sra_stat, work_dir, files=None, file_state=None):
    cached_ext = sra_stat.get('current_ext', None)
    if cached_ext is not None:
        return cached_ext
    if isinstance(file_state, RunFileState):
        files = file_state.files
    detected_ext = get_newest_intermediate_file_extension(sra_stat, work_dir=work_dir, files=files)
    sra_stat['current_ext'] = detected_ext
    return detected_ext

def set_current_intermediate_extension(sra_stat, ext):
    sra_stat['current_ext'] = ext

def append_file_binary(src_path, dst_path, chunk_size=8 * 1024 * 1024):
    if chunk_size is None:
        chunk_size = 8 * 1024 * 1024
    try:
        chunk_size = int(chunk_size)
    except (TypeError, ValueError):
        chunk_size = 8 * 1024 * 1024
    chunk_size = max(1, chunk_size)

    with open(src_path, 'rb') as src_handle, open(dst_path, 'ab') as dst_handle:
        src_fd = src_handle.fileno()
        dst_fd = dst_handle.fileno()
        offset = 0
        if hasattr(os, 'sendfile'):
            try:
                src_size = os.fstat(src_fd).st_size
                while offset < src_size:
                    sent = os.sendfile(dst_fd, src_fd, offset, src_size - offset)
                    if sent is None:
                        break
                    sent = int(sent)
                    if sent <= 0:
                        break
                    offset += sent
                if offset >= src_size:
                    return
            except OSError:
                pass
        if offset > 0:
            src_handle.seek(offset)
        shutil.copyfileobj(src_handle, dst_handle, length=chunk_size)

def list_run_dir_files(work_dir):
    try:
        with os.scandir(work_dir) as entries:
            return {
                entry.name
                for entry in entries
                if entry.is_file()
            }
    except FileNotFoundError:
        return set()

class RunFileState:
    def __init__(self, work_dir, files=None):
        self.work_dir = work_dir
        if files is None:
            self.files = list_run_dir_files(work_dir)
        else:
            self.files = set(files)

    def has(self, filename):
        return filename in self.files

    def path(self, filename):
        return os.path.join(self.work_dir, filename)

    def add(self, filename):
        self.files.add(filename)

    def discard(self, filename):
        self.files.discard(filename)

    def to_set(self):
        return set(self.files)

def _resolve_run_file_state(work_dir, files=None, file_state=None):
    if isinstance(file_state, RunFileState):
        return file_state
    return RunFileState(work_dir=work_dir, files=files)

def getfastq_search_term(ncbi_id, additional_search_term=None):
    # https://www.ncbi.nlm.nih.gov/books/NBK49540/
    if additional_search_term is None:
        search_term = ncbi_id
    else:
        search_term = ncbi_id + ' AND ' + additional_search_term
    return search_term

def _read_sra_id_list(path_id_list):
    with open(path_id_list) as fin:
        return [
            line.strip()
            for line in fin
            if (line.strip() != '') and (not line.lstrip().startswith('#'))
        ]

def _fetch_single_sra_metadata_frame(args, sra_id, log_details=True):
    search_term = getfastq_search_term(sra_id, args.entrez_additional_search_term)
    if log_details:
        print('Entrez search term:', search_term)
        xml_root = getfastq_getxml(search_term)
    else:
        xml_root = getfastq_getxml(search_term, verbose=False)
    metadata_dict_tmp = Metadata.from_xml(xml_root)
    if metadata_dict_tmp.df.shape[0] == 0:
        if log_details:
            print('No associated SRA. Skipping {}'.format(sra_id))
        return None
    if log_details:
        print('Filtering SRA entry with --layout:', args.layout)
    layout = get_layout(args, metadata_dict_tmp)
    layout_series = metadata_dict_tmp.df['lib_layout'].fillna('').astype(str).str.strip().str.lower()
    metadata_dict_tmp.df = metadata_dict_tmp.df.loc[(layout_series == layout), :]
    if args.sci_name is not None:
        if log_details:
            print('Filtering SRA entry with --sci_name:', args.sci_name)
        sci_series = metadata_dict_tmp.df['scientific_name'].fillna('').astype(str).str.strip()
        metadata_dict_tmp.df = metadata_dict_tmp.df.loc[(sci_series == str(args.sci_name).strip()), :]
    return metadata_dict_tmp.df

def _fetch_id_list_metadata_frames(args, sra_id_list):
    if len(sra_id_list) == 0:
        return []
    unique_sra_ids = list(dict.fromkeys(sra_id_list))
    requested_workers = getattr(args, 'internal_jobs', 'auto')
    requested_threads = getattr(args, 'threads', 'auto')
    internal_cpu_budget = getattr(args, 'internal_cpu_budget', 'auto')
    effective_workers, _ = resolve_worker_allocation(
        requested_workers=requested_workers,
        requested_threads=requested_threads,
        internal_cpu_budget=internal_cpu_budget,
        worker_option_name='internal_jobs',
        context='getfastq metadata fetch:',
    )
    max_workers = min(ID_LIST_METADATA_MAX_WORKERS, effective_workers, len(unique_sra_ids))
    log_details = len(unique_sra_ids) <= ID_LIST_VERBOSE_LIMIT
    if max_workers > 1:
        print(
            'Fetching SRA metadata for {:,} unique IDs ({:,} requested) with {:,} worker(s) (cap {}).'.format(
                len(unique_sra_ids),
                len(sra_id_list),
                max_workers,
                ID_LIST_METADATA_MAX_WORKERS,
            ),
            flush=True,
        )
    if not log_details:
        print(
            'Per-ID metadata retrieval logs suppressed for {:,} IDs (limit {}).'.format(
                len(unique_sra_ids),
                ID_LIST_VERBOSE_LIMIT,
            ),
            flush=True,
        )
    frames_by_sra, failures = run_tasks_with_optional_threads(
        task_items=unique_sra_ids,
        task_fn=lambda sra_id: _fetch_single_sra_metadata_frame(args, sra_id, log_details=log_details),
        max_workers=max_workers,
    )
    for sra_id, exc in failures:
        print('Failed to retrieve metadata for {}. Skipping. {}'.format(sra_id, exc), flush=True)
    metadata_frames = []
    for sra_id in sra_id_list:
        metadata_frame = frames_by_sra.get(sra_id)
        if metadata_frame is None:
            continue
        metadata_frames.append(metadata_frame)
    return metadata_frames

def getfastq_getxml(search_term, retmax=1000, verbose=True):
    def merge_xml_chunk(root, chunk):
        # Merge package-set chunks directly to avoid nested container nodes.
        if (chunk.tag == root.tag) and root.tag.endswith('_SET'):
            root.extend(list(chunk))
            return root
        if root.tag.endswith('_SET'):
            root.append(chunk)
            return root
        # If roots are non-container records, wrap them to preserve all entries.
        container_tag = root.tag + '_SET'
        wrapped = ET.Element(container_tag)
        wrapped.append(root)
        if (chunk.tag == root.tag) and (not chunk.tag.endswith('_SET')):
            wrapped.append(chunk)
        elif chunk.tag == container_tag:
            wrapped.extend(list(chunk))
        else:
            wrapped.append(chunk)
        return wrapped

    def fetch_xml_chunk(record_ids, start, end, retmax, max_retry=10):
        for _ in range(max_retry):
            try:
                handle = Entrez.efetch(
                    db="sra",
                    id=record_ids[start:end],
                    rettype="full",
                    retmode="xml",
                    retmax=retmax,
                )
            except (HTTPError, URLError) as e:
                if verbose:
                    print(e, '- Trying Entrez.efetch() again...')
                continue
            try:
                return ET.parse(handle).getroot()
            except ET.ParseError:
                if verbose:
                    print('XML may be truncated. Retrying...', flush=True)
                continue
        raise RuntimeError(
            'Failed to parse Entrez XML chunk after {} retries (records {}-{}).'.format(
                max_retry,
                start,
                end - 1,
            )
        )

    entrez_db = 'sra'
    try:
        sra_handle = Entrez.esearch(db=entrez_db, term=search_term, retmax=10000000)
    except (HTTPError, URLError) as e:
        print(e, '- Trying Entrez.esearch() again...')
        sra_handle = Entrez.esearch(db=entrez_db, term=search_term, retmax=10000000)
    sra_record = Entrez.read(sra_handle)
    record_ids = sra_record["IdList"]
    num_record = len(record_ids)
    if verbose:
        print('Number of SRA records:', num_record)
    if num_record == 0:
        return ET.Element('EXPERIMENT_PACKAGE_SET')
    root = None
    for start in range(0, num_record, retmax):
        end = min(start + retmax, num_record)
        if verbose:
            print('processing SRA records:', start, '-', end - 1, flush=True)
        chunk = fetch_xml_chunk(record_ids=record_ids, start=start, end=end, retmax=retmax, max_retry=10)
        if root is None:
            root = chunk
        else:
            root = merge_xml_chunk(root, chunk)
    error_node = root.find('.//Error')
    if error_node is not None:
        error_text = ''.join(error_node.itertext()).strip()
        if error_text != '':
            print(error_text)
        raise RuntimeError('Error found in Entrez XML response. Search term: ' + search_term)
    return root

def get_range(sra_stat, offset, total_sra_bp, max_bp):
    if (total_sra_bp <= max_bp):
        start = 1
        end = sra_stat['total_spot']
    else:
        if (sra_stat['total_spot'] > (sra_stat['num_read_per_sra'] + offset)):
            start = offset
            end = offset + sra_stat['num_read_per_sra']
        elif (sra_stat['total_spot'] > sra_stat['num_read_per_sra']):
            start = sra_stat['total_spot'] - sra_stat['num_read_per_sra']
            end = sra_stat['total_spot']
        elif (sra_stat['total_spot'] <= sra_stat['num_read_per_sra']):
            start = 1
            end = sra_stat['total_spot']
    return start, end

def detect_concat_input_files(output_files, run_ids, inext):
    expected_names = set()
    for run_id in run_ids:
        expected_names.add(run_id + inext)
        expected_names.add(run_id + '_1' + inext)
        expected_names.add(run_id + '_2' + inext)
    return sorted([
        f for f in output_files
        if f in expected_names
    ])


def get_concat_output_basename(args, infile, paired=False):
    if args.id is not None:
        id_prefix = str(args.id)
        if paired:
            return id_prefix + re.sub('.*(_[1-2])', r'\g<1>', infile)
        # Avoid duplicated prefixes when --id is the same run identifier.
        if infile.startswith(id_prefix + '.') or infile.startswith(id_prefix + '_'):
            return infile
        return id_prefix + infile
    if args.id_list is not None:
        if paired:
            return os.path.basename(args.id_list) + re.sub('.*(_[1-2])', r'\g<1>', infile)
        return os.path.basename(args.id_list) + infile
    return infile


def _validate_rename_source_and_destination(infile_path, outfile_path):
    if not os.path.exists(infile_path):
        raise FileNotFoundError('Concatenation input file not found for rename shortcut: {}'.format(infile_path))
    if not os.path.isfile(infile_path):
        raise IsADirectoryError('Concatenation input path exists but is not a file: {}'.format(infile_path))
    if os.path.exists(outfile_path) and (not os.path.isfile(outfile_path)):
        raise IsADirectoryError('Concatenation output path exists but is not a file: {}'.format(outfile_path))


def _has_single_complete_pair(infiles, inext):
    suffix1 = '_1' + inext
    suffix2 = '_2' + inext
    pair1 = [infile for infile in infiles if infile.endswith(suffix1)]
    pair2 = [infile for infile in infiles if infile.endswith(suffix2)]
    if (len(pair1) != 1) or (len(pair2) != 1):
        return False
    run1 = pair1[0][:-len(suffix1)]
    run2 = pair2[0][:-len(suffix2)]
    return run1 == run2


def maybe_rename_without_concat(layout, infiles, inext, args, output_dir, run_ids):
    num_inext_files = len(infiles)
    is_single_run = len(run_ids) == 1
    if (layout == 'single') and is_single_run and (num_inext_files == 1):
        print('Only 1', inext, 'file was detected. No concatenation will happen.', flush=True)
        infile = infiles[0]
        outfile = get_concat_output_basename(args=args, infile=infile, paired=False)
        infile_path = os.path.join(output_dir, infile)
        outfile_path = os.path.join(output_dir, outfile)
        _validate_rename_source_and_destination(infile_path=infile_path, outfile_path=outfile_path)
        if infile != outfile:
            print('Replacing ID in the output file name:', infile, outfile)
            os.replace(infile_path, outfile_path)
        return True
    if (layout == 'paired') and is_single_run and (num_inext_files == 2) and _has_single_complete_pair(infiles, inext):
        print('Only 1 pair of', inext, 'files were detected. No concatenation will happen.', flush=True)
        for infile in infiles:
            outfile = get_concat_output_basename(args=args, infile=infile, paired=True)
            infile_path = os.path.join(output_dir, infile)
            outfile_path = os.path.join(output_dir, outfile)
            _validate_rename_source_and_destination(infile_path=infile_path, outfile_path=outfile_path)
            if infile != outfile:
                print('Replacing ID in the output file name:', infile, outfile)
                os.replace(infile_path, outfile_path)
        return True
    return False


def resolve_concat_subexts(layout):
    if layout == 'single':
        return ['']
    if layout == 'paired':
        return ['_1', '_2']
    return []


def build_concat_output_path(args, output_dir, subext, outext):
    if args.id is not None:
        return os.path.join(output_dir, args.id + subext + outext)
    return os.path.join(output_dir, os.path.basename(args.id_list) + subext + outext)


def concatenate_files_with_system_cat(infile_paths, outfile_path):
    cat_exe = shutil.which('cat')
    if cat_exe is None:
        return False
    with open(outfile_path, 'wb') as out_handle:
        cat_out = subprocess.run([cat_exe] + infile_paths, stdout=out_handle, stderr=subprocess.PIPE)
    if cat_out.returncode == 0:
        return True
    if os.path.exists(outfile_path):
        os.remove(outfile_path)
    sys.stderr.write('System concat with cat failed for {}. Falling back to Python concat.\n'.format(outfile_path))
    return False


def concat_fastq_files_for_subext(run_ids, subext, inext, output_dir, outfile_path):
    infiles = [run_id + subext + inext for run_id in run_ids]
    infile_paths = []
    if os.path.exists(outfile_path):
        if not os.path.isfile(outfile_path):
            raise IsADirectoryError('Concatenation output path exists but is not a file: {}'.format(outfile_path))
        os.remove(outfile_path)
    for infile in infiles:
        infile_path = os.path.join(output_dir, infile)
        if os.path.exists(infile_path) and (not os.path.isfile(infile_path)):
            raise IsADirectoryError('Concatenation input path exists but is not a file: {}'.format(infile_path))
        if not os.path.exists(infile_path):
            raise FileNotFoundError('Dumped fastq not found: ' + infile_path)
        print('Concatenated file:', infile_path, flush=True)
        infile_paths.append(infile_path)
    if not concatenate_files_with_system_cat(infile_paths, outfile_path):
        for infile_path in infile_paths:
            append_file_binary(infile_path, outfile_path)
    print('')


def cleanup_concat_tmp_files(args, metadata, g, output_dir, run_ids, output_files):
    if not args.remove_tmp:
        return
    for sra_id in run_ids:
        sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
        ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir, files=output_files)
        remove_intermediate_files(sra_stat, ext=ext, work_dir=output_dir)


def concat_fastq(args, metadata, output_dir, g):
    layout = get_layout(args, metadata)
    inext = '.amalgkit.fastq.gz'
    run_ids = collect_valid_run_ids(metadata.df.loc[:, 'run'].tolist(), unique=True)
    if len(run_ids) == 0:
        raise ValueError('No valid Run IDs were found for concatenation.')
    output_files = list_run_dir_files(output_dir)
    infiles = detect_concat_input_files(output_files=output_files, run_ids=run_ids, inext=inext)
    if maybe_rename_without_concat(
        layout=layout,
        infiles=infiles,
        inext=inext,
        args=args,
        output_dir=output_dir,
        run_ids=run_ids,
    ):
        return None
    print('Concatenating files with the extension:', inext)
    outext = '.amalgkit.fastq.gz'
    subexts = resolve_concat_subexts(layout)
    try:
        requested_workers = int(getattr(args, 'threads', 1))
    except (TypeError, ValueError):
        requested_workers = 1
    requested_workers = max(1, requested_workers)
    worker_count = min(len(subexts), requested_workers)
    if worker_count <= 1:
        for subext in subexts:
            outfile_path = build_concat_output_path(args=args, output_dir=output_dir, subext=subext, outext=outext)
            concat_fastq_files_for_subext(
                run_ids=run_ids,
                subext=subext,
                inext=inext,
                output_dir=output_dir,
                outfile_path=outfile_path,
            )
    else:
        print(
            'Concatenating {:,} read groups with {:,} parallel worker(s).'.format(
                len(subexts),
                worker_count,
            ),
            flush=True,
        )
        def concat_single_subext(subext):
            outfile_path = build_concat_output_path(args=args, output_dir=output_dir, subext=subext, outext=outext)
            concat_fastq_files_for_subext(
                run_ids=run_ids,
                subext=subext,
                inext=inext,
                output_dir=output_dir,
                outfile_path=outfile_path,
            )
            return subext
        _, failures = run_tasks_with_optional_threads(
            task_items=subexts,
            task_fn=concat_single_subext,
            max_workers=worker_count,
        )
        if failures:
            details = '; '.join(['{}: {}'.format(subext, exc) for subext, exc in failures])
            raise RuntimeError('Concatenation failed for {}/{} read group(s). {}'.format(
                len(failures),
                len(subexts),
                details,
            ))
    cleanup_concat_tmp_files(
        args=args,
        metadata=metadata,
        g=g,
        output_dir=output_dir,
        run_ids=run_ids,
        output_files=output_files,
    )
    return None

def remove_sra_files(metadata, amalgkit_out_dir):
    def is_sra_artifact_name(entry_name, sra_id):
        sra_base = sra_id + '.sra'
        return (entry_name == sra_base) or entry_name.startswith(sra_base + '.')

    print('Starting SRA file removal.', flush=True)
    getfastq_root = os.path.join(os.path.realpath(amalgkit_out_dir), 'getfastq')
    for sra_id in collect_valid_run_ids(metadata.df['run'].tolist(), unique=True):
        sra_dir = os.path.join(getfastq_root, sra_id)
        sra_pattern = os.path.join(sra_dir, sra_id + '.sra*')
        path_downloaded_sras = []
        try:
            with os.scandir(sra_dir) as entries:
                for entry in entries:
                    if (not entry.is_file()) or (not is_sra_artifact_name(entry.name, sra_id)):
                        continue
                    path_downloaded_sras.append(entry.path)
        except (FileNotFoundError, NotADirectoryError):
            path_downloaded_sras = []
        if len(path_downloaded_sras) > 0:
            for path_downloaded_sra in sorted(path_downloaded_sras):
                print('Deleting SRA file: {}'.format(path_downloaded_sra))
                os.remove(path_downloaded_sra)
        else:
            print('SRA file not found. Pattern searched: {}'.format(sra_pattern))
    print('')

def get_layout(args, metadata):
    layout_arg = str(args.layout).strip().lower()
    if layout_arg == 'auto':
        if 'lib_layout' not in metadata.df.columns:
            raise ValueError('Column "lib_layout" is required in metadata when --layout auto.')
        normalized_layouts = (
            metadata.df['lib_layout']
            .fillna('')
            .astype(str)
            .str.strip()
            .str.lower()
        )
        layouts = sorted(set([layout for layout in normalized_layouts.tolist() if layout in ['single', 'paired']]))
        if len(layouts) == 0:
            raise ValueError('No valid lib_layout value was found in metadata. Expected "single" or "paired".')
        if len(layouts) != 1:
            print('Detected multiple layouts in the metadata:', layouts)
        layout = 'paired' if 'paired' in layouts else 'single'
    else:
        if layout_arg not in ['single', 'paired']:
            raise ValueError('--layout must be one of: auto, single, paired.')
        layout = layout_arg
    return layout

def remove_old_intermediate_files(sra_id, work_dir, files=None):
    if files is None:
        files = list_run_dir_files(work_dir)
    files = [
        f for f in find_run_prefixed_entries(files, sra_id)
        if (not f.endswith('.sra')) and (os.path.isfile(os.path.join(work_dir, f)))
    ]
    for f in files:
        f_path = os.path.join(work_dir, f)
        print('Deleting old intermediate file:', f_path)
        os.remove(f_path)


def remove_intermediate_files(sra_stat, ext, work_dir):
    file_paths = list()
    if sra_stat['layout'] == 'single':
        file_paths.append(os.path.join(work_dir, sra_stat['sra_id'] + ext))
    elif sra_stat['layout'] == 'paired':
        for i in [1, 2]:
            file_paths.append(os.path.join(work_dir, sra_stat['sra_id'] + '_' + str(i) + ext))
    for file_path in file_paths:
        if os.path.exists(file_path):
            if not os.path.isfile(file_path):
                raise IsADirectoryError('Intermediate path exists but is not a file: {}'.format(file_path))
            print('Deleting intermediate file:', file_path)
            os.remove(file_path)
        else:
            print('Tried to delete but file not found:', file_path)

def normalize_url_for_urllib(source_name, source_url, gcp_project=''):
    source_url = str(source_url).strip()
    parsed = urllib.parse.urlparse(source_url)
    scheme = parsed.scheme.lower()
    if (source_name == 'GCP') and (scheme == 'gs'):
        bucket = parsed.netloc
        object_path = parsed.path.lstrip('/')
        if (bucket == '') or (object_path == ''):
            return source_url
        normalized = 'https://storage.googleapis.com/{}/{}'.format(bucket, object_path)
        project = str(gcp_project).strip()
        if project != '':
            query_pairs = urllib.parse.parse_qsl(parsed.query, keep_blank_values=True)
            if not any((k == 'userProject') for k, _ in query_pairs):
                query_pairs.append(('userProject', project))
            query = urllib.parse.urlencode(query_pairs)
            if query != '':
                normalized = normalized + '?' + query
        return normalized
    return source_url

def resolve_sra_download_sources(metadata, sra_stat, args):
    sra_sources = dict()
    sra_id = sra_stat['sra_id']
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_id))
    if args.aws:
        aws_link = metadata.df.at[ind_sra, 'AWS_Link']
        if aws_link == '':
            sys.stderr.write('AWS_Link is empty and will be skipped.\n')
        else:
            sra_sources['AWS'] = aws_link
    if args.gcp:
        gcp_link = metadata.df.at[ind_sra, 'GCP_Link']
        if gcp_link == '':
            sys.stderr.write('GCP_Link is empty and will be skipped.\n')
        else:
            sra_sources['GCP'] = gcp_link
    if args.ncbi:
        ncbi_link = metadata.df.at[ind_sra, 'NCBI_Link']
        if ncbi_link == '':
            sys.stderr.write('NCBI_Link is empty and will be skipped.\n')
        else:
            sra_sources['NCBI'] = ncbi_link
    return sra_sources

def normalize_sra_download_method(args):
    method_raw = str(getattr(args, 'sra_download_method', 'auto')).strip().lower()
    if method_raw in ['auto', 'urllib', 'curl']:
        return method_raw
    sys.stderr.write('Unknown --sra_download_method "{}". Falling back to "auto".\n'.format(method_raw))
    return 'auto'

def download_with_curl(source_url, path_downloaded_sra, args, sra_source_name):
    curl_exe = shutil.which('curl')
    if curl_exe is None:
        return False
    tmp_path = path_downloaded_sra + '.curltmp.{}'.format(time.time_ns())
    command = [
        curl_exe,
        '-L',
        '--fail',
        '--retry', '3',
        '--retry-delay', '2',
        '--connect-timeout', '20',
        '-o', tmp_path,
        source_url,
    ]
    print('Command:', ' '.join(command))
    out = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if should_print_getfastq_command_output(args):
        print('curl stdout:')
        print(out.stdout.decode('utf8', errors='replace'))
        print('curl stderr:')
        print(out.stderr.decode('utf8', errors='replace'))
    if out.returncode != 0:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        sys.stderr.write('curl failed SRA download from {}.\n'.format(sra_source_name))
        return False
    if not os.path.exists(tmp_path):
        sys.stderr.write('curl download did not create output file for {}.\n'.format(sra_source_name))
        return False
    os.replace(tmp_path, path_downloaded_sra)
    print('SRA file was downloaded with curl from {}'.format(sra_source_name), flush=True)
    return True


def download_sra_from_source(sra_id, sra_source_name, source_url_original, path_downloaded_sra, args):
    source_url = normalize_url_for_urllib(
        source_name=sra_source_name,
        source_url=source_url_original,
        gcp_project=getattr(args, 'gcp_project', ''),
    )
    print("Trying to fetch {} from {}: {}".format(sra_id, sra_source_name, source_url_original))
    if source_url != source_url_original:
        print("Converted {} URL for urllib: {}".format(sra_source_name, source_url))
    if source_url == 'nan':
        sys.stderr.write("Skipping. No URL for {}.\n".format(sra_source_name))
        return False
    scheme = urllib.parse.urlparse(source_url).scheme.lower()
    if scheme not in ('http', 'https', 'ftp'):
        msg = 'Skipping {} download source due to unsupported URL scheme for urllib: {}\n'
        sys.stderr.write(msg.format(sra_source_name, scheme if scheme else '(none)'))
        return False
    method = normalize_sra_download_method(args)
    if method in ['auto', 'curl']:
        if download_with_curl(
            source_url=source_url,
            path_downloaded_sra=path_downloaded_sra,
            args=args,
            sra_source_name=sra_source_name,
        ):
            return True
        if method == 'curl':
            sys.stderr.write('Falling back to urllib.request after curl failure.\n')
    try:
        urllib.request.urlretrieve(source_url, path_downloaded_sra)
        if os.path.exists(path_downloaded_sra):
            print('SRA file was downloaded with urllib.request from {}'.format(sra_source_name), flush=True)
            return True
    except urllib.error.HTTPError as e:
        if (sra_source_name == 'GCP') and (e.code == 400):
            details = ''
            try:
                details = e.read().decode('utf-8', errors='ignore')
            except Exception:
                details = ''
            if ('UserProjectMissing' in details) and (str(getattr(args, 'gcp_project', '')).strip() == ''):
                txt = 'GCP requester-pays bucket requires --gcp_project for billing context. '
                txt += 'Continuing with other download sources.\n'
                sys.stderr.write(txt)
        sys.stderr.write("urllib.request failed SRA download from {}.\n".format(sra_source_name))
    except urllib.error.URLError:
        sys.stderr.write("urllib.request failed SRA download from {}.\n".format(sra_source_name))
    return False


def download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
    path_downloaded_sra = os.path.join(work_dir, sra_stat['sra_id'] + '.sra')
    if os.path.exists(path_downloaded_sra):
        if not os.path.isfile(path_downloaded_sra):
            raise IsADirectoryError('SRA path exists but is not a file: {}'.format(path_downloaded_sra))
        print('Previously-downloaded sra file was detected at: {}'.format(path_downloaded_sra))
        if (overwrite):
            print('Removing', path_downloaded_sra)
            print('New sra file will be downloaded.')
            os.remove(path_downloaded_sra)
        else:
            return None
    else:
        print('Previously-downloaded sra file was not detected. New sra file will be downloaded.')

    if (args.aws) or (args.ncbi) or (args.gcp):
        sra_id = sra_stat['sra_id']
        sra_sources = resolve_sra_download_sources(metadata=metadata, sra_stat=sra_stat, args=args)
        if len(sra_sources)==0:
            print('No source URL is available. Check whether --aws, --gcp, and --ncbi are properly set.')
        is_sra_download_completed = False
        for sra_source_name in sra_sources.keys():
            source_url_original = str(sra_sources[sra_source_name])
            is_sra_download_completed = download_sra_from_source(
                sra_id=sra_id,
                sra_source_name=sra_source_name,
                source_url_original=source_url_original,
                path_downloaded_sra=path_downloaded_sra,
                args=args,
            )
            if is_sra_download_completed:
                break
        if not is_sra_download_completed:
            sys.stderr.write("Exhausted all sources of download.\n")
        else:
            if not os.path.exists(path_downloaded_sra):
                raise FileNotFoundError('SRA file download failed: ' + sra_stat['sra_id'])
            return
    err_txt = 'SRA file download failed for {}. Expected PATH: {}. '
    err_txt += 'Cloud URL download sources were exhausted and prefetch fallback is obsolete.'
    if not os.path.exists(path_downloaded_sra):
        raise FileNotFoundError(err_txt.format(sra_stat['sra_id'], path_downloaded_sra))

def check_getfastq_dependency(args):
    obsolete_pfd = getattr(args, 'obsolete_pfd', getattr(args, 'pfd', None))
    obsolete_pfd_exe = getattr(args, 'obsolete_pfd_exe', getattr(args, 'pfd_exe', None))
    obsolete_fastq_dump_exe = getattr(args, 'obsolete_fastq_dump_exe', getattr(args, 'fastq_dump_exe', None))
    if obsolete_pfd is not None:
        sys.stderr.write('--pfd is obsolete and ignored. fasterq-dump is always used for extraction.\n')
    if obsolete_pfd_exe:
        sys.stderr.write('--pfd_exe is obsolete and ignored.\n')
    if obsolete_fastq_dump_exe:
        sys.stderr.write('--fastq_dump_exe is obsolete and ignored.\n')
    obsolete_prefetch_exe = getattr(args, 'obsolete_prefetch_exe', getattr(args, 'prefetch_exe', None))
    if obsolete_prefetch_exe:
        sys.stderr.write('--prefetch_exe is obsolete and ignored.\n')

    def probe_command(command, label):
        try:
            out = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError as exc:
            raise FileNotFoundError('{} executable not found: {}'.format(label, command[0])) from exc
        if out.returncode != 0:
            raise RuntimeError(
                '{} dependency probe failed with exit code {}: {}'.format(
                    label,
                    out.returncode,
                    ' '.join(command),
                )
            )
        return out

    fasterq_dump_exe = getattr(args, 'fasterq_dump_exe', 'fasterq-dump')
    probe_command([fasterq_dump_exe, '-h'], 'fasterq-dump')
    seqkit_exe = resolve_seqkit_exe(args)
    probe_command([seqkit_exe, '--help'], 'seqkit')
    if bool(getattr(args, 'fastp', False)):
        fastp_exe = getattr(args, 'fastp_exe', 'fastp')
        probe_command([fastp_exe, '--help'], 'fastp')
    if is_rrna_filter_enabled(args):
        ensure_sortmerna_download_dir_exists(resolve_sortmerna_download_dir(args))
        sortmerna_exe = getattr(args, 'sortmerna_exe', 'sortmerna')
        probe_command([sortmerna_exe, '--version'], 'sortmerna')
    return None

def remove_sra_path(path_downloaded_sra):
    if not os.path.lexists(path_downloaded_sra):
        return
    print('Removing cached SRA object before retry: {}'.format(path_downloaded_sra))
    if os.path.isdir(path_downloaded_sra) and (not os.path.islink(path_downloaded_sra)):
        shutil.rmtree(path_downloaded_sra)
    else:
        os.remove(path_downloaded_sra)

def should_print_getfastq_command_output(args):
    dump_print = getattr(args, 'dump_print', None)
    if dump_print is not None:
        return dump_print
    legacy_dump_print = getattr(args, 'pfd_print', None)
    if legacy_dump_print is not None:
        return legacy_dump_print
    return True

def resolve_seqkit_exe(args):
    seqkit_exe = getattr(args, 'seqkit_exe', 'seqkit')
    if seqkit_exe is None:
        return 'seqkit'
    seqkit_exe = str(seqkit_exe).strip()
    if seqkit_exe == '':
        return 'seqkit'
    return seqkit_exe

def resolve_seqkit_threads(args):
    try:
        seqkit_threads = int(getattr(args, 'threads', 1))
    except (TypeError, ValueError):
        seqkit_threads = 1
    return max(1, seqkit_threads)

def normalize_seqkit_threads(seqkit_threads):
    try:
        normalized_threads = int(seqkit_threads)
    except (TypeError, ValueError):
        normalized_threads = 1
    return max(1, normalized_threads)

def run_seqkit_seq_command(input_paths, output_path, args, command_label, seqkit_threads=None):
    if len(input_paths) == 0:
        raise ValueError('No input FASTQ path was provided for {}.'.format(command_label))
    seqkit_exe = resolve_seqkit_exe(args)
    if seqkit_threads is None:
        seqkit_threads = resolve_seqkit_threads(args)
    seqkit_threads = normalize_seqkit_threads(seqkit_threads)
    command = [seqkit_exe, 'seq', '-j', str(seqkit_threads), '-w', '0', '-o', output_path] + list(input_paths)
    print('Command:', ' '.join(command))
    out = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if should_print_getfastq_command_output(args):
        print('{} stdout:'.format(command_label))
        print(out.stdout.decode('utf8', errors='replace'))
        print('{} stderr:'.format(command_label))
        print(out.stderr.decode('utf8', errors='replace'))
    if out.returncode != 0:
        raise RuntimeError('{} failed.'.format(command_label))
    return out

def run_seqkit_range_command(input_path, output_path, start, end, args, command_label, seqkit_threads=None):
    if not os.path.exists(input_path):
        raise FileNotFoundError('FASTQ input not found: {}'.format(input_path))
    if not os.path.isfile(input_path):
        raise IsADirectoryError('FASTQ input path exists but is not a file: {}'.format(input_path))
    if os.path.exists(output_path) and (not os.path.isfile(output_path)):
        raise IsADirectoryError('FASTQ output path exists but is not a file: {}'.format(output_path))
    if os.path.exists(output_path):
        os.remove(output_path)
    seqkit_exe = resolve_seqkit_exe(args)
    if seqkit_threads is None:
        seqkit_threads = resolve_seqkit_threads(args)
    seqkit_threads = normalize_seqkit_threads(seqkit_threads)
    command = [
        seqkit_exe, 'range',
        '-j', str(seqkit_threads),
        '-r', '{}:{}'.format(int(start), int(end)),
        '-w', '0',
        '-o', output_path,
        input_path,
    ]
    print('Command:', ' '.join(command))
    out = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if should_print_getfastq_command_output(args):
        print('{} stdout:'.format(command_label))
        print(out.stdout.decode('utf8', errors='replace'))
        print('{} stderr:'.format(command_label))
        print(out.stderr.decode('utf8', errors='replace'))
    if out.returncode != 0:
        raise RuntimeError('{} failed.'.format(command_label))
    return out

def run_seqkit_replace_command(input_path, output_path, suffix, args, command_label, seqkit_threads=None):
    if not os.path.exists(input_path):
        raise FileNotFoundError('FASTQ input not found: {}'.format(input_path))
    if not os.path.isfile(input_path):
        raise IsADirectoryError('FASTQ input path exists but is not a file: {}'.format(input_path))
    if os.path.exists(output_path) and (not os.path.isfile(output_path)):
        raise IsADirectoryError('FASTQ output path exists but is not a file: {}'.format(output_path))
    if os.path.exists(output_path):
        os.remove(output_path)
    seqkit_exe = resolve_seqkit_exe(args)
    if seqkit_threads is None:
        seqkit_threads = resolve_seqkit_threads(args)
    seqkit_threads = normalize_seqkit_threads(seqkit_threads)
    command = [
        seqkit_exe, 'replace',
        '-j', str(seqkit_threads),
        '-w', '0',
        '-p', r'^(\S+).*',
        '-r', '$1' + str(suffix),
        '-o', output_path,
        input_path,
    ]
    print('Command:', ' '.join(command))
    out = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if should_print_getfastq_command_output(args):
        print('{} stdout:'.format(command_label))
        print(out.stdout.decode('utf8', errors='replace'))
        print('{} stderr:'.format(command_label))
        print(out.stderr.decode('utf8', errors='replace'))
    if out.returncode != 0:
        raise RuntimeError('{} failed.'.format(command_label))
    return out

def _parse_seqkit_stats_rows(stdout_txt):
    lines = [line.rstrip('\n') for line in (stdout_txt or '').splitlines() if line.strip() != '']
    if len(lines) < 2:
        raise RuntimeError('seqkit stats output was empty.')
    header_fields = lines[0].split('\t')
    rows = []
    for line in lines[1:]:
        value_fields = line.split('\t')
        row = {}
        for idx, field_name in enumerate(header_fields):
            row[field_name] = value_fields[idx] if idx < len(value_fields) else ''
        rows.append(row)
    if len(rows) == 0:
        raise RuntimeError('seqkit stats output did not include data rows.')
    return rows

def _parse_seqkit_stats_numeric(row, key, path_fastq):
    if key not in row:
        raise RuntimeError('seqkit stats output is missing "{}" for {}'.format(key, path_fastq))
    raw = str(row.get(key, '')).strip().replace(',', '')
    if raw == '':
        raise RuntimeError('seqkit stats output has an empty "{}" value for {}'.format(key, path_fastq))
    try:
        return int(float(raw))
    except (TypeError, ValueError) as exc:
        raise RuntimeError('Failed to parse seqkit stats "{}" for {}'.format(key, path_fastq)) from exc

def _parse_seqkit_stats_row_records_and_bases(row, path_fastq):
    num_records = _parse_seqkit_stats_numeric(row=row, key='num_seqs', path_fastq=path_fastq)
    if ('sum_len' in row) and (str(row.get('sum_len', '')).strip() not in ['', 'NA', 'nan', 'NaN']):
        num_bases = _parse_seqkit_stats_numeric(row=row, key='sum_len', path_fastq=path_fastq)
    else:
        avg_len = _parse_seqkit_stats_numeric(row=row, key='avg_len', path_fastq=path_fastq)
        num_bases = int(num_records * avg_len)
    return num_records, num_bases

def scan_fastq_records_and_bases_with_seqkit_batch(path_fastq_paths, args, seqkit_threads=None):
    if path_fastq_paths is None:
        return {}
    deduplicated_paths = []
    seen = set()
    for path_fastq in path_fastq_paths:
        if path_fastq is None:
            continue
        abs_path = os.path.abspath(path_fastq)
        if abs_path in seen:
            continue
        seen.add(abs_path)
        deduplicated_paths.append(path_fastq)
    if len(deduplicated_paths) == 0:
        return {}
    seqkit_exe = resolve_seqkit_exe(args)
    if seqkit_threads is None:
        seqkit_threads = resolve_seqkit_threads(args)
    seqkit_threads = normalize_seqkit_threads(seqkit_threads)
    command = [
        seqkit_exe, 'stats',
        '-T',
        '-a',
        '-j', str(seqkit_threads),
    ] + list(deduplicated_paths)
    print('Command:', ' '.join(command))
    out = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if should_print_getfastq_command_output(args):
        print('FASTQ stats scan with seqkit stdout:')
        print(out.stdout.decode('utf8', errors='replace'))
        print('FASTQ stats scan with seqkit stderr:')
        print(out.stderr.decode('utf8', errors='replace'))
    if out.returncode != 0:
        stderr_txt = out.stderr.decode('utf8', errors='replace')
        raise RuntimeError('FASTQ stats scan with seqkit failed (exit code {}): {}'.format(out.returncode, stderr_txt))
    rows = _parse_seqkit_stats_rows(out.stdout.decode('utf8', errors='replace'))
    stats_by_path = {}
    fallback_idx = 0
    fallback_paths = [os.path.abspath(path_fastq) for path_fastq in deduplicated_paths]
    for row in rows:
        row_path = None
        for key in ['file', 'path', 'filename']:
            raw_path = str(row.get(key, '')).strip()
            if raw_path != '':
                row_path = os.path.abspath(raw_path)
                break
        if (row_path is None) and (fallback_idx < len(fallback_paths)):
            row_path = fallback_paths[fallback_idx]
            fallback_idx += 1
        if row_path is None:
            continue
        stats_by_path[row_path] = _parse_seqkit_stats_row_records_and_bases(row=row, path_fastq=row_path)
    missing_paths = [path_fastq for path_fastq in fallback_paths if path_fastq not in stats_by_path]
    if len(missing_paths) > 0:
        raise RuntimeError('seqkit stats output did not include all requested FASTQ files: {}'.format(missing_paths))
    return stats_by_path

def compress_fastq_with_seqkit(
    path_fastq,
    path_fastq_gz,
    args,
    command_label='FASTQ compression with seqkit',
    seqkit_threads=None,
):
    if not os.path.exists(path_fastq):
        raise FileNotFoundError('FASTQ input not found: {}'.format(path_fastq))
    if not os.path.isfile(path_fastq):
        raise IsADirectoryError('FASTQ input path exists but is not a file: {}'.format(path_fastq))
    if os.path.exists(path_fastq_gz) and (not os.path.isfile(path_fastq_gz)):
        raise IsADirectoryError('FASTQ output path exists but is not a file: {}'.format(path_fastq_gz))
    if os.path.exists(path_fastq_gz):
        os.remove(path_fastq_gz)
    run_seqkit_seq_command(
        input_paths=[path_fastq],
        output_path=path_fastq_gz,
        args=args,
        command_label=command_label,
        seqkit_threads=seqkit_threads,
    )
    os.remove(path_fastq)

def count_fastq_records(path_fastq):
    num_lines = 0
    open_func = gzip.open if path_fastq.endswith('.gz') else open
    with open_func(path_fastq, 'rb') as f:
        while True:
            chunk = f.read(FASTQ_RECORD_COUNT_CHUNK_BYTES)
            if not chunk:
                break
            num_lines += chunk.count(b'\n')
    if (num_lines % 4) != 0:
        txt = 'FASTQ line count is not divisible by 4 and may be truncated: {}\n'
        sys.stderr.write(txt.format(path_fastq))
    return num_lines // 4

def count_fastq_records_and_bases(path_fastq):
    num_records = 0
    num_bases = 0
    open_func = gzip.open if path_fastq.endswith('.gz') else open
    with open_func(path_fastq, 'rb') as fin:
        while True:
            line1 = fin.readline()
            if line1 == b'':
                break
            line2 = fin.readline()
            line3 = fin.readline()
            line4 = fin.readline()
            if (line2 == b'') or (line3 == b'') or (line4 == b''):
                txt = 'FASTQ record seems truncated while counting bases: {}\n'
                sys.stderr.write(txt.format(path_fastq))
                break
            num_records += 1
            num_bases += len(line2.rstrip(b'\r\n'))
    return num_records, num_bases

def summarize_layout_fastq_records_and_bases(layout, single_path=None, pair1_path=None, pair2_path=None, args=None):
    num_single, bp_single = (0, 0)
    num_pair1, bp_pair1 = (0, 0)
    num_pair2, bp_pair2 = (0, 0)
    if args is not None:
        try:
            paths = [single_path, pair1_path, pair2_path]
            stats_by_path = scan_fastq_records_and_bases_with_seqkit_batch(
                path_fastq_paths=paths,
                args=args,
            )
            if single_path is not None:
                num_single, bp_single = stats_by_path.get(os.path.abspath(single_path), (0, 0))
            if pair1_path is not None:
                num_pair1, bp_pair1 = stats_by_path.get(os.path.abspath(pair1_path), (0, 0))
            if pair2_path is not None:
                num_pair2, bp_pair2 = stats_by_path.get(os.path.abspath(pair2_path), (0, 0))
        except Exception as exc:
            sys.stderr.write('seqkit stats scan failed. Falling back to Python FASTQ parser. {}\n'.format(exc))
            if single_path is not None:
                num_single, bp_single = count_fastq_records_and_bases(single_path)
            if pair1_path is not None:
                num_pair1, bp_pair1 = count_fastq_records_and_bases(pair1_path)
            if pair2_path is not None:
                num_pair2, bp_pair2 = count_fastq_records_and_bases(pair2_path)
    else:
        if single_path is not None:
            num_single, bp_single = count_fastq_records_and_bases(single_path)
        if pair1_path is not None:
            num_pair1, bp_pair1 = count_fastq_records_and_bases(pair1_path)
        if pair2_path is not None:
            num_pair2, bp_pair2 = count_fastq_records_and_bases(pair2_path)
    if layout == 'paired':
        num_spots = max(num_pair1, num_pair2)
        bp_total = bp_pair1 + bp_pair2
    else:
        num_spots = num_single
        bp_total = bp_single
    return num_spots, bp_total

def is_rrna_filter_enabled(args):
    raw = getattr(args, 'rrna_filter', False)
    if isinstance(raw, str):
        return bool(strtobool(raw))
    return bool(raw)

def resolve_filter_order(args):
    raw_order = getattr(args, 'filter_order', 'fastp_first')
    if raw_order is None:
        return 'fastp_first'
    order = str(raw_order).strip().lower()
    if order == '':
        return 'fastp_first'
    if order in ['fastp_first', 'rrna_first']:
        return order
    raise ValueError('--filter_order must be "fastp_first" or "rrna_first".')

def get_filter_execution_order(args):
    fastp_enabled = bool(getattr(args, 'fastp', False))
    rrna_enabled = is_rrna_filter_enabled(args)
    if fastp_enabled and rrna_enabled:
        if resolve_filter_order(args) == 'rrna_first':
            return ['rrna', 'fastp']
        return ['fastp', 'rrna']
    if rrna_enabled:
        return ['rrna']
    if fastp_enabled:
        return ['fastp']
    return []

def resolve_bp_amalgkit_source_column(args):
    filter_order = get_filter_execution_order(args)
    if len(filter_order) == 0:
        return 'bp_written'
    last_filter = filter_order[-1]
    if last_filter == 'fastp':
        return 'bp_fastp_out'
    if last_filter == 'rrna':
        return 'bp_rrna_out'
    raise ValueError('Unknown filter name in execution order: {}'.format(last_filter))

def parse_sortmerna_option_args(sortmerna_option):
    if not sortmerna_option:
        return []
    try:
        return shlex.split(sortmerna_option)
    except ValueError as e:
        raise ValueError('Invalid --sortmerna_option string: {}'.format(sortmerna_option)) from e

def resolve_shared_download_dir(args):
    out_dir = getattr(args, 'out_dir', './')
    return os.path.join(os.path.realpath(out_dir), 'downloads')

def resolve_sortmerna_download_dir(args):
    raw_dir = getattr(args, 'download_dir', 'inferred')
    if raw_dir is None:
        return resolve_shared_download_dir(args)
    normalized = str(raw_dir).strip()
    if normalized.lower() in ['', 'inferred']:
        return resolve_shared_download_dir(args)
    return os.path.realpath(normalized)

def ensure_sortmerna_download_dir_exists(download_dir):
    if os.path.exists(download_dir) and (not os.path.isdir(download_dir)):
        raise NotADirectoryError(
            'SortMeRNA download path exists but is not a directory: {}'.format(download_dir)
        )
    os.makedirs(download_dir, exist_ok=True)
    return download_dir

def download_sortmerna_reference_gz(url, gz_path, label):
    tmp_path = gz_path + '.tmp.{}'.format(time.time_ns())
    if os.path.exists(tmp_path):
        if os.path.isdir(tmp_path):
            raise IsADirectoryError('Temporary SortMeRNA download path is a directory: {}'.format(tmp_path))
        os.remove(tmp_path)
    print('Downloading {} reference: {}'.format(label, url), flush=True)
    try:
        urllib.request.urlretrieve(url, tmp_path)
    except Exception as e:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise RuntimeError('Failed to download {} reference: {}'.format(label, url)) from e
    os.replace(tmp_path, gz_path)

def expand_sortmerna_reference_gz(gz_path, fasta_path):
    tmp_path = fasta_path + '.tmp.{}'.format(time.time_ns())
    if os.path.exists(tmp_path):
        if os.path.isdir(tmp_path):
            raise IsADirectoryError('Temporary SortMeRNA expansion path is a directory: {}'.format(tmp_path))
        os.remove(tmp_path)
    with gzip.open(gz_path, 'rb') as fin, open(tmp_path, 'wb') as fout:
        shutil.copyfileobj(fin, fout)
    os.replace(tmp_path, fasta_path)

def ensure_sortmerna_reference_files_exist(args):
    download_dir = ensure_sortmerna_download_dir_exists(resolve_sortmerna_download_dir(args))
    refs = []
    for ref_spec in SORTMERNA_REFERENCE_SPECS:
        fasta_path = os.path.join(download_dir, ref_spec['fasta_filename'])
        gz_path = os.path.join(download_dir, ref_spec['gz_filename'])
        if os.path.exists(fasta_path):
            if not os.path.isfile(fasta_path):
                raise IsADirectoryError('SortMeRNA reference path exists but is not a file: {}'.format(fasta_path))
            refs.append(fasta_path)
            continue
        if os.path.exists(gz_path) and (not os.path.isfile(gz_path)):
            raise IsADirectoryError('SortMeRNA compressed reference path is not a file: {}'.format(gz_path))
        if not os.path.exists(gz_path):
            download_sortmerna_reference_gz(
                url=ref_spec['url'],
                gz_path=gz_path,
                label=ref_spec['label'],
            )
        print('Expanding SortMeRNA reference: {}'.format(gz_path), flush=True)
        expand_sortmerna_reference_gz(gz_path=gz_path, fasta_path=fasta_path)
        refs.append(fasta_path)
    return refs

def resolve_sortmerna_refs(args):
    return ensure_sortmerna_reference_files_exist(args)

def estimate_num_written_spots_from_fastq(sra_stat, files=None, file_state=None):
    work_dir = sra_stat['getfastq_sra_dir']
    sra_id = sra_stat['sra_id']
    run_file_state = _resolve_run_file_state(work_dir=work_dir, files=files, file_state=file_state)
    def detect_fastq_path(prefix):
        for ext in ['.fastq.gz', '.fastq']:
            filename = prefix + ext
            if run_file_state.has(filename):
                return run_file_state.path(filename)
        return None
    single_path = detect_fastq_path(sra_id)
    pair1_path = detect_fastq_path(sra_id + '_1')
    pair2_path = detect_fastq_path(sra_id + '_2')
    num_single = count_fastq_records(single_path) if single_path is not None else 0
    num_pair1 = count_fastq_records(pair1_path) if pair1_path is not None else 0
    num_pair2 = count_fastq_records(pair2_path) if pair2_path is not None else 0
    return max(num_pair1, num_pair2) + num_single

def trim_fastq_by_spot_range(path_fastq, start, end, args=None, compress_to_gz=False):
    if not os.path.exists(path_fastq):
        return 0
    start = max(1, int(start))
    end = max(0, int(end))
    if compress_to_gz:
        if args is None:
            raise ValueError('args is required when compress_to_gz=True.')
        path_fastq_gz = path_fastq + '.gz'
        run_seqkit_range_command(
            input_path=path_fastq,
            output_path=path_fastq_gz,
            start=start,
            end=end,
            args=args,
            command_label='FASTQ spot-range trim/compress with seqkit',
        )
        kept = count_fastq_records(path_fastq_gz)
        os.remove(path_fastq)
        return kept
    tmp_path = path_fastq + '.trimtmp'
    kept = 0
    spot_index = 0
    with open(path_fastq, 'rb') as fin, open(tmp_path, 'wb') as fout:
        while True:
            line1 = fin.readline()
            if line1 == b'':
                break
            line2 = fin.readline()
            line3 = fin.readline()
            line4 = fin.readline()
            if (line2 == b'') or (line3 == b'') or (line4 == b''):
                raise ValueError('Malformed FASTQ (record truncated): {}'.format(path_fastq))
            spot_index += 1
            if spot_index < start:
                continue
            if spot_index > end:
                break
            fout.writelines((line1, line2, line3, line4))
            kept += 1
    os.replace(tmp_path, path_fastq)
    return kept

def trim_fasterq_output_files(
    sra_stat,
    start,
    end,
    files=None,
    file_state=None,
    args=None,
    compress_to_gz=False,
    return_file_state=False,
):
    work_dir = sra_stat['getfastq_sra_dir']
    sra_id = sra_stat['sra_id']
    run_file_state = _resolve_run_file_state(work_dir=work_dir, files=files, file_state=file_state)
    kept_counts = {'': 0, '_1': 0, '_2': 0}
    for suffix in ['', '_1', '_2']:
        filename = sra_id + suffix + '.fastq'
        if not run_file_state.has(filename):
            continue
        path_fastq = run_file_state.path(filename)
        kept_counts[suffix] = trim_fastq_by_spot_range(
            path_fastq=path_fastq,
            start=start,
            end=end,
            args=args,
            compress_to_gz=compress_to_gz,
        )
        if compress_to_gz:
            run_file_state.discard(filename)
            run_file_state.add(filename + '.gz')
    if return_file_state:
        return kept_counts, run_file_state
    return kept_counts

def calculate_written_spots_from_trim_counts(trimmed_counts):
    if not isinstance(trimmed_counts, dict):
        return None
    try:
        num_single = int(trimmed_counts.get('', 0))
        num_pair1 = int(trimmed_counts.get('_1', 0))
        num_pair2 = int(trimmed_counts.get('_2', 0))
    except (TypeError, ValueError):
        return None
    return max(num_pair1, num_pair2) + num_single

def parse_fasterq_dump_written_spots(stdout_txt, stderr_txt):
    combined = '\n'.join([stdout_txt or '', stderr_txt or ''])
    matched = re.findall(r'^\s*spots\s+written\s*:\s*([0-9][0-9,]*)\s*$', combined, flags=re.IGNORECASE | re.MULTILINE)
    if len(matched) == 0:
        return None
    return int(matched[-1].replace(',', ''))

def parse_fasterq_dump_written_reads(stdout_txt, stderr_txt):
    combined = '\n'.join([stdout_txt or '', stderr_txt or ''])
    matched = re.findall(r'^\s*reads\s+written\s*:\s*([0-9][0-9,]*)\s*$', combined, flags=re.IGNORECASE | re.MULTILINE)
    if len(matched) == 0:
        return None
    return int(matched[-1].replace(',', ''))

def infer_written_spots_from_written_reads(sra_stat, written_reads, run_file_state):
    try:
        written_reads = int(written_reads)
    except (TypeError, ValueError):
        return None
    if written_reads < 0:
        return None
    layout = str(sra_stat.get('layout', '')).strip().lower()
    if layout == 'single':
        return written_reads
    if layout != 'paired':
        return None
    # split-3 paired output can include singleton reads in <run>.fastq, which prevents exact inference.
    singleton_name = sra_stat['sra_id'] + '.fastq'
    if run_file_state.has(singleton_name):
        return None
    if (written_reads % 2) != 0:
        return None
    return int(written_reads / 2)

def is_full_requested_spot_range(sra_stat, start, end):
    total_spot = sra_stat.get('total_spot', None)
    if total_spot is None:
        return False
    return (int(start) <= 1) and (int(end) >= int(total_spot))

def compress_fasterq_output_files(sra_stat, args, files=None, file_state=None, return_file_state=False):
    work_dir = sra_stat['getfastq_sra_dir']
    sra_id = sra_stat['sra_id']
    run_file_state = _resolve_run_file_state(work_dir=work_dir, files=files, file_state=file_state)
    fastq_paths = list()
    fastq_filenames = list()
    for suffix in ['', '_1', '_2']:
        filename = sra_id + suffix + '.fastq'
        if run_file_state.has(filename):
            fastq_paths.append(run_file_state.path(filename))
            fastq_filenames.append(filename)
    if len(fastq_paths) == 0:
        if return_file_state:
            return run_file_state
        return run_file_state.to_set()
    total_seqkit_threads = resolve_seqkit_threads(args)
    worker_count = min(len(fastq_paths), total_seqkit_threads)
    worker_count = max(1, worker_count)
    threads_per_worker = max(1, int(total_seqkit_threads / worker_count))
    if worker_count <= 1:
        for path_fastq in fastq_paths:
            print('Compressing with seqkit: {} -> {}'.format(path_fastq, path_fastq + '.gz'))
            compress_fastq_with_seqkit(
                path_fastq=path_fastq,
                path_fastq_gz=path_fastq + '.gz',
                args=args,
                command_label='FASTQ compression with seqkit',
                seqkit_threads=threads_per_worker,
            )
    else:
        print(
            'Compressing {:,} FASTQ files with seqkit using {:,} parallel worker(s), {:,} thread(s) each.'.format(
                len(fastq_paths),
                worker_count,
                threads_per_worker,
            ),
            flush=True,
        )
        task_indices = list(range(len(fastq_paths)))
        def run_single_compression(task_idx):
            path_fastq = fastq_paths[task_idx]
            print('Compressing with seqkit: {} -> {}'.format(path_fastq, path_fastq + '.gz'))
            compress_fastq_with_seqkit(
                path_fastq=path_fastq,
                path_fastq_gz=path_fastq + '.gz',
                args=args,
                command_label='FASTQ compression with seqkit',
                seqkit_threads=threads_per_worker,
            )
            return path_fastq
        _, failures = run_tasks_with_optional_threads(
            task_items=task_indices,
            task_fn=run_single_compression,
            max_workers=worker_count,
        )
        if failures:
            details = '; '.join(['{}: {}'.format(fastq_paths[idx], exc) for idx, exc in failures])
            raise RuntimeError('FASTQ compression with seqkit failed for {}/{} file(s). {}'.format(
                len(failures),
                len(fastq_paths),
                details,
            ))
    for filename in fastq_filenames:
        run_file_state.discard(filename)
        run_file_state.add(filename + '.gz')
    if return_file_state:
        return run_file_state
    return run_file_state.to_set()

def normalize_requested_spot_range(start, end):
    start = max(1, int(start))
    end = int(end)
    if end < start:
        end = start
    return start, end


def build_fasterq_dump_command(args, sra_stat, path_downloaded_sra, size_check, start=None, end=None):
    fasterq_dump_exe = getattr(args, 'fasterq_dump_exe', 'fasterq-dump')
    command = [
        fasterq_dump_exe,
        '--split-3',
        '--skip-technical',
        '--min-read-len', str(args.min_read_length),
        '--size-check', str(size_check),
        '-e', str(max(1, args.threads)),
        '-O', sra_stat['getfastq_sra_dir'],
        '-t', sra_stat['getfastq_sra_dir'],
    ]
    disk_limit = getattr(args, 'fasterq_disk_limit', None)
    disk_limit_tmp = getattr(args, 'fasterq_disk_limit_tmp', None)
    if disk_limit:
        command.extend(['--disk-limit', str(disk_limit)])
    if disk_limit_tmp:
        command.extend(['--disk-limit-tmp', str(disk_limit_tmp)])
    if (start is not None) and (end is not None) and (not is_full_requested_spot_range(sra_stat, start, end)):
        normalized_start, normalized_end = normalize_requested_spot_range(start=start, end=end)
        command.extend(['-N', str(normalized_start), '-X', str(normalized_end)])
    command.append(path_downloaded_sra)
    return command


def execute_fasterq_dump_command(fasterq_dump_command, args, prefix='Command'):
    print('{}: {}'.format(prefix, ' '.join(fasterq_dump_command)))
    fqd_out = subprocess.run(fasterq_dump_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if should_print_getfastq_command_output(args):
        print('fasterq-dump stdout:')
        print(fqd_out.stdout.decode('utf8', errors='replace'))
        print('fasterq-dump stderr:')
        print(fqd_out.stderr.decode('utf8', errors='replace'))
    return fqd_out


def run_fasterq_dump_with_retry(fasterq_dump_command, path_downloaded_sra, metadata, sra_stat, args):
    fqd_out = execute_fasterq_dump_command(fasterq_dump_command, args, prefix='Command')
    if fqd_out.returncode == 0:
        return fqd_out

    sys.stderr.write("fasterq-dump did not finish safely. Removing the cached SRA file and retrying once.\n")
    remove_sra_path(path_downloaded_sra)
    download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=sra_stat['getfastq_sra_dir'], overwrite=True)
    fqd_out = execute_fasterq_dump_command(fasterq_dump_command, args, prefix='Retry command')
    if fqd_out.returncode != 0:
        sys.stderr.write("fasterq-dump did not finish safely after re-download.\n")
        sys.exit(1)
    return fqd_out


def ensure_fasterq_output_files_exist(sra_stat):
    sra_id = sra_stat['sra_id']
    work_dir = sra_stat['getfastq_sra_dir']
    candidate_paths = [
        os.path.join(work_dir, sra_id + '.fastq'),
        os.path.join(work_dir, sra_id + '_1.fastq'),
        os.path.join(work_dir, sra_id + '_2.fastq'),
    ]
    detected = []
    for path in candidate_paths:
        if not os.path.exists(path):
            continue
        if not os.path.isfile(path):
            raise IsADirectoryError('fasterq-dump output path exists but is not a file: {}'.format(path))
        detected.append(path)
    if len(detected) == 0:
        raise FileNotFoundError(
            'fasterq-dump did not generate FASTQ files for {} under {}'.format(sra_id, work_dir)
        )


def resolve_written_spots_from_fasterq_output(
    sra_stat,
    start,
    end,
    fqd_out,
    run_file_state,
    args=None,
    compress_trimmed_to_gz=False,
):
    stdout_txt = fqd_out.stdout.decode('utf8', errors='replace')
    stderr_txt = fqd_out.stderr.decode('utf8', errors='replace')
    written_spots = parse_fasterq_dump_written_spots(
        stdout_txt=stdout_txt,
        stderr_txt=stderr_txt,
    )
    if written_spots is not None:
        print('Using fasterq-dump reported spot count: {:,}'.format(written_spots))
        return written_spots, run_file_state
    written_reads = parse_fasterq_dump_written_reads(
        stdout_txt=stdout_txt,
        stderr_txt=stderr_txt,
    )
    if written_reads is not None:
        written_spots = infer_written_spots_from_written_reads(
            sra_stat=sra_stat,
            written_reads=written_reads,
            run_file_state=run_file_state,
        )
        if written_spots is not None:
            txt = 'Using fasterq-dump reported read count to infer spot count: {:,} (reads written: {:,})'
            print(txt.format(written_spots, written_reads))
            return written_spots, run_file_state
        txt = 'fasterq-dump reported reads written ({:,}), but exact spot inference was not possible. '
        txt += 'Falling back to FASTQ record counting.'
        print(txt.format(written_reads))
    else:
        print('fasterq-dump did not report written spots. Falling back to FASTQ record counting.')
    return written_spots, run_file_state


def calculate_requested_spots(start, end):
    normalized_start, normalized_end = normalize_requested_spot_range(start=start, end=end)
    return max(0, (normalized_end - normalized_start + 1))

def should_compress_fasterq_output_before_filters(args):
    filter_order = get_filter_execution_order(args)
    return len(filter_order) == 0


def update_extraction_counts(metadata, ind_sra, written_spots, dumped_spots, spot_length):
    written_spots = max(0, int(written_spots))
    dumped_spots = max(0, int(dumped_spots))
    if written_spots > dumped_spots:
        dumped_spots = written_spots
    rejected_spots = dumped_spots - written_spots
    metadata.df.at[ind_sra, 'num_dumped'] += dumped_spots
    metadata.df.at[ind_sra, 'num_rejected'] += rejected_spots
    metadata.df.at[ind_sra, 'num_written'] += written_spots
    metadata.df.at[ind_sra, 'bp_dumped'] += dumped_spots * spot_length
    metadata.df.at[ind_sra, 'bp_rejected'] += rejected_spots * spot_length
    metadata.df.at[ind_sra, 'bp_written'] += written_spots * spot_length


def run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
    path_downloaded_sra = os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + '.sra')
    start, end = normalize_requested_spot_range(start=start, end=end)
    raw_size_check = getattr(args, 'fasterq_size_check', True)
    size_check = normalize_fasterq_size_check(raw_size_check)
    fasterq_dump_command = build_fasterq_dump_command(
        args=args,
        sra_stat=sra_stat,
        path_downloaded_sra=path_downloaded_sra,
        size_check=size_check,
        start=start,
        end=end,
    )
    print('Total sampled bases:', "{:,}".format(sra_stat['spot_length'] * (end - start + 1)), 'bp')
    fqd_out = run_fasterq_dump_with_retry(
        fasterq_dump_command=fasterq_dump_command,
        path_downloaded_sra=path_downloaded_sra,
        metadata=metadata,
        sra_stat=sra_stat,
        args=args,
    )
    ensure_fasterq_output_files_exist(sra_stat=sra_stat)
    run_file_state = RunFileState(work_dir=sra_stat['getfastq_sra_dir'])
    written_spots, run_file_state = resolve_written_spots_from_fasterq_output(
        sra_stat=sra_stat,
        start=start,
        end=end,
        fqd_out=fqd_out,
        run_file_state=run_file_state,
    )
    should_compress = should_compress_fasterq_output_before_filters(args)
    if should_compress:
        updated_file_state = compress_fasterq_output_files(
            sra_stat=sra_stat,
            args=args,
            file_state=run_file_state,
            return_file_state=True,
        )
        if isinstance(updated_file_state, RunFileState):
            run_file_state = updated_file_state
        else:
            run_file_state = RunFileState(work_dir=sra_stat['getfastq_sra_dir'])
        set_current_intermediate_extension(sra_stat, '.fastq.gz')
    else:
        print('Skipping seqkit compression after fasterq-dump because downstream filtering is enabled.')
        set_current_intermediate_extension(sra_stat, '.fastq')
    if written_spots is None:
        written_spots = estimate_num_written_spots_from_fastq(sra_stat, file_state=run_file_state)
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_stat['sra_id']))
    requested_spots = calculate_requested_spots(start=start, end=end)
    update_extraction_counts(
        metadata=metadata,
        ind_sra=ind_sra,
        written_spots=written_spots,
        dumped_spots=requested_spots,
        spot_length=sra_stat['spot_length'],
    )
    sra_stat = detect_layout_from_file(sra_stat, files=run_file_state.files)
    updated_file_state = remove_unpaired_files(sra_stat, file_state=run_file_state, return_file_state=True)
    if isinstance(updated_file_state, RunFileState):
        run_file_state = updated_file_state
    metadata.df.at[ind_sra,'layout_amalgkit'] = sra_stat['layout']
    if return_file_state:
        return metadata, sra_stat, run_file_state
    if return_files:
        return metadata, sra_stat, run_file_state.to_set()
    return metadata,sra_stat

def normalize_fasterq_size_check(raw_size_check):
    if isinstance(raw_size_check, str):
        normalized_size_check = raw_size_check.strip().lower()
        if normalized_size_check in {'on', 'off', 'only'}:
            return normalized_size_check
        if normalized_size_check in {'yes', 'y', 'true', '1'}:
            return 'on'
        if normalized_size_check in {'no', 'n', 'false', '0'}:
            return 'off'
        return 'on'
    return 'on' if bool(raw_size_check) else 'off'

def remove_unpaired_files(sra_stat, files=None, file_state=None, return_file_state=False):
    run_file_state = _resolve_run_file_state(
        work_dir=sra_stat['getfastq_sra_dir'],
        files=files,
        file_state=file_state,
    )
    if (sra_stat['layout']=='paired'):
        # Order is important in this list. More downstream should come first.
        extensions = ['.amalgkit.fastq.gz', '.rename.fastq.gz', '.fastp.fastq.gz', '.fastq.gz', '.fastq']
        for ext in extensions:
            single_filename = sra_stat['sra_id'] + ext
            if run_file_state.has(single_filename):
                single_fastq_file = run_file_state.path(single_filename)
                print('Removing 3rd fastq file: {}'.format(single_fastq_file), flush=True)
                os.remove(single_fastq_file)
                run_file_state.discard(single_filename)
    if return_file_state:
        return run_file_state
    return run_file_state.to_set()

def get_identical_paired_ratio(read1_path, read2_path, num_checked_reads=IDENTICAL_PAIRED_CHECKED_READS):
    num_identical = 0
    num_checked = 0
    read_length = 0
    open_read1 = gzip.open if str(read1_path).endswith('.gz') else open
    open_read2 = gzip.open if str(read2_path).endswith('.gz') else open
    with open_read1(read1_path, 'rt') as read1, open_read2(read2_path, 'rt') as read2:
        while num_checked < num_checked_reads:
            block1 = [read1.readline() for _ in range(4)]
            block2 = [read2.readline() for _ in range(4)]
            if any([line == '' for line in block1 + block2]):
                break
            seq1 = block1[1].strip()
            seq2 = block2[1].strip()
            if num_checked == 0:
                read_length = len(seq1)
            if seq1 == seq2:
                num_identical += 1
            num_checked += 1
    identical_ratio = num_identical / num_checked if num_checked else 0
    return identical_ratio, num_checked, read_length


def finalize_maybe_treat_paired_result(metadata, sra_stat, run_file_state, return_files=False, return_file_state=False):
    if return_file_state:
        return metadata, sra_stat, run_file_state
    if return_files:
        return metadata, sra_stat, run_file_state.to_set()
    return metadata, sra_stat


def resolve_paired_read_paths_for_identical_check(sra_stat, work_dir, run_file_state):
    if sra_stat['layout'] != 'paired':
        return None
    inext = get_or_detect_intermediate_extension(sra_stat, work_dir=work_dir, file_state=run_file_state)
    if inext == 'no_extension_found':
        return None
    read1_name = sra_stat['sra_id'] + '_1' + inext
    read2_name = sra_stat['sra_id'] + '_2' + inext
    if (not run_file_state.has(read1_name)) or (not run_file_state.has(read2_name)):
        return None
    return {
        'inext': inext,
        'read1_name': read1_name,
        'read2_name': read2_name,
        'read1_path': run_file_state.path(read1_name),
        'read2_path': run_file_state.path(read2_name),
    }


def convert_paired_files_to_single(sra_stat, metadata, run_file_state, read_info, identical_ratio, num_checked, read_length):
    txt = 'Read1 and Read2 are nearly identical ({:.2%} of {:,} pairs): {}. '
    txt += 'Treating as single-end reads and removing redundant read2 file.\n'
    sys.stderr.write(txt.format(identical_ratio, num_checked, sra_stat['sra_id']))
    inext = read_info['inext']
    read1_name = read_info['read1_name']
    read2_name = read_info['read2_name']
    read1_path = read_info['read1_path']
    read2_path = read_info['read2_path']
    single_name = sra_stat['sra_id'] + inext
    single_path = run_file_state.path(single_name)
    if run_file_state.has(single_name):
        os.remove(single_path)
        run_file_state.discard(single_name)
    os.rename(read1_path, single_path)
    os.remove(read2_path)
    run_file_state.discard(read1_name)
    run_file_state.discard(read2_name)
    run_file_state.add(single_name)
    sra_stat['layout'] = 'single'
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_stat['sra_id']))
    if 'layout_amalgkit' in metadata.df.columns:
        metadata.df.at[ind_sra, 'layout_amalgkit'] = 'single'
    if (read_length > 0) and ('spot_length' in metadata.df.columns):
        sra_stat['spot_length'] = read_length
        metadata.df.at[ind_sra, 'spot_length'] = read_length
        if 'spot_length_amalgkit' in metadata.df.columns:
            metadata.df.at[ind_sra, 'spot_length_amalgkit'] = read_length
    set_current_intermediate_extension(sra_stat, inext)


def maybe_treat_paired_as_single(
    sra_stat,
    metadata,
    work_dir,
    threshold=IDENTICAL_PAIRED_RATIO_THRESHOLD,
    num_checked_reads=IDENTICAL_PAIRED_CHECKED_READS,
    files=None,
    file_state=None,
    return_files=False,
    return_file_state=False,
):
    run_file_state = _resolve_run_file_state(work_dir=work_dir, files=files, file_state=file_state)
    read_info = resolve_paired_read_paths_for_identical_check(
        sra_stat=sra_stat,
        work_dir=work_dir,
        run_file_state=run_file_state,
    )
    if read_info is None:
        return finalize_maybe_treat_paired_result(
            metadata=metadata,
            sra_stat=sra_stat,
            run_file_state=run_file_state,
            return_files=return_files,
            return_file_state=return_file_state,
        )
    identical_ratio, num_checked, read_length = get_identical_paired_ratio(
        read1_path=read_info['read1_path'],
        read2_path=read_info['read2_path'],
        num_checked_reads=num_checked_reads,
    )
    if (num_checked > 0) and (identical_ratio >= threshold):
        convert_paired_files_to_single(
            sra_stat=sra_stat,
            metadata=metadata,
            run_file_state=run_file_state,
            read_info=read_info,
            identical_ratio=identical_ratio,
            num_checked=num_checked,
            read_length=read_length,
        )
    return finalize_maybe_treat_paired_result(
        metadata=metadata,
        sra_stat=sra_stat,
        run_file_state=run_file_state,
        return_files=return_files,
        return_file_state=return_file_state,
    )

def parse_fastp_metrics(stderr_txt):
    duplication_rate = numpy.nan
    insert_size_peak = numpy.nan
    for line in stderr_txt.split('\n'):
        if line.startswith('Duplication rate:'):
            matched = re.search(r'Duplication rate:\s*([0-9.]+)%', line)
            if matched:
                duplication_rate = float(matched.group(1))
        if line.startswith('Insert size peak'):
            matched = re.search(r'Insert size peak.*:\s*([0-9.]+)', line)
            if matched:
                insert_size_peak = float(matched.group(1))
    return duplication_rate, insert_size_peak

def parse_fastp_summary_counts(stderr_txt):
    lines = stderr_txt.split('\n')
    num_in = []
    num_out = []
    bp_in = []
    bp_out = []

    def _parse_total_line(lines, idx, expected_prefix, section_name):
        if idx >= len(lines):
            raise RuntimeError('Unexpected fastp stderr format: missing "{}" in {} section.'.format(expected_prefix, section_name))
        raw = lines[idx].strip()
        if not raw.startswith(expected_prefix):
            raise RuntimeError('Unexpected fastp stderr format: expected "{}" in {} section, got "{}".'.format(
                expected_prefix, section_name, raw
            ))
        value_txt = raw.replace(expected_prefix, '').replace(',', '').strip()
        try:
            return int(value_txt)
        except ValueError as e:
            raise RuntimeError('Unexpected fastp stderr format: non-integer value "{}" for {}.'.format(value_txt, expected_prefix)) from e

    for i, line in enumerate(lines):
        if ' before filtering:' in line:
            num_in.append(_parse_total_line(lines, i + 1, 'total reads:', 'before filtering'))
            bp_in.append(_parse_total_line(lines, i + 2, 'total bases:', 'before filtering'))
        if (' after filtering:' in line) or (' aftering filtering:' in line):
            num_out.append(_parse_total_line(lines, i + 1, 'total reads:', 'after filtering'))
            bp_out.append(_parse_total_line(lines, i + 2, 'total bases:', 'after filtering'))

    return num_in, num_out, bp_in, bp_out

def update_fastp_metrics(metadata, ind_sra, current_num_in, duplication_rate, insert_size_peak):
    previous_num_in = metadata.df.at[ind_sra, 'num_fastp_in']
    total_num_in = previous_num_in + current_num_in
    metric_pairs = [
        ('fastp_duplication_rate', duplication_rate),
        ('fastp_insert_size_peak', insert_size_peak),
    ]
    for metric_key, current_value in metric_pairs:
        if numpy.isnan(current_value):
            continue
        previous_value = metadata.df.at[ind_sra, metric_key]
        if (current_num_in <= 0) or numpy.isnan(previous_value) or (previous_num_in <= 0):
            metadata.df.at[ind_sra, metric_key] = current_value
        else:
            weighted_mean = ((previous_value * previous_num_in) + (current_value * current_num_in)) / total_num_in
            metadata.df.at[ind_sra, metric_key] = weighted_mean

def write_fastp_stats(sra_stat, metadata, output_dir):
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_stat['sra_id']))
    out_df = pandas.DataFrame([{
        'run': sra_stat['sra_id'],
        'fastp_duplication_rate': metadata.df.at[ind_sra, 'fastp_duplication_rate'],
        'fastp_insert_size_peak': metadata.df.at[ind_sra, 'fastp_insert_size_peak'],
        'num_fastp_in': metadata.df.at[ind_sra, 'num_fastp_in'],
        'num_fastp_out': metadata.df.at[ind_sra, 'num_fastp_out'],
        'bp_fastp_in': metadata.df.at[ind_sra, 'bp_fastp_in'],
        'bp_fastp_out': metadata.df.at[ind_sra, 'bp_fastp_out'],
    }])
    out_path = os.path.join(output_dir, 'fastp_stats.tsv')
    out_df.to_csv(out_path, sep='\t', index=False)


def write_getfastq_stats(sra_stat, metadata, output_dir):
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_stat['sra_id']))
    row = {'run': sra_stat['sra_id']}
    for col in GETFASTQ_STATS_COLUMNS:
        if col in metadata.df.columns:
            row[col] = metadata.df.at[ind_sra, col]
        else:
            row[col] = numpy.nan
    out_df = pandas.DataFrame([row])
    out_path = os.path.join(output_dir, 'getfastq_stats.tsv')
    out_df.to_csv(out_path, sep='\t', index=False)


def resolve_fastp_threads(num_threads):
    if num_threads > 16:
        print('Too many threads for fastp (--threads {}). Only 16 threads will be used.'.format(num_threads))
        return 16
    return num_threads


def parse_fastp_option_args(fastp_option):
    if not fastp_option:
        return []
    try:
        return shlex.split(fastp_option)
    except ValueError as e:
        raise ValueError('Invalid --fastp_option string: {}'.format(fastp_option)) from e


def build_fastp_io_arguments(sra_stat, output_dir, inext, outext):
    input_names = []
    output_names = []
    io_args = []
    if sra_stat['layout'] == 'single':
        infile = os.path.join(output_dir, sra_stat['sra_id'])
        input_names = [sra_stat['sra_id'] + inext]
        output_names = [sra_stat['sra_id'] + outext]
        io_args = ['--in1', infile + inext, '--out1', infile + outext]
    elif sra_stat['layout'] == 'paired':
        infile1 = os.path.join(output_dir, sra_stat['sra_id'] + '_1')
        infile2 = os.path.join(output_dir, sra_stat['sra_id'] + '_2')
        input_names = [sra_stat['sra_id'] + '_1' + inext, sra_stat['sra_id'] + '_2' + inext]
        output_names = [sra_stat['sra_id'] + '_1' + outext, sra_stat['sra_id'] + '_2' + outext]
        io_args = [
            '--in1', infile1 + inext, '--out1', infile1 + outext,
            '--in2', infile2 + inext, '--out2', infile2 + outext,
        ]
    return io_args, input_names, output_names


def finalize_run_fastp_return(metadata, run_file_state, return_files=False, return_file_state=False):
    if return_file_state:
        return metadata, run_file_state
    if return_files:
        if run_file_state is None:
            return metadata, None
        return metadata, run_file_state.to_set()
    return metadata


def execute_fastp_command(fp_command, fastp_print=False):
    print('Command:', ' '.join(fp_command))
    fp_out = subprocess.run(fp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if fastp_print:
        print('fastp stdout:')
        print(fp_out.stdout.decode('utf8', errors='replace'))
        print('fastp stderr:')
        print(fp_out.stderr.decode('utf8', errors='replace'))
    if fp_out.returncode != 0:
        raise RuntimeError(
            'fastp did not finish safely (exit code {}). Command: {}\n{}'.format(
                fp_out.returncode,
                ' '.join(fp_command),
                fp_out.stderr.decode('utf8', errors='replace'),
            )
        )
    return fp_out


def ensure_fastp_output_files(output_dir, output_names):
    for output_name in output_names:
        output_path = os.path.join(output_dir, output_name)
        if not os.path.exists(output_path):
            raise FileNotFoundError('fastp output file was not generated: {}'.format(output_path))
        if not os.path.isfile(output_path):
            raise IsADirectoryError('fastp output path exists but is not a file: {}'.format(output_path))


def update_metadata_after_fastp(metadata, sra_stat, output_dir, fp_stderr):
    num_in, num_out, bp_in, bp_out = parse_fastp_summary_counts(fp_stderr)
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_stat['sra_id']))
    current_num_in = sum(num_in)
    duplication_rate, insert_size_peak = parse_fastp_metrics(fp_stderr)
    update_fastp_metrics(metadata, ind_sra, current_num_in, duplication_rate, insert_size_peak)
    metadata.df.at[ind_sra, 'num_fastp_in'] += sum(num_in)
    metadata.df.at[ind_sra, 'num_fastp_out'] += sum(num_out)
    metadata.df.at[ind_sra, 'bp_fastp_in'] += sum(bp_in)
    metadata.df.at[ind_sra, 'bp_fastp_out'] += sum(bp_out)
    write_fastp_stats(sra_stat, metadata, output_dir)
    return metadata


def run_fastp(
    sra_stat,
    args,
    output_dir,
    metadata,
    files=None,
    file_state=None,
    return_files=False,
    return_file_state=False,
):
    run_file_state = None
    if (files is not None) or isinstance(file_state, RunFileState):
        run_file_state = _resolve_run_file_state(work_dir=output_dir, files=files, file_state=file_state)
    inext = get_or_detect_intermediate_extension(sra_stat, work_dir=output_dir, file_state=run_file_state)
    outext = '.fastp.fastq.gz'
    fastp_thread = resolve_fastp_threads(args.threads)
    fastp_exe = getattr(args, 'fastp_exe', 'fastp')
    fastp_option_args = parse_fastp_option_args(args.fastp_option)
    fp_command = [fastp_exe, '--thread', str(fastp_thread), '--length_required', str(args.min_read_length)] + fastp_option_args
    io_args, input_names, output_names = build_fastp_io_arguments(
        sra_stat=sra_stat,
        output_dir=output_dir,
        inext=inext,
        outext=outext,
    )
    fp_command = fp_command + io_args
    fp_command = [fc for fc in fp_command if fc != '']
    fp_out = execute_fastp_command(fp_command, fastp_print=args.fastp_print)
    ensure_fastp_output_files(output_dir=output_dir, output_names=output_names)
    if run_file_state is not None:
        for output_name in output_names:
            run_file_state.add(output_name)
    if args.remove_tmp:
        remove_intermediate_files(sra_stat, ext=inext, work_dir=output_dir)
        if run_file_state is not None:
            for input_name in input_names:
                run_file_state.discard(input_name)
    fp_stderr = fp_out.stderr.decode('utf8', errors='replace')
    metadata = update_metadata_after_fastp(metadata, sra_stat, output_dir, fp_stderr)
    set_current_intermediate_extension(sra_stat, outext)
    return finalize_run_fastp_return(
        metadata=metadata,
        run_file_state=run_file_state,
        return_files=return_files,
        return_file_state=return_file_state,
    )

def prepare_sortmerna_workdir(work_dir):
    os.makedirs(work_dir, exist_ok=True)
    for sub_dir in ['kvdb', 'out']:
        sub_path = os.path.join(work_dir, sub_dir)
        if not os.path.exists(sub_path):
            continue
        if os.path.isdir(sub_path):
            shutil.rmtree(sub_path)
        else:
            os.remove(sub_path)
    os.makedirs(os.path.join(work_dir, 'out'), exist_ok=True)

def build_sortmerna_command(args, sortmerna_refs, input_paths, layout, work_dir, other_prefix):
    sortmerna_exe = getattr(args, 'sortmerna_exe', 'sortmerna')
    command = [sortmerna_exe]
    for sortmerna_ref in sortmerna_refs:
        command.extend(['--ref', sortmerna_ref])
    for input_path in input_paths:
        command.extend(['--reads', input_path])
    command.extend([
        '--workdir', work_dir,
        '--threads', str(max(1, int(getattr(args, 'threads', 1)))),
        '--fastx',
        '--other', other_prefix,
    ])
    if layout == 'paired':
        command.extend(['--paired_in', '--out2'])
    command.extend(parse_sortmerna_option_args(getattr(args, 'sortmerna_option', '')))
    return command

def execute_sortmerna_command(sortmerna_command, args):
    print('Command:', ' '.join(sortmerna_command))
    sm_out = subprocess.run(sortmerna_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if should_print_getfastq_command_output(args):
        print('sortmerna stdout:')
        print(sm_out.stdout.decode('utf8', errors='replace'))
        print('sortmerna stderr:')
        print(sm_out.stderr.decode('utf8', errors='replace'))
    if sm_out.returncode != 0:
        raise RuntimeError(
            'SortMeRNA did not finish safely (exit code {}). Command: {}\n{}'.format(
                sm_out.returncode,
                ' '.join(sortmerna_command),
                sm_out.stderr.decode('utf8', errors='replace'),
            )
        )
    return sm_out

def _parse_int_matches_from_text(pattern, text):
    if text is None:
        return []
    matched_values = re.findall(pattern, text, flags=re.IGNORECASE)
    parsed = []
    for raw in matched_values:
        try:
            parsed.append(int(str(raw).replace(',', '')))
        except (TypeError, ValueError):
            continue
    return parsed

def infer_sortmerna_spot_divisor(layout, reads_total, num_input_files=1, known_input_counts=None):
    if layout != 'paired':
        return 1
    default_divisor = max(1, int(num_input_files))
    parsed_known_input = None
    if known_input_counts is not None:
        try:
            parsed_known_input = parse_known_fastq_counts(known_input_counts)
        except Exception:
            parsed_known_input = None
    if parsed_known_input is None:
        return default_divisor
    known_input_spots = parsed_known_input[0]
    if known_input_spots <= 0:
        return default_divisor
    ratio = float(reads_total) / float(known_input_spots)
    candidate_divisors = sorted(set([1, 2, default_divisor]))
    best_divisor = min(candidate_divisors, key=lambda candidate: abs(ratio - float(candidate)))
    return max(1, int(best_divisor))

def parse_sortmerna_output_counts(stdout_txt, stderr_txt, layout, num_input_files=1, known_input_counts=None):
    combined_txt = '\n'.join([stdout_txt or '', stderr_txt or ''])
    reads_added_values = _parse_int_matches_from_text(
        pattern=r'Reads\s+added:\s*([0-9][0-9,]*)',
        text=combined_txt,
    )
    aligned_values = _parse_int_matches_from_text(
        pattern=r'Num\s+aligned\s+reads\s*\(passing\s+E-value\)\s*:\s*([0-9][0-9,]*)',
        text=combined_txt,
    )
    all_reads_count_values = _parse_int_matches_from_text(
        pattern=r'all_reads_count\s*=\s*([0-9][0-9,]*)',
        text=combined_txt,
    )
    total_reads_mapped_values = _parse_int_matches_from_text(
        pattern=r'total_reads_mapped\s*=\s*([0-9][0-9,]*)',
        text=combined_txt,
    )
    all_reads_len_values = _parse_int_matches_from_text(
        pattern=r'all_reads_len\s*=\s*([0-9][0-9,]*)',
        text=combined_txt,
    )
    reads_total = max(reads_added_values) if len(reads_added_values) > 0 else None
    if reads_total is None:
        reads_total = max(all_reads_count_values) if len(all_reads_count_values) > 0 else None
    mapped_total = max(aligned_values) if len(aligned_values) > 0 else None
    if mapped_total is None:
        mapped_total = max(total_reads_mapped_values) if len(total_reads_mapped_values) > 0 else None
    if (reads_total is None) or (mapped_total is None):
        return None
    reads_non_rrna_total = max(0, int(reads_total) - int(mapped_total))
    divisor = infer_sortmerna_spot_divisor(
        layout=layout,
        reads_total=reads_total,
        num_input_files=num_input_files,
        known_input_counts=known_input_counts,
    )
    num_out_spots = int(round(float(reads_non_rrna_total) / float(divisor)))
    num_out_spots = max(0, num_out_spots)
    bp_out = None
    parsed_known_input = None
    if known_input_counts is not None:
        try:
            parsed_known_input = parse_known_fastq_counts(known_input_counts)
        except Exception:
            parsed_known_input = None
    if (parsed_known_input is not None) and (parsed_known_input[0] > 0):
        bp_per_spot = float(parsed_known_input[1]) / float(parsed_known_input[0])
        bp_out = int(round(float(num_out_spots) * bp_per_spot))
    elif (len(all_reads_len_values) > 0) and (reads_total > 0):
        all_reads_len_total = max(all_reads_len_values)
        bp_per_spot = (float(all_reads_len_total) * float(divisor)) / float(reads_total)
        bp_out = int(round(float(num_out_spots) * bp_per_spot))
    if bp_out is None:
        return {'num_spots': num_out_spots}
    return {
        'num_spots': num_out_spots,
        'bp_total': max(0, int(bp_out)),
    }

def collect_sortmerna_other_fastq_candidates(other_prefix):
    out_dir = os.path.dirname(other_prefix)
    out_basename = os.path.basename(other_prefix)
    if not os.path.isdir(out_dir):
        return []
    candidates = []
    valid_exts = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
    with os.scandir(out_dir) as entries:
        for entry in entries:
            if (not entry.is_file()) or (not entry.name.startswith(out_basename)):
                continue
            if not entry.name.lower().endswith(valid_exts):
                continue
            candidates.append(entry.path)
    return sorted(candidates)

def resolve_sortmerna_single_output_path(other_prefix, candidates):
    valid_exts = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
    for ext in valid_exts:
        preferred = other_prefix + ext
        if preferred in candidates:
            return preferred
    if len(candidates) == 1:
        return candidates[0]
    return None

def resolve_sortmerna_paired_output_paths(other_prefix, candidates):
    valid_exts = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
    pair_suffixes = [
        ('_fwd', '_rev'),
        ('-fwd', '-rev'),
        ('.fwd', '.rev'),
        ('_paired_fwd', '_paired_rev'),
        ('-paired-fwd', '-paired-rev'),
        ('_R1', '_R2'),
        ('_r1', '_r2'),
        ('_1', '_2'),
        ('.1', '.2'),
        ('-1', '-2'),
    ]
    candidate_set = set(candidates)
    for suffix1, suffix2 in pair_suffixes:
        for ext in valid_exts:
            path1 = other_prefix + suffix1 + ext
            path2 = other_prefix + suffix2 + ext
            if (path1 in candidate_set) and (path2 in candidate_set):
                return path1, path2
    fwd_tokens = ['fwd', 'forward', '_r1', '.r1', '-r1', '_1.', '.1.', '-1.']
    rev_tokens = ['rev', 'reverse', '_r2', '.r2', '-r2', '_2.', '.2.', '-2.']
    fwd_candidates = []
    rev_candidates = []
    for candidate in candidates:
        lower_name = os.path.basename(candidate).lower()
        if any([token in lower_name for token in fwd_tokens]):
            fwd_candidates.append(candidate)
            continue
        if any([token in lower_name for token in rev_tokens]):
            rev_candidates.append(candidate)
            continue
    if (len(fwd_candidates) == 1) and (len(rev_candidates) == 1):
        return fwd_candidates[0], rev_candidates[0]
    return None, None

def create_empty_fastq_gz(path_fastq_gz):
    if os.path.exists(path_fastq_gz):
        if not os.path.isfile(path_fastq_gz):
            raise IsADirectoryError('Output path exists but is not a file: {}'.format(path_fastq_gz))
        os.remove(path_fastq_gz)
    with gzip.open(path_fastq_gz, 'wb'):
        pass

def materialize_sortmerna_output(source_path, dest_path, args, seqkit_threads=None):
    if source_path is None:
        create_empty_fastq_gz(dest_path)
        return
    if not os.path.exists(source_path):
        raise FileNotFoundError('SortMeRNA output file not found: {}'.format(source_path))
    if not os.path.isfile(source_path):
        raise IsADirectoryError('SortMeRNA output path exists but is not a file: {}'.format(source_path))
    if os.path.exists(dest_path):
        if not os.path.isfile(dest_path):
            raise IsADirectoryError('SortMeRNA destination path exists but is not a file: {}'.format(dest_path))
        os.remove(dest_path)
    if source_path.endswith('.gz'):
        shutil.move(source_path, dest_path)
    else:
        compress_fastq_with_seqkit(
            path_fastq=source_path,
            path_fastq_gz=dest_path,
            args=args,
            command_label='SortMeRNA output compression with seqkit',
            seqkit_threads=seqkit_threads,
        )

def parse_known_fastq_counts(known_counts):
    if known_counts is None:
        return None
    if isinstance(known_counts, dict):
        num_value = known_counts.get('num_spots', known_counts.get('num_records', None))
        bp_value = known_counts.get('bp_total', known_counts.get('num_bases', known_counts.get('bp', None)))
    elif isinstance(known_counts, (tuple, list)) and (len(known_counts) >= 2):
        num_value = known_counts[0]
        bp_value = known_counts[1]
    else:
        raise ValueError('known_input_counts must be a dict or a 2-item tuple/list.')
    if (num_value is None) or (bp_value is None):
        raise ValueError('known_input_counts must include read-count and base-count values.')
    try:
        num_spots = max(0, int(num_value))
        bp_total = max(0, int(bp_value))
    except (TypeError, ValueError) as exc:
        raise ValueError('known_input_counts values must be integer-like.') from exc
    return num_spots, bp_total

def update_metadata_after_rrna_filter(
    metadata,
    sra_stat,
    input_paths,
    output_paths,
    args=None,
    known_input_counts=None,
    known_output_counts=None,
):
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_stat['sra_id']))
    layout = sra_stat['layout']
    parsed_known_counts = None
    if known_input_counts is not None:
        try:
            parsed_known_counts = parse_known_fastq_counts(known_input_counts)
        except Exception as exc:
            sys.stderr.write(
                'Invalid known_input_counts for rRNA metadata update. Falling back to FASTQ scan. {}\n'.format(exc)
            )
    if parsed_known_counts is not None:
        num_in, bp_in = parsed_known_counts
    else:
        num_in, bp_in = summarize_layout_fastq_records_and_bases(
            layout=layout,
            single_path=input_paths.get(''),
            pair1_path=input_paths.get('_1'),
            pair2_path=input_paths.get('_2'),
            args=args,
        )
    parsed_known_output = None
    if known_output_counts is not None:
        try:
            parsed_known_output = parse_known_fastq_counts(known_output_counts)
        except Exception as exc:
            sys.stderr.write(
                'Invalid known_output_counts for rRNA metadata update. Falling back to FASTQ scan. {}\n'.format(exc)
            )
    if parsed_known_output is not None:
        num_out, bp_out = parsed_known_output
    else:
        num_out, bp_out = summarize_layout_fastq_records_and_bases(
            layout=layout,
            single_path=output_paths.get(''),
            pair1_path=output_paths.get('_1'),
            pair2_path=output_paths.get('_2'),
            args=args,
        )
    metadata.df.at[ind_sra, 'num_rrna_in'] += num_in
    metadata.df.at[ind_sra, 'num_rrna_out'] += num_out
    metadata.df.at[ind_sra, 'bp_rrna_in'] += bp_in
    metadata.df.at[ind_sra, 'bp_rrna_out'] += bp_out
    return metadata

def run_sortmerna_rrna_filter(
    sra_stat,
    args,
    output_dir,
    metadata,
    files=None,
    file_state=None,
    known_input_counts=None,
    return_files=False,
    return_file_state=False,
):
    if not is_rrna_filter_enabled(args):
        run_file_state = _resolve_run_file_state(work_dir=output_dir, files=files, file_state=file_state)
        return finalize_run_fastp_return(
            metadata=metadata,
            run_file_state=run_file_state,
            return_files=return_files,
            return_file_state=return_file_state,
        )
    run_file_state = _resolve_run_file_state(work_dir=output_dir, files=files, file_state=file_state)
    inext = get_or_detect_intermediate_extension(sra_stat, work_dir=output_dir, file_state=run_file_state)
    sra_id = sra_stat['sra_id']
    layout = sra_stat['layout']
    if layout == 'single':
        input_names = [sra_id + inext]
        output_name_by_suffix = {'': sra_id + '.rrna.fastq.gz'}
    elif layout == 'paired':
        input_names = [sra_id + '_1' + inext, sra_id + '_2' + inext]
        output_name_by_suffix = {
            '_1': sra_id + '_1.rrna.fastq.gz',
            '_2': sra_id + '_2.rrna.fastq.gz',
        }
    else:
        raise ValueError('Unsupported library layout for SortMeRNA filtering: {}'.format(layout))
    input_paths = {}
    for input_name in input_names:
        if not run_file_state.has(input_name):
            raise FileNotFoundError('SortMeRNA input file not found: {}'.format(run_file_state.path(input_name)))
        input_paths[input_name] = run_file_state.path(input_name)
    sortmerna_refs = resolve_sortmerna_refs(args)
    sortmerna_work_dir = os.path.join(output_dir, 'sortmerna_work')
    prepare_sortmerna_workdir(sortmerna_work_dir)
    other_prefix = os.path.join(sortmerna_work_dir, 'out', 'non_rrna')
    sortmerna_command = build_sortmerna_command(
        args=args,
        sortmerna_refs=sortmerna_refs,
        input_paths=[input_paths[input_name] for input_name in input_names],
        layout=layout,
        work_dir=sortmerna_work_dir,
        other_prefix=other_prefix,
    )
    sm_out = execute_sortmerna_command(sortmerna_command=sortmerna_command, args=args)
    sm_stdout_txt = sm_out.stdout.decode('utf8', errors='replace')
    sm_stderr_txt = sm_out.stderr.decode('utf8', errors='replace')
    known_output_counts = parse_sortmerna_output_counts(
        stdout_txt=sm_stdout_txt,
        stderr_txt=sm_stderr_txt,
        layout=layout,
        num_input_files=len(input_names),
        known_input_counts=known_input_counts,
    )
    if isinstance(known_output_counts, dict):
        if ('num_spots' in known_output_counts) and ('bp_total' in known_output_counts):
            print(
                'Using SortMeRNA log-derived output counts: {:,} reads, {:,} bp.'.format(
                    int(known_output_counts['num_spots']),
                    int(known_output_counts['bp_total']),
                ),
                flush=True,
            )
        elif 'num_spots' in known_output_counts:
            print(
                'Using SortMeRNA log-derived output read count: {:,} reads.'.format(
                    int(known_output_counts['num_spots'])
                ),
                flush=True,
            )
    candidates = collect_sortmerna_other_fastq_candidates(other_prefix=other_prefix)
    output_paths = {}
    if layout == 'single':
        source_path = resolve_sortmerna_single_output_path(other_prefix=other_prefix, candidates=candidates)
        if (source_path is None) and (len(candidates) > 0):
            txt = 'Could not resolve non-rRNA output file from SortMeRNA. Candidates: {}'
            raise FileNotFoundError(txt.format(', '.join(candidates)))
        output_name = output_name_by_suffix['']
        output_path = os.path.join(output_dir, output_name)
        materialize_sortmerna_output(
            source_path=source_path,
            dest_path=output_path,
            args=args,
            seqkit_threads=resolve_seqkit_threads(args),
        )
        run_file_state.add(output_name)
        output_paths[''] = output_path
    else:
        source_path1, source_path2 = resolve_sortmerna_paired_output_paths(
            other_prefix=other_prefix,
            candidates=candidates,
        )
        if (source_path1 is None) and (source_path2 is None) and (len(candidates) > 0):
            txt = 'Could not resolve paired non-rRNA output files from SortMeRNA. Candidates: {}'
            raise FileNotFoundError(txt.format(', '.join(candidates)))
        if (source_path1 is None) != (source_path2 is None):
            txt = 'Could not resolve paired non-rRNA output files from SortMeRNA. Candidates: {}'
            raise FileNotFoundError(txt.format(', '.join(candidates)))
        output_name1 = output_name_by_suffix['_1']
        output_name2 = output_name_by_suffix['_2']
        output_path1 = os.path.join(output_dir, output_name1)
        output_path2 = os.path.join(output_dir, output_name2)
        total_seqkit_threads = resolve_seqkit_threads(args)
        materialize_jobs = [
            (source_path1, output_path1, output_name1, '_1'),
            (source_path2, output_path2, output_name2, '_2'),
        ]
        worker_count = min(len(materialize_jobs), total_seqkit_threads)
        worker_count = max(1, worker_count)
        threads_per_worker = max(1, int(total_seqkit_threads / worker_count))
        if worker_count <= 1:
            for source_path, output_path, _output_name, _suffix in materialize_jobs:
                materialize_sortmerna_output(
                    source_path=source_path,
                    dest_path=output_path,
                    args=args,
                    seqkit_threads=threads_per_worker,
                )
        else:
            print(
                'Materializing SortMeRNA paired outputs with {:,} parallel worker(s), {:,} thread(s) each.'.format(
                    worker_count,
                    threads_per_worker,
                ),
                flush=True,
            )
            def materialize_single_job(job_idx):
                source_path, output_path, _output_name, _suffix = materialize_jobs[job_idx]
                materialize_sortmerna_output(
                    source_path=source_path,
                    dest_path=output_path,
                    args=args,
                    seqkit_threads=threads_per_worker,
                )
                return job_idx
            _, failures = run_tasks_with_optional_threads(
                task_items=list(range(len(materialize_jobs))),
                task_fn=materialize_single_job,
                max_workers=worker_count,
            )
            if failures:
                details = '; '.join([
                    '{}: {}'.format(materialize_jobs[job_idx][3], exc)
                    for job_idx, exc in failures
                ])
                raise RuntimeError('SortMeRNA output materialization failed for paired FASTQ files. {}'.format(details))
        run_file_state.add(output_name1)
        run_file_state.add(output_name2)
        output_paths['_1'] = output_path1
        output_paths['_2'] = output_path2
    input_paths_by_suffix = {}
    if layout == 'single':
        input_paths_by_suffix[''] = input_paths[sra_id + inext]
    else:
        input_paths_by_suffix['_1'] = input_paths[sra_id + '_1' + inext]
        input_paths_by_suffix['_2'] = input_paths[sra_id + '_2' + inext]
    metadata = update_metadata_after_rrna_filter(
        metadata=metadata,
        sra_stat=sra_stat,
        input_paths=input_paths_by_suffix,
        output_paths=output_paths,
        args=args,
        known_input_counts=known_input_counts,
        known_output_counts=known_output_counts,
    )
    if bool(getattr(args, 'remove_tmp', False)):
        remove_intermediate_files(sra_stat, ext=inext, work_dir=output_dir)
        for input_name in input_names:
            run_file_state.discard(input_name)
    set_current_intermediate_extension(sra_stat, '.rrna.fastq.gz')
    return finalize_run_fastp_return(
        metadata=metadata,
        run_file_state=run_file_state,
        return_files=return_files,
        return_file_state=return_file_state,
    )

def rename_reads(sra_stat, args, output_dir, files=None, file_state=None, return_file_state=False):
    run_file_state = _resolve_run_file_state(work_dir=output_dir, files=files, file_state=file_state)
    inext = get_or_detect_intermediate_extension(sra_stat, work_dir=output_dir, file_state=run_file_state)
    outext = '.rename.fastq.gz'
    def rewrite_headers(infile, outfile, suffix, seqkit_threads):
        run_seqkit_replace_command(
            input_path=infile,
            output_path=outfile,
            suffix=suffix,
            args=args,
            command_label='Trinity read-header rewrite with seqkit',
            seqkit_threads=seqkit_threads,
        )

    total_seqkit_threads = resolve_seqkit_threads(args)
    if sra_stat['layout'] == 'single':
        inbase = os.path.join(output_dir, sra_stat['sra_id'])
        if run_file_state.has(sra_stat['sra_id'] + inext):
            infile = inbase + inext
            outfile = inbase + outext
            print('Rewriting read headers for Trinity format: {} -> {}'.format(infile, outfile))
            rewrite_headers(infile=infile, outfile=outfile, suffix='/1', seqkit_threads=total_seqkit_threads)
            run_file_state.add(sra_stat['sra_id'] + outext)
    elif sra_stat['layout'] == 'paired':
        inbase1 = os.path.join(output_dir, sra_stat['sra_id'] + '_1')
        inbase2 = os.path.join(output_dir, sra_stat['sra_id'] + '_2')
        if run_file_state.has(sra_stat['sra_id'] + '_1' + inext):
            infile1 = inbase1 + inext
            infile2 = inbase2 + inext
            outfile1 = inbase1 + outext
            outfile2 = inbase2 + outext
            rewrite_jobs = [
                (infile1, outfile1, '/1'),
                (infile2, outfile2, '/2'),
            ]
            worker_count = min(len(rewrite_jobs), total_seqkit_threads)
            worker_count = max(1, worker_count)
            threads_per_worker = max(1, int(total_seqkit_threads / worker_count))
            if worker_count <= 1:
                for infile, outfile, suffix in rewrite_jobs:
                    print('Rewriting read headers for Trinity format: {} -> {}'.format(infile, outfile))
                    rewrite_headers(
                        infile=infile,
                        outfile=outfile,
                        suffix=suffix,
                        seqkit_threads=threads_per_worker,
                    )
            else:
                print(
                    'Rewriting paired read headers with {:,} parallel worker(s), {:,} thread(s) each.'.format(
                        worker_count,
                        threads_per_worker,
                    ),
                    flush=True,
                )
                def rewrite_single_job(job_idx):
                    infile, outfile, suffix = rewrite_jobs[job_idx]
                    print('Rewriting read headers for Trinity format: {} -> {}'.format(infile, outfile))
                    rewrite_headers(
                        infile=infile,
                        outfile=outfile,
                        suffix=suffix,
                        seqkit_threads=threads_per_worker,
                    )
                    return job_idx
                _, failures = run_tasks_with_optional_threads(
                    task_items=list(range(len(rewrite_jobs))),
                    task_fn=rewrite_single_job,
                    max_workers=worker_count,
                )
                if failures:
                    details = '; '.join([
                        '{}: {}'.format(rewrite_jobs[job_idx][0], exc)
                        for job_idx, exc in failures
                    ])
                    raise RuntimeError('Trinity read-header rewrite with seqkit failed for paired FASTQ(s). {}'.format(details))
            run_file_state.add(sra_stat['sra_id'] + '_1' + outext)
            run_file_state.add(sra_stat['sra_id'] + '_2' + outext)
    if args.remove_tmp:
        remove_intermediate_files(sra_stat, ext=inext, work_dir=output_dir)
        if sra_stat['layout'] == 'single':
            run_file_state.discard(sra_stat['sra_id'] + inext)
        elif sra_stat['layout'] == 'paired':
            run_file_state.discard(sra_stat['sra_id'] + '_1' + inext)
            run_file_state.discard(sra_stat['sra_id'] + '_2' + inext)
    set_current_intermediate_extension(sra_stat, outext)
    if return_file_state:
        return run_file_state
    return run_file_state.to_set()

def rename_fastq(sra_stat, output_dir, inext, outext):
    def rename_single_file(in_path, out_path):
        if not os.path.exists(in_path):
            raise FileNotFoundError('Intermediate fastq file not found for renaming: {}'.format(in_path))
        if not os.path.isfile(in_path):
            raise IsADirectoryError('Intermediate path exists but is not a file: {}'.format(in_path))
        if os.path.exists(out_path) and (not os.path.isfile(out_path)):
            raise IsADirectoryError('Renaming destination exists but is not a file: {}'.format(out_path))
        os.replace(in_path, out_path)

    if sra_stat['layout'] == 'single':
        inbase = os.path.join(output_dir, sra_stat['sra_id'])
        rename_single_file(inbase + inext, inbase + outext)
    elif sra_stat['layout'] == 'paired':
        inbase1 = os.path.join(output_dir, sra_stat['sra_id'] + '_1')
        inbase2 = os.path.join(output_dir, sra_stat['sra_id'] + '_2')
        rename_single_file(inbase1 + inext, inbase1 + outext)
        rename_single_file(inbase2 + inext, inbase2 + outext)
    set_current_intermediate_extension(sra_stat, outext)

def calc_2nd_ranges(metadata):
    df = metadata.df
    sra_target_bp = pandas.to_numeric(df.loc[:, 'bp_until_target_size'], errors='coerce').to_numpy(dtype=float)
    rate_obtained = pandas.to_numeric(df.loc[:, 'rate_obtained'], errors='coerce').to_numpy(dtype=float)
    if 'spot_length_amalgkit' in df.columns:
        spot_length_series = df.loc[:, 'spot_length_amalgkit']
    elif 'spot_length' in df.columns:
        spot_length_series = df.loc[:, 'spot_length']
    else:
        spot_length_series = pandas.Series(numpy.nan, index=df.index)
    spot_lengths = pandas.to_numeric(spot_length_series, errors='coerce').to_numpy(dtype=float)
    total_spots = pandas.to_numeric(df.loc[:, 'total_spots'], errors='coerce').to_numpy(dtype=float)
    start_2nds = pandas.to_numeric(df.loc[:, 'spot_end_1st'], errors='coerce').to_numpy(dtype=float) + 1.0

    sra_target_bp = numpy.where(numpy.isfinite(sra_target_bp) & (sra_target_bp > 0), sra_target_bp, 0.0)
    start_2nds = numpy.where(numpy.isfinite(start_2nds) & (start_2nds > 0), start_2nds, 1.0)
    valid_spot_length = numpy.isfinite(spot_lengths) & (spot_lengths > 0)
    safe_spot_lengths = numpy.where(valid_spot_length, spot_lengths, 1.0)
    safe_total_spots = numpy.where(
        numpy.isfinite(total_spots) & (total_spots >= start_2nds),
        total_spots,
        start_2nds,
    )

    target_reads_base = numpy.divide(
        sra_target_bp,
        safe_spot_lengths,
        out=numpy.zeros_like(sra_target_bp, dtype=float),
        where=valid_spot_length,
    )
    use_rate = numpy.isfinite(rate_obtained) & (rate_obtained > 0)
    compensated_target_reads = numpy.divide(
        target_reads_base,
        rate_obtained,
        out=target_reads_base.copy(),
        where=use_rate,
    )
    compensated_target_reads = numpy.where(
        numpy.isfinite(compensated_target_reads) & (compensated_target_reads > 0),
        compensated_target_reads,
        0.0,
    )
    sra_target_reads = compensated_target_reads.astype(int) + 1
    desired_end_2nds = start_2nds + sra_target_reads
    end_2nds = numpy.minimum(desired_end_2nds, safe_total_spots)
    overflow_reads = numpy.maximum(desired_end_2nds - safe_total_spots, 0.0)
    pooled_missing_bp = float((overflow_reads * safe_spot_lengths).sum())

    for _ in range(1000):
        if pooled_missing_bp <= 0:
            print('Enough read numbers were assigned for the 2nd round sequence extraction.', flush=True)
            break
        made_progress = False
        for i in range(end_2nds.shape[0]):
            if pooled_missing_bp <= 0:
                break
            if not valid_spot_length[i]:
                continue
            remaining_reads = int(max(0.0, safe_total_spots[i] - end_2nds[i]))
            if remaining_reads <= 0:
                continue
            alloc_reads = int(pooled_missing_bp / safe_spot_lengths[i])
            if alloc_reads <= 0:
                # If a full read cannot be represented for this SRA's spot length, try next SRA.
                continue
            alloc_reads = min(alloc_reads, remaining_reads)
            if alloc_reads <= 0:
                continue
            end_2nds[i] += alloc_reads
            pooled_missing_bp -= alloc_reads * safe_spot_lengths[i]
            made_progress = True
        if not made_progress:
            break
    if pooled_missing_bp > 0:
        print('Reached total spots in all SRAs.', flush=True)
    metadata.df.loc[:, 'spot_start_2nd'] = start_2nds.astype(int)
    metadata.df.loc[:, 'spot_end_2nd'] = end_2nds.astype(int)
    return metadata

def has_remaining_spots_after_first_round(metadata):
    df = metadata.df
    if ('total_spots' not in df.columns) or ('spot_end_1st' not in df.columns):
        return True
    total_spots = pandas.to_numeric(df.loc[:, 'total_spots'], errors='coerce').to_numpy(dtype=float)
    spot_end_1st = pandas.to_numeric(df.loc[:, 'spot_end_1st'], errors='coerce').to_numpy(dtype=float)
    valid_mask = numpy.isfinite(total_spots) & numpy.isfinite(spot_end_1st)
    if not valid_mask.any():
        return True
    remaining = numpy.maximum(total_spots[valid_mask] - spot_end_1st[valid_mask], 0.0)
    return bool(float(remaining.sum()) > 0.0)

def is_2nd_round_needed(rate_obtained_1st, tol):
    # tol is acceptable percentage loss relative to --max_bp.
    required_rate = 1 - (tol * 0.01)
    return rate_obtained_1st < required_rate

def print_read_stats(args, metadata, g, sra_stat=None, individual=False):
    rrna_filter = is_rrna_filter_enabled(args)
    if sra_stat is None:
        df = metadata.df
        print('Target size (--max_bp): {:,} bp'.format(g['max_bp']))
    else:
        df = metadata.df.loc[(metadata.df['run']==sra_stat['sra_id']),:]
        print('Individual target size: {:,} bp'.format(g['num_bp_per_sra']))
    print('Sum of extracted dumped reads: {:,} bp'.format(df['bp_dumped'].sum()))
    print('Sum of extracted rejected reads: {:,} bp'.format(df['bp_rejected'].sum()))
    print('Sum of extracted written reads: {:,} bp'.format(df['bp_written'].sum()))
    if args.fastp:
        print('Sum of fastp input reads: {:,} bp'.format(df['bp_fastp_in'].sum()))
        print('Sum of fastp output reads: {:,} bp'.format(df['bp_fastp_out'].sum()))
    if rrna_filter:
        print('Sum of SortMeRNA input reads: {:,} bp'.format(df['bp_rrna_in'].sum()))
        print('Sum of SortMeRNA non-rRNA reads: {:,} bp'.format(df['bp_rrna_out'].sum()))
    if individual:
        print('Individual SRA IDs:', ' '.join(df['run'].values))
        read_types = list()
        keys = list()
        read_types = read_types + ['extracted dumped reads', 'extracted rejected reads',
                                   'extracted written reads']
        keys = keys + ['bp_dumped', 'bp_rejected', 'bp_written']
        if args.fastp:
            read_types = read_types + ['fastp input reads', 'fastp output reads']
            keys = keys + ['bp_fastp_in', 'bp_fastp_out']
        if rrna_filter:
            read_types = read_types + ['SortMeRNA input reads', 'SortMeRNA non-rRNA reads']
            keys = keys + ['bp_rrna_in', 'bp_rrna_out']
        if len(read_types) > 0:
            for rt, key in zip(read_types, keys):
                values = ['{:,}'.format(s) for s in df[key].values]
                txt = ' '.join(values)
                print('Individual {} (bp): {}'.format(rt, txt))
    print('')


def _apply_batch_to_entrez_metadata(metadata, args):
    batch_value = getattr(args, 'batch', None)
    if batch_value is None:
        return metadata
    batch = int(batch_value)
    if batch <= 0:
        raise ValueError('--batch must be >= 1.')
    total_rows = int(metadata.df.shape[0])
    if total_rows == 0:
        return metadata
    print('--batch is specified with --id/--id_list. Processing one SRA per job.', flush=True)
    txt = 'This is {:,}th job. In total, {:,} jobs will be necessary for this metadata table.'
    print(txt.format(batch, total_rows), flush=True)
    if batch > total_rows:
        sys.stderr.write('--batch {} is too large. Exiting.\n'.format(batch_value))
        sys.exit(0)
    metadata.df = metadata.df.reset_index(drop=True).loc[[batch - 1], :].copy()
    return metadata


def getfastq_metadata(args):
    if (args.id is not None) and (args.id_list is not None):
        raise ValueError('--id and --id_list are mutually exclusive. Specify only one.')
    if args.id is not None:
        print('--id is specified. Downloading SRA metadata from Entrez.')
        Entrez.email = args.entrez_email
        sra_id = args.id
        search_term = getfastq_search_term(sra_id, args.entrez_additional_search_term)
        print('Entrez search term:', search_term)
        xml_root = getfastq_getxml(search_term=search_term)
        metadata = Metadata.from_xml(xml_root=xml_root)
        print('Filtering SRA entry with --layout:', args.layout)
        layout = get_layout(args, metadata)
        layout_series = metadata.df['lib_layout'].fillna('').astype(str).str.strip().str.lower()
        metadata.df = metadata.df.loc[(layout_series == layout), :]
        if args.sci_name is not None:
            print('Filtering SRA entry with --sci_name:', args.sci_name)
            sci_series = metadata.df['scientific_name'].fillna('').astype(str).str.strip()
            metadata.df = metadata.df.loc[(sci_series == str(args.sci_name).strip()), :]
    if args.id_list is not None:
        print('--id_list is specified. Downloading SRA metadata from Entrez.')
        Entrez.email = args.entrez_email
        id_list_path = os.path.realpath(args.id_list)
        if not os.path.exists(id_list_path):
            raise FileNotFoundError('SRA ID list file not found: {}'.format(id_list_path))
        if not os.path.isfile(id_list_path):
            raise IsADirectoryError('SRA ID list path exists but is not a file: {}'.format(id_list_path))
        sra_id_list = _read_sra_id_list(id_list_path)
        metadata_frames = _fetch_id_list_metadata_frames(args, sra_id_list)
        if len(metadata_frames)==0:
            print('No associated SRA is found with --id_list. Exiting.')
            sys.exit(1)
        metadata = Metadata.from_DataFrame(pandas.concat(metadata_frames, ignore_index=True))
        print('Filtering SRA entries with --layout:', args.layout)
        layout = get_layout(args, metadata)
        layout_series = metadata.df['lib_layout'].fillna('').astype(str).str.strip().str.lower()
        metadata.df = metadata.df.loc[(layout_series == layout), :]
        if args.sci_name is not None:
            print('Filtering SRA entries with --sci_name:', args.sci_name)
            sci_series = metadata.df['scientific_name'].fillna('').astype(str).str.strip()
            metadata.df = metadata.df.loc[(sci_series == str(args.sci_name).strip()), :]
    if (args.id is None)&(args.id_list is None):
        metadata = load_metadata(args)
    else:
        metadata = _apply_batch_to_entrez_metadata(metadata, args)
    metadata.df['total_bases'] = pandas.to_numeric(metadata.df['total_bases'], errors='coerce')
    metadata.df['spot_length'] = pandas.to_numeric(metadata.df['spot_length'], errors='coerce')
    return metadata

def is_getfastq_output_present(sra_stat, files=None):
    if files is None:
        files = list_run_dir_files(sra_stat['getfastq_sra_dir'])
    sra_stat = detect_layout_from_file(sra_stat, files=files)
    prefixes = [sra_stat['sra_id'], ]
    if sra_stat['layout'] == 'single':
        sub_exts = ['', ]
    elif sra_stat['layout'] == 'paired':
        sub_exts = ['_1', '_2']
    exts = ['.amalgkit.fastq.gz', ]
    is_output_present = True
    for prefix, sub_ext, ext in itertools.product(prefixes, sub_exts, exts):
        out_name1 = prefix + sub_ext + ext
        out_name2 = prefix + sub_ext + ext + '.safely_removed'
        out_path1 = os.path.join(sra_stat['getfastq_sra_dir'], out_name1)
        out_path2 = os.path.join(sra_stat['getfastq_sra_dir'], out_name2)
        is_out1 = out_name1 in files
        is_out2 = out_name2 in files
        if is_out1:
            print('getfastq output detected: {}'.format(out_path1))
        if is_out2:
            print('getfastq output detected: {}'.format(out_path2))
        is_output_present = bool(is_output_present and (is_out1 or is_out2))
    return is_output_present

def remove_experiment_without_run(metadata):
    num_all_run = metadata.df.shape[0]
    run_ids = metadata.df.loc[:, 'run'].fillna('').astype(str).str.strip()
    metadata.df['run'] = run_ids
    is_missing_run = (run_ids == '')
    num_missing_run = is_missing_run.sum()
    if (num_missing_run > 0):
        print('There are {} out of {} Experiments without Run ID. Removing.'.format(num_missing_run, num_all_run))
        metadata.df = metadata.df.loc[~is_missing_run, :]
    return metadata


def collect_valid_run_ids(run_values, unique=False):
    run_ids = []
    seen = set()
    for run_value in run_values:
        if pandas.isna(run_value):
            continue
        run_id = str(run_value).strip()
        if run_id == '':
            continue
        if unique and (run_id in seen):
            continue
        seen.add(run_id)
        run_ids.append(run_id)
    return run_ids

def initialize_columns(metadata, g):
    time_keys = ['time_start_1st', 'time_end_1st', 'time_start_2nd', 'time_end_2nd',]
    spot_keys = ['spot_start_1st', 'spot_end_1st', 'spot_start_2nd', 'spot_end_2nd',]
    keys = (['num_dumped', 'num_rejected', 'num_written', 'num_fastp_in', 'num_fastp_out',
             'num_rrna_in', 'num_rrna_out', 'bp_amalgkit',
             'bp_dumped', 'bp_rejected', 'bp_written', 'bp_fastp_in', 'bp_fastp_out',
             'bp_rrna_in', 'bp_rrna_out', 'bp_discarded', 'spot_length_amalgkit',
             'bp_still_available', 'bp_specified_for_extraction','rate_obtained','layout_amalgkit',
             'fastp_duplication_rate', 'fastp_insert_size_peak',]
            + time_keys + spot_keys)
    for key in keys:
        if key=='layout_amalgkit':
            metadata.df.loc[:,key] = ''
        elif key in ['rate_obtained', 'fastp_duplication_rate', 'fastp_insert_size_peak']:
            metadata.df.loc[:, key] = numpy.nan
        else:
            metadata.df.loc[:,key] = 0
    metadata.df.loc[:, 'bp_until_target_size'] = g['num_bp_per_sra']
    cols = ['total_spots','total_bases','size','nominal_length','nominal_sdev','spot_length']
    for col in cols:
        if any([ dtype in str(metadata.df[col].dtype) for dtype in ['str','object'] ]):
            metadata.df[col] = metadata.df.loc[:,col].astype(str).str.replace('^$', 'nan', regex=True).astype(float)
    for col in time_keys:
        metadata.df[col] = metadata.df[col].astype(float)
    return metadata

def sequence_extraction(args, sra_stat, metadata, g, start, end):
    sra_id = sra_stat['sra_id']
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_id))
    prev_num_written = metadata.df.at[ind_sra, 'num_written']
    prev_bp_dumped = metadata.df.at[ind_sra, 'bp_dumped']
    prev_bp_written = metadata.df.at[ind_sra, 'bp_written']
    metadata, sra_stat, run_file_state = run_fasterq_dump(
        sra_stat,
        args,
        metadata,
        start,
        end,
        return_file_state=True,
    )
    metadata, sra_stat, run_file_state = maybe_treat_paired_as_single(
        sra_stat=sra_stat,
        metadata=metadata,
        work_dir=sra_stat['getfastq_sra_dir'],
        file_state=run_file_state,
        return_file_state=True,
    )
    delta_num_written = metadata.df.at[ind_sra, 'num_written'] - prev_num_written
    delta_bp_dumped = metadata.df.at[ind_sra, 'bp_dumped'] - prev_bp_dumped
    delta_bp_written = metadata.df.at[ind_sra, 'bp_written'] - prev_bp_written
    latest_stage_counts = {
        'num_spots': max(0, int(delta_num_written)),
        'bp_total': max(0, int(delta_bp_written)),
    }
    latest_stage_source = 'fasterq'
    metadata.df.at[ind_sra, 'bp_discarded'] += max(0, delta_bp_dumped - delta_bp_written)
    metadata.df.at[ind_sra,'layout_amalgkit'] = sra_stat['layout']
    no_read_written = (metadata.df.at[ind_sra, 'num_written'] == 0)
    if no_read_written:
        write_getfastq_stats(sra_stat=sra_stat, metadata=metadata, output_dir=sra_stat['getfastq_sra_dir'])
        return metadata
    filter_order = get_filter_execution_order(args)
    for filter_name in filter_order:
        if filter_name == 'fastp':
            prev_num_fastp_out = metadata.df.at[ind_sra, 'num_fastp_out']
            prev_bp_fastp_in = metadata.df.at[ind_sra, 'bp_fastp_in']
            prev_bp_fastp_out = metadata.df.at[ind_sra, 'bp_fastp_out']
            metadata, run_file_state = run_fastp(
                sra_stat,
                args,
                sra_stat['getfastq_sra_dir'],
                metadata,
                file_state=run_file_state,
                return_file_state=True,
            )
            delta_num_out = metadata.df.at[ind_sra, 'num_fastp_out'] - prev_num_fastp_out
            delta_in = metadata.df.at[ind_sra, 'bp_fastp_in'] - prev_bp_fastp_in
            delta_out = metadata.df.at[ind_sra, 'bp_fastp_out'] - prev_bp_fastp_out
            metadata.df.at[ind_sra, 'bp_discarded'] += max(0, delta_in - delta_out)
            latest_stage_counts = {
                'num_spots': max(0, int(delta_num_out)),
                'bp_total': max(0, int(delta_out)),
            }
            latest_stage_source = 'fastp'
            continue
        if filter_name == 'rrna':
            prev_rrna_in = metadata.df.at[ind_sra, 'bp_rrna_in']
            prev_rrna_out = metadata.df.at[ind_sra, 'bp_rrna_out']
            prev_num_rrna_out = metadata.df.at[ind_sra, 'num_rrna_out']
            rrna_known_input_counts = latest_stage_counts if (latest_stage_source == 'fastp') else None
            metadata, run_file_state = run_sortmerna_rrna_filter(
                sra_stat=sra_stat,
                args=args,
                output_dir=sra_stat['getfastq_sra_dir'],
                metadata=metadata,
                file_state=run_file_state,
                known_input_counts=rrna_known_input_counts,
                return_file_state=True,
            )
            delta_num_out = metadata.df.at[ind_sra, 'num_rrna_out'] - prev_num_rrna_out
            delta_in = metadata.df.at[ind_sra, 'bp_rrna_in'] - prev_rrna_in
            delta_out = metadata.df.at[ind_sra, 'bp_rrna_out'] - prev_rrna_out
            metadata.df.at[ind_sra, 'bp_discarded'] += max(0, delta_in - delta_out)
            latest_stage_counts = {
                'num_spots': max(0, int(delta_num_out)),
                'bp_total': max(0, int(delta_out)),
            }
            latest_stage_source = 'rrna'
            continue
        raise ValueError('Unsupported filter name in execution order: {}'.format(filter_name))
    if args.read_name == 'trinity':
        run_file_state = rename_reads(
            sra_stat,
            args,
            sra_stat['getfastq_sra_dir'],
            file_state=run_file_state,
            return_file_state=True,
        )
    inext = get_or_detect_intermediate_extension(sra_stat, work_dir=sra_stat['getfastq_sra_dir'])
    outext = '.amalgkit.fastq.gz'
    rename_fastq(sra_stat, sra_stat['getfastq_sra_dir'], inext, outext)
    metadata.df.at[ind_sra,'bp_still_available'] = sra_stat['spot_length'] * (sra_stat['total_spot'] - end)
    bp_specified_for_extraction = sra_stat['spot_length'] * (end - start)
    metadata.df.at[ind_sra, 'bp_specified_for_extraction'] += bp_specified_for_extraction
    bp_source_col = resolve_bp_amalgkit_source_column(args)
    metadata.df.at[ind_sra, 'bp_amalgkit'] = metadata.df.at[ind_sra, bp_source_col]
    metadata.df.at[ind_sra, 'rate_obtained'] = metadata.df.at[ind_sra, 'bp_amalgkit'] / g['num_bp_per_sra']
    metadata.df.at[ind_sra, 'bp_until_target_size'] -= metadata.df.at[ind_sra, 'bp_amalgkit']
    write_getfastq_stats(sra_stat=sra_stat, metadata=metadata, output_dir=sra_stat['getfastq_sra_dir'])
    return metadata

def sequence_extraction_1st_round(args, sra_stat, metadata, g):
    offset = 10000
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_stat['sra_id']))
    metadata.df.at[ind_sra, 'time_start_1st'] = time.time()
    start, end = get_range(sra_stat, offset, g['total_sra_bp'], g['max_bp'])
    metadata.df.at[ind_sra, 'spot_length_amalgkit'] = sra_stat['spot_length']
    metadata.df.at[ind_sra,'spot_start_1st'] = start
    metadata.df.at[ind_sra,'spot_end_1st'] = end
    metadata = sequence_extraction(args, sra_stat, metadata, g, start, end)
    txt = 'Time elapsed for 1st-round sequence extraction: {}, {:,.1f} sec'
    print(txt.format(sra_stat['sra_id'], int(time.time() - g['start_time'])))
    print('\n--- getfastq 1st-round sequence generation report ---')
    print_read_stats(args, metadata, g, sra_stat, individual=False)
    txt = '{:.2f}% of reads were obtained in the 1st-round sequence generation: {:,} bp out of the individual target amount of {:,} bp'
    percent_obtained = metadata.df.at[ind_sra,'rate_obtained']*100
    bp_amalgkit = metadata.df.at[ind_sra, 'bp_amalgkit']
    print(txt.format(percent_obtained, bp_amalgkit, g['num_bp_per_sra']), flush=True)
    metadata.df.at[ind_sra, 'time_end_1st'] = time.time()
    elapsed_time = metadata.df.at[ind_sra, 'time_end_1st'] - metadata.df.at[ind_sra, 'time_start_1st']
    txt = 'Time elapsed for 1st-round sequence extraction: {}, {:,.1f} sec'
    print(txt.format(sra_stat['sra_id'], elapsed_time))
    print('')
    return metadata

def sequence_extraction_2nd_round(args, sra_stat, metadata, g):
    print('Starting the 2nd-round sequence extraction.')
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_stat['sra_id']))
    metadata.df.at[ind_sra,'time_start_2nd'] = time.time()
    ext_main = '.amalgkit.fastq.gz'
    ext_1st_tmp = '.amalgkit_1st.fastq.gz'
    print('')
    sra_id = sra_stat['sra_id']
    layout = sra_stat['layout']
    start = metadata.df.at[ind_sra,'spot_start_2nd']
    end = metadata.df.at[ind_sra,'spot_end_2nd']
    if (start >= end):
        txt = '{}: All spots have been extracted in the 1st trial. Cancelling the 2nd trial. start={:,}, end={:,}'
        print(txt.format(sra_id, start, end))
        return metadata
    no_read_in_1st = (metadata.df.at[ind_sra, 'bp_written'] == 0)
    if no_read_in_1st:
        print('No read was extracted in 1st round. Skipping 2nd round: {}'.format(sra_id))
        return metadata
    else:
        rename_fastq(sra_stat, sra_stat['getfastq_sra_dir'], inext=ext_main, outext=ext_1st_tmp)
    metadata = sequence_extraction(args, sra_stat, metadata, g, start, end)
    if (layout == 'single'):
        subexts = ['']
    elif (layout == 'paired'):
        subexts = ['_1', '_2']
    try:
        requested_workers = int(getattr(args, 'threads', 1))
    except (TypeError, ValueError):
        requested_workers = 1
    requested_workers = max(1, requested_workers)
    worker_count = min(len(subexts), requested_workers)

    def merge_single_subext(subext):
        added_path = os.path.join(sra_stat['getfastq_sra_dir'], sra_id + subext + ext_1st_tmp)
        adding_path = os.path.join(sra_stat['getfastq_sra_dir'], sra_id + subext + ext_main)
        if os.path.exists(added_path) and (not os.path.isfile(added_path)):
            raise IsADirectoryError('Dumped fastq path exists but is not a file: {}'.format(added_path))
        if os.path.exists(adding_path) and (not os.path.isfile(adding_path)):
            raise IsADirectoryError('Dumped fastq path exists but is not a file: {}'.format(adding_path))
        if not os.path.exists(added_path):
            raise FileNotFoundError('Dumped fastq not found: ' + added_path)
        if not os.path.exists(adding_path):
            raise FileNotFoundError('Dumped fastq not found: ' + adding_path)
        append_file_binary(adding_path, added_path)
        os.remove(adding_path)
        os.rename(added_path, adding_path)
        return subext

    if worker_count <= 1:
        for subext in subexts:
            merge_single_subext(subext)
    else:
        print(
            'Merging 2nd-round FASTQ chunks for {} using {:,} parallel worker(s).'.format(
                sra_id,
                worker_count,
            ),
            flush=True,
        )
        _, failures = run_tasks_with_optional_threads(
            task_items=subexts,
            task_fn=merge_single_subext,
            max_workers=worker_count,
        )
        if failures:
            details = '; '.join(['{}: {}'.format(subext, exc) for subext, exc in failures])
            raise RuntimeError('2nd-round FASTQ merge failed for {}. {}'.format(sra_id, details))
    metadata.df.at[ind_sra, 'time_end_2nd'] = time.time()
    elapsed_time = metadata.df.at[ind_sra, 'time_end_2nd'] - metadata.df.at[ind_sra, 'time_start_2nd']
    txt = 'Time elapsed for 2nd-round sequence extraction: {}, {:,} sec'
    print(txt.format(sra_stat['sra_id'], elapsed_time))
    print('')
    return metadata

def sequence_extraction_private(metadata, sra_stat, args):
    ind_sra = sra_stat.get('metadata_idx', get_metadata_row_index_by_run(metadata, sra_stat['sra_id']))
    for col in ['read1_path','read2_path']:
        path_from_raw = metadata.df.at[ind_sra, col]
        if pandas.isna(path_from_raw):
            sys.stderr.write('Private fastq file path is missing in metadata column "{}".\n'.format(col))
            continue
        path_from = str(path_from_raw).strip()
        if path_from == '' or path_from.lower() == 'nan':
            sys.stderr.write('Private fastq file path is missing in metadata column "{}".\n'.format(col))
            continue
        path_to = os.path.join(sra_stat['getfastq_sra_dir'], os.path.basename(path_from))
        path_to = path_to.replace('.fq', '.fastq')
        if not path_to.endswith('.gz'):
            path_to = path_to+'.gz' # .gz is necessary even if the original file is not compressed.
        if os.path.isfile(path_from):
            if os.path.lexists(path_to):
                if os.path.isdir(path_to) and (not os.path.islink(path_to)):
                    raise IsADirectoryError(
                        'Private output path exists but is not a file/symlink: {}'.format(path_to)
                    )
                os.remove(path_to)
            os.symlink(src=path_from, dst=path_to)
        elif os.path.exists(path_from):
            sys.stderr.write('Private fastq path exists but is not a file: {}\n'.format(path_from))
        else:
            sys.stderr.write('Private fastq file not found: {}\n'.format(path_from))
    set_current_intermediate_extension(sra_stat, '.fastq.gz')
    latest_stage_counts = None
    latest_stage_source = None
    for filter_name in get_filter_execution_order(args):
        if filter_name == 'fastp':
            prev_num_fastp_out = metadata.df.at[ind_sra, 'num_fastp_out']
            prev_bp_fastp_out = metadata.df.at[ind_sra, 'bp_fastp_out']
            metadata = run_fastp(sra_stat, args, sra_stat['getfastq_sra_dir'], metadata)
            delta_num_out = metadata.df.at[ind_sra, 'num_fastp_out'] - prev_num_fastp_out
            delta_bp_out = metadata.df.at[ind_sra, 'bp_fastp_out'] - prev_bp_fastp_out
            latest_stage_counts = {
                'num_spots': max(0, int(delta_num_out)),
                'bp_total': max(0, int(delta_bp_out)),
            }
            latest_stage_source = 'fastp'
            continue
        if filter_name == 'rrna':
            rrna_known_input_counts = latest_stage_counts if (latest_stage_source == 'fastp') else None
            metadata = run_sortmerna_rrna_filter(
                sra_stat=sra_stat,
                args=args,
                output_dir=sra_stat['getfastq_sra_dir'],
                metadata=metadata,
                known_input_counts=rrna_known_input_counts,
            )
            latest_stage_source = 'rrna'
            continue
        raise ValueError('Unsupported filter name in execution order: {}'.format(filter_name))
    inext = get_or_detect_intermediate_extension(sra_stat, work_dir=sra_stat['getfastq_sra_dir'])
    if (inext=='no_extension_found')&(sra_stat['layout']=='paired'):
        raise ValueError('Paired-end file names may be invalid. They should contain _1 and _2 to indicate a pair: {}'.format(sra_stat['sra_id']))
    outext = '.amalgkit.fastq.gz'
    rename_fastq(sra_stat, sra_stat['getfastq_sra_dir'], inext, outext)
    write_getfastq_stats(sra_stat=sra_stat, metadata=metadata, output_dir=sra_stat['getfastq_sra_dir'])
    return metadata

def check_metadata_validity(metadata):
    required_columns = ['run', 'total_bases', 'total_spots', 'spot_length']
    missing_columns = [col for col in required_columns if col not in metadata.df.columns]
    if len(missing_columns) > 0:
        raise ValueError(
            'Missing required metadata column(s) for getfastq: {}'.format(', '.join(missing_columns))
        )
    if metadata.df.shape[0] == 0:
        raise ValueError('No SRA entry found. Make sure whether --id or --id_list is compatible with --sci_name and --layout.')
    run_ids = metadata.df.loc[:, 'run'].fillna('').astype(str).str.strip()
    metadata.df['run'] = run_ids
    is_missing_run = (run_ids == '')
    if is_missing_run.any():
        raise ValueError('Missing Run ID(s) were detected in metadata.')
    duplicate_mask = run_ids.duplicated(keep=False)
    if duplicate_mask.any():
        duplicated_runs = run_ids.loc[duplicate_mask].drop_duplicates().tolist()
        raise ValueError('Duplicate Run ID(s) were detected in metadata: {}'.format(', '.join(duplicated_runs)))
    total_bases = pandas.to_numeric(metadata.df.loc[:, 'total_bases'], errors='coerce')
    is_total_bases_invalid = (~numpy.isfinite(total_bases.to_numpy(dtype=float))) | (total_bases <= 0)
    if is_total_bases_invalid.any():
        txt = 'Invalid value(s) of total_bases were detected in {}. Filling a placeholder value 999,999,999,999\n'
        run_ids = metadata.df.loc[is_total_bases_invalid, 'run'].astype(str).tolist()
        sys.stderr.write(txt.format(', '.join(run_ids)))
        total_bases.loc[is_total_bases_invalid] = 999999999999
    metadata.df['total_bases'] = total_bases.astype(int)

    total_spots = pandas.to_numeric(metadata.df.loc[:, 'total_spots'], errors='coerce')
    is_total_spots_invalid = (~numpy.isfinite(total_spots.to_numpy(dtype=float))) | (total_spots <= 0)
    if is_total_spots_invalid.any():
        spot_length = pandas.to_numeric(metadata.df.loc[is_total_spots_invalid, 'spot_length'], errors='coerce')
        new_values = metadata.df.loc[is_total_spots_invalid, 'total_bases'] / spot_length
        txt = 'Invalid value(s) of total_spots were detected in {}. Filling a placeholder value 999,999,999,999\n'
        run_ids = metadata.df.loc[is_total_spots_invalid, 'run'].astype(str).tolist()
        sys.stderr.write(txt.format(', '.join(run_ids)))
        invalid_new_values = (~numpy.isfinite(new_values.to_numpy(dtype=float))) | (new_values <= 0)
        new_values.loc[invalid_new_values] = 999999999999  # https://github.com/kfuku52/amalgkit/issues/110
        total_spots.loc[is_total_spots_invalid] = new_values.astype(int)
    metadata.df['total_spots'] = total_spots.astype(int)
    for run_id, total_bases in zip(metadata.df['run'].tolist(), metadata.df['total_bases'].tolist()):
        txt = 'Individual SRA size of {}: {:,} bp'
        print(txt.format(run_id, int(total_bases)))
    return metadata

def initialize_global_params(args, metadata):
    g = dict()
    g['start_time'] = time.time()
    max_bp_raw = str(args.max_bp).replace(',', '').strip()
    try:
        g['max_bp'] = int(max_bp_raw)
    except ValueError as exc:
        raise ValueError('--max_bp must be a positive integer: {}'.format(args.max_bp)) from exc
    if g['max_bp'] <= 0:
        raise ValueError('--max_bp must be > 0.')
    g['num_sra'] = metadata.df.shape[0]
    if g['num_sra'] <= 0:
        raise ValueError('No SRA entries were found in metadata.')
    g['num_bp_per_sra'] = int(g['max_bp'] / g['num_sra'])
    if g['num_bp_per_sra'] <= 0:
        raise ValueError(
            '--max_bp ({}) is too small for {:,} SRA runs. '
            'Increase --max_bp to at least {:,}.'.format(
                g['max_bp'],
                g['num_sra'],
                g['num_sra'],
            )
        )
    g['total_sra_bp'] = metadata.df.loc[:,'total_bases'].sum()
    print('Number of SRAs to be processed: {:,}'.format(g['num_sra']))
    print('Total target size (--max_bp): {:,} bp'.format(g['max_bp']))
    print('The sum of SRA sizes: {:,} bp'.format(g['total_sra_bp']))
    print('Target size per SRA: {:,} bp'.format(g['num_bp_per_sra']))
    return g

def process_getfastq_run(args, row_index, sra_id, run_row_df, g):
    print('')
    print('Processing SRA ID: {}'.format(sra_id))
    run_metadata = Metadata.from_DataFrame(run_row_df)
    sra_stat = get_sra_stat(sra_id, run_metadata, g['num_bp_per_sra'])
    sra_stat['getfastq_sra_dir'] = get_getfastq_run_dir(args, sra_id)
    run_dir_files = list_run_dir_files(sra_stat['getfastq_sra_dir'])
    is_output_present = is_getfastq_output_present(sra_stat, files=run_dir_files)
    if is_output_present and (not args.redo):
        txt = 'Output file(s) detected. Skipping {}. Set "--redo yes" for reanalysis.'
        print(txt.format(sra_id), flush=True)
        return {
            'row_index': row_index,
            'sra_id': sra_id,
            'row': run_metadata.df.iloc[0].copy(),
            'flag_any_output_file_present': True,
            'flag_private_file': False,
            'getfastq_sra_dir': sra_stat['getfastq_sra_dir'],
        }
    remove_old_intermediate_files(sra_id=sra_id, work_dir=sra_stat['getfastq_sra_dir'], files=run_dir_files)
    print('Library layout:', sra_stat['layout'])
    print('Number of reads:', "{:,}".format(sra_stat['total_spot']))
    print('Single/Paired read length:', sra_stat['spot_length'], 'bp')
    print('Total bases:', "{:,}".format(int(run_metadata.df.at[0, 'total_bases'])), 'bp')
    flag_private_file = False
    if 'private_file' in run_metadata.df.columns:
        if run_metadata.df.at[0, 'private_file'] == 'yes':
            print('Processing {} as private data. --max_bp is disabled.'.format(sra_id), flush=True)
            flag_private_file = True
            sequence_extraction_private(run_metadata, sra_stat, args)
    if not flag_private_file:
        print('Processing {} as publicly available data from SRA.'.format(sra_id), flush=True)
        download_sra(run_metadata, sra_stat, args, sra_stat['getfastq_sra_dir'], overwrite=False)
        run_metadata = sequence_extraction_1st_round(args, sra_stat, run_metadata, g)
    return {
        'row_index': row_index,
        'sra_id': sra_id,
        'row': run_metadata.df.iloc[0].copy(),
        'flag_any_output_file_present': bool(is_output_present),
        'flag_private_file': bool(flag_private_file),
        'getfastq_sra_dir': sra_stat['getfastq_sra_dir'],
    }


def run_first_round_getfastq(args, metadata, run_rows, g, jobs):
    run_results_by_id = dict()
    if (jobs == 1) or (len(run_rows) <= 1):
        for row_index, sra_id in run_rows:
            run_results_by_id[sra_id] = process_getfastq_run(
                args=args,
                row_index=row_index,
                sra_id=sra_id,
                run_row_df=metadata.df.loc[[row_index], :].copy(),
                g=g,
            )
        return run_results_by_id

    max_workers = min(jobs, len(run_rows))
    print('Running getfastq for {:,} SRA runs with {:,} parallel jobs.'.format(len(run_rows), max_workers), flush=True)
    results_by_row, failures = run_tasks_with_optional_threads(
        task_items=run_rows,
        task_fn=lambda run_row: process_getfastq_run(
            args,
            run_row[0],
            run_row[1],
            metadata.df.loc[[run_row[0]], :].copy(),
            g,
        ),
        max_workers=max_workers,
    )
    for run_row, run_result in results_by_row.items():
        run_results_by_id[run_row[1]] = run_result
    if failures:
        details = '; '.join(['{}: {}'.format(run_row[1], err) for run_row, err in failures])
        raise RuntimeError('getfastq failed for {}/{} SRA runs. {}'.format(len(failures), len(run_rows), details))
    return run_results_by_id


def apply_first_round_getfastq_results(metadata, run_rows, run_results_by_id):
    flag_private_file = False
    flag_any_output_file_present = False
    last_getfastq_sra_dir = None
    if len(run_rows) == 0:
        return metadata, flag_private_file, flag_any_output_file_present, last_getfastq_sra_dir
    first_row_series = run_results_by_id[run_rows[0][1]]['row']
    common_cols = [col for col in metadata.df.columns if col in first_row_series.index]
    for row_index, sra_id in run_rows:
        run_result = run_results_by_id[sra_id]
        row_series = run_result['row']
        metadata.df.loc[row_index, common_cols] = row_series.loc[common_cols].to_numpy()
        flag_any_output_file_present |= run_result['flag_any_output_file_present']
        flag_private_file |= run_result['flag_private_file']
        last_getfastq_sra_dir = run_result['getfastq_sra_dir']
    return metadata, flag_private_file, flag_any_output_file_present, last_getfastq_sra_dir


def maybe_run_getfastq_second_round(args, metadata, run_rows, g, flag_private_file, flag_any_output_file_present):
    if flag_private_file or flag_any_output_file_present:
        return metadata
    g['rate_obtained_1st'] = metadata.df.loc[:, 'bp_amalgkit'].sum() / g['max_bp']
    if is_2nd_round_needed(g['rate_obtained_1st'], args.tol):
        txt = 'Only {:,.2f}% ({:,}/{:,}) of the target size (--max_bp) was obtained in the 1st round. Proceeding to the 2nd round read extraction.'
        print(txt.format(g['rate_obtained_1st'] * 100, metadata.df.loc[:, 'bp_amalgkit'].sum(), g['max_bp']), flush=True)
        if not has_remaining_spots_after_first_round(metadata):
            print('All spots were already extracted in the 1st round. Skipping the 2nd round.', flush=True)
        else:
            metadata = calc_2nd_ranges(metadata)
            for _, sra_id in run_rows:
                sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
                sra_stat['getfastq_sra_dir'] = get_getfastq_run_dir(args, sra_id)
                metadata = sequence_extraction_2nd_round(args, sra_stat, metadata, g)
    else:
        print('Sufficient data were obtained in the 1st-round sequence extraction. Proceeding without the 2nd round.')
    g['rate_obtained_2nd'] = metadata.df.loc[:, 'bp_amalgkit'].sum() / g['max_bp']
    txt = '2nd round read extraction improved % bp from {:,.2f}% to {:,.2f}%'
    print(txt.format(g['rate_obtained_1st'] * 100, g['rate_obtained_2nd'] * 100), flush=True)
    return metadata


def run_getfastq_postprocessing(args, metadata, last_getfastq_sra_dir, flag_any_output_file_present, g):
    if args.remove_sra:
        remove_sra_files(metadata=metadata, amalgkit_out_dir=args.out_dir)
    else:
        if last_getfastq_sra_dir is not None:
            print('SRA files not removed: {}'.format(last_getfastq_sra_dir))
        else:
            print('SRA files not removed.')
    if not flag_any_output_file_present:
        print('')
        print('\n--- getfastq final report ---')
        print_read_stats(args, metadata, g, sra_stat=None, individual=True)


def getfastq_main(args):
    threads, jobs, _ = resolve_thread_worker_allocation(
        requested_threads=getattr(args, 'threads', 'auto'),
        requested_workers=getattr(args, 'internal_jobs', 'auto'),
        internal_cpu_budget=getattr(args, 'internal_cpu_budget', 'auto'),
        worker_option_name='internal_jobs',
        context='getfastq:',
        disable_workers=(getattr(args, 'batch', None) is not None),
    )
    args.threads = threads
    args.internal_jobs = jobs
    check_getfastq_dependency(args)
    metadata = getfastq_metadata(args)
    metadata = remove_experiment_without_run(metadata)
    metadata = check_metadata_validity(metadata)
    g = initialize_global_params(args, metadata)
    metadata = initialize_columns(metadata, g)
    run_rows = list(zip(metadata.df.index.tolist(), metadata.df['run'].tolist()))
    run_results_by_id = run_first_round_getfastq(args=args, metadata=metadata, run_rows=run_rows, g=g, jobs=jobs)
    metadata, flag_private_file, flag_any_output_file_present, last_getfastq_sra_dir = apply_first_round_getfastq_results(
        metadata=metadata,
        run_rows=run_rows,
        run_results_by_id=run_results_by_id,
    )
    metadata = maybe_run_getfastq_second_round(
        args=args,
        metadata=metadata,
        run_rows=run_rows,
        g=g,
        flag_private_file=flag_private_file,
        flag_any_output_file_present=flag_any_output_file_present,
    )
    run_getfastq_postprocessing(
        args=args,
        metadata=metadata,
        last_getfastq_sra_dir=last_getfastq_sra_dir,
        flag_any_output_file_present=flag_any_output_file_present,
        g=g,
    )
