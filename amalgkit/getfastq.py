from Bio import Entrez
import itertools
import numpy
import gzip
import pandas
import shlex

from amalgkit.util import *

import glob
import xml.etree.ElementTree as ET
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from urllib.error import HTTPError

IDENTICAL_PAIRED_RATIO_THRESHOLD = 0.99
IDENTICAL_PAIRED_CHECKED_READS = 2000

def append_file_binary(src_path, dst_path, chunk_size=1024 * 1024):
    with open(src_path, 'rb') as src_handle, open(dst_path, 'ab') as dst_handle:
        shutil.copyfileobj(src_handle, dst_handle, length=chunk_size)

def getfastq_search_term(ncbi_id, additional_search_term=None):
    # https://www.ncbi.nlm.nih.gov/books/NBK49540/
    if additional_search_term is None:
        search_term = ncbi_id
    else:
        search_term = ncbi_id + ' AND ' + additional_search_term
    return search_term

def getfastq_getxml(search_term, retmax=1000):
    entrez_db = 'sra'
    try:
        sra_handle = Entrez.esearch(db=entrez_db, term=search_term, retmax=10000000)
    except HTTPError as e:
        print(e, '- Trying Entrez.esearch() again...')
        sra_handle = Entrez.esearch(db=entrez_db, term=search_term, retmax=10000000)
    sra_record = Entrez.read(sra_handle)
    record_ids = sra_record["IdList"]
    num_record = len(record_ids)
    print('Number of SRA records:', num_record)
    if num_record == 0:
        return ET.Element('EXPERIMENT_PACKAGE_SET')
    root = None
    for start in range(0, num_record, retmax):
        end = min(start + retmax, num_record)
        print('processing SRA records:', start, '-', end - 1, flush=True)
        try:
            handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
        except HTTPError as e:
            print(e, '- Trying Entrez.efetch() again...')
            handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
        chunk = ET.parse(handle).getroot()
        if root is None:
            root = chunk
        else:
            root.append(chunk)
    xml_string = ET.tostring(root, encoding='unicode')
    for line in xml_string.split('\n'):
        if '<Error>' in line:
            print(line)
            raise Exception('<Error> found in the xml. Search term: ' + search_term)
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

def concat_fastq(args, metadata, output_dir, g):
    layout = get_layout(args, metadata)
    inext = '.amalgkit.fastq.gz'
    infiles = list()
    for sra_id in metadata.df.loc[:, 'run']:
        infiles.append([f for f in os.listdir(output_dir) if (f.endswith(inext)) & (f.startswith(sra_id))])
    infiles = [item for sublist in infiles for item in sublist]
    num_inext_files = len(infiles)
    if (layout == 'single') & (num_inext_files == 1):
        print('Only 1', inext, 'file was detected. No concatenation will happen.', flush=True)
        if args.id is not None:
            outfile = args.id + infiles[0]
        elif args.id_list is not None:
            outfile = os.path.basename(args.id_list) + infiles[0]
        if infiles[0] != outfile:
            print('Replacing ID in the output file name:', infiles[0], outfile)
            infile_path = os.path.join(output_dir, infiles[0])
            outfile_path = os.path.join(output_dir, outfile)
            os.rename(infile_path, outfile_path)
        return None
    elif (layout == 'paired') & (num_inext_files == 2):
        print('Only 1 pair of', inext, 'files were detected. No concatenation will happen.', flush=True)
        for infile in infiles:
            if args.id is not None:
                outfile = args.id + re.sub('.*(_[1-2])', r'\g<1>', infile)
            elif args.id_list is not None:
                outfile = os.path.basename(args.id_list) + re.sub('.*(_[1-2])', r'\g<1>', infile)
            if infile != outfile:
                print('Replacing ID in the output file name:', infile, outfile)
                infile_path = os.path.join(output_dir, infile)
                outfile_path = os.path.join(output_dir, outfile)
                os.rename(infile_path, outfile_path)
        return None
    else:
        print('Concatenating files with the extension:', inext)
        outext = '.amalgkit.fastq.gz'
        if layout == 'single':
            subexts = ['', ]
        elif layout == 'paired':
            subexts = ['_1', '_2', ]
        for subext in subexts:
            infiles = metadata.df['run'].replace('$', subext + inext, regex=True)
            if args.id is not None:
                outfile_path = os.path.join(output_dir, args.id + subext + outext)
            elif args.id_list is not None:
                outfile_path = os.path.join(output_dir, os.path.basename(args.id_list) + subext + outext)
            if os.path.exists(outfile_path):
                os.remove(outfile_path)
            for infile in infiles:
                infile_path = os.path.join(output_dir, infile)
                assert os.path.exists(infile_path), 'Dumped fastq not found: ' + infile_path
                print('Concatenated file:', infile_path, flush=True)
                append_file_binary(infile_path, outfile_path)
            print('')
        if args.remove_tmp:
            for i in metadata.df.index:
                sra_id = metadata.df.loc[i, 'run']
                sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
                ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
                remove_intermediate_files(sra_stat, ext=ext, work_dir=output_dir)
        return None

def remove_sra_files(metadata, amalgkit_out_dir):
    print('Starting SRA file removal.', flush=True)
    for sra_id in metadata.df['run']:
        sra_pattern = os.path.join(os.path.realpath(amalgkit_out_dir), 'getfastq', sra_id, sra_id + '.sra*')
        path_downloaded_sras = glob.glob(sra_pattern)
        if len(path_downloaded_sras) > 0:
            for path_downloaded_sra in path_downloaded_sras:
                print('Deleting SRA file: {}'.format(path_downloaded_sra))
                os.remove(path_downloaded_sra)
        else:
            print('SRA file not found. Pattern searched: {}'.format(sra_pattern))
    print('')

def get_layout(args, metadata):
    if args.layout == 'auto':
        layouts = metadata.df['lib_layout'].unique().tolist()
        if (len(layouts) != 1):
            print('Detected multiple layouts in the metadata:', layouts)
        layout = 'paired' if 'paired' in layouts else 'single'
    else:
        layout = args.layout
    return layout

def remove_old_intermediate_files(sra_id, work_dir):
    old_files = os.listdir(work_dir)
    files = [f for f in old_files if
             (f.startswith(sra_id)) & (not f.endswith('.sra')) & (os.path.isfile(os.path.join(work_dir, f)))]
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

def download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
    path_downloaded_sra = os.path.join(work_dir, sra_stat['sra_id'] + '.sra')
    if os.path.exists(path_downloaded_sra):
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
        sra_sources = dict()
        sra_id = sra_stat['sra_id']
        is_sra = (metadata.df['run']==sra_stat['sra_id'])
        if args.aws:
            aws_link = metadata.df.loc[is_sra,'AWS_Link'].values[0]
            if aws_link=='':
                sys.stderr.write('AWS_Link is empty and will be skipped.\n')
            else:
                sra_sources['AWS'] = aws_link
        if args.gcp:
            gcp_link = metadata.df.loc[is_sra,'GCP_Link'].values[0]
            if gcp_link=='':
                sys.stderr.write('GCP_Link is empty and will be skipped.\n')
            else:
                sra_sources['GCP'] = gcp_link
        if args.ncbi:
            ncbi_link = metadata.df.loc[is_sra,'NCBI_Link'].values[0]
            if ncbi_link=='':
                sys.stderr.write('NCBI_Link is empty and will be skipped.\n')
            else:
                sra_sources['NCBI'] = ncbi_link
        if len(sra_sources)==0:
            print('No source URL is available. Check whether --aws, --gcp, and --ncbi are properly set.')
        is_sra_download_completed = False
        for sra_source_name in sra_sources.keys():
            source_url_original = str(sra_sources[sra_source_name])
            source_url = normalize_url_for_urllib(
                source_name=sra_source_name,
                source_url=source_url_original,
                gcp_project=getattr(args, 'gcp_project', ''),
            )
            print("Trying to fetch {} from {}: {}".format(sra_id, sra_source_name, source_url_original))
            if source_url != source_url_original:
                print("Converted {} URL for urllib: {}".format(sra_source_name, source_url))
            if source_url=='nan':
                sys.stderr.write("Skipping. No URL for {}.\n".format(sra_source_name))
                continue
            scheme = urllib.parse.urlparse(source_url).scheme.lower()
            if scheme not in ('http', 'https', 'ftp'):
                msg = 'Skipping {} download source due to unsupported URL scheme for urllib: {}\n'
                sys.stderr.write(msg.format(sra_source_name, scheme if scheme else '(none)'))
                continue
            try:
                urllib.request.urlretrieve(source_url, path_downloaded_sra)
                if os.path.exists(path_downloaded_sra):
                    is_sra_download_completed = True
                    print('SRA file was downloaded with urllib.request from {}'.format(sra_source_name), flush=True)
                    break
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
        if not is_sra_download_completed:
            sys.stderr.write("Exhausted all sources of download.\n")
        else:
            assert os.path.exists(path_downloaded_sra), 'SRA file download failed: ' + sra_stat['sra_id']
            return
    err_txt = 'SRA file download failed for {}. Expected PATH: {}. '
    err_txt += 'Cloud URL download sources were exhausted and prefetch fallback is obsolete.'
    assert os.path.exists(path_downloaded_sra), err_txt.format(sra_stat['sra_id'], path_downloaded_sra)

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
    fasterq_dump_exe = getattr(args, 'fasterq_dump_exe', 'fasterq-dump')
    test_fqd = subprocess.run([fasterq_dump_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert (test_fqd.returncode == 0), "fasterq-dump PATH cannot be found: " + fasterq_dump_exe
    if args.fastp:
        test_fp = subprocess.run([args.fastp_exe, '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_fp.returncode == 0), "fastp PATH cannot be found: " + args.fastp_exe
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

def count_fastq_records(path_fastq):
    num_lines = 0
    open_func = gzip.open if path_fastq.endswith('.gz') else open
    with open_func(path_fastq, 'rt') as f:
        for _ in f:
            num_lines += 1
    if (num_lines % 4) != 0:
        txt = 'FASTQ line count is not divisible by 4 and may be truncated: {}\n'
        sys.stderr.write(txt.format(path_fastq))
    return int(num_lines / 4)

def estimate_num_written_spots_from_fastq(sra_stat):
    work_dir = sra_stat['getfastq_sra_dir']
    sra_id = sra_stat['sra_id']
    def detect_fastq_path(prefix):
        for ext in ['.fastq.gz', '.fastq']:
            path = os.path.join(work_dir, prefix + ext)
            if os.path.exists(path):
                return path
        return None
    single_path = detect_fastq_path(sra_id)
    pair1_path = detect_fastq_path(sra_id + '_1')
    pair2_path = detect_fastq_path(sra_id + '_2')
    num_single = count_fastq_records(single_path) if single_path is not None else 0
    num_pair1 = count_fastq_records(pair1_path) if pair1_path is not None else 0
    num_pair2 = count_fastq_records(pair2_path) if pair2_path is not None else 0
    return max(num_pair1, num_pair2) + num_single

def trim_fastq_by_spot_range(path_fastq, start, end):
    if not os.path.exists(path_fastq):
        return 0
    start = max(1, int(start))
    end = max(0, int(end))
    tmp_path = path_fastq + '.trimtmp'
    kept = 0
    spot_index = 0
    with open(path_fastq, 'rt') as fin, open(tmp_path, 'wt') as fout:
        while True:
            line1 = fin.readline()
            if line1 == '':
                break
            line2 = fin.readline()
            line3 = fin.readline()
            line4 = fin.readline()
            if any([line == '' for line in [line2, line3, line4]]):
                raise Exception('Malformed FASTQ (record truncated): {}'.format(path_fastq))
            spot_index += 1
            if spot_index < start:
                continue
            if spot_index > end:
                break
            fout.write(line1)
            fout.write(line2)
            fout.write(line3)
            fout.write(line4)
            kept += 1
    os.replace(tmp_path, path_fastq)
    return kept

def trim_fasterq_output_files(sra_stat, start, end):
    work_dir = sra_stat['getfastq_sra_dir']
    sra_id = sra_stat['sra_id']
    for suffix in ['', '_1', '_2']:
        path_fastq = os.path.join(work_dir, sra_id + suffix + '.fastq')
        trim_fastq_by_spot_range(path_fastq=path_fastq, start=start, end=end)

def compress_fasterq_output_files(sra_stat, args):
    work_dir = sra_stat['getfastq_sra_dir']
    sra_id = sra_stat['sra_id']
    files = list()
    for suffix in ['', '_1', '_2']:
        path_fastq = os.path.join(work_dir, sra_id + suffix + '.fastq')
        if os.path.exists(path_fastq):
            files.append(path_fastq)
    if len(files) == 0:
        return
    pigz_exe = shutil.which('pigz')
    def compress_fastq_with_python(path_fastq):
        path_fastq_gz = path_fastq + '.gz'
        with open(path_fastq, 'rb') as f_in, gzip.open(path_fastq_gz, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(path_fastq)
    if (pigz_exe is not None) and (args.threads > 1):
        command = [pigz_exe, '-p', str(args.threads)] + files
        print('Command:', ' '.join(command))
        out = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if should_print_getfastq_command_output(args):
            print('Compression stdout:')
            print(out.stdout.decode('utf8'))
            print('Compression stderr:')
            print(out.stderr.decode('utf8'))
        if out.returncode != 0:
            raise Exception('Compression failed.')
    else:
        print('Using Python gzip fallback for compression.')
        for path_fastq in files:
            print('Compressing: {} -> {}'.format(path_fastq, path_fastq + '.gz'))
            compress_fastq_with_python(path_fastq)

def run_fasterq_dump(sra_stat, args, metadata, start, end):
    path_downloaded_sra = os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + '.sra')
    fasterq_dump_exe = getattr(args, 'fasterq_dump_exe', 'fasterq-dump')
    fasterq_dump_command = [
        fasterq_dump_exe,
        '--split-3',
        '--skip-technical',
        '--min-read-len', str(args.min_read_length),
        '-e', str(max(1, args.threads)),
        '-O', sra_stat['getfastq_sra_dir'],
        '-t', sra_stat['getfastq_sra_dir'],
        path_downloaded_sra,
    ]
    print('Total sampled bases:', "{:,}".format(sra_stat['spot_length'] * (end - start + 1)), 'bp')

    def run_fasterq_dump_command(prefix='Command'):
        print('{}: {}'.format(prefix, ' '.join(fasterq_dump_command)))
        fqd_out = subprocess.run(fasterq_dump_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if should_print_getfastq_command_output(args):
            print('fasterq-dump stdout:')
            print(fqd_out.stdout.decode('utf8'))
            print('fasterq-dump stderr:')
            print(fqd_out.stderr.decode('utf8'))
        return fqd_out

    fqd_out = run_fasterq_dump_command(prefix='Command')
    if (fqd_out.returncode != 0):
        sys.stderr.write("fasterq-dump did not finish safely. Removing the cached SRA file and retrying once.\n")
        remove_sra_path(path_downloaded_sra)
        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=sra_stat['getfastq_sra_dir'], overwrite=True)
        fqd_out = run_fasterq_dump_command(prefix='Retry command')
        if (fqd_out.returncode != 0):
            sys.stderr.write("fasterq-dump did not finish safely after re-download.\n")
            sys.exit(1)
    trim_fasterq_output_files(sra_stat=sra_stat, start=start, end=end)
    compress_fasterq_output_files(sra_stat=sra_stat, args=args)
    nw = [estimate_num_written_spots_from_fastq(sra_stat), ]
    nd = [sum(nw), ]
    nr = [0, ]
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_stat['sra_id']].values[0]
    metadata.df.at[ind_sra,'num_dumped'] += sum(nd)
    metadata.df.at[ind_sra,'num_rejected'] += sum(nr)
    metadata.df.at[ind_sra,'num_written'] += sum(nw)
    metadata.df.at[ind_sra,'bp_dumped'] += sum(nd) * sra_stat['spot_length']
    metadata.df.at[ind_sra,'bp_rejected'] += sum(nr) * sra_stat['spot_length']
    metadata.df.at[ind_sra,'bp_written'] += sum(nw) * sra_stat['spot_length']
    sra_stat = detect_layout_from_file(sra_stat)
    remove_unpaired_files(sra_stat)
    metadata.df.at[ind_sra,'layout_amalgkit'] = sra_stat['layout']
    return metadata,sra_stat

def remove_unpaired_files(sra_stat):
    if (sra_stat['layout']=='paired'):
        # Order is important in this list. More downstream should come first.
        extensions = ['.amalgkit.fastq.gz', '.rename.fastq.gz', '.fastp.fastq.gz', '.fastq.gz']
        for ext in extensions:
            single_fastq_file = os.path.join(sra_stat['getfastq_sra_dir'], sra_stat['sra_id'] + ext)
            if os.path.exists(single_fastq_file):
                print('Removing 3rd fastq file: {}'.format(single_fastq_file), flush=True)
                os.remove(single_fastq_file)

def get_identical_paired_ratio(read1_path, read2_path, num_checked_reads=IDENTICAL_PAIRED_CHECKED_READS):
    num_identical = 0
    num_checked = 0
    read_length = 0
    with gzip.open(read1_path, 'rt') as read1, gzip.open(read2_path, 'rt') as read2:
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

def maybe_treat_paired_as_single(sra_stat, metadata, work_dir,
                                 threshold=IDENTICAL_PAIRED_RATIO_THRESHOLD,
                                 num_checked_reads=IDENTICAL_PAIRED_CHECKED_READS):
    if sra_stat['layout'] != 'paired':
        return metadata, sra_stat
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=work_dir)
    if inext == 'no_extension_found':
        return metadata, sra_stat
    read1_path = os.path.join(work_dir, sra_stat['sra_id'] + '_1' + inext)
    read2_path = os.path.join(work_dir, sra_stat['sra_id'] + '_2' + inext)
    if not (os.path.exists(read1_path) and os.path.exists(read2_path)):
        return metadata, sra_stat
    identical_ratio, num_checked, read_length = get_identical_paired_ratio(
        read1_path=read1_path,
        read2_path=read2_path,
        num_checked_reads=num_checked_reads,
    )
    if (num_checked > 0) and (identical_ratio >= threshold):
        txt = 'Read1 and Read2 are nearly identical ({:.2%} of {:,} pairs): {}. '
        txt += 'Treating as single-end reads and removing redundant read2 file.\n'
        sys.stderr.write(txt.format(identical_ratio, num_checked, sra_stat['sra_id']))
        single_path = os.path.join(work_dir, sra_stat['sra_id'] + inext)
        if os.path.exists(single_path):
            os.remove(single_path)
        os.rename(read1_path, single_path)
        os.remove(read2_path)
        sra_stat['layout'] = 'single'
        ind_sra = metadata.df.index[metadata.df.loc[:, 'run'] == sra_stat['sra_id']].values[0]
        if 'layout_amalgkit' in metadata.df.columns:
            metadata.df.at[ind_sra, 'layout_amalgkit'] = 'single'
        if (read_length > 0) and ('spot_length' in metadata.df.columns):
            sra_stat['spot_length'] = read_length
            metadata.df.at[ind_sra, 'spot_length'] = read_length
            if 'spot_length_amalgkit' in metadata.df.columns:
                metadata.df.at[ind_sra, 'spot_length_amalgkit'] = read_length
    return metadata, sra_stat

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
    ind_sra = metadata.df.index[metadata.df.loc[:, 'run'] == sra_stat['sra_id']].values[0]
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

def run_fastp(sra_stat, args, output_dir, metadata):
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
    outext = '.fastp.fastq.gz'
    if args.threads > 16:
        print('Too many threads for fastp (--threads {}). Only 16 threads will be used.'.format(args.threads))
        fastp_thread = 16
    else:
        fastp_thread = args.threads
    fastp_exe = getattr(args, 'fastp_exe', 'fastp')
    if args.fastp_option:
        try:
            fastp_option_args = shlex.split(args.fastp_option)
        except ValueError as e:
            raise ValueError('Invalid --fastp_option string: {}'.format(args.fastp_option)) from e
    else:
        fastp_option_args = []
    fp_command = [fastp_exe, '--thread', str(fastp_thread), '--length_required', str(args.min_read_length)] + fastp_option_args
    if sra_stat['layout'] == 'single':
        infile = os.path.join(output_dir, sra_stat['sra_id'])
        fp_command = fp_command + ['--in1', infile + inext, '--out1', infile + outext]
    elif sra_stat['layout'] == 'paired':
        infile1 = os.path.join(output_dir, sra_stat['sra_id'] + '_1')
        infile2 = os.path.join(output_dir, sra_stat['sra_id'] + '_2')
        fp_command = fp_command + ['--in1', infile1 + inext, '--out1', infile1 + outext, '--in2', infile2 + inext,
                                   '--out2', infile2 + outext]
    fp_command = [fc for fc in fp_command if fc != '']
    print('Command:', ' '.join(fp_command))
    fp_out = subprocess.run(fp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if args.fastp_print:
        print('fastp stdout:')
        print(fp_out.stdout.decode('utf8'))
        print('fastp stderr:')
        print(fp_out.stderr.decode('utf8'))
    if fp_out.returncode != 0:
        raise RuntimeError(
            'fastp did not finish safely (exit code {}). Command: {}\n{}'.format(
                fp_out.returncode,
                ' '.join(fp_command),
                fp_out.stderr.decode('utf8', errors='replace'),
            )
        )
    if args.remove_tmp:
        remove_intermediate_files(sra_stat, ext=inext, work_dir=output_dir)
    fp_stderr = fp_out.stderr.decode('utf8')
    num_in, num_out, bp_in, bp_out = parse_fastp_summary_counts(fp_stderr)
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_stat['sra_id']].values[0]
    current_num_in = sum(num_in)
    duplication_rate, insert_size_peak = parse_fastp_metrics(fp_stderr)
    update_fastp_metrics(metadata, ind_sra, current_num_in, duplication_rate, insert_size_peak)
    metadata.df.at[ind_sra,'num_fastp_in'] += sum(num_in)
    metadata.df.at[ind_sra,'num_fastp_out'] += sum(num_out)
    metadata.df.at[ind_sra,'bp_fastp_in'] += sum(bp_in)
    metadata.df.at[ind_sra,'bp_fastp_out'] += sum(bp_out)
    write_fastp_stats(sra_stat, metadata, output_dir)
    return metadata

def rename_reads(sra_stat, args, output_dir):
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
    outext = '.rename.fastq.gz'
    def open_fastq(path, mode):
        if path.endswith('.gz'):
            return gzip.open(path, mode)
        return open(path, mode)

    def rewrite_headers(infile, outfile, suffix):
        with open_fastq(infile, 'rt') as fin, open_fastq(outfile, 'wt') as fout:
            while True:
                line1 = fin.readline()
                if line1 == '':
                    break
                line2 = fin.readline()
                line3 = fin.readline()
                line4 = fin.readline()
                if any([line == '' for line in [line2, line3, line4]]):
                    raise Exception('Malformed FASTQ (record truncated): {}'.format(infile))
                header = line1.rstrip('\n').split()[0]
                fout.write(header + suffix + '\n')
                fout.write(line2)
                fout.write(line3)
                fout.write(line4)

    if sra_stat['layout'] == 'single':
        inbase = os.path.join(output_dir, sra_stat['sra_id'])
        if os.path.exists(inbase + inext):
            infile = inbase + inext
            outfile = inbase + outext
            print('Rewriting read headers for Trinity format: {} -> {}'.format(infile, outfile))
            rewrite_headers(infile=infile, outfile=outfile, suffix='/1')
    elif sra_stat['layout'] == 'paired':
        inbase1 = os.path.join(output_dir, sra_stat['sra_id'] + '_1')
        inbase2 = os.path.join(output_dir, sra_stat['sra_id'] + '_2')
        if os.path.exists(inbase1 + inext):
            infile1 = inbase1 + inext
            infile2 = inbase2 + inext
            outfile1 = inbase1 + outext
            outfile2 = inbase2 + outext
            print('Rewriting read headers for Trinity format: {} -> {}'.format(infile1, outfile1))
            rewrite_headers(infile=infile1, outfile=outfile1, suffix='/1')
            print('Rewriting read headers for Trinity format: {} -> {}'.format(infile2, outfile2))
            rewrite_headers(infile=infile2, outfile=outfile2, suffix='/2')
    if args.remove_tmp:
        remove_intermediate_files(sra_stat, ext=inext, work_dir=output_dir)

def rename_fastq(sra_stat, output_dir, inext, outext):
    if sra_stat['layout'] == 'single':
        inbase = os.path.join(output_dir, sra_stat['sra_id'])
        os.rename(inbase + inext, inbase + outext)
    elif sra_stat['layout'] == 'paired':
        inbase1 = os.path.join(output_dir, sra_stat['sra_id'] + '_1')
        inbase2 = os.path.join(output_dir, sra_stat['sra_id'] + '_2')
        os.rename(inbase1 + inext, inbase1 + outext)
        os.rename(inbase2 + inext, inbase2 + outext)

def calc_2nd_ranges(metadata):
    sra_target_bp = metadata.df.loc[:,'bp_until_target_size']
    rate_obtained = metadata.df.loc[:,'rate_obtained']
    spot_lengths = metadata.df.loc[:,'spot_length_amalgkit']
    total_spots = metadata.df.loc[:,'total_spots']
    sra_target_reads = numpy.zeros_like(sra_target_bp)
    for i, ind in enumerate(metadata.df.index):
        if numpy.isnan(rate_obtained.loc[ind]):
            sra_target_reads[i] = (sra_target_bp.loc[ind]/spot_lengths.loc[ind]).astype(int)+1 # If no read was extracted in 1st.
        else:
            sra_target_reads[i] = ((sra_target_bp.loc[ind]/spot_lengths.loc[ind])/rate_obtained.loc[ind]).astype(int)+1
    start_2nds = metadata.df.loc[:,'spot_end_1st'] + 1
    end_2nds = start_2nds + sra_target_reads
    pooled_missing_bp = metadata.df.loc[:,'bp_until_target_size'].sum()
    for dummy in range(1000):
        current_total_bp = 0
        for ind in end_2nds.index:
            pooled_missing_read = (pooled_missing_bp / spot_lengths.loc[ind]).astype(int)
            if ((end_2nds.loc[ind] + pooled_missing_read) < total_spots.loc[ind]):
                pooled_missing_bp = 0
                end_2nds.loc[ind] = end_2nds.loc[ind] + pooled_missing_bp
            elif (end_2nds.loc[ind] + pooled_missing_read > total_spots.loc[ind]):
                pooled_missing_bp = (end_2nds.loc[ind] + pooled_missing_read - total_spots.loc[ind]) * \
                                    spot_lengths.loc[ind]
                end_2nds.loc[ind] = total_spots.loc[ind]
            current_total_bp += end_2nds.loc[ind] * spot_lengths.loc[ind]
        all_equal_total_spots = all([e2 == ts for e2, ts in zip(end_2nds, total_spots)])
        is_enough_read = (current_total_bp >= metadata.df.loc[:,'bp_until_target_size'].sum())
        if all_equal_total_spots:
            print('Reached total spots in all SRAs.', flush=True)
            break
        if is_enough_read:
            print('Enough read numbers were assigned for the 2nd round sequence extraction.', flush=True)
            break
    metadata.df.loc[:,'spot_start_2nd'] = start_2nds
    metadata.df.loc[:,'spot_end_2nd'] = end_2nds
    return metadata

def is_2nd_round_needed(rate_obtained_1st, tol):
    # tol is acceptable percentage loss relative to --max_bp.
    required_rate = 1 - (tol * 0.01)
    return rate_obtained_1st < required_rate

def print_read_stats(args, metadata, g, sra_stat=None, individual=False):
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
        if len(read_types) > 0:
            for rt, key in zip(read_types, keys):
                values = ['{:,}'.format(s) for s in df[key].values]
                txt = ' '.join(values)
                print('Individual {} (bp): {}'.format(rt, txt))
    print('')

def getfastq_metadata(args):
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
        metadata.df = metadata.df.loc[(metadata.df['lib_layout'] == layout), :]
        if args.sci_name is not None:
            print('Filtering SRA entry with --sci_name:', args.sci_name)
            metadata.df = metadata.df.loc[(metadata.df['scientific_name'] == args.sci_name), :]
    if args.id_list is not None:
        print('--id_list is specified. Downloading SRA metadata from Entrez.')
        Entrez.email = args.entrez_email
        with open(args.id_list) as fin:
            sra_id_list = [
                line.strip()
                for line in fin
                if (line.strip() != '') and (not line.lstrip().startswith('#'))
            ]
        metadata_dict = dict()
        for sra_id in sra_id_list:
            search_term = getfastq_search_term(sra_id, args.entrez_additional_search_term)
            print('Entrez search term:', search_term)
            xml_root = getfastq_getxml(search_term)
            metadata_dict_tmp = Metadata.from_xml(xml_root)
            if metadata_dict_tmp.df.shape[0]==0:
                print('No associated SRA. Skipping {}'.format(sra_id))
                continue
            metadata_dict[sra_id] = metadata_dict_tmp
            print('Filtering SRA entry with --layout:', args.layout)
            layout = get_layout(args, metadata_dict[sra_id])
            metadata_dict[sra_id].df = metadata_dict[sra_id].df.loc[(metadata_dict[sra_id].df['lib_layout'] == layout), :]
            if args.sci_name is not None:
                print('Filtering SRA entry with --sci_name:', args.sci_name)
                metadata_dict[sra_id].df = metadata_dict[sra_id].df.loc[(metadata_dict[sra_id].df['scientific_name'] == args.sci_name), :]
        if len(metadata_dict)==0:
            print('No associated SRA is found with --id_list. Exiting.')
            sys.exit(1)
        metadata = list(metadata_dict.values())[0]
        metadata.df = pandas.concat([ v.df for v in metadata_dict.values() ], ignore_index=True)
    if (args.id is None)&(args.id_list is None):
        metadata = load_metadata(args)
    metadata.df['total_bases'] = metadata.df.loc[:,'total_bases'].replace('', numpy.nan).astype(float)
    metadata.df['spot_length'] = metadata.df.loc[:, 'spot_length'].replace('', numpy.nan).astype(float)
    return metadata

def is_getfastq_output_present(sra_stat):
    sra_stat = detect_layout_from_file(sra_stat)
    prefixes = [sra_stat['sra_id'], ]
    if sra_stat['layout'] == 'single':
        sub_exts = ['', ]
    elif sra_stat['layout'] == 'paired':
        sub_exts = ['_1', '_2']
    exts = ['.amalgkit.fastq.gz', ]
    is_output_present = True
    for prefix, sub_ext, ext in itertools.product(prefixes, sub_exts, exts):
        out_path1 = os.path.join(sra_stat['getfastq_sra_dir'], prefix + sub_ext + ext)
        out_path2 = os.path.join(sra_stat['getfastq_sra_dir'], prefix + sub_ext + ext + '.safely_removed')
        is_out1 = os.path.exists(out_path1)
        is_out2 = os.path.exists(out_path2)
        if is_out1:
            print('getfastq output detected: {}'.format(out_path1))
        if is_out2:
            print('getfastq output detected: {}'.format(out_path2))
        is_output_present *= (is_out1 | is_out2)
    return is_output_present

def remove_experiment_without_run(metadata):
    num_all_run = metadata.df.shape[0]
    is_missing_run = (metadata.df.loc[:, 'run'] == '')
    num_missing_run = is_missing_run.sum()
    if (num_missing_run > 0):
        print('There are {} out of {} Experiments without Run ID. Removing.'.format(num_missing_run, num_all_run))
        metadata.df = metadata.df.loc[~is_missing_run, :]
    return metadata

def initialize_columns(metadata, g):
    time_keys = ['time_start_1st', 'time_end_1st', 'time_start_2nd', 'time_end_2nd',]
    spot_keys = ['spot_start_1st', 'spot_end_1st', 'spot_start_2nd', 'spot_end_2nd',]
    keys = (['num_dumped', 'num_rejected', 'num_written', 'num_fastp_in', 'num_fastp_out','bp_amalgkit',
            'bp_dumped', 'bp_rejected', 'bp_written', 'bp_fastp_in', 'bp_fastp_out', 'bp_discarded',
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
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_id].values[0]
    metadata,sra_stat = run_fasterq_dump(sra_stat, args, metadata, start, end)
    metadata, sra_stat = maybe_treat_paired_as_single(
        sra_stat=sra_stat,
        metadata=metadata,
        work_dir=sra_stat['getfastq_sra_dir'],
    )
    bp_discarded = metadata.df.at[ind_sra,'bp_dumped'] - metadata.df.at[ind_sra,'bp_written']
    metadata.df.at[ind_sra,'bp_discarded'] += bp_discarded
    metadata.df.at[ind_sra,'layout_amalgkit'] = sra_stat['layout']
    no_read_written = (metadata.df.loc[(metadata.df.loc[:,'run']==sra_id),'num_written'].values[0]==0)
    if no_read_written:
        return metadata
    if args.fastp:
        metadata = run_fastp(sra_stat, args, sra_stat['getfastq_sra_dir'], metadata)
        bp_discarded = metadata.df.at[ind_sra,'bp_dumped'] - metadata.df.at[ind_sra,'bp_fastp_out']
        metadata.df.at[ind_sra,'bp_discarded'] += bp_discarded
    if args.read_name == 'trinity':
        rename_reads(sra_stat, args, sra_stat['getfastq_sra_dir'])
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=sra_stat['getfastq_sra_dir'])
    outext = '.amalgkit.fastq.gz'
    rename_fastq(sra_stat, sra_stat['getfastq_sra_dir'], inext, outext)
    metadata.df.at[ind_sra,'bp_still_available'] = sra_stat['spot_length'] * (sra_stat['total_spot'] - end)
    bp_specified_for_extraction = sra_stat['spot_length'] * (end - start)
    metadata.df.at[ind_sra, 'bp_specified_for_extraction'] += bp_specified_for_extraction
    if args.fastp:
        metadata.df.at[ind_sra, 'bp_amalgkit'] = metadata.df.at[ind_sra,'bp_fastp_out']
    else:
        metadata.df.at[ind_sra, 'bp_amalgkit'] = metadata.df.at[ind_sra,'bp_written']
    metadata.df.at[ind_sra, 'rate_obtained'] = metadata.df.at[ind_sra, 'bp_amalgkit'] / g['num_bp_per_sra']
    metadata.df.at[ind_sra, 'bp_until_target_size'] -= metadata.df.at[ind_sra, 'bp_amalgkit']
    return metadata

def sequence_extraction_1st_round(args, sra_stat, metadata, g):
    offset = 10000
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_stat['sra_id']].values[0]
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
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_stat['sra_id']].values[0]
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
    no_read_in_1st = (metadata.df.loc[(metadata.df.loc[:,'run']==sra_id),'bp_written'].values[0]==0)
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
    for subext in subexts:
        added_path = os.path.join(sra_stat['getfastq_sra_dir'], sra_id + subext + ext_1st_tmp)
        adding_path = os.path.join(sra_stat['getfastq_sra_dir'], sra_id + subext + ext_main)
        assert os.path.exists(added_path), 'Dumped fastq not found: ' + added_path
        assert os.path.exists(adding_path), 'Dumped fastq not found: ' + adding_path
        append_file_binary(adding_path, added_path)
        os.remove(adding_path)
        os.rename(added_path, adding_path)
    metadata.df.at[ind_sra, 'time_end_2nd'] = time.time()
    elapsed_time = metadata.df.at[ind_sra, 'time_end_2nd'] - metadata.df.at[ind_sra, 'time_start_2nd']
    txt = 'Time elapsed for 2nd-round sequence extraction: {}, {:,} sec'
    print(txt.format(sra_stat['sra_id'], elapsed_time))
    print('')
    return metadata

def sequence_extraction_private(i, metadata, sra_stat, args):
    for col in ['read1_path','read2_path']:
        path_from = metadata.df.at[i,col]
        path_to = os.path.join(sra_stat['getfastq_sra_dir'], os.path.basename(path_from))
        path_to = path_to.replace('.fq', '.fastq')
        if not path_to.endswith('.gz'):
            path_to = path_to+'.gz' # .gz is necessary even if the original file is not compressed.
        if os.path.exists(path_from):
            if os.path.lexists(path_to):
                os.remove(path_to)
            os.symlink(src=path_from, dst=path_to)
        else:
            sys.stderr.write('Private fastq file not found: {}\n'.format(path_from))
    if args.fastp:
        metadata = run_fastp(sra_stat, args, sra_stat['getfastq_sra_dir'], metadata)
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=sra_stat['getfastq_sra_dir'])
    if (inext=='no_extension_found')&(sra_stat['layout']=='paired'):
        raise Exception('Paired-end file names may be invalid. They should contain _1 and _2 to indicate a pair: {}'.format(sra_stat['sra_id']))
    outext = '.amalgkit.fastq.gz'
    rename_fastq(sra_stat, sra_stat['getfastq_sra_dir'], inext, outext)
    return metadata

def check_metadata_validity(metadata):
    assert metadata.df.shape[0] > 0, 'No SRA entry found. Make sure whether --id or --id_list is compatible with --sci_name and --layout.'
    is_total_bases_na = metadata.df.loc[:,'total_bases'].isnull()
    is_total_bases_na |= (metadata.df.loc[:, 'total_bases']==0)
    is_total_bases_na |= (metadata.df.loc[:, 'total_bases']=='')
    if is_total_bases_na.any():
        txt = 'Empty value(s) of total_bases were detected in {}. Filling a placeholder value 999,999,999,999\n'
        sys.stderr.write(txt.format(', '.join(metadata.df.loc[is_total_bases_na, 'run'])))
        metadata.df.loc[is_total_bases_na,'total_bases'] = 999999999999
        metadata.df['total_bases'] = metadata.df.loc[:, 'total_bases'].astype(int)
    is_total_spots_na = metadata.df.loc[:, 'total_spots'].isnull()
    is_total_spots_na |=  (metadata.df.loc[:, 'total_spots']==0)
    is_total_spots_na |=  (metadata.df.loc[:, 'total_spots']=='')
    if is_total_spots_na.any():
        new_values = metadata.df.loc[is_total_spots_na,'total_bases'] / metadata.df.loc[is_total_spots_na,'spot_length']
        if is_total_spots_na.any():
            txt = 'Empty value(s) of total_spots were detected in {}. Filling a placeholder value 999,999,999,999\n'
            sys.stderr.write(txt.format(', '.join(metadata.df.loc[is_total_spots_na,'run'])))
            new_values.loc[new_values.isnull()] = 999999999999 # https://github.com/kfuku52/amalgkit/issues/110
        new_values = new_values.astype(int)
        metadata.df.loc[is_total_spots_na, 'total_spots'] = new_values
    for i in metadata.df.index:
        txt = 'Individual SRA size of {}: {:,} bp'
        print(txt.format(metadata.df.at[i, 'run'], metadata.df.at[i, 'total_bases']))
    return metadata

def initialize_global_params(args, metadata):
    g = dict()
    g['start_time'] = time.time()
    g['max_bp'] = int(args.max_bp.replace(',', ''))
    g['num_sra'] = metadata.df.shape[0]
    g['num_bp_per_sra'] = int(g['max_bp'] / g['num_sra'])
    g['total_sra_bp'] = metadata.df.loc[:,'total_bases'].sum()
    print('Number of SRAs to be processed: {:,}'.format(g['num_sra']))
    print('Total target size (--max_bp): {:,} bp'.format(g['max_bp']))
    print('The sum of SRA sizes: {:,} bp'.format(g['total_sra_bp']))
    print('Target size per SRA: {:,} bp'.format(g['num_bp_per_sra']))
    return g

def getfastq_main(args):
    check_getfastq_dependency(args)
    metadata = getfastq_metadata(args)
    metadata = remove_experiment_without_run(metadata)
    metadata = check_metadata_validity(metadata)
    g = initialize_global_params(args, metadata)
    metadata = initialize_columns(metadata, g)
    flag_private_file = False
    flag_any_output_file_present = False
    # 1st round sequence extraction
    for i in metadata.df.index:
        print('')
        sra_id = metadata.df.at[i, 'run']
        print('Processing SRA ID: {}'.format(sra_id))
        sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
        sra_stat['getfastq_sra_dir'] = get_getfastq_run_dir(args, sra_id)
        if (is_getfastq_output_present(sra_stat)):
            flag_any_output_file_present =True
            if not args.redo:
                txt = 'Output file(s) detected. Skipping {}. Set "--redo yes" for reanalysis.'
                print(txt.format(sra_id), flush=True)
                continue
        remove_old_intermediate_files(sra_id=sra_id, work_dir=sra_stat['getfastq_sra_dir'])
        print('Library layout:', sra_stat['layout'])
        print('Number of reads:', "{:,}".format(sra_stat['total_spot']))
        print('Single/Paired read length:', sra_stat['spot_length'], 'bp')
        print('Total bases:', "{:,}".format(int(metadata.df.loc[i, 'total_bases'])), 'bp')
        flag_private_file = False
        if 'private_file' in metadata.df.columns:
            if metadata.df.at[i,'private_file']=='yes':
                print('Processing {} as private data. --max_bp is disabled.'.format(sra_id), flush=True)
                flag_private_file = True
                sequence_extraction_private(i, metadata, sra_stat, args)
        if not flag_private_file:
            print('Processing {} as publicly available data from SRA.'.format(sra_id), flush=True)
            download_sra(metadata, sra_stat, args, sra_stat['getfastq_sra_dir'], overwrite=False)
            metadata = sequence_extraction_1st_round(args, sra_stat, metadata, g)
    # 2nd round sequence extraction
    if (not flag_private_file) & (not flag_any_output_file_present):
        g['rate_obtained_1st'] = metadata.df.loc[:,'bp_amalgkit'].sum() / g['max_bp']
        if is_2nd_round_needed(g['rate_obtained_1st'], args.tol):
            txt = 'Only {:,.2f}% ({:,}/{:,}) of the target size (--max_bp) was obtained in the 1st round. Proceeding to the 2nd round read extraction.'
            print(txt.format(g['rate_obtained_1st']*100, metadata.df.loc[:,'bp_amalgkit'].sum(), g['max_bp']), flush=True)
            metadata = calc_2nd_ranges(metadata)
            for i in metadata.df.index:
                sra_id = metadata.df.at[i, 'run']
                sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
                sra_stat['getfastq_sra_dir'] = get_getfastq_run_dir(args, sra_id)
                metadata = sequence_extraction_2nd_round(args, sra_stat, metadata, g)
        else:
            print('Sufficient data were obtained in the 1st-round sequence extraction. Proceeding without the 2nd round.')
        g['rate_obtained_2nd'] = metadata.df.loc[:, 'bp_amalgkit'].sum() / g['max_bp']
        txt = '2nd round read extraction improved % bp from {:,.2f}% to {:,.2f}%'
        print(txt.format(g['rate_obtained_1st']*100, g['rate_obtained_2nd']*100), flush=True)
    # Postprocessing
    if args.remove_sra:
        remove_sra_files(metadata=metadata, amalgkit_out_dir=args.out_dir)
    else:
        print('SRA files not removed: {}'.format(sra_stat['getfastq_sra_dir']))
    if (not flag_any_output_file_present):
        print('')
        print('\n--- getfastq final report ---')
        print_read_stats(args, metadata, g, sra_stat=None, individual=True)
