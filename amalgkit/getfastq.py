from Bio import Entrez
import itertools
import numpy

from amalgkit.util import *

import glob
import lxml
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from urllib.error import HTTPError

def getfastq_search_term(ncbi_id, additional_search_term=None):
    # https://www.ncbi.nlm.nih.gov/books/NBK49540/
    if additional_search_term is None:
        search_term = ncbi_id
    else:
        search_term = ncbi_id + ' AND ' + additional_search_term
    #   search_term = '"'+ncbi_id+'"'+'[Accession]'
    return search_term


def getfastq_getxml(search_term, save_xml=True, retmax=1000):
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
    root = None
    for i in numpy.arange(numpy.ceil(num_record // retmax) + 1):
        start = int(i * retmax)
        end = int(((i + 1) * retmax) - 1) if num_record >= int(((i + 1) * retmax) - 1) else num_record
        print('processing SRA records:', start, '-', end, flush=True)
        try:
            handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
        except HTTPError as e:
            print(e, '- Trying Entrez.efetch() again...')
            handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
        chunk = lxml.etree.parse(handle).getroot()
        if root is None:
            root = chunk
        else:
            root.append(chunk)
    xml_string = lxml.etree.tostring(root, pretty_print=True)
    for line in str(xml_string).split('\n'):
        if '<Error>' in line:
            print(line)
            raise Exception(species_name, ': <Error> found in the xml.')
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


def concat_fastq(args, metadata, output_dir, num_bp_per_sra):
    layout = get_layout(args, metadata)
    inext = '.amalgkit.fastq.gz'
    infiles = list()
    for sra_id in metadata.df.loc[:, 'run']:
        infiles.append([f for f in os.listdir(output_dir) if (f.endswith(inext)) & (f.startswith(sra_id))])
    infiles = [item for sublist in infiles for item in sublist]
    num_inext_files = len(infiles)
    if (layout == 'single') & (num_inext_files == 1):
        print('Only 1', inext, 'file was detected. No concatenation will happen.', flush=True)
        outfile = args.id + inext
        if infiles[0] != outfile:
            print('Replacing ID in the output file name:', infiles[0], outfile)
            infile_path = os.path.join(output_dir, infiles[0])
            outfile_path = os.path.join(output_dir, outfile)
            os.rename(infile_path, outfile_path)
        return None
    elif (layout == 'paired') & (num_inext_files == 2):
        print('Only 1 pair of', inext, 'files were detected. No concatenation will happen.', flush=True)
        for infile in infiles:
            outfile = args.id + re.sub('.*(_[1-2])', '\g<1>', infile)
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
            outfile_path = os.path.join(output_dir, args.id + subext + outext)
            if os.path.exists(outfile_path):
                os.remove(outfile_path)
            # with gzip.open(outfile_path, 'wb') as outfile:
            #    for each_infile in infiles:
            #        with gzip.open(os.path.join(args.out_dir, each_infile), 'rb') as infile:
            #            shutil.copyfileobj(infile, outfile) # unacceptably slow
            if os.path.exists(outfile_path):
                os.remove(outfile_path)
            for infile in infiles:
                infile_path = os.path.join(output_dir, infile)
                assert os.path.exists(infile_path), 'Dumped fastq not found: ' + infile_path
                print('Concatenated file:', infile_path, flush=True)
                os.system('cat "' + infile_path + '" >> "' + outfile_path + '"')
            print('')
        if args.remove_tmp == 'yes':
            for i in metadata.df.index:
                sra_id = metadata.df.loc[i, 'run']
                sra_stat = get_sra_stat(sra_id, metadata, num_bp_per_sra)
                ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
                remove_intermediate_files(sra_stat, ext=ext, work_dir=output_dir)
        return None


def remove_sra_files(metadata, sra_dir):
    for sra_id in metadata.df['run']:
        sra_pattern = os.path.join(sra_dir, sra_id + '.sra*')
        sra_paths = glob.glob(sra_pattern)
        if len(sra_paths) > 0:
            for sra_path in sra_paths:
                print('Deleting:', sra_path)
                os.remove(sra_path)
        else:
            print('SRA file not found. Pattern searched:', sra_pattern)
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


def download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
    sra_path = os.path.join(work_dir, sra_stat['sra_id'] + '.sra')
    individual_sra_tmp_dir = os.path.join(work_dir, sra_stat['sra_id'] + '/')

    if os.path.exists(sra_path):
        print('Previously-downloaded sra file was detected.')
        if (overwrite):
            print('Removing', sra_path)
            print('New sra file will be downloaded.')
            os.remove(sra_path)
        else:
            return None
    else:
        print('Previously-downloaded sra file was not detected. New sra file will be downloaded.')

    if (args.aws) or (args.ncbi) or (args.gcp):
        source = []
        sra_source_list = []
        sra_id = sra_stat['sra_id']

        if args.gcp:
            source.append('GCP')
            sra_source_list.append(metadata.df.loc[metadata.df['run'] == sra_id]['GCP_Link'].values[0])

        if args.aws:
            source.append('AWS')
            sra_source_list.append(metadata.df.loc[metadata.df['run'] == sra_id]['AWS_Link'].values[0])

        if args.ncbi:
            source.append('NCBI')
            sra_source_list.append(metadata.df.loc[metadata.df['run'] == sra_id]['NCBI_Link'].values[0])

        if len(sra_source_list) > 1:
            print("Multiple sources set. Trying one by one.")

        dl_status = 'failed'

        for sra_source in sra_source_list:
            print("trying to fetch {} from {}".format(sra_id, sra_source))
            try:
                urllib.request.urlretrieve(str(sra_source), os.path.join(work_dir, (str(sra_id + '.sra'))))
                dl_status = 'success'
                break
            except urllib.error.URLError:
                print("ERROR: urllib.request did not work. Trying wget")
                try:
                    import wget
                    wget.download(str(sra_source), os.path.join(work_dir, (str(sra_id + '.sra'))))
                except ModuleNotFoundError:
                    print("ERROR: Could not find wget")
                except urllib.error.URLError:
                    print(
                        "ERROR: Could not download from " + sra_source + ".")
                    dl_status = 'failed'
                    continue

                print(
                    "ERROR: Could not download from " + sra_source + ".")
                dl_status = 'failed'
                continue

        if dl_status == 'failed':
            print("Exhausted all Sources, trying prefetch.")
        else:
            print("done!")
            assert os.path.exists(sra_path), 'SRA file download failed: ' + sra_stat['sra_id']
            return

    if not os.path.exists(sra_path):
        print('Trying to download the SRA file using prefetch.')
        if os.path.exists(individual_sra_tmp_dir):
            shutil.rmtree(individual_sra_tmp_dir)
        prefetch_command = [args.prefetch_exe, '--force', 'no', '--max-size', '100G',
                            '--output-directory', './', str(sra_stat['sra_id'])]
        print('Command:', ' '.join(prefetch_command))
        prefetch_out = subprocess.run(prefetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('prefetch stdout:')
        print(prefetch_out.stdout.decode('utf8'))
        print('prefetch stderr:')
        print(prefetch_out.stderr.decode('utf8'))
        if (prefetch_out.returncode !=0):
            sys.stderr.write("prefetch did not finish safely.\n")
        if (prefetch_out.returncode):
            sys.stderr.write('Trying prefetch again...\n')
            prefetch_command = [args.prefetch_exe, '--force', 'no', '--max-size', '100G',
                                sra_stat['sra_id']]
            prefetch_out = subprocess.run(prefetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print('prefetch stdout:')
            print(prefetch_out.stdout.decode('utf8'))
            print('prefetch stderr:')
            print(prefetch_out.stderr.decode('utf8'))
            if (prefetch_out.returncode !=0):
                sys.stderr.write("prefetch did not finish safely.\n")
    # Move files downloaded by prefetch. This is necessary because absolute path didn't work for prefetch --output-directory
    if os.path.exists(os.path.join('./', sra_stat['sra_id'] + '/', sra_stat['sra_id'] + '.sra')):
        subprocess.run(['mv', os.path.join('./', sra_stat['sra_id'] + '/', sra_stat['sra_id'] + '.sra'), sra_path],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        shutil.rmtree(os.path.join('./', sra_stat['sra_id'] + '/'))
    elif os.path.exists(os.path.expanduser(os.path.join('~/ncbi/public/sra/', sra_stat['sra_id'] + '.sra'))):
        subprocess.run(
            ['mv', os.path.expanduser(os.path.join('~/ncbi/public/sra/', sra_stat['sra_id'] + '.sra')), sra_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Move files downloaded by ascp
    if os.path.exists(os.path.join(individual_sra_tmp_dir, sra_stat['sra_id'] + '.sra')):
        subprocess.run(['mv', os.path.join(individual_sra_tmp_dir, sra_stat['sra_id'] + '.sra'), sra_path],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        shutil.rmtree(individual_sra_tmp_dir)
    assert os.path.exists(sra_path), 'SRA file download failed: ' + sra_stat['sra_id']


def check_getfastq_dependency(args):
    if args.pfd == 'yes':
        test_pfd = subprocess.run([args.pfd_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_pfd.returncode == 0), "parallel-fastq-dump PATH cannot be found: " + args.pfd_exe
        test_prefetch = subprocess.run([args.prefetch_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_prefetch.returncode == 0), "prefetch (SRA toolkit) PATH cannot be found: " + args.prefetch_exe
    if args.fastp:
        test_fp = subprocess.run([args.fastp_exe, '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_fp.returncode == 0), "fastp PATH cannot be found: " + args.fastp_exe
    test_pigz = subprocess.run(['pigz', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if test_pigz.returncode == 0:
        print('pigz found. It will be used for compression/decompression in read name formatting.')
        gz_exe = 'pigz -p ' + str(args.threads)
        ungz_exe = 'unpigz -p' + str(args.threads)
    else:
        print('pigz not found. gzip/gunzip will be used for compression/decompression in read name formatting.')
        gz_exe = 'gzip'
        ungz_exe = 'gunzip'
    return gz_exe, ungz_exe


def set_getfastq_directories(args, sra_id):
    if args.out_dir.startswith('./'):
        args.out_dir = args.out_dir.replace('.', os.getcwd())
    if args.id is not None:
        output_dir = os.path.join(args.out_dir, 'getfastq', args.id)
    elif args.metadata is not None:
        output_dir = os.path.join(args.out_dir, 'getfastq', sra_id)
    elif args.id_list is not None:
        output_dir = os.path.join(args.out_dir, 'getfastq', sra_id)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir


def run_pfd(sra_stat, args, output_dir, seq_summary, start, end):
    sra_path = os.path.join(output_dir, sra_stat['sra_id'] + '.sra')
    pfd_command = ['parallel-fastq-dump', '-t', str(args.threads), '--minReadLen', str(args.min_read_length),
                   '--qual-filter-1',
                   '--skip-technical', '--split-3', '--clip', '--gzip', '--outdir', output_dir,
                   '--tmpdir', output_dir]
    print('Total sampled bases:', "{:,}".format(sra_stat['spot_length'] * (end - start + 1)), 'bp')
    pfd_command = pfd_command + ['--minSpotId', str(start), '--maxSpotId', str(end)]
    # If sra_stat['sra_id'], not sra_path, is provided, pfd couldn't find pre-downloaded .sra files
    # and start downloading it to $HOME/ncbi/public/sra/
    pfd_command = pfd_command + ['-s', sra_path]
    # pfd_command = pfd_command + ['-s', sra_stat['sra_id']]
    print('Command:', ' '.join(pfd_command))
    pfd_out = subprocess.run(pfd_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if args.pfd_print == 'yes':
        print('parallel-fastq-dump stdout:')
        print(pfd_out.stdout.decode('utf8'))
        print('parallel-fastq-dump stderr:')
        print(pfd_out.stderr.decode('utf8'))
    if (pfd_out.returncode != 0):
        sys.stderr.write("pfd did not finish safely.\n")
        sys.exit(1)
    stdout = pfd_out.stdout.decode('utf8')
    nd = [int(line.replace('Read ', '').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Read')]
    nr = [int(line.replace('Rejected ', '').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Rejected')]
    nw = [int(line.replace('Written ', '').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Written')]
    seq_summary['num_dumped'].loc[sra_stat['sra_id']] += sum(nd)
    seq_summary['num_rejected'].loc[sra_stat['sra_id']] += sum(nr)
    seq_summary['num_written'].loc[sra_stat['sra_id']] += sum(nw)
    seq_summary['bp_dumped'].loc[sra_stat['sra_id']] += sum(nd) * sra_stat['spot_length']
    seq_summary['bp_rejected'].loc[sra_stat['sra_id']] += sum(nr) * sra_stat['spot_length']
    seq_summary['bp_written'].loc[sra_stat['sra_id']] += sum(nw) * sra_stat['spot_length']
    if (sra_stat['layout'] == 'paired'):
        fastq_files = glob.glob(os.path.join(output_dir, sra_stat['sra_id']+'*.fastq*'))
        unpaired_file = os.path.join(output_dir, sra_stat['sra_id'] + '.fastq.gz')
        if (os.path.exists(unpaired_file))&(len(fastq_files)==3):
            print('layout = {}; Deleting unpaired file: {}'.format(sra_stat['layout'], unpaired_file))
            os.remove(unpaired_file)
        elif (os.path.exists(unpaired_file))&(len(fastq_files)==1):
            sys.stderr.write('Single-end fastq was generated even though layout = {}\n'.format(sra_stat['layout']))
            sys.stderr.write('This sample will be treated as single-end sequencing.\n')
            sra_stat['layout'] = 'single'
        else:
            print('Unpaired file not found:', unpaired_file)
    seq_summary['layout'].loc[sra_stat['sra_id']] = sra_stat['layout']
    return seq_summary,sra_stat


def run_fastp(sra_stat, args, output_dir, seq_summary):
    print('Running fastp.')
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
    outext = '.fastp.fastq.gz'
    if args.threads > 16:
        print('Too many threads for fastp (--threads {}). Only 16 threads will be used.'.format(args.threads))
        fastp_thread = 16
    else:
        fastp_thread = args.threads
    fp_command = ['fastp', '--thread', str(fastp_thread), '--length_required',
                  str(args.min_read_length)] + args.fastp_option.split(' ')
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
    try:
        fp_out = subprocess.run(fp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' returned with error (code {}): {}".format(e.cmd, e.returncode, e.output))
    if args.fastp_print == 'yes':
        print('fastp stdout:')
        print(fp_out.stdout.decode('utf8'))
        print('fastp stderr:')
        print(fp_out.stderr.decode('utf8'))
    if args.remove_tmp == 'yes':
        remove_intermediate_files(sra_stat, ext=inext, work_dir=output_dir)
    bps = fp_out.stderr.decode('utf8').split('\n')
    num_in = list()
    num_out = list()
    bp_in = list()
    bp_out = list()
    for i in range(len(bps)):
        if (' before filtering:' in bps[i]):
            num_in.append(int(bps[i + 1].replace('total reads: ', '')))
            bp_in.append(int(bps[i + 2].replace('total bases: ', '')))
        if (' after filtering:' in bps[i]) | (' aftering filtering:' in bps[i]):
            num_out.append(int(bps[i + 1].replace('total reads: ', '')))
            bp_out.append(int(bps[i + 2].replace('total bases: ', '')))
    seq_summary['num_fastp_in'].loc[sra_stat['sra_id']] += sum(num_in)
    seq_summary['num_fastp_out'].loc[sra_stat['sra_id']] += sum(num_out)
    seq_summary['bp_fastp_in'].loc[sra_stat['sra_id']] += sum(bp_in)
    seq_summary['bp_fastp_out'].loc[sra_stat['sra_id']] += sum(bp_out)
    return seq_summary


def rename_reads(sra_stat, args, output_dir, gz_exe, ungz_exe):
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
    outext = '.rename.fastq.gz'
    if sra_stat['layout'] == 'single':
        inbase = os.path.join(output_dir, sra_stat['sra_id'])
        if os.path.exists(inbase + inext):
            infile = inbase + inext
            os.system(
                ungz_exe + ' -c "' + infile + '" | sed -e "s|[[:space:]].*|/1|" | ' + gz_exe + ' -c > "' + inbase + outext + '"')
    elif sra_stat['layout'] == 'paired':
        inbase1 = os.path.join(output_dir, sra_stat['sra_id'] + '_1')
        inbase2 = os.path.join(output_dir, sra_stat['sra_id'] + '_2')
        if os.path.exists(inbase1 + inext):
            infile1 = inbase1 + inext
            infile2 = inbase2 + inext
            os.system(
                ungz_exe + ' -c "' + infile1 + '" | sed -e "s|[[:space:]].*|/1|" | ' + gz_exe + ' -c > "' + inbase1 + outext + '"')
            os.system(
                ungz_exe + ' -c "' + infile2 + '" | sed -e "s|[[:space:]].*|/2|" | ' + gz_exe + ' -c > "' + inbase2 + outext + '"')
    if args.remove_tmp == 'yes':
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


def sequence_extraction(args, sra_stat, start, end, output_dir, seq_summary, gz_exe, ungz_exe, total_sra_bp):
    sra_id = sra_stat['sra_id']
    if args.pfd == 'yes':
        seq_summary,sra_stat = run_pfd(sra_stat, args, output_dir, seq_summary, start, end)
        bp_remaining = seq_summary['bp_dumped'].loc[sra_id] - seq_summary['bp_written'].loc[sra_id]
        seq_summary['bp_remaining'].loc[sra_id] = bp_remaining
    if args.fastp:
        seq_summary = run_fastp(sra_stat, args, output_dir, seq_summary)
        bp_remaining = seq_summary['bp_dumped'].loc[sra_id] - seq_summary['bp_fastp_out'].loc[sra_id]
        seq_summary['bp_remaining'].loc[sra_id] = bp_remaining
    if args.read_name == 'trinity':
        rename_reads(sra_stat, args, output_dir, gz_exe, ungz_exe)
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
    outext = '.amalgkit.fastq.gz'
    rename_fastq(sra_stat, output_dir, inext, outext)
    seq_summary['bp_still_available'].loc[sra_id] = sra_stat['spot_length'] * (sra_stat['total_spot'] - end)
    seq_summary['rate_passed'].loc[sra_id] = (total_sra_bp - seq_summary['bp_remaining'].loc[sra_id]) / total_sra_bp
    seq_summary['spot_length'].loc[sra_id] = sra_stat['spot_length']
    seq_summary['total_spot'].loc[sra_id] = sra_stat['total_spot']
    seq_summary['start_1st'].loc[sra_id] = start
    seq_summary['end_1st'].loc[sra_id] = end
    return seq_summary


def calc_2nd_ranges(seq_summary):
    sra_target_bp = seq_summary['total_bp_remaining'] / len(seq_summary['bp_remaining'])
    sra_target_reads = (sra_target_bp / seq_summary['spot_length']).astype(int)
    start_2nds = seq_summary['end_1st'] + 1
    end_2nds = start_2nds + (sra_target_reads / seq_summary['rate_passed']).astype(int)
    total_spots = seq_summary['total_spot']
    spot_lengths = seq_summary['spot_length']
    pooled_missing_bp = seq_summary['total_bp_remaining']
    for dummy in range(1000):
        current_total_bp = 0
        for sra_id in end_2nds.index:
            pooled_missing_read = (pooled_missing_bp / spot_lengths.loc[sra_id]).astype(int)
            if (end_2nds.loc[sra_id] + pooled_missing_read < total_spots.loc[sra_id]):
                pooled_missing_bp = 0
                end_2nds.loc[sra_id] = end_2nds.loc[sra_id] + pooled_missing_bp
            elif (end_2nds.loc[sra_id] + pooled_missing_read > total_spots.loc[sra_id]):
                pooled_missing_bp = (end_2nds.loc[sra_id] + pooled_missing_read - total_spots.loc[sra_id]) * \
                                    spot_lengths.loc[sra_id]
                end_2nds.loc[sra_id] = total_spots.loc[sra_id]
            current_total_bp += end_2nds.loc[sra_id] * spot_lengths.loc[sra_id]
        all_equal_total_spots = all([e2 == ts for e2, ts in zip(end_2nds, total_spots)])
        is_enough_read = (current_total_bp >= seq_summary['total_bp_remaining'])
        if all_equal_total_spots:
            print('Reached total spots in all SRAs.', flush=True)
            break
        if is_enough_read:
            print('Enough read numbers were assigned for the 2nd round sequence extraction.', flush=True)
            break
    seq_summary['start_2nd'] = start_2nds
    seq_summary['end_2nd'] = end_2nds
    return seq_summary


def print_read_stats(args, seq_summary, max_bp, individual=True):
    print('Target size (--max_bp): {:,} bp'.format(max_bp))
    if args.pfd == 'yes':
        print('Sum of fastq_dump dumped reads: {:,} bp'.format(seq_summary['bp_dumped'].sum()))
        print('Sum of fastq_dump rejected reads: {:,} bp'.format(seq_summary['bp_rejected'].sum()))
        print('Sum of fastq_dump written reads: {:,} bp'.format(seq_summary['bp_written'].sum()))
    if args.fastp:
        print('Sum of fastp input reads: {:,} bp'.format(seq_summary['bp_fastp_in'].sum()))
        print('Sum of fastp output reads: {:,} bp'.format(seq_summary['bp_fastp_out'].sum()))
    if individual:
        print('Individual SRA IDs:', ' '.join(seq_summary['bp_dumped'].index.values))
        read_types = list()
        keys = list()
        if args.pfd == 'yes':
            read_types = read_types + ['fastq_dump dumped reads', 'fastq_dump rejected reads',
                                       'fastq_dump written reads']
            keys = keys + ['bp_dumped', 'bp_rejected', 'bp_written']
        if args.fastp:
            read_types = read_types + ['fastp input reads', 'fastp output reads']
            keys = keys + ['bp_fastp_in', 'bp_fastp_out']
        if len(read_types) > 0:
            for rt, key in zip(read_types, keys):
                values = ['{:,}'.format(s) for s in seq_summary[key].values]
                txt = ' '.join(values)
                print('Individual {} (bp): {}'.format(rt, txt))
    print('')


def write_updated_metadata(args, metadata, seq_summary, sra_id):
    is_sra = (metadata.df['run'] == sra_id)
    if metadata.df.loc[is_sra, 'lib_layout'].values[0] == 'single':
        denom = 1
    elif metadata.df.loc[is_sra, 'lib_layout'].values[0] == 'paired':
        denom = 2
    if args.pfd == 'yes':
        metadata.df.loc[is_sra, 'num_read_fastq_dumped'] = seq_summary['num_dumped'].at[sra_id] / denom
        metadata.df.loc[is_sra, 'num_read_fastq_rejected'] = seq_summary['num_rejected'].at[sra_id] / denom
        metadata.df.loc[is_sra, 'num_read_fastq_written'] = seq_summary['num_written'].at[sra_id] / denom
    if args.fastp:
        metadata.df.loc[is_sra, 'num_read_fastp_input'] = seq_summary['num_fastp_in'].at[sra_id] / denom
        metadata.df.loc[is_sra, 'num_read_fastp'] = seq_summary['num_fastp_out'].at[sra_id] / denom
    metadata.df.loc[is_sra, 'spot_length'] = seq_summary['spot_length'].at[sra_id]
    metadata_output_dir = os.path.join(args.out_dir, 'metadata', 'updated_metadata')
    if not os.path.exists(metadata_output_dir):
        os.makedirs(metadata_output_dir)
    outpath = os.path.join(metadata_output_dir, 'metadata_' + sra_id + '.tsv')
    print('Writing updated metadata: {}'.format(outpath))
    metadata.df.to_csv(outpath, sep='\t', index=False)

def getfastq_metadata(args):
   # assert (args.id is None and args.id_list is None) != (
    #            args.metadata is None), 'Either --id, --id_list or --metadata should be specified.'
    if args.id is not None:
        print('--id is specified. Downloading SRA metadata from Entrez.')
        assert (args.entrez_email != 'aaa@bbb.com'), "Provide your email address. No worry, you won't get spam emails."
        Entrez.email = args.entrez_email
        sra_id = args.id
        search_term = getfastq_search_term(sra_id, args.entrez_additional_search_term)
        print('Entrez search term:', search_term)
        xml_root = getfastq_getxml(search_term)
        metadata = Metadata.from_xml(xml_root)
        print('Filtering SRA entry with --layout:', args.layout)
        layout = get_layout(args, metadata)
        metadata.df = metadata.df.loc[(metadata.df['lib_layout'] == layout), :]
        if args.sci_name is not None:
            print('Filtering SRA entry with --sci_name:', args.sci_name)
            metadata.df = metadata.df.loc[(metadata.df['scientific_name'] == args.sci_name), :]
    if args.id_list is not None:
        print('--id is specified. Downloading SRA metadata from Entrez.')
        assert (args.entrez_email != 'aaa@bbb.com'), "Provide your email address. No worry, you won't get spam emails."
        Entrez.email = args.entrez_email
        sra_id_list = [line.rstrip('\n') for line in open(args.id_list)]
        for sra_id in sra_id_list:
            search_term = getfastq_search_term(sra_id, args.entrez_additional_search_term)
            print('Entrez search term:', search_term)
            xml_root = getfastq_getxml(search_term)
            metadata = Metadata.from_xml(xml_root)
            print('Filtering SRA entry with --layout:', args.layout)
            layout = get_layout(args, metadata)
            metadata.df = metadata.df.loc[(metadata.df['lib_layout'] == layout), :]
            if args.sci_name is not None:
                print('Filtering SRA entry with --sci_name:', args.sci_name)
                metadata.df = metadata.df.loc[(metadata.df['scientific_name'] == args.sci_name), :]

    if args.metadata is not None:
        print('--metadata is specified. Reading existing metadata table.')
        assert args.concat == 'no', '--concat should be set "no" when --metadata is specified.'
        metadata = load_metadata(args)
    return metadata


def is_getfastq_output_present(args, sra_stat, output_dir):
    if args.concat == 'yes':
        prefixes = [args.id, ]
    else:
        prefixes = [sra_stat['sra_id'], ]
    if sra_stat['layout'] == 'single':
        sub_exts = ['', ]
    elif sra_stat['layout'] == 'paired':
        sub_exts = ['_1', '_2']
    exts = ['.amalgkit.fastq.gz', ]
    is_output_present = True
    for prefix, sub_ext, ext in itertools.product(prefixes, sub_exts, exts):
        out_path1 = os.path.join(output_dir, prefix + sub_ext + ext)
        out_path2 = os.path.join(output_dir, prefix + sub_ext + ext + '.safely_removed')
        is_out1 = os.path.exists(out_path1)
        is_out2 = os.path.exists(out_path2)
        if is_out1:
            print('getfastq output detected: {}'.format(out_path1))
        if is_out2:
            print('getfastq output detected: {}'.format(out_path2))
        if (sra_stat['layout'] == 'paired')&(not(is_out1 | is_out2)):
            sub_exts = ['', ]
            for prefix, sub_ext, ext in itertools.product(prefixes, sub_exts, exts):
                out_path_single_like1 = os.path.join(output_dir, prefixes[0] + sub_ext + exts[0])
                out_path_single_like2 = os.path.join(output_dir, prefixes[0] + sub_ext + exts[0] + '.safely_removed')
                is_out1 = os.path.exists(out_path_single_like1)
                is_out2 = os.path.exists(out_path_single_like2)
                txt = 'Single-end getfastq output was generated even though layout = paired: {}'
                if is_out1:
                    print(txt.format(out_path_single_like1))
                if is_out2:
                    print(txt.format(out_path_single_like2))
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

def initialize_seq_summary(metadata):
    sra_ids = metadata.df.loc[:, 'run'].values
    seq_summary = dict()
    keys = ['num_dumped', 'num_rejected', 'num_written', 'num_fastp_in', 'num_fastp_out',
            'bp_dumped', 'bp_rejected', 'bp_written', 'bp_fastp_in', 'bp_fastp_out', 'bp_remaining',
            'bp_still_available', 'layout',
            'spot_length', 'total_spot', 'rate_passed', 'start_1st', 'end_1st', 'start_2nd', 'end_2nd', ]
    for key in keys:
        if key=='layout':
            seq_summary[key] = pandas.Series('', index=sra_ids)
        else:
            seq_summary[key] = pandas.Series(0, index=sra_ids)
    return seq_summary

def sequence_extraction_1st_round(args, sra_stat, output_dir, seq_summary, gz_exe, ungz_exe, total_sra_bp, max_bp):
    offset = 10000  # https://edwards.sdsu.edu/research/fastq-dump/
    start, end = get_range(sra_stat, offset, total_sra_bp, max_bp)
    seq_summary = sequence_extraction(args, sra_stat, start, end, output_dir, seq_summary,
                                      gz_exe, ungz_exe, total_sra_bp)
    txt = 'Time elapsed for 1st-round sequence extraction: {}, {:,} sec'
    print(txt.format(sra_stat['sra_id'], int(time.time() - seq_summary['start_time'])))

    print('\n--- getfastq 1st-round sequence generation report ---')
    print_read_stats(args, seq_summary, max_bp)
    seq_summary['total_bp_remaining'] = seq_summary['bp_remaining'].sum()
    seq_summary['percent_remaining'] = seq_summary['total_bp_remaining'] / max_bp * 100
    if args.fastp:
        seq_summary['total_bp_out'] = seq_summary['bp_fastp_out'].sum()
    else:
        seq_summary['total_bp_out'] = seq_summary['bp_written'].sum()
    txt = '{:.2f}% of reads were obtained in the 1st-round sequence generation.'
    print(txt.format(100-seq_summary['percent_remaining']), flush=True)
    txt = 'The amount of generated reads were {:.2f}% ({:,}/{:,}) smaller than the target size (tol={}%).'
    print(txt.format(seq_summary['percent_remaining'], seq_summary['total_bp_out'], max_bp, args.tol), flush=True)
    return seq_summary

def sequence_extraction_2st_round(args, sra_stat, output_dir, seq_summary, gz_exe, ungz_exe, total_sra_bp, max_bp):
    print('Starting the 2nd-round sequence extraction to compensate it.')
    seq_summary = calc_2nd_ranges(seq_summary)
    ext_main = '.amalgkit.fastq.gz'
    ext_1st_tmp = '.amalgkit_1st.fastq.gz'

    print('')
    seq_summary['start_time_2nd'] = time.time()
    sra_id = sra_stat['sra_id']
    start = seq_summary['start_2nd'].loc[sra_id]
    end = seq_summary['end_2nd'].loc[sra_id]
    layout = sra_stat['layout']
    if (start >= end):
        txt = '{}: All spots have been extracted in the 1st trial. Cancelling the 2nd trial. start={:,}, end={:,}'
        print(txt.format(sra_id, start, end))
        return seq_summary
    rename_fastq(sra_stat, output_dir, inext=ext_main, outext=ext_1st_tmp)
    seq_summary = sequence_extraction(args, sra_stat, start, end, output_dir, seq_summary,
                                      gz_exe, ungz_exe, total_sra_bp)
    if (layout == 'single'):
        subexts = ['']
    elif (layout == 'paired'):
        subexts = ['_1', '_2']
    for subext in subexts:
        added_path = os.path.join(output_dir, sra_id + subext + ext_1st_tmp)
        adding_path = os.path.join(output_dir, sra_id + subext + ext_main)
        assert os.path.exists(added_path), 'Dumped fastq not found: ' + added_path
        assert os.path.exists(adding_path), 'Dumped fastq not found: ' + adding_path
        os.system('cat "' + adding_path + '" >> "' + added_path + '"')
        os.remove(adding_path)
        os.rename(added_path, adding_path)
    txt = 'Time elapsed for 2nd-round sequence extraction: {}, {:,} sec'
    print(txt.format(sra_stat['sra_id'], int(time.time() - seq_summary['start_time_2nd'])))
    return seq_summary

def sequence_extraction_private(i, metadata, sra_stat, seq_summary, output_dir, args):
    for col in ['read1_path','read2_path']:
        path_from = metadata.df.at[i,col]
        path_to = os.path.join(output_dir, os.path.basename(path_from))
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
        run_fastp(sra_stat, args, output_dir, seq_summary)
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
    outext = '.amalgkit.fastq.gz'
    rename_fastq(sra_stat, output_dir, inext, outext)

def getfastq_main(args):
    # sra_dir = os.path.join(os.path.expanduser("~"), 'ncbi/public/sra')
    gz_exe, ungz_exe = check_getfastq_dependency(args)
    metadata = getfastq_metadata(args)
    assert metadata.df.shape[0] > 0, 'No SRA entry found. Make sure if --id is compatible with --sci_name and --layout.'
    metadata = remove_experiment_without_run(metadata)
    print('SRA IDs:', ' '.join(metadata.df['run'].tolist()))
    max_bp = int(args.max_bp.replace(',', ''))
    num_sra = metadata.df.shape[0]
    num_bp_per_sra = int(max_bp / num_sra)
    if metadata.df.loc[:,'total_bases'].isna().any():
        raise Exception('Empty value(s) of total_bases were detected in the metadata table.')
    total_sra_bp = metadata.df.loc[:,'total_bases'].sum()
    print('Number of SRAs:', num_sra)
    print('Total target size (--max_bp):', "{:,}".format(max_bp), 'bp')
    print('Total SRA size:', "{:,}".format(total_sra_bp), 'bp')
    print('Target size per SRA:', "{:,}".format(num_bp_per_sra), 'bp')
    for i in metadata.df.index:
        txt = 'Individual SRA size of {}: {:,}'
        print(txt.format(metadata.df.at[i, 'run'], metadata.df.at[i, 'total_bases']))
    for i in metadata.df.index:
        print('')
        seq_summary = initialize_seq_summary(metadata)
        seq_summary['start_time'] = time.time()
        sra_id = metadata.df.at[i, 'run']
        sra_stat = get_sra_stat(sra_id, metadata, num_bp_per_sra)
        output_dir = set_getfastq_directories(args, sra_id)
        if (is_getfastq_output_present(args, sra_stat, output_dir)) & (args.redo == 'no'):
            if metadata.df.shape[0] == 1:
                print('Output file(s) detected. Exiting.  Set "--redo yes" for reanalysis.')
                sys.exit()
            else:
                print('Output file(s) detected. Skipping {}.  Set "--redo yes" for reanalysis.'.format(sra_id))
                continue
        remove_old_intermediate_files(sra_id=sra_id, work_dir=output_dir)
        print('SRA ID:', sra_stat['sra_id'])
        print('Library layout:', sra_stat['layout'])
        print('Number of reads:', "{:,}".format(sra_stat['total_spot']))
        print('Single/Paired read length:', sra_stat['spot_length'], 'bp')
        print('Total bases:', "{:,}".format(int(metadata.df.loc[i, 'total_bases'])), 'bp')
        flag_private_file = False
        if 'private_file' in metadata.df.columns:
            if metadata.df.at[i,'private_file']=='yes':
                print('Processing {} as private data. --max_bp is disabled.'.format(sra_id), flush=True)
                flag_private_file = True
                sequence_extraction_private(i, metadata, sra_stat, seq_summary, output_dir, args)
        if not flag_private_file:
            print('Processing {} as publicly available data from SRA.'.format(sra_id), flush=True)
            download_sra(metadata, sra_stat, args, output_dir, overwrite=False)
            seq_summary = sequence_extraction_1st_round(args, sra_stat, output_dir, seq_summary,
                                                        gz_exe, ungz_exe, total_sra_bp, max_bp)
            if (seq_summary['percent_remaining'] < args.tol):
                print('Enough data were obtained in the 1st-round sequence extraction. Proceeding without the 2nd round.')
            else:
                seq_summary = sequence_extraction_2st_round(args, sra_stat, output_dir, seq_summary, gz_exe, ungz_exe,
                                                            total_sra_bp, max_bp)
    for sra_id in metadata.df.loc[:,'run'].values:
        write_updated_metadata(args, metadata, seq_summary, sra_id)
    print('')
    if args.concat == 'yes':
        concat_fastq(args, metadata, output_dir, num_bp_per_sra)
    if args.remove_sra == 'yes':
        remove_sra_files(metadata, sra_dir=output_dir)
    else:
        if args.pfd == 'yes':
            print('SRA files not removed:', output_dir)
    print('\n--- getfastq final report ---')
    print_read_stats(args, seq_summary, max_bp)