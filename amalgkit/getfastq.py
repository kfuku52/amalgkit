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
            raise Exception('\<Error> found in the xml. Search term: '+search_term)
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
                outfile = args.id + re.sub('.*(_[1-2])', '\g<1>', infile)
            elif args.id_list is not None:
                outfile = os.path.basename(args.id_list) + re.sub('.*(_[1-2])', '\g<1>', infile)
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
        if args.remove_tmp:
            for i in metadata.df.index:
                sra_id = metadata.df.loc[i, 'run']
                sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
                ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
                remove_intermediate_files(sra_stat, ext=ext, work_dir=output_dir)
        return None


def remove_sra_files(metadata, sra_id, amalgkit_out_dir):
    for sra_id in metadata.df['run']:
        sra_pattern = os.path.join(os.path.realpath(amalgkit_out_dir), 'getfastq', sra_id, sra_id + '.sra*')
        path_downloaded_sras = glob.glob(sra_pattern)
        if len(path_downloaded_sras) > 0:
            for path_downloaded_sra in path_downloaded_sras:
                print('Deleting: {}'.format(path_downloaded_sra))
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


def download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
    path_downloaded_sra = os.path.join(work_dir, sra_stat['sra_id'] + '.sra')
    individual_sra_tmp_dir = os.path.join(work_dir, sra_stat['sra_id'] + '/')

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
            print("Trying to fetch {} from {}: {}".format(sra_id, sra_source_name, sra_sources[sra_source_name]))
            if str(sra_sources[sra_source_name])=='nan':
                sys.stderr.write("Skipping. No URL for {}.\n".format(sra_source_name))
                continue
            try:
                urllib.request.urlretrieve(str(sra_sources[sra_source_name]), path_downloaded_sra)
                if os.path.exists(path_downloaded_sra):
                    is_sra_download_completed = True
                    print('SRA file was downloaded with urllib.request from {}'.format(sra_source_name), flush=True)
                    break
            except urllib.error.URLError:
                sys.stderr.write("urllib.request failed SRA download from {}.\n".format(sra_source_name))
        if not is_sra_download_completed:
            sys.stderr.write("Exhausted all sources of download.\n")
        else:
            assert os.path.exists(path_downloaded_sra), 'SRA file download failed: ' + sra_stat['sra_id']
            return

    if not os.path.exists(path_downloaded_sra):
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
        if (prefetch_out.returncode):
            sys.stderr.write("prefetch did not finish safely. Trying prefetch again.\n")
            prefetch_command = [args.prefetch_exe, '--force', 'no', '--max-size', '100G',
                                sra_stat['sra_id']]
            prefetch_out = subprocess.run(prefetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print('prefetch stdout:')
            print(prefetch_out.stdout.decode('utf8'))
            print('prefetch stderr:')
            print(prefetch_out.stderr.decode('utf8'))
            if (prefetch_out.returncode !=0):
                sys.stderr.write("Again, prefetch did not finish safely.\n")
    # Move files downloaded by prefetch. This is necessary because absolute path didn't work for prefetch --output-directory
    if os.path.exists(os.path.join('./', sra_stat['sra_id'] + '/', sra_stat['sra_id'] + '.sra')):
        subprocess.run(['mv', os.path.join('./', sra_stat['sra_id'] + '/', sra_stat['sra_id'] + '.sra'), path_downloaded_sra],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        shutil.rmtree(os.path.join('./', sra_stat['sra_id'] + '/'))
    elif os.path.exists(os.path.expanduser(os.path.join('~/ncbi/public/sra/', sra_stat['sra_id'] + '.sra'))):
        subprocess.run(
            ['mv', os.path.expanduser(os.path.join('~/ncbi/public/sra/', sra_stat['sra_id'] + '.sra')), path_downloaded_sra],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Move files downloaded by ascp
    if os.path.exists(os.path.join(individual_sra_tmp_dir, sra_stat['sra_id'] + '.sra')):
        subprocess.run(['mv', os.path.join(individual_sra_tmp_dir, sra_stat['sra_id'] + '.sra'), path_downloaded_sra],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        shutil.rmtree(individual_sra_tmp_dir)
    assert os.path.exists(path_downloaded_sra), 'SRA file download failed: ' + sra_stat['sra_id']


def check_getfastq_dependency(args):
    if args.pfd:
        test_pfd = subprocess.run([args.pfd_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_pfd.returncode == 0), "parallel-fastq-dump PATH cannot be found: " + args.pfd_exe
        # commented out because prefetch is often not activatable in containers and no longer strictly required for getfastq.
        #test_prefetch = subprocess.run([args.prefetch_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #assert (test_prefetch.returncode == 0), "prefetch (SRA toolkit) PATH cannot be found: " + args.prefetch_exe
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


def get_getfastq_run_dir(args, sra_id):
    amalgkit_out_dir = os.path.realpath(args.out_dir)
    if args.id is not None:
        run_output_dir = os.path.join(amalgkit_out_dir, 'getfastq', args.id)
    elif (args.id_list is not None) & args.concat:
        run_output_dir = os.path.join(amalgkit_out_dir, 'getfastq', os.path.basename(args.id_list))
    elif (args.id_list is not None) & (not args.concat):
        run_output_dir = os.path.join(amalgkit_out_dir, 'getfastq', sra_id)
    else:
        run_output_dir = os.path.join(amalgkit_out_dir, 'getfastq', sra_id)
    if not os.path.exists(run_output_dir):
        os.makedirs(run_output_dir)
    return run_output_dir

def run_pfd(sra_stat, args, metadata, start, end):
    path_downloaded_sra = os.path.join(sra_stat['output_dir'], sra_stat['sra_id'] + '.sra')
    pfd_command = ['parallel-fastq-dump', '-t', str(args.threads), '--minReadLen', str(args.min_read_length),
                   '--qual-filter-1',
                   '--skip-technical', '--split-3', '--clip', '--gzip', '--outdir', sra_stat['output_dir'],
                   '--tmpdir', sra_stat['output_dir']]
    print('Total sampled bases:', "{:,}".format(sra_stat['spot_length'] * (end - start + 1)), 'bp')
    pfd_command = pfd_command + ['--minSpotId', str(int(start)), '--maxSpotId', str(int(end))]
    # If sra_stat['sra_id'], not path_downloaded_sra, is provided, pfd couldn't find pre-downloaded .sra files
    # and start downloading it to $HOME/ncbi/public/sra/
    pfd_command = pfd_command + ['-s', path_downloaded_sra]
    # pfd_command = pfd_command + ['-s', sra_stat['sra_id']]
    print('Command:', ' '.join(pfd_command))
    pfd_out = subprocess.run(pfd_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if args.pfd_print:
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
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_stat['sra_id']].values[0]
    metadata.df.at[ind_sra,'num_dumped'] += sum(nd)
    metadata.df.at[ind_sra,'num_rejected'] += sum(nr)
    metadata.df.at[ind_sra,'num_written'] += sum(nw)
    metadata.df.at[ind_sra,'bp_dumped'] += sum(nd) * sra_stat['spot_length']
    metadata.df.at[ind_sra,'bp_rejected'] += sum(nr) * sra_stat['spot_length']
    metadata.df.at[ind_sra,'bp_written'] += sum(nw) * sra_stat['spot_length']
    paired_fastq_files = [
        os.path.join(sra_stat['output_dir'], sra_stat['sra_id'] + '_1.fastq.gz'),
        os.path.join(sra_stat['output_dir'], sra_stat['sra_id'] + '_2.fastq.gz'),
    ]
    single_fastq_file = os.path.join(sra_stat['output_dir'], sra_stat['sra_id'] + '.fastq.gz')
    is_paired_end = all([os.path.exists(f) for f in paired_fastq_files])
    is_unpaird_file = os.path.exists(single_fastq_file)
    if (not is_paired_end) & is_unpaird_file:
        is_single_end = True
    else:
        is_single_end = False
    assert is_paired_end | is_single_end, 'No fastq file was generated.'
    assert (not is_paired_end) | (not is_unpaird_file), 'Paired-end/single-end layout cannot be determined.'
    if is_paired_end & is_unpaird_file:
        print('layout = {}; Deleting unpaired file: {}'.format(sra_stat['layout'], unpaired_file))
        os.remove(unpaired_file)
    if is_single_end & (sra_stat['layout'] == 'paired'):
        txt = 'Single-end fastq was generated even though layout in the metadata = {}. '
        txt += 'This sample will be treated as single-end reads: {}\n'
        txt = txt.format(sra_stat['layout'], sra_stat['sra_id'])
        sys.stderr.write(txt)
        sra_stat['layout'] = 'single'
    if is_paired_end & (sra_stat['layout'] == 'single'):
        txt = 'Paired-end fastq was generated even though layout in the metadata = {}. '
        txt += 'This sample will be treated as paired-end reads: {}\n'
        txt = txt.format(sra_stat['layout'], sra_stat['sra_id'])
        sys.stderr.write(txt)
        sra_stat['layout'] = 'paired'
    metadata.df.at[ind_sra,'layout_amalgkit'] = sra_stat['layout']
    return metadata,sra_stat

def run_fastp(sra_stat, args, output_dir, metadata):
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
    if args.fastp_print:
        print('fastp stdout:')
        print(fp_out.stdout.decode('utf8'))
        print('fastp stderr:')
        print(fp_out.stderr.decode('utf8'))
    if args.remove_tmp:
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
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_stat['sra_id']].values[0]
    metadata.df.at[ind_sra,'num_fastp_in'] += sum(num_in)
    metadata.df.at[ind_sra,'num_fastp_out'] += sum(num_out)
    metadata.df.at[ind_sra,'bp_fastp_in'] += sum(bp_in)
    metadata.df.at[ind_sra,'bp_fastp_out'] += sum(bp_out)
    return metadata


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
    for i in numpy.arange(sra_target_reads.shape[0]):
        if numpy.isnan(rate_obtained[i]):
            sra_target_reads[i] = (sra_target_bp[i]/spot_lengths[i]).astype(int)+1 # If no read was extracted in 1st.
        else:
            sra_target_reads[i] = ((sra_target_bp[i]/spot_lengths[i])/rate_obtained[i]).astype(int)+1
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


def print_read_stats(args, metadata, g, sra_stat=None, individual=False):
    if sra_stat is None:
        df = metadata.df
        print('Target size (--max_bp): {:,} bp'.format(g['max_bp']))
    else:
        df = metadata.df.loc[(metadata.df['run']==sra_stat['sra_id']),:]
        print('Individual target size: {:,} bp'.format(g['num_bp_per_sra']))
    if args.pfd:
        print('Sum of fastq_dump dumped reads: {:,} bp'.format(df['bp_dumped'].sum()))
        print('Sum of fastq_dump rejected reads: {:,} bp'.format(df['bp_rejected'].sum()))
        print('Sum of fastq_dump written reads: {:,} bp'.format(df['bp_written'].sum()))
    if args.fastp:
        print('Sum of fastp input reads: {:,} bp'.format(df['bp_fastp_in'].sum()))
        print('Sum of fastp output reads: {:,} bp'.format(df['bp_fastp_out'].sum()))
    if individual:
        print('Individual SRA IDs:', ' '.join(df['run'].values))
        read_types = list()
        keys = list()
        if args.pfd:
            read_types = read_types + ['fastq_dump dumped reads', 'fastq_dump rejected reads',
                                       'fastq_dump written reads']
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
        xml_root = getfastq_getxml(search_term)
        metadata = Metadata.from_xml(xml_root)
        print('Filtering SRA entry with --layout:', args.layout)
        layout = get_layout(args, metadata)
        metadata.df = metadata.df.loc[(metadata.df['lib_layout'] == layout), :]
        if args.sci_name is not None:
            print('Filtering SRA entry with --sci_name:', args.sci_name)
            metadata.df = metadata.df.loc[(metadata.df['scientific_name'] == args.sci_name), :]
    if args.id_list is not None:
        print('--id_list is specified. Downloading SRA metadata from Entrez.')
        Entrez.email = args.entrez_email
        sra_id_list = [line.rstrip('\n') for line in open(args.id_list) if not line.startswith('#')]
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
        assert args.concat == False, '--concat should be set "no" with the input from --metadata.'
        metadata = load_metadata(args)
    metadata.df['total_bases'] = metadata.df.loc[:,'total_bases'].replace('', numpy.nan).astype(float)
    metadata.df['spot_length'] = metadata.df.loc[:, 'spot_length'].replace('', numpy.nan).astype(float)
    return metadata


def is_getfastq_output_present(args, sra_stat):
    if args.concat:
        if args.id is not None:
            prefixes = [args.id, ]
        elif args.id_list is not None:
            prefixes = [os.path.basename(args.id_list), ]
    else:
        prefixes = [sra_stat['sra_id'], ]
    if sra_stat['layout'] == 'single':
        sub_exts = ['', ]
    elif sra_stat['layout'] == 'paired':
        sub_exts = ['_1', '_2']
    exts = ['.amalgkit.fastq.gz', ]
    is_output_present = True
    for prefix, sub_ext, ext in itertools.product(prefixes, sub_exts, exts):
        out_path1 = os.path.join(sra_stat['output_dir'], prefix + sub_ext + ext)
        out_path2 = os.path.join(sra_stat['output_dir'], prefix + sub_ext + ext + '.safely_removed')
        is_out1 = os.path.exists(out_path1)
        is_out2 = os.path.exists(out_path2)
        if is_out1:
            print('getfastq output detected: {}'.format(out_path1))
        if is_out2:
            print('getfastq output detected: {}'.format(out_path2))
        if (sra_stat['layout'] == 'paired')&(not(is_out1 | is_out2)):
            sub_exts = ['', ]
            for prefix, sub_ext, ext in itertools.product(prefixes, sub_exts, exts):
                out_path_single_like1 = os.path.join(sra_stat['output_dir'], prefixes[0] + sub_ext + exts[0])
                out_path_single_like2 = os.path.join(sra_stat['output_dir'], prefixes[0] + sub_ext + exts[0] + '.safely_removed')
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

def initialize_columns(metadata, g):
    keys = ['num_dumped', 'num_rejected', 'num_written', 'num_fastp_in', 'num_fastp_out','bp_amalgkit',
            'bp_dumped', 'bp_rejected', 'bp_written', 'bp_fastp_in', 'bp_fastp_out', 'bp_discarded',
            'bp_still_available', 'bp_specified_for_extraction','rate_obtained','layout_amalgkit',
            'time_start_1st', 'time_end_1st', 'time_start_2nd', 'time_end_2nd',
            'spot_start_1st', 'spot_end_1st', 'spot_start_2nd', 'spot_end_2nd', ]
    for key in keys:
        if key=='layout_amalgkit':
            metadata.df.loc[:,key] = ''
        elif key=='rate_obtained':
            metadata.df.loc[:, key] = numpy.nan
        else:
            metadata.df.loc[:,key] = 0
    metadata.df.loc[:, 'bp_until_target_size'] = g['num_bp_per_sra']
    cols = ['total_spots','total_bases','size','nominal_length','nominal_sdev','spot_length']
    for col in cols:
        if any([ dtype in str(metadata.df[col].dtype) for dtype in ['str','object'] ]):
            metadata.df[col] = metadata.df.loc[:,col].str.replace('^$', 'nan', regex=True).astype(float)
    return metadata

def sequence_extraction(args, sra_stat, metadata, g, start, end):
    sra_id = sra_stat['sra_id']
    ind_sra = metadata.df.index[metadata.df.loc[:,'run'] == sra_id].values[0]
    if args.pfd:
        metadata,sra_stat = run_pfd(sra_stat, args, metadata, start, end)
        bp_discarded = metadata.df.at[ind_sra,'bp_dumped'] - metadata.df.at[ind_sra,'bp_written']
        metadata.df.at[ind_sra,'bp_discarded'] += bp_discarded
    no_read_written = (metadata.df.loc[(metadata.df.loc[:,'run']==sra_id),'num_written'].values[0]==0)
    if no_read_written:
        return metadata
    if args.fastp:
        metadata = run_fastp(sra_stat, args, sra_stat['output_dir'], metadata)
        bp_discarded = metadata.df.at[ind_sra,'bp_dumped'] - metadata.df.at[ind_sra,'bp_fastp_out']
        metadata.df.at[ind_sra,'bp_discarded'] += bp_discarded
    if args.read_name == 'trinity':
        rename_reads(sra_stat, args, sra_stat['output_dir'], g['gz_exe'], g['ungz_exe'])
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=sra_stat['output_dir'])
    outext = '.amalgkit.fastq.gz'
    rename_fastq(sra_stat, sra_stat['output_dir'], inext, outext)
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
    offset = 10000  # https://edwards.sdsu.edu/research/fastq-dump/
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
        rename_fastq(sra_stat, sra_stat['output_dir'], inext=ext_main, outext=ext_1st_tmp)
    metadata = sequence_extraction(args, sra_stat, metadata, g, start, end)
    if (layout == 'single'):
        subexts = ['']
    elif (layout == 'paired'):
        subexts = ['_1', '_2']
    for subext in subexts:
        added_path = os.path.join(sra_stat['output_dir'], sra_id + subext + ext_1st_tmp)
        adding_path = os.path.join(sra_stat['output_dir'], sra_id + subext + ext_main)
        assert os.path.exists(added_path), 'Dumped fastq not found: ' + added_path
        assert os.path.exists(adding_path), 'Dumped fastq not found: ' + adding_path
        os.system('cat "' + adding_path + '" >> "' + added_path + '"')
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
        path_to = os.path.join(sra_stat['output_dir'], os.path.basename(path_from))
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
        metadata = run_fastp(sra_stat, args, sra_stat['output_dir'], metadata)
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=sra_stat['output_dir'])
    outext = '.amalgkit.fastq.gz'
    rename_fastq(sra_stat, sra_stat['output_dir'], inext, outext)
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

def initialize_global_params(args, metadata, gz_exe, ungz_exe):
    g = dict()
    g['start_time'] = time.time()
    g['max_bp'] = int(args.max_bp.replace(',', ''))
    g['num_sra'] = metadata.df.shape[0]
    g['num_bp_per_sra'] = int(g['max_bp'] / g['num_sra'])
    g['total_sra_bp'] = metadata.df.loc[:,'total_bases'].sum()
    g['gz_exe'] = gz_exe
    g['ungz_exe'] = ungz_exe
    print('Number of SRAs to be processed: {:,}'.format(g['num_sra']))
    print('Total target size (--max_bp): {:,} bp'.format(g['max_bp']))
    print('The sum of SRA sizes: {:,} bp'.format(g['total_sra_bp']))
    print('Target size per SRA: {:,} bp'.format(g['num_bp_per_sra']))
    return g

def getfastq_main(args):
    gz_exe, ungz_exe = check_getfastq_dependency(args)
    metadata = getfastq_metadata(args)
    metadata = remove_experiment_without_run(metadata)
    metadata = check_metadata_validity(metadata)
    g = initialize_global_params(args, metadata, gz_exe, ungz_exe)
    metadata = initialize_columns(metadata, g)
    flag_private_file = False
    flag_any_output_file_present = False
    # 1st round sequence extraction
    for i in metadata.df.index:
        print('')
        sra_id = metadata.df.at[i, 'run']
        print('Processing SRA ID: {}'.format(sra_id))
        sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
        sra_stat['output_dir'] = get_getfastq_run_dir(args, sra_id)
        if (is_getfastq_output_present(args, sra_stat)) & (args.redo == False):
            print('Output file(s) detected. Skipping {}.  Set "--redo yes" for reanalysis.'.format(sra_id))
            flag_any_output_file_present =True
            continue
        remove_old_intermediate_files(sra_id=sra_id, work_dir=sra_stat['output_dir'])
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
            download_sra(metadata, sra_stat, args, sra_stat['output_dir'], overwrite=False)
            metadata = sequence_extraction_1st_round(args, sra_stat, metadata, g)
    # 2nd round sequence extraction
    if (not flag_private_file) & (not flag_any_output_file_present):
        g['rate_obtained_1st'] = metadata.df.loc[:,'bp_amalgkit'].sum() / g['max_bp']
        if (g['rate_obtained_1st'] < (args.tol*0.01)):
            print('Enough data were obtained in the 1st-round sequence extraction. Proceeding without the 2nd round.')
        else:
            txt = 'Only {:,.2f}% ({:,}/{:,}) of the target size (--max_bp) was obtained in the 1st round. Proceeding to the 2nd round read extraction.'
            print(txt.format(g['rate_obtained_1st']*100, metadata.df.loc[:,'bp_amalgkit'].sum(), g['max_bp']), flush=True)
            metadata = calc_2nd_ranges(metadata)
            for i in metadata.df.index:
                sra_id = metadata.df.at[i, 'run']
                sra_stat = get_sra_stat(sra_id, metadata, g['num_bp_per_sra'])
                sra_stat['output_dir'] = get_getfastq_run_dir(args, sra_id)
                metadata = sequence_extraction_2nd_round(args, sra_stat, metadata, g)
        g['rate_obtained_2nd'] = metadata.df.loc[:, 'bp_amalgkit'].sum() / g['max_bp']
        txt = '2nd round read extraction improved % bp from {:,.2f}% to {:,.2f}%'
        print(txt.format(g['rate_obtained_1st']*100, g['rate_obtained_2nd']*100), flush=True)
    # Postprocessing
    if (not flag_any_output_file_present):
        print('')
        if args.concat:
            concat_fastq(args, metadata, sra_stat['output_dir'], g)
        if args.remove_sra:
            for i in metadata.df.index:
                sra_id = metadata.df.at[i, 'run']
                remove_sra_files(metadata, sra_id=sra_id, amalgkit_out_dir=args.out_dir)
        else:
            print('SRA files not removed: {}'.format(sra_stat['output_dir']))
        print('\n--- getfastq final report ---')
        print_read_stats(args, metadata, g, sra_stat=None, individual=True)