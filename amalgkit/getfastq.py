from Bio import Entrez
from urllib.error import HTTPError
import itertools
import numpy, pandas
import time, datetime, lxml, subprocess, os, shutil, gzip, glob, sys
from amalgkit.metadata import Metadata
from amalgkit.util import *

def getfastq_search_term(ncbi_id, additional_search_term=None):
    # https://www.ncbi.nlm.nih.gov/books/NBK49540/
    if additional_search_term is None:
        search_term = ncbi_id
    else:
        search_term = ncbi_id+' AND '+additional_search_term
    #search_term = '"'+ncbi_id+'"'+'[Accession]'
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
    for i in numpy.arange(numpy.ceil(num_record//retmax)+1):
        start = int(i*retmax)
        end = int(((i+1)*retmax)-1) if num_record >= int(((i+1)*retmax)-1) else num_record
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
    if (total_sra_bp<=max_bp):
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
    return start,end

def concat_fastq(args, metadata, output_dir, num_bp_per_sra):
    import re
    layout = get_layout(args, metadata)
    inext = '.amalgkit.fastq.gz'
    infiles = list()
    for sra_id in metadata.df.loc[:,'run']:
        infiles.append([ f for f in os.listdir(output_dir) if (f.endswith(inext))&(f.startswith(sra_id)) ])
    infiles = [ item for sublist in infiles for item in sublist ]
    num_inext_files = len(infiles)
    if (layout=='single')&(num_inext_files==1):
        print('Only 1', inext, 'file was detected. No concatenation will happen.')
        outfile = sra_id+inext
        if infiles[0]!=outfile:
            print('Replacing ID in the output file name:', infiles[0], outfile)
            infile_path = os.path.join(output_dir, infiles[0])
            outfile_path = os.path.join(output_dir, outfile)
            os.rename(infile_path, outfile_path)
        return None
    elif (layout=='paired')&(num_inext_files==2):
        print('Only 1 pair of', inext, 'files were detected. No concatenation will happen.')
        for infile in infiles:
            outfile = sra_id+re.sub('.*(_[1-2])', '\g<1>', infile)
            if infile!=outfile:
                print('Replacing ID in the output file name:', infile, outfile)
                infile_path = os.path.join(output_dir, infile)
                outfile_path = os.path.join(output_dir, outfile)
                os.rename(infile_path, outfile_path)
        return None
    else:
        print('Concatenating files with the extension:', inext)
        outext = '.amalgkit.fastq.gz'
        if layout=='single':
            subexts = ['',]
        elif layout=='paired':
            subexts = ['_1','_2',]
        for subext in subexts:
            infiles = metadata.df['run'].replace('$', subext+inext, regex=True)
            outfile_path = os.path.join(output_dir, sra_id+subext+outext)
            if os.path.exists(outfile_path):
                os.remove(outfile_path)
            #with gzip.open(outfile_path, 'wb') as outfile:
            #    for each_infile in infiles:
            #        with gzip.open(os.path.join(args.work_dir, each_infile), 'rb') as infile:
            #            shutil.copyfileobj(infile, outfile) # unacceptably slow
            if os.path.exists(outfile_path):
                os.remove(outfile_path)
            for infile in infiles:
                infile_path = os.path.join(output_dir, infile)
                assert os.path.exists(infile_path), 'Dumped fastq not found: '+infile_path
                os.system('cat "'+infile_path+'" >> "'+outfile_path+'"')
            print('')
        if args.remove_tmp=='yes':
            for i in metadata.df.index:
                sra_id = metadata.df.loc[i,'run']
                sra_stat = get_sra_stat(sra_id, metadata, num_bp_per_sra)
                ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
                remove_intermediate_files(sra_stat, ext=ext, work_dir=output_dir)
        return None

def remove_sra_files(metadata, sra_dir):
    for sra_id in metadata.df['run']:
        sra_pattern = os.path.join(sra_dir, sra_id+'.sra*')
        sra_paths = glob.glob(sra_pattern)
        if len(sra_paths)>0:
            for sra_path in sra_paths:
                print('Deleting:', sra_path)
                os.remove(sra_path)
        else:
            print('SRA file not found. Pattern searched:', sra_pattern)
    print('')

def get_layout(args, metadata):
    if args.layout=='auto':
        layouts = metadata.df['lib_layout'].unique().tolist()
        if (len(layouts)!=1):
            print('Detected multiple layouts in the metadata:', layouts)
        layout = 'paired' if 'paired' in layouts else 'single'
    else:
        layout = args.layout
    return layout

def remove_old_intermediate_files(metadata, work_dir):
    for i in metadata.df.index:
        sra_id = metadata.df.loc[i,'run']
        old_files = os.listdir(work_dir)
        files = [ f for f in old_files if (f.startswith(sra_id))&(not f.endswith('.sra'))&(os.path.isfile(os.path.join(work_dir,f))) ]
        for f in files:
            f_path = os.path.join(work_dir, f)
            print('Deleting old intermediate file:', f_path)
            os.remove(f_path)

def remove_intermediate_files(sra_stat, ext, work_dir):
    file_paths = list()
    if sra_stat['layout']=='single':
        file_paths.append(os.path.join(work_dir,sra_stat['sra_id']+ext))
    elif sra_stat['layout']=='paired':
        for i in [1,2]:
            file_paths.append(os.path.join(work_dir,sra_stat['sra_id']+'_'+str(i)+ext))
    for file_path in file_paths:
        if os.path.exists(file_path):
            print('Deleting intermediate file:', file_path)
            os.remove(file_path)
        else:
            print('Tried to delete but file not found:', file_path)

def get_newest_intermediate_file_extension(sra_stat, work_dir):
    # Order is important in this list. More downstream should come first.
    extensions = ['.amalgkit.fastq.gz','.rename.fastq.gz','.fastp.fastq.gz','.fastq.gz']
    if sra_stat['layout']=='single':
        subext = ''
    elif sra_stat['layout']=='paired':
        subext = '_1'
    files = os.listdir(work_dir)
    for ext in extensions:
        if any([ f==sra_stat['sra_id']+subext+ext for f in files ]):
            ext_out = ext
            break
    assert 'ext_out' in locals(), 'None of expected extensions ('+' '.join(extensions)+') found in '+work_dir
    return ext_out

def download_sra(sra_stat, args, work_dir, overwrite=False):
    sra_path = os.path.join(work_dir, sra_stat['sra_id']+'.sra')
    individual_sra_tmp_dir = os.path.join(work_dir, sra_stat['sra_id']+'/')
    if os.path.exists(sra_path):
        print('Previously-downloaded sra file was detected.')
        if (overwrite):
            print('Removing', sra_path)
            print('New sra file will be downloaded.')
            os.remove(sra_path)
    else:
        print('Previously-downloaded sra file was not detected. New sra file will be downloaded.')
    if (args.ascp=='yes')&(not os.path.exists(sra_path)):
        print('Trying to download the SRA file using ascp.')
        sra_site = 'anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/'+sra_stat['sra_id'][0:3]+'/'+sra_stat['sra_id'][0:6]+'/'+sra_stat['sra_id']+'/'+sra_stat['sra_id']+'.sra'
        ascp_command = [args.ascp_exe, '-v', '-i', '"'+args.ascp_key+'"', '-k', '1', '-T', '-l', '300m', sra_site, '"'+work_dir+'"']
        print('Command:', ' '.join(ascp_command))
        ascp_out = subprocess.run(ascp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('ascp stdout:')
        print(ascp_out.stdout.decode('utf8'))
        print('ascp stderr:')
        print(ascp_out.stderr.decode('utf8'))
    if not os.path.exists(sra_path):
        print('Trying to download the SRA file using prefetch (fasp protocol).')
        if os.path.exists(individual_sra_tmp_dir):
            shutil.rmtree(individual_sra_tmp_dir)
        prefetch_command = [args.prefetch_exe, '--force', 'no', '--transport', 'fasp', '--max-size', '100G',
                            '--output-directory', './', sra_stat['sra_id']]
        print('Command:', ' '.join(prefetch_command))
        prefetch_out = subprocess.run(prefetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('prefetch stdout:')
        print(prefetch_out.stdout.decode('utf8'))
        print('prefetch stderr:')
        print(prefetch_out.stderr.decode('utf8'))
    if (not os.path.exists(sra_path))&(not os.path.exists(os.path.join(work_dir, sra_stat['sra_id']+'/', sra_stat['sra_id']+'.sra'))):
        print('Trying to download the SRA file using prefetch (http protocol).')
        if os.path.exists(individual_sra_tmp_dir):
            shutil.rmtree(individual_sra_tmp_dir)
        prefetch_command = [args.prefetch_exe, '--force', 'no', '--transport', 'http', '--max-size', '100G',
                            '--output-directory', './', sra_stat['sra_id']]
        print('Command:', ' '.join(prefetch_command))
        prefetch_out = subprocess.run(prefetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('prefetch stdout:')
        print(prefetch_out.stdout.decode('utf8'))
        print('prefetch stderr:')
        print(prefetch_out.stderr.decode('utf8'))
    # Move files downloaded by prefetch. This is necessary because absolute path didn't work for prefetch --output-directory
    if os.path.exists(os.path.join('./', sra_stat['sra_id']+'/', sra_stat['sra_id']+'.sra')):
        subprocess.run(['mv', os.path.join('./', sra_stat['sra_id']+'/', sra_stat['sra_id']+'.sra'), sra_path],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        shutil.rmtree(os.path.join('./', sra_stat['sra_id']+'/'))
    # Move files downloaded by ascp
    if os.path.exists(os.path.join(individual_sra_tmp_dir, sra_stat['sra_id']+'.sra')):
        subprocess.run(['mv', os.path.join(individual_sra_tmp_dir, sra_stat['sra_id']+'.sra'), sra_path],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        shutil.rmtree(individual_sra_tmp_dir)
    assert os.path.exists(sra_path), 'SRA file download failed: '+sra_stat['sra_id']

def check_getfastq_dependency(args):
    if args.pfd=='yes':
        test_pfd = subprocess.run([args.pfd_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_pfd.returncode==0), "parallel-fastq-dump PATH cannot be found: "+args.pfd_exe
        test_prefetch = subprocess.run([args.prefetch_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_prefetch.returncode==0), "prefetch (SRA toolkit) PATH cannot be found: "+args.prefetch_exe
    if args.ascp=='yes':
        test_ascp = subprocess.run([args.ascp_exe, '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_ascp.returncode==0), "ascp (Aspera Connect) PATH cannot be found: "+args.ascp_exe
    if args.fastp=='yes':
        test_fp = subprocess.run([args.fastp_exe, '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_fp.returncode==0), "fastp PATH cannot be found: "+args.fastp_exe
    test_pigz = subprocess.run(['pigz', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if test_pigz.returncode==0:
        print('pigz found. It will be used for compression/decompression in read name formatting.')
        gz_exe = 'pigz -p '+str(args.threads)
        ungz_exe = 'unpigz -p'+str(args.threads)
    else:
        print('pigz not found. gzip/gunzip will be used for compression/decompression in read name formatting.')
        gz_exe = 'gzip'
        ungz_exe = 'gunzip'
    return gz_exe,ungz_exe

def set_getfastq_directories(args, sra_id):

    if args.work_dir.startswith('./'):
        args.work_dir = args.work_dir.replace('.', os.getcwd())
    if args.id is not None:
        output_dir = os.path.join(args.work_dir, 'getfastq', args.id)
    elif args.metadata is not None:
        output_dir = os.path.join(args.work_dir, 'getfastq', sra_id)
    elif args.id_list is not None:
        output_dir = os.path.join(args.work_dir, 'getfastq', sra_id)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir

def run_pfd(sra_stat, args, output_dir, seq_summary, start, end):
    sra_path = os.path.join(output_dir, sra_stat['sra_id']+'.sra')
    pfd_command = ['parallel-fastq-dump', '-t', str(args.threads), '--minReadLen', str(args.min_read_length), '--qual-filter-1',
                   '--skip-technical', '--split-3', '--clip', '--gzip', '--outdir', output_dir,
                   '--tmpdir', output_dir]
    print('Total sampled bases:', "{:,}".format(sra_stat['spot_length']*(end-start+1)), 'bp')
    pfd_command = pfd_command + ['--minSpotId', str(start), '--maxSpotId', str(end)]
    # If sra_stat['sra_id'], not sra_path, is provided, pfd couldn't find pre-downloaded .sra files
    # and start downloading it to $HOME/ncbi/public/sra/
    pfd_command = pfd_command + ['-s', sra_path]
    #pfd_command = pfd_command + ['-s', sra_stat['sra_id']]
    print('Command:', ' '.join(pfd_command))
    pfd_out = subprocess.run(pfd_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if args.pfd_print=='yes':
        print('parallel-fastq-dump stdout:')
        print(pfd_out.stdout.decode('utf8'))
        print('parallel-fastq-dump stderr:')
        print(pfd_out.stderr.decode('utf8'))
    stdout = pfd_out.stdout.decode('utf8')
    nd = [ int(line.replace('Read ','').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Read') ]
    nr = [ int(line.replace('Rejected ','').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Rejected') ]
    nw = [ int(line.replace('Written ','').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Written') ]
    seq_summary['bp_dumped'].loc[sra_stat['sra_id']] += sum(nd) * sra_stat['spot_length']
    seq_summary['bp_rejected'].loc[sra_stat['sra_id']] += sum(nr) * sra_stat['spot_length']
    seq_summary['bp_written'].loc[sra_stat['sra_id']] += sum(nw) * sra_stat['spot_length']
    if (sra_stat['layout']=='paired'):
        unpaired_file = os.path.join(output_dir, sra_stat['sra_id']+'.fastq.gz')
        if os.path.exists(unpaired_file):
            print('layout =', sra_stat['layout'], '; Deleting unpaired file:', unpaired_file)
            os.remove(unpaired_file)
        else:
            print('Unpaired file not found:', unpaired_file)
    return seq_summary

def run_fastp(sra_stat, args, output_dir, seq_summary):
    print('Running fastp.')
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
    outext = '.fastp.fastq.gz'
    if args.threads>16:
        print('Too many threads for fastp (--threads {}). Only 16 threads will be used.'.format(args.threads))
        fastp_thread = 16
    else:
        fastp_thread = args.threads
    fp_command = ['fastp', '--thread', str(fastp_thread), '--length_required', str(args.min_read_length)] + args.fastp_option.split(' ')
    if sra_stat['layout']=='single':
        infile = os.path.join(output_dir,sra_stat['sra_id'])
        fp_command = fp_command + ['--in1',infile+inext,'--out1',infile+outext]
    elif sra_stat['layout']=='paired':
        infile1 = os.path.join(output_dir,sra_stat['sra_id']+'_1')
        infile2 = os.path.join(output_dir,sra_stat['sra_id']+'_2')
        fp_command = fp_command + ['--in1',infile1+inext,'--out1',infile1+outext,'--in2',infile2+inext,'--out2',infile2+outext]
    fp_command = [ fc for fc in fp_command if fc!='' ]
    print('Command:', ' '.join(fp_command))
    fp_out = subprocess.run(fp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if args.fastp_print=='yes':
        print('fastp stdout:')
        print(fp_out.stdout.decode('utf8'))
        print('fastp stderr:')
        print(fp_out.stderr.decode('utf8'))
    if args.remove_tmp=='yes':
        remove_intermediate_files(sra_stat, ext=inext, work_dir=output_dir)
    bps = fp_out.stderr.decode('utf8').split('\n')
    bp_in = list()
    bp_out = list()
    for i in range(len(bps)):
        if (' before filtering:' in bps[i]):
            bp_in.append(int(bps[i+2].replace('total bases: ', '')))
        if (' after filtering:' in bps[i])|(' aftering filtering:' in bps[i]):
            bp_out.append(int(bps[i+2].replace('total bases: ', '')))
    seq_summary['bp_fastp_in'].loc[sra_stat['sra_id']] += sum(bp_in)
    seq_summary['bp_fastp_out'].loc[sra_stat['sra_id']] += sum(bp_out)
    return seq_summary

def rename_reads(sra_stat, args, output_dir, gz_exe, ungz_exe):
    inext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir)
    outext = '.rename.fastq.gz'
    if sra_stat['layout']=='single':
        inbase = os.path.join(output_dir,sra_stat['sra_id'])
        if os.path.exists(inbase+inext):
            infile = inbase+inext
            os.system(ungz_exe+' -c "'+infile+'" | sed -e "s|[[:space:]].*|/1|" | '+gz_exe+' -c > "'+inbase+outext+'"')
    elif sra_stat['layout']=='paired':
        inbase1 = os.path.join(output_dir,sra_stat['sra_id']+'_1')
        inbase2 = os.path.join(output_dir,sra_stat['sra_id']+'_2')
        if os.path.exists(inbase1+inext):
            infile1 = inbase1+inext
            infile2 = inbase2+inext
            os.system(ungz_exe+' -c "'+infile1+'" | sed -e "s|[[:space:]].*|/1|" | '+gz_exe+' -c > "'+inbase1+outext+'"')
            os.system(ungz_exe+' -c "'+infile2+'" | sed -e "s|[[:space:]].*|/2|" | '+gz_exe+' -c > "'+inbase2+outext+'"')
    if args.remove_tmp=='yes':
        remove_intermediate_files(sra_stat, ext=inext, work_dir=output_dir)

def rename_fastq(sra_stat, output_dir, inext, outext):
    if sra_stat['layout']=='single':
        inbase = os.path.join(output_dir, sra_stat['sra_id'])
        os.rename(inbase+inext, inbase+outext)
    elif sra_stat['layout']=='paired':
        inbase1 = os.path.join(output_dir, sra_stat['sra_id']+'_1')
        inbase2 = os.path.join(output_dir, sra_stat['sra_id']+'_2')
        os.rename(inbase1+inext, inbase1+outext)
        os.rename(inbase2+inext, inbase2+outext)

def get_sra_stat(sra_id, metadata, num_bp_per_sra):
    sra_stat = dict()
    sra_stat['sra_id'] = sra_id
    is_sra = (metadata.df.loc[:,'run']==sra_id)
    assert is_sra.sum()==1, 'There are multiple metadata rows with the same SRA ID: '+sra_id
    sra_stat['layout'] = metadata.df.loc[is_sra,'lib_layout'].values[0]
    sra_stat['total_spot'] = int(metadata.df.loc[is_sra,'total_spots'].values[0])
    try:
        sra_stat['spot_length'] = int(metadata.df.loc[is_sra,'spot_length'].values[0])
    except ValueError as e:
        sra_stat['spot_length'] = int(int(metadata.df.loc[is_sra,'total_bases'].values[0])/int(sra_stat['total_spot']))
        print('spot_length cannot be obtained directly from the metadata.')
        print('Using total_bases/total_spots ( =', str(sra_stat['spot_length']), ') instead.')
    sra_stat['num_read_per_sra'] = int(num_bp_per_sra/sra_stat['spot_length'])
    return sra_stat

def sequence_extraction(args, sra_stat, start, end, output_dir, seq_summary, gz_exe, ungz_exe, total_sra_bp):
    sra_id = sra_stat['sra_id']
    if args.pfd=='yes':
        seq_summary = run_pfd(sra_stat, args, output_dir, seq_summary, start, end)
        seq_summary['bp_remaining'].loc[sra_id] = seq_summary['bp_dumped'].loc[sra_id] - seq_summary['bp_written'].loc[sra_id]
    if args.fastp=='yes':
        seq_summary = run_fastp(sra_stat, args, output_dir, seq_summary)
        seq_summary['bp_remaining'].loc[sra_id] = seq_summary['bp_dumped'].loc[sra_id] - seq_summary['bp_fastp_out'].loc[sra_id]
    if args.read_name=='trinity':
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

def calc_2nd_ranges(args, total_bp_remaining, seq_summary):
    sra_target_bp = total_bp_remaining / len(seq_summary['bp_remaining'])
    sra_target_reads = (sra_target_bp / seq_summary['spot_length']).astype(int)
    start_2nds = seq_summary['end_1st'] + 1
    end_2nds = start_2nds + (sra_target_reads / seq_summary['rate_passed']).astype(int)
    total_spots = seq_summary['total_spot']
    spot_lengths = seq_summary['spot_length']
    pooled_missing_bp = total_bp_remaining
    for dummy in range(1000):
        current_total_bp = 0
        for sra_id in end_2nds.index:
            pooled_missing_read = (pooled_missing_bp/spot_lengths.loc[sra_id]).astype(int)
            if (end_2nds.loc[sra_id] + pooled_missing_read < total_spots.loc[sra_id]):
                pooled_missing_bp = 0
                end_2nds.loc[sra_id] = end_2nds.loc[sra_id] + pooled_missing_bp
            elif (end_2nds.loc[sra_id] + pooled_missing_read > total_spots.loc[sra_id]):
                pooled_missing_bp = (end_2nds.loc[sra_id] + pooled_missing_read - total_spots.loc[sra_id]) * spot_lengths.loc[sra_id]
                end_2nds.loc[sra_id] = total_spots.loc[sra_id]
            current_total_bp += end_2nds.loc[sra_id] * spot_lengths.loc[sra_id]
        all_equal_total_spots = all([ e2==ts for e2,ts in zip(end_2nds,total_spots) ])
        is_enough_read = (current_total_bp>=total_bp_remaining)
        if all_equal_total_spots:
            print('Reached total spots in all SRAs.')
            break
        if is_enough_read:
            print('Enough read numbers were assigned for the 2nd round sequence extraction.')
            break
    seq_summary['start_2nd'] = start_2nds
    seq_summary['end_2nd'] = end_2nds
    return seq_summary

def print_read_stats(args, seq_summary, max_bp, individual=True):
    print('Target size (--max_bp):', "{:,}".format(max_bp), 'bp')
    if args.pfd=='yes':
        print('Sum of fastq_dump dumped reads:', "{:,}".format(seq_summary['bp_dumped'].sum()), 'bp')
        print('Sum of fastq_dump rejected reads:', "{:,}".format(seq_summary['bp_rejected'].sum()), 'bp')
        print('Sum of fastq_dump written reads:', "{:,}".format(seq_summary['bp_written'].sum()), 'bp')
    if args.fastp=='yes':
        print('Sum of fastp input reads:', "{:,}".format(seq_summary['bp_fastp_in'].sum()), 'bp')
        print('Sum of fastp output reads:', "{:,}".format(seq_summary['bp_fastp_out'].sum()), 'bp')
    if individual:
        print('Individual SRA IDs:', ' '.join(seq_summary['bp_dumped'].index.values))
        read_types = list()
        keys = list()
        if args.pfd=='yes':
            read_types = read_types + ['fastq_dump dumped reads', 'fastq_dump rejected reads', 'fastq_dump written reads']
            keys = keys + ['bp_dumped', 'bp_rejected', 'bp_written']
        if args.fastp=='yes':
            read_types = read_types + ['fastp input reads', 'fastp output reads']
            keys = keys + ['bp_fastp_in', 'bp_fastp_out']
        if len(read_types)>0:
            for rt,key in zip(read_types, keys):
                values = [ '{:,}'.format(s) for s in seq_summary[key].values ]
                txt = ' '.join(values)
                print('Individual {} (bp): {}'.format(rt, txt))
    print('')

def getfastq_metadata(args):
    assert (args.id is None and args.id_list is None)!=(args.metadata is None), 'Either --id, --id_list or --metadata should be specified.'
    if args.id is not None:
        print('--id is specified. Downloading SRA metadata from Entrez.')
        assert (args.entrez_email!='aaa@bbb.com'), "Provide your email address. No worry, you won't get spam emails."
        Entrez.email = args.entrez_email
        sra_id = args.id
        output_dir = set_getfastq_directories(args, sra_id)
        search_term = getfastq_search_term(sra_id, args.entrez_additional_search_term)
        print('Entrez search term:', search_term)
        xml_root = getfastq_getxml(search_term)
        metadata = Metadata.from_xml(xml_root)
       # if args.save_metadata:
            #metadata.df.to_csv(os.path.join(output_dir,'metadata_all.tsv'), sep='\t', index=False)
        print('Filtering SRA entry with --layout:', args.layout)
        layout = get_layout(args, metadata)
        metadata.df = metadata.df.loc[(metadata.df['lib_layout']==layout),:]
        if args.sci_name is not None:
            print('Filtering SRA entry with --sci_name:', args.sci_name)
            metadata.df = metadata.df.loc[(metadata.df['scientific_name']==args.sci_name),:]
        if args.save_metadata:
            metadata.df.to_csv(os.path.join(output_dir,'metadata_',sra_id,'.tsv'), sep='\t', index=False)
    if args.id_list is not None:
        print('--id is specified. Downloading SRA metadata from Entrez.')
        assert (args.entrez_email!='aaa@bbb.com'), "Provide your email address. No worry, you won't get spam emails."
        Entrez.email = args.entrez_email
        sra_id_list = [line.rstrip('\n') for line in open(args.id_list)]
        for sra_id in sra_id_list:
            output_dir = set_getfastq_directories(args, sra_id)
            search_term = getfastq_search_term(sra_id, args.entrez_additional_search_term)
            print('Entrez search term:', search_term)
            xml_root = getfastq_getxml(search_term)
            metadata = Metadata.from_xml(xml_root)
            #if args.save_metadata:
             #   metadata.df.to_csv(os.path.join(output_dir,'metadata_all.tsv'), sep='\t', index=False)
            print('Filtering SRA entry with --layout:', args.layout)
            layout = get_layout(args, metadata)
            metadata.df = metadata.df.loc[(metadata.df['lib_layout']==layout),:]
            if args.sci_name is not None:
                print('Filtering SRA entry with --sci_name:', args.sci_name)
                metadata.df = metadata.df.loc[(metadata.df['scientific_name']==args.sci_name),:]
            if args.save_metadata:
                metadata.df.to_csv(os.path.join(output_dir,'metadata_',sra_id,'.tsv'), sep='\t', index=False)

    if args.metadata is not None:
        print('--metadata is specified. Reading existing metadata table.')
        assert args.concat=='no', '--concat should be set "no" when --metadata is specified.'
        metadata = load_metadata(args)
    return metadata

def is_getfastq_output_present(args, sra_stat, output_dir):
    if args.concat=='yes':
        prefixes = [args.id,]
    else:
        prefixes = [sra_stat['sra_id'],]
    if sra_stat['layout']=='single':
        sub_exts = ['',]
    elif sra_stat['layout']=='paired':
        sub_exts = ['_1', '_2']
    exts = ['.amalgkit.fastq.gz',]
    is_output_present = True
    for prefix,sub_ext,ext in itertools.product(prefixes, sub_exts, exts):
        out_path1 = os.path.join(output_dir, prefix+sub_ext+ext)
        out_path2 = os.path.join(output_dir, prefix+sub_ext+ext+'.safely_removed')
        print(out_path1)
        is_output_present *= (os.path.exists(out_path1)|os.path.exists(out_path2))
    return is_output_present

def getfastq_main(args):
    #sra_dir = os.path.join(os.path.expanduser("~"), 'ncbi/public/sra')
    gz_exe,ungz_exe = check_getfastq_dependency(args)
    metadata = getfastq_metadata(args)
    assert metadata.df.shape[0] > 0, 'No SRA entry found. Make sure if --id is compatible with --sci_name and --layout.'
    print('SRA IDs:', ' '.join(metadata.df['run'].tolist()))
    max_bp = int(args.max_bp.replace(',',''))
    num_sra = metadata.df.shape[0]
    num_bp_per_sra = int(max_bp/num_sra)
    total_sra_bp = int(metadata.df['total_bases'].astype(int).sum())
    offset = 10000 # https://edwards.sdsu.edu/research/fastq-dump/
    print('Number of SRAs:', num_sra)
    print('Total target size (--max_bp):', "{:,}".format(max_bp), 'bp')
    print('Total SRA size:', "{:,}".format(total_sra_bp), 'bp')
    print('Target size per SRA:', "{:,}".format(num_bp_per_sra), 'bp')
    for i in metadata.df.index:
        print('Individual SRA size :', metadata.df.loc[i,'run'], ':', "{:,}".format(int(metadata.df.loc[i,'total_bases'])), 'bp')
    sra_ids = metadata.df.loc[:,'run'].values
    seq_summary = dict()
    keys = ['bp_dumped','bp_rejected','bp_written','bp_fastp_in','bp_fastp_out','bp_remaining','bp_still_available',
            'spot_length','total_spot','rate_passed','start_1st','end_1st','start_2nd','end_2nd',]
    for key in keys:
        seq_summary[key] = pandas.Series(0, index=sra_ids)
    for i in metadata.df.index:
        print('')
        start_time = time.time()
        sra_id = metadata.df.loc[i,'run']
        sra_stat = get_sra_stat(sra_id, metadata, num_bp_per_sra)
        output_dir = set_getfastq_directories(args, sra_id)
        if (is_getfastq_output_present(args, sra_stat, output_dir))&(args.redo=='no'):
            if metadata.df.shape[0]==1:
                print('Output file(s) detected. Exiting.  Set "--redo yes" for reanalysis.')
                sys.exit()
            else:
                print('Output file(s) detected. Skipping {}.  Set "--redo yes" for reanalysis.'.format(sra_id))
                continue
        remove_old_intermediate_files(metadata, work_dir=output_dir)
        print('SRA ID:', sra_stat['sra_id'])
        print('Library layout:', sra_stat['layout'])
        print('Number of reads:', "{:,}".format(sra_stat['total_spot']))
        print('Single/Paired read length:', sra_stat['spot_length'], 'bp')
        print('Total bases:', "{:,}".format(int(metadata.df.loc[i,'total_bases'])), 'bp')
        download_sra(sra_stat, args, output_dir, overwrite=False)
        start,end = get_range(sra_stat, offset, total_sra_bp, max_bp)
        seq_summary = sequence_extraction(args, sra_stat, start, end, output_dir, seq_summary,
                                          gz_exe, ungz_exe, total_sra_bp)
        print('Time elapsed for 1st-round sequence extraction: {}, {:,} sec'.format(sra_stat['sra_id'], int(time.time()-start_time)))

    print('\n--- getfastq 1st-round sequence generation report ---')
    print_read_stats(args, seq_summary, max_bp)
    total_bp_remaining = seq_summary['bp_remaining'].sum()
    percent_remaining = total_bp_remaining/max_bp*100
    if args.pfd=='yes':
        total_bp_dumped = seq_summary['bp_dumped'].sum()
        total_bp_out = seq_summary['bp_written'].sum()
    if args.fastp=='yes':
        total_bp_out = seq_summary['bp_fastp_out'].sum()
    percent_dropped = total_bp_out/total_bp_dumped*100
    txt = '{:.2f}% of reads were dropped in the 1st-round sequence generation.'
    print(txt.format(percent_remaining, args.tol))
    txt = 'The amount of generated reads were {:.2f}% ({:,}/{:,}) smaller than the target size (tol={}%).'
    print(txt.format(percent_remaining, total_bp_out, max_bp, args.tol))
    if (percent_remaining<args.tol):
        print('Proceeding without 2nd-round sequence extraction.')
    else:
        print('Starting the 2nd-round sequence extraction to compensate it.')
        seq_summary = calc_2nd_ranges(args, total_bp_remaining, seq_summary)
        ext_main = '.amalgkit.fastq.gz'
        ext_1st_tmp = '.amalgkit_1st.fastq.gz'
        for i in metadata.df.index:
            print('')
            start_time = time.time()
            sra_id = metadata.df.loc[i,'run']
            sra_stat = get_sra_stat(sra_id, metadata, num_bp_per_sra)
            start = seq_summary['start_2nd'].loc[sra_id]
            end = seq_summary['end_2nd'].loc[sra_id]
            if (start>=end):
                txt = '{}: All spots were used in the 1st trial. Cancelling the 2nd trial. start={:,}, end={:,}'
                print(txt.format(sra_id, start, end))
                continue
            rename_fastq(sra_stat, output_dir, inext=ext_main, outext=ext_1st_tmp)
            seq_summary = sequence_extraction(args, sra_stat, start, end, output_dir, output_dir, seq_summary,
                                              gz_exe, ungz_exe, total_sra_bp)
            if (layout=='single'):
                subexts = ['']
            elif (layout=='paired'):
                subexts = ['_1','_2']
            for subext in subexts:
                added_path = os.path.join(output_dir, sra_id+subext+ext_1st_tmp)
                adding_path = os.path.join(output_dir, sra_id+subext+ext_main)
                assert os.path.exists(added_path), 'Dumped fastq not found: '+added_path
                assert os.path.exists(adding_path), 'Dumped fastq not found: '+adding_path
                os.system('cat "'+adding_path+'" >> "'+added_path+'"')
                os.remove(adding_path)
                os.rename(added_path, adding_path)
            print('Time elapsed for 2nd-round sequence extraction: {}, {:,} sec'.format(sra_stat['sra_id'], int(time.time()-start_time)))

    print('')
    if args.concat=='yes':
        concat_fastq(args, metadata, output_dir, num_bp_per_sra)
    if args.remove_sra=='yes':
        remove_sra_files(metadata, sra_dir=output_dir)
    else:
        if args.pfd=='yes':
            print('SRA files not removed:', output_dir)
    print('\n--- getfastq final report ---')
    print_read_stats(args, seq_summary, max_bp)

