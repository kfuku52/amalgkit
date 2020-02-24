from Bio import Entrez
from amalgkit.metadata import Metadata
from urllib.error import HTTPError
import numpy, pandas
import time, datetime, lxml, subprocess, os, shutil, gzip, glob

def getfastq_search_term(ncbi_id, additional_search_term):
    # https://www.ncbi.nlm.nih.gov/books/NBK49540/
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

def get_range(max_bp, total_sra_bp, total_spot, num_read_per_sra, offset):
    if (total_sra_bp<=max_bp):
        start = 1
        end = total_spot
    else:
        if (total_spot > (num_read_per_sra + offset)):
            start = offset
            end = offset + num_read_per_sra
        elif (total_spot > num_read_per_sra):
            start = total_spot - num_read_per_sra
            end = total_spot
        elif (total_spot <= num_read_per_sra):
            start = 1
            end = total_spot
    return start,end

def concat_fastq(args, metadata, sra_output_dir):
    import re
    layout = get_layout(args, metadata)
    inext = '.amalgkit.fastq.gz'
    infiles = [ f for f in os.listdir(sra_output_dir) if f.endswith(inext) ]
    num_inext_files = len(infiles)
    if (layout=='single')&(num_inext_files==1):
        print('Only 1', inext, 'file was detected. No concatenation will happen. Just replacing ID in the output file name.')
        outfile = args.id+inext
        os.rename(infiles[0], outfile)
        return None
    elif (layout=='paired')&(num_inext_files==2):
        print('Only 1 pair of', inext, 'files were detected. No concatenation will happen. Just replacing IDs in the output file names.')
        for infile in infiles:
            outfile = args.id+re.sub('.*(_[1-2])', '\g<1>', infile)
            os.rename(infile, outfile)
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
            outfile_path = os.path.join(sra_output_dir, args.id+subext+outext)
            if os.path.exists(outfile_path):
                os.remove(outfile_path)
            #with gzip.open(outfile_path, 'wb') as outfile:
            #    for each_infile in infiles:
            #        with gzip.open(os.path.join(args.work_dir, each_infile), 'rb') as infile:
            #            shutil.copyfileobj(infile, outfile) # unacceptably slow
            if os.path.exists(outfile_path):
                os.remove(outfile_path)
            for infile in infiles:
                infile_path = os.path.join(sra_output_dir, infile)
                assert os.path.exists(infile_path), 'Dumped fastq not found: '+infile_path
                os.system('cat "'+infile_path+'" >> "'+outfile_path+'"')
            print('')
        if args.remove_tmp=='yes':
            for i in metadata.df.index:
                sra_id = metadata.df.loc[i,'run']
                ext = get_newest_intermediate_file_extension(sra_id=sra_id, work_dir=sra_output_dir)
                remove_intermediate_files(sra_id=sra_id, layout=layout, ext=ext, work_dir=sra_output_dir)

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

def remove_intermediate_files(sra_id, layout, ext, work_dir):
    file_paths = list()
    if layout=='single':
        file_paths.append(os.path.join(work_dir,sra_id+ext))
    elif layout=='paired':
        for i in [1,2]:
            file_paths.append(os.path.join(work_dir,sra_id+'_'+str(i)+ext))
    for file_path in file_paths:
        if os.path.exists(file_path):
            print('Deleting intermediate file:', file_path)
            os.remove(file_path)
        else:
            print('Tried to delete but file not found:', file_path)

def get_newest_intermediate_file_extension(sra_id, work_dir):
    # Order is important in this list. More downstream should come first.
    extensions = ['.amalgkit.fastq.gz','.rename.fastq.gz','.fastp.fastq.gz','.fastq.gz']
    files = os.listdir(work_dir)
    for ext in extensions:
        if any([ (f.startswith(sra_id))&(f.endswith(ext)) for f in files ]):
            ext_out = ext
            break
    assert 'ext_out' in locals(), 'Any of expected extensions ('+' '.join(extensions)+') found in '+work_dir
    return ext_out

def download_sra(sra_id, args, work_dir, overwrite=False):
    sra_path = os.path.join(work_dir, sra_id+'.sra')
    individual_sra_tmp_dir = os.path.join(work_dir, sra_id+'/')
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
        sra_site = 'anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/'+sra_id[0:3]+'/'+sra_id[0:6]+'/'+sra_id+'/'+sra_id+'.sra'
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
                            '--output-directory', './', sra_id]
        print('Command:', ' '.join(prefetch_command))
        prefetch_out = subprocess.run(prefetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('prefetch stdout:')
        print(prefetch_out.stdout.decode('utf8'))
        print('prefetch stderr:')
        print(prefetch_out.stderr.decode('utf8'))
    if (not os.path.exists(sra_path))&(not os.path.exists(os.path.join(work_dir, sra_id+'/', sra_id+'.sra'))):
        print('Trying to download the SRA file using prefetch (http protocol).')
        if os.path.exists(individual_sra_tmp_dir):
            shutil.rmtree(individual_sra_tmp_dir)
        prefetch_command = [args.prefetch_exe, '--force', 'no', '--transport', 'http', '--max-size', '100G',
                            '--output-directory', './', sra_id]
        print('Command:', ' '.join(prefetch_command))
        prefetch_out = subprocess.run(prefetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print('prefetch stdout:')
        print(prefetch_out.stdout.decode('utf8'))
        print('prefetch stderr:')
        print(prefetch_out.stderr.decode('utf8'))
    # Move files downloaded by prefetch. This is necessary because absolute path didn't work for prefetch --output-directory
    if os.path.exists(os.path.join('./', sra_id+'/', sra_id+'.sra')):
        subprocess.run(['mv', os.path.join('./', sra_id+'/', sra_id+'.sra'), sra_path],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        shutil.rmtree(os.path.join('./', sra_id+'/'))
    # Move files downloaded by ascp
    if os.path.exists(os.path.join(individual_sra_tmp_dir, sra_id+'.sra')):
        subprocess.run(['mv', os.path.join(individual_sra_tmp_dir, sra_id+'.sra'), sra_path],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        shutil.rmtree(individual_sra_tmp_dir)
    assert os.path.exists(sra_path), 'SRA file download failed: '+sra_id

def getfastq_main(args):
    #sra_dir = os.path.join(os.path.expanduser("~"), 'ncbi/public/sra')
    assert (args.entrez_email!='aaa@bbb.com'), "Provide your email address. No worry, you won't get spam emails."
    if args.pfd=='yes':
        test_pfd = subprocess.run([args.pfd_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_pfd.returncode==0), "parallel-fastq-dump PATH cannot be found."
        test_prefetch = subprocess.run([args.prefetch_exe, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_prefetch.returncode==0), "prefetch (SRA toolkit) PATH cannot be found."
    if args.ascp=='yes':
        test_ascp = subprocess.run([args.ascp_exe, '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_ascp.returncode==0), "ascp (Aspera Connect) PATH cannot be found."
    if args.fastp=='yes':
        test_fp = subprocess.run([args.fastp_exe, '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_fp.returncode==0), "fastp PATH cannot be found."
    if args.read_name!='default':
        test_pigz = subprocess.run(['pigz', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if test_pigz.returncode==0:
            print('pigz found. It will be used for compression/decompression in read name formatting.')
            gz_exe = 'pigz -p '+str(args.threads)
            ungz_exe = 'unpigz -p'+str(args.threads)
        else:
            print('pigz not found. gzip/gunzip will be used for compression/decompression in read name formatting.')
            gz_exe = 'gzip'
            ungz_exe = 'gunzip'
    if args.work_dir.startswith('./'):
        args.work_dir = args.work_dir.replace('.', os.getcwd())
    if not os.path.exists(args.work_dir):
        os.makedirs(args.work_dir)
    if args.auto_dir == 'yes':
        output_dir = os.path.join(args.work_dir, 'getfastq_'+args.id)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        output_dir = args.work_dir

    sra_output_dir = os.path.join(output_dir, 'sra_files')
    if not os.path.exists(sra_output_dir):
        os.mkdir(sra_output_dir)

    sra_temp_dir = os.path.join(output_dir, 'temp')
    if not os.path.exists(sra_temp_dir):
        os.mkdir(sra_temp_dir)

    Entrez.email = args.entrez_email
    search_term = getfastq_search_term(args.id, args.entrez_additional_search_term)
    print('Entrez search term:', search_term)
    xml_root = getfastq_getxml(search_term)
    metadata = Metadata.from_xml(xml_root)
    if args.save_metadata:
        metadata.df.to_csv(os.path.join(output_dir,'metadata_all.tsv'), sep='\t', index=False)
    print('Filtering SRA entry with --layout:', args.layout)
    layout = get_layout(args, metadata)
    metadata.df = metadata.df.loc[(metadata.df['lib_layout']==layout),:]
    if args.sci_name is not None:
        print('Filtering SRA entry with --sci_name:', args.sci_name)
        metadata.df = metadata.df.loc[(metadata.df['scientific_name']==args.sci_name),:]
    remove_old_intermediate_files(metadata, work_dir=sra_output_dir)
    if args.save_metadata:
        metadata.df.to_csv(os.path.join(output_dir,'metadata_target.tsv'), sep='\t', index=False)
    assert metadata.df.shape[0] > 0, 'No SRA entry found. Make sure if --id is compatible with --sci_name and --layout.'
    print('SRA IDs:', ' '.join(metadata.df['run'].tolist()))
    max_bp = int(args.max_bp.replace(',',''))
    num_sra = metadata.df.shape[0]
    num_bp_per_sra = int(max_bp/num_sra)
    total_sra_bp = int(metadata.df['total_bases'].astype(int).sum())
    offset = 10000 # https://edwards.sdsu.edu/research/fastq-dump/
    print('Number of SRA:', num_sra)
    print('max_bp:', "{:,}".format(max_bp), 'bp')
    print('Total SRA size:', total_sra_bp, 'bp')
    print('Max size per SRA:', "{:,}".format(num_bp_per_sra), 'bp')
    for i in metadata.df.index:
        print('Individual SRA size :', metadata.df.loc[i,'run'], ':', "{:,}".format(int(metadata.df.loc[i,'total_bases'])), 'bp')
    bp_dumped = 0
    bp_rejected = 0
    bp_written = 0
    bp_fastp_in = 0
    bp_fastp_out = 0
    for i in metadata.df.index:
        print('')
        start_time = time.time()
        sra_id = metadata.df.loc[i,'run']
        total_spot = int(metadata.df.loc[i,'total_spots'])
        try:
            spot_length = int(metadata.df.loc[i,'spot_length'])
        except ValueError as e:
            print('INFO: spot_length cannot be obtained directly from the metadata. Using total_bases/total_spots instead.')
            spot_length = int(int(metadata.df.loc[i,'total_bases'])/int(metadata.df.loc[i,'total_spots']))
        num_read_per_sra = int(num_bp_per_sra/spot_length)
        print('SRA ID:', sra_id)
        print('Library layout:', layout)
        print('Number of reads:', "{:,}".format(total_spot))
        print('Single/Paired read length:', spot_length, 'bp')
        print('Total bases:', "{:,}".format(int(metadata.df.loc[i,'total_bases'])), 'bp')
        if args.pfd=='yes':
            sra_path = os.path.join(sra_output_dir, sra_id+'.sra')
            download_sra(sra_id, args, sra_output_dir, overwrite=False)
            pfd_command = ['parallel-fastq-dump', '-t', str(args.threads), '--minReadLen', str(args.min_read_length), '--qual-filter-1',
                           '--skip-technical', '--split-3', '--clip', '--gzip', '--outdir', sra_output_dir,
                           '--tmpdir', sra_temp_dir]
            start,end = get_range(max_bp, total_sra_bp, total_spot, num_read_per_sra, offset)
            print('Total sampled bases:', "{:,}".format(spot_length*(end-start+1)), 'bp')
            pfd_command = pfd_command + ['--minSpotId', str(start), '--maxSpotId', str(end)]
            # If sra_id, not sra_path, is provided, pfd couldn't find pre-downloaded .sra files
            # and start downloading it to $HOME/ncbi/public/sra/
            pfd_command = pfd_command + ['-s', sra_path]
            #pfd_command = pfd_command + ['-s', sra_id]
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
            bp_dumped += sum(nd) * spot_length
            bp_rejected += sum(nr) * spot_length
            bp_written += sum(nw) * spot_length
            if (layout=='paired'):
                unpaired_file = os.path.join(sra_output_dir, sra_id+'.fastq.gz')
                if os.path.exists(unpaired_file):
                    print('Layout=',layout,'; Deleting unpaired file:', unpaired_file)
                    os.remove(unpaired_file)
                else:
                    print('Unpaired file not found:', unpaired_file)
        if args.fastp=='yes':
            print('Running fastp.')
            inext = get_newest_intermediate_file_extension(sra_id=sra_id, work_dir=sra_output_dir)
            outext = '.fastp.fastq.gz'
            if args.threads>16:
                print('Too many threads for fastp (--threads {}). Use 16 threads only.'.format(args.threads))
                fastp_thread = 16
            else:
                fastp_thread = args.threads
            fp_command = ['fastp', '--thread', str(fastp_thread), '--length_required', str(args.min_read_length)] + args.fastp_option.split(' ')
            if layout=='single':
                infile = os.path.join(sra_output_dir,sra_id)
                fp_command = fp_command + ['--in1',infile+inext,'--out1',infile+outext]
            elif layout=='paired':
                infile1 = os.path.join(sra_output_dir,sra_id+'_1')
                infile2 = os.path.join(sra_output_dir,sra_id+'_2')
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
                remove_intermediate_files(sra_id=sra_id, layout=layout, ext=inext, work_dir=sra_output_dir)
            bps = fp_out.stderr.decode('utf8').split('\n')
            bp_in = [ int(line.replace('total bases: ','').split(' ')[0]) for line in bps if line.startswith('total bases') ][0::2]
            bp_out = [ int(line.replace('total bases: ','').split(' ')[0]) for line in bps if line.startswith('total bases') ][1::2]
            bp_fastp_in += sum(bp_in)
            bp_fastp_out += sum(bp_out)
        if args.read_name=='trinity':
            inext = get_newest_intermediate_file_extension(sra_id=sra_id, work_dir=sra_output_dir)
            outext = '.rename.fastq.gz'
            if layout=='single':
                inbase = os.path.join(sra_output_dir,sra_id)
                if os.path.exists(inbase+inext):
                    infile = inbase+inext
                os.system(ungz_exe+' -c "'+infile+'" | sed -e "s|[[:space:]].*|/1|" | '+gz_exe+' -c > "'+inbase+outext+'"')
            elif layout=='paired':
                inbase1 = os.path.join(sra_output_dir,sra_id+'_1')
                inbase2 = os.path.join(sra_output_dir,sra_id+'_2')
                if os.path.exists(inbase1+inext):
                    infile1 = inbase1+inext
                    infile2 = inbase2+inext
                os.system(ungz_exe+' -c "'+infile1+'" | sed -e "s|[[:space:]].*|/1|" | '+gz_exe+' -c > "'+inbase1+outext+'"')
                os.system(ungz_exe+' -c "'+infile2+'" | sed -e "s|[[:space:]].*|/2|" | '+gz_exe+' -c > "'+inbase2+outext+'"')
            if args.remove_tmp=='yes':
                remove_intermediate_files(sra_id=sra_id, layout=layout, ext=inext, work_dir=sra_output_dir)
        if args.pfd=='yes':
            inext = get_newest_intermediate_file_extension(sra_id=sra_id, work_dir=sra_output_dir)
            outext = '.amalgkit.fastq.gz'
            if layout=='single':
                inbase = os.path.join(sra_output_dir, sra_id)
                os.rename(inbase+inext, inbase+outext)
            elif layout=='paired':
                inbase1 = os.path.join(sra_output_dir, sra_id+'_1')
                inbase2 = os.path.join(sra_output_dir, sra_id+'_2')
                os.rename(inbase1+inext, inbase1+outext)
                os.rename(inbase2+inext, inbase2+outext)
        print('Time elapsed:', sra_id, int(time.time()-start_time), '[sec]')
    print('')
    if args.concat=='yes':
        print('max_bp:', "{:,}".format(max_bp), 'bp')
        print('Sum of dumped reads:', "{:,}".format(bp_dumped), 'bp')
        print('Sum of rejected reads:', "{:,}".format(bp_rejected), 'bp')
        print('Sum of written reads:', "{:,}".format(bp_written), 'bp')
        if args.fastp=='yes':
            print('Sum of fastp input reads:', "{:,}".format(bp_fastp_in), 'bp')
            print('Sum of fastp output reads:', "{:,}".format(bp_fastp_out), 'bp')
        print('')
        concat_fastq(args, metadata, sra_output_dir)
    if args.remove_sra=='yes':
        remove_sra_files(metadata, sra_dir=sra_output_dir)
    else:
        if args.pfd=='yes':
            print('SRA files not removed:', sra_output_dir)

