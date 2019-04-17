from Bio import Entrez
from amalgkit.metadata import Metadata
from urllib.error import HTTPError
import numpy, pandas
import time, datetime, lxml, subprocess, os, shutil, gzip, glob

def getfastq_search_term(ncbi_id):
    # https://www.ncbi.nlm.nih.gov/books/NBK49540/
    search_term = ncbi_id+' AND "platform illumina"[Properties] AND "type rnaseq"[Filter] AND "sra biosample"[Filter]'
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

def concat_fastq(args, metadata, clean=True):
    is_output_exist = True
    ext = '.fastp.fastq.gz' if args.fastp=='yes' else '.fastq.gz'
    if args.layout=='single':
        subexts = ['',]
    elif args.layout=='paired':
        subexts = ['_1','_2',]
    for subext in subexts:
        infiles = metadata.df['run'].replace('$',subext+ext, regex=True)
        outfile_path = os.path.join(args.work_dir,args.id+subext+ext)
        #with gzip.open(outfile_path, 'wb') as outfile:
        #    for each_infile in infiles:
        #        with gzip.open(os.path.join(args.work_dir, each_infile), 'rb') as infile:
        #            shutil.copyfileobj(infile, outfile) # unacceptably slow
        if os.path.exists(outfile_path):
            os.remove(outfile_path)
        for infile in infiles:
            infile_path = os.path.join(args.work_dir, infile)
            if os.path.exists(infile_path):
                os.system('cat '+infile_path+' >> '+outfile_path)
            else:
                print('Dumped fastq not found:', infile_path)
        is_current_output_exist = os.path.exists(outfile_path)
        is_output_exist = is_output_exist & is_current_output_exist
        if is_current_output_exist:
            print('Output written to:', outfile_path)
            if clean:
                for infile in infiles:
                    infile_path = os.path.join(args.work_dir, infile)
                    print('Deleting:', infile_path)
                    os.remove(infile_path)
        print('')
    if (clean)&(args.layout=='paired'):
        for sra_id in metadata.df['run']:
            unpaired_file = os.path.join(args.work_dir, sra_id+'.fastq.gz')
            if os.path.exists(unpaired_file):
                print('Deleting:', unpaired_file)
                os.remove(unpaired_file)
            else:
                print('Unpaired file not found:', unpaired_file)
        print('')

def remove_sra_files(metadata, sra_dir):
    for sra_id in metadata.df['run']:
        sra_pattern = os.path.join(sra_dir, sra_id+'.*')
        sra_paths = glob.glob(sra_pattern)
        if len(sra_paths)>0:
            for sra_path in sra_paths:
                print('Deleting:', sra_path)
                os.remove(sra_path)
        else:
            print('SRA file not found:', sra_pattern)
    print('')

def getfastq_main(args):
    sra_dir = os.path.join(os.path.expanduser("~"), 'ncbi/public/sra')
    assert (args.entrez_email!='aaa@bbb.com'), "Provide your email address. No worry, you won't get spam emails."
    test_pfd = subprocess.run(['parallel-fastq-dump', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert (test_pfd.returncode==0), "parallel-fastq-dump is not installed."
    if args.fastp=='yes':
        test_fp = subprocess.run(['fastp', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert (test_fp.returncode==0), "fastp is not installed."
    if not os.path.exists(args.work_dir):
        os.makedirs(args.work_dir)
    Entrez.email = args.entrez_email
    search_term = getfastq_search_term(args.id)
    print('Entrez search term:', search_term)
    xml_root = getfastq_getxml(search_term)
    metadata = Metadata.from_xml(xml_root)
    if args.save_metadata:
        metadata.df.to_csv(os.path.join(args.work_dir,'getfastq_sra_metadata.tsv'), sep='\t', index=False)
    metadata.df = metadata.df.loc[(metadata.df['lib_layout']==args.layout),:]
    print('SRA IDs:', ' '.join(metadata.df['run'].tolist()))
    max_bp = int(args.max_bp.replace(',',''))
    num_sra = metadata.df.shape[0]
    num_bp_per_sra = int(max_bp/num_sra)
    total_sra_bp = int(metadata.df['total_bases'].astype(int).sum())
    offset = 10000 # https://edwards.sdsu.edu/research/fastq-dump/
    print('Number of SRA', num_sra)
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
        sra_id = metadata.df.loc[i,'run']
        layout = metadata.df.loc[i,'lib_layout']
        total_spot = int(metadata.df.loc[i,'total_spots'])
        spot_length = int(metadata.df.loc[i,'spot_length'])
        num_read_per_sra = int(num_bp_per_sra/spot_length)
        print('SRA ID:', sra_id)
        print('Library layout:', layout)
        print('Number of reads:', "{:,}".format(total_spot))
        print('Single/Paired read length:', spot_length, 'bp')
        print('Total bases:', "{:,}".format(int(metadata.df.loc[i,'total_bases'])), 'bp')
        pfd_command = ['parallel-fastq-dump', '-t', str(args.threads), '--minReadLen', '25', '--qual-filter-1',
                       '--skip-technical', '--split-3', '--clip', '--gzip', '--outdir', args.work_dir,
                       '--tmpdir', args.work_dir]
        start,end = get_range(max_bp, total_sra_bp, total_spot, num_read_per_sra, offset)
        print('Total sampled bases:', "{:,}".format(spot_length*(end-start+1)), 'bp')
        pfd_command = pfd_command + ['--minSpotId', str(start), '--maxSpotId', str(end)]
        pfd_command = pfd_command + ['--sra-id', sra_id]
        print('Command:', ' '.join(pfd_command))
        if args.pfd=='yes':
            pfd_out = subprocess.run(pfd_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout = pfd_out.stdout.decode('utf8')
            nd = [ int(line.replace('Read ','').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Read') ]
            nr = [ int(line.replace('Rejected ','').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Rejected') ]
            nw = [ int(line.replace('Written ','').split(' ')[0]) for line in stdout.split('\n') if line.startswith('Written') ]
            bp_dumped += sum(nd) * spot_length
            bp_rejected += sum(nr) * spot_length
            bp_written += sum(nw) * spot_length
        if args.fastp=='yes':
            print('Running fastp.')
            inext = '.fastq.gz'
            outext = '.fastp.fastq.gz'
            fp_command = ['fastp', '--thread', str(args.threads)] + args.fastp_option.split(' ')
            if layout=='single':
                infile = os.path.join(args.work_dir,sra_id)
                fp_command = fp_command + ['--in1',infile+inext,'--out1',infile+outext]
            elif layout=='paired':
                infile1 = os.path.join(args.work_dir,sra_id+'_1')
                infile2 = os.path.join(args.work_dir,sra_id+'_2')
                fp_command = fp_command + ['--in1',infile1+inext,'--out1',infile1+outext,'--in2',infile2+inext,'--out2',infile2+outext]
            fp_command = [ fc for fc in fp_command if fc!='' ]
            print('Command:', ' '.join(fp_command))
            fp_out = subprocess.run(fp_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if layout=='single':
                if os.path.exists(infile+outext):
                    os.remove(infile+inext)
            if layout=='paired':
                if os.path.exists(infile1+outext)&os.path.exists(infile2+outext):
                    os.remove(infile1+inext)
                    os.remove(infile2+inext)
            bps = fp_out.stderr.decode('utf8').split('\n')
            bp_in = [ int(line.replace('total bases: ','').split(' ')[0]) for line in bps if line.startswith('total bases') ][0::2]
            bp_out = [ int(line.replace('total bases: ','').split(' ')[0]) for line in bps if line.startswith('total bases') ][1::2]
            bp_fastp_in += sum(bp_in)
            bp_fastp_out += sum(bp_out)
    print('')
    if args.pfd=='yes':
        print('max_bp:', "{:,}".format(max_bp), 'bp')
        print('Sum of dumped reads:', "{:,}".format(bp_dumped), 'bp')
        print('Sum of rejected reads:', "{:,}".format(bp_rejected), 'bp')
        print('Sum of written reads:', "{:,}".format(bp_written), 'bp')
        if args.fastp=='yes':
            print('Sum of fastp input reads:', "{:,}".format(bp_fastp_in), 'bp')
            print('Sum of fastp output reads:', "{:,}".format(bp_fastp_out), 'bp')
        print('')
        concat_fastq(args, metadata, clean=True)
    if args.remove_sra=='yes':
        remove_sra_files(metadata, sra_dir)
    else:
        if args.pfd=='yes':
            print('SRA files not removed:', sra_dir)
