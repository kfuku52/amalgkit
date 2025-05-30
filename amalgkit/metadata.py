from Bio import Entrez
import lxml.etree
import numpy
import pandas

import datetime
import os
import re
import sys
import time
import warnings

from amalgkit.util import *

from urllib.error import HTTPError

def fetch_sra_xml(search_term, retmax=1000):
    try:
        sra_handle = Entrez.esearch(db="sra", term=search_term, retmax=10000000)
    except HTTPError as e:
        print(e, '- Trying Entrez.esearch() again...')
        sra_handle = Entrez.esearch(db="sra", term=search_term, retmax=10000000)
    sra_record = Entrez.read(sra_handle)
    record_ids = sra_record["IdList"]
    num_record = len(record_ids)
    print('Number of SRA records: {:,}'.format(num_record))
    start_time = time.time()
    print('{}: SRA XML retrieval started.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    root = None
    max_retry = 10
    for i in numpy.arange(numpy.ceil(num_record//retmax)+1):
        start = int(i*retmax)
        end = int(((i+1)*retmax)-1) if num_record >= int(((i+1)*retmax)-1) else num_record
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print('{}: Retrieving SRA XML: {:,}-{:,} of {:,} records'.format(now, start, end, num_record), flush=True)
        for j in range(max_retry):
            try:
                handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
            except HTTPError as e:
                sleep_second = 60
                print('{} - Trying Entrez.efetch() again after {:,} seconds...'.format(e, sleep_second), flush=True)
                time.sleep(sleep_second)
                continue
            try:
                chunk = lxml.etree.parse(handle).getroot()
            except:
                print('XML may be truncated. Retrying...', flush=True)
                continue
            break
        if root is None:
            root = chunk
        else:
            root.append(chunk)
    elapsed_time = int(time.time() - start_time)
    print('{}: SRA XML retrieval ended.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print('SRA XML retrieval time: {:,.1f} sec'.format(elapsed_time), flush=True)
    xml_string = lxml.etree.tostring(root, encoding='UTF-8', pretty_print=True)
    for line in str(xml_string).split('\n'):
        if '<Error>' in line:
            print(line)
            raise Exception(',Error. found in the xml.')
    return root

def metadata_main(args):
    Entrez.email = args.entrez_email
    metadata_dir = os.path.join(args.out_dir, 'metadata')
    metadata_outfile_path = os.path.join(metadata_dir, "metadata.tsv")
    if os.path.exists(metadata_outfile_path):
        if args.redo:
            os.remove(metadata_outfile_path)
        else:
            print('Exiting. --redo is specified and the output file already exists at: {}'.format(metadata_outfile_path))
            sys.exit(0)
    for path_dir in [args.out_dir, metadata_dir]:
        if not os.path.exists(path_dir):
            print('Creating directory: {}'.format(path_dir))
            os.mkdir(path_dir)
    search_term = args.search_string
    print('Entrez search term:', search_term)
    root = fetch_sra_xml(search_term=search_term)
    metadata = Metadata.from_xml(xml_root=root)
    metadata.df['tissue'] = metadata.df['tissue'].astype(str)
    metadata.df.loc[(metadata.df['tissue']=='nan'), 'tissue'] = ''
    metadata.df.loc[:, 'sample_group'] = metadata.df.loc[:, 'tissue'].str.lower()
    metadata.reorder(omit_misc=False)
    metadata.df.to_csv(metadata_outfile_path, sep="\t", index=False)
    if metadata.df.shape[0]==0:
        txt = 'No entry was found/survived in the metadata processing. Please reconsider the --search_term specification.\n'
        sys.stderr.write(txt)
