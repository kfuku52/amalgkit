from Bio import Entrez
import xml.etree.ElementTree as ET
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

def merge_xml_chunk(root, chunk):
    # Merge package-set chunks directly to avoid nested container nodes.
    if (chunk.tag == root.tag) and (len(chunk) > 0):
        root.extend(list(chunk))
    else:
        root.append(chunk)


def esearch_sra_with_retry(search_term):
    try:
        sra_handle = Entrez.esearch(db="sra", term=search_term, retmax=10000000)
    except HTTPError as e:
        print(e, '- Trying Entrez.esearch() again...')
        sra_handle = Entrez.esearch(db="sra", term=search_term, retmax=10000000)
    return Entrez.read(sra_handle)


def fetch_sra_xml_chunk(record_ids, start, end, retmax, max_retry=10):
    for _ in range(max_retry):
        try:
            handle = Entrez.efetch(db="sra", id=record_ids[start:end], rettype="full", retmode="xml", retmax=retmax)
        except HTTPError as e:
            sleep_second = 60
            print('{} - Trying Entrez.efetch() again after {:,} seconds...'.format(e, sleep_second), flush=True)
            time.sleep(sleep_second)
            continue
        try:
            return ET.parse(handle).getroot()
        except Exception:
            print('XML may be truncated. Retrying...', flush=True)
            continue
    raise RuntimeError('Failed to parse Entrez XML chunk after {} retries (records {}-{}).'.format(
        max_retry, start, end - 1
    ))


def raise_if_xml_has_error(root):
    error_node = root.find('.//Error')
    if error_node is None:
        return
    error_text = ''.join(error_node.itertext()).strip()
    if error_text != '':
        print(error_text)
    raise Exception(',Error. found in the xml.')


def fetch_sra_xml(search_term, retmax=1000):
    sra_record = esearch_sra_with_retry(search_term)
    record_ids = sra_record["IdList"]
    num_record = len(record_ids)
    print('Number of SRA records: {:,}'.format(num_record))
    if num_record == 0:
        return ET.Element('EXPERIMENT_PACKAGE_SET')
    start_time = time.time()
    print('{}: SRA XML retrieval started.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    root = None
    for start in range(0, num_record, retmax):
        end = min(start + retmax, num_record)
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print('{}: Retrieving SRA XML: {:,}-{:,} of {:,} records'.format(now, start, end - 1, num_record), flush=True)
        chunk = fetch_sra_xml_chunk(record_ids, start, end, retmax, max_retry=10)
        if root is None:
            root = chunk
        else:
            merge_xml_chunk(root, chunk)
    elapsed_time = int(time.time() - start_time)
    print('{}: SRA XML retrieval ended.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print('SRA XML retrieval time: {:,.1f} sec'.format(elapsed_time), flush=True)
    raise_if_xml_has_error(root)
    return root

def metadata_main(args):
    Entrez.email = args.entrez_email
    metadata_dir = os.path.join(args.out_dir, 'metadata')
    metadata_outfile_path = os.path.join(metadata_dir, "metadata.tsv")
    if os.path.exists(metadata_outfile_path):
        if args.redo:
            os.remove(metadata_outfile_path)
        else:
            print('Exiting. Output file already exists (set --redo yes to overwrite): {}'.format(metadata_outfile_path))
            sys.exit(0)
    for path_dir in [args.out_dir, metadata_dir]:
        if not os.path.exists(path_dir):
            print('Creating directory: {}'.format(path_dir))
        os.makedirs(path_dir, exist_ok=True)
    search_term = args.search_string
    print('Entrez search term:', search_term)
    root = fetch_sra_xml(search_term=search_term)
    metadata = Metadata.from_xml(xml_root=root)
    metadata.df['tissue'] = metadata.df['tissue'].astype(str)
    metadata.df.loc[(metadata.df['tissue']=='nan'), 'tissue'] = ''
    metadata.df.loc[:, 'sample_group'] = metadata.df.loc[:, 'tissue'].str.lower()
    metadata.df['taxid'] = pandas.to_numeric(metadata.df['taxid'], errors='coerce').astype('Int64')
    metadata.add_standard_rank_taxids()
    if args.resolve_names:
        metadata.resolve_scientific_names()
    metadata.reorder(omit_misc=False)
    metadata.df.to_csv(metadata_outfile_path, sep="\t", index=False)
    if metadata.df.shape[0]==0:
        txt = 'No entry was found/survived in the metadata processing. Please reconsider the --search_term specification.\n'
        sys.stderr.write(txt)
