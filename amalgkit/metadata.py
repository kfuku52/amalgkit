import pandas

import os
import sys
import time  # noqa: F401
import xml.etree.ElementTree as ET  # noqa: F401

from Bio import Entrez

from amalgkit.exceptions import AmalgkitExit
from amalgkit.metadata_utils import Metadata
from amalgkit.output_utils import atomic_write_dataframe
from amalgkit.sra import (
    esearch_sra_with_retry as _esearch_sra_with_retry,
    fetch_sra_xml as _fetch_sra_xml,
    fetch_sra_xml_chunk as _fetch_sra_xml_chunk,
    iter_sra_xml_chunks as _iter_sra_xml_chunks,
    merge_xml_chunk as _merge_xml_chunk,
    raise_if_xml_has_error as _raise_if_xml_has_error,
    search_sra_record_ids as _search_sra_record_ids,
)

def merge_xml_chunk(root, chunk):
    return _merge_xml_chunk(root, chunk)


def esearch_sra_with_retry(search_term):
    return _esearch_sra_with_retry(search_term)


def fetch_sra_xml_chunk(record_ids, start, end, retmax, max_retry=10):
    return _fetch_sra_xml_chunk(record_ids, start, end, retmax, max_retry=max_retry, verbose=True)


def raise_if_xml_has_error(root):
    return _raise_if_xml_has_error(root)


def search_sra_record_ids(search_term):
    return _search_sra_record_ids(search_term, verbose=True)


def iter_sra_xml_chunks(record_ids, retmax=1000):
    for chunk in _iter_sra_xml_chunks(
        record_ids=record_ids,
        retmax=retmax,
        verbose=True,
        timestamp_logs=True,
        progress_label='Retrieving SRA XML',
    ):
        raise_if_xml_has_error(chunk)
        yield chunk


def fetch_sra_xml(search_term, retmax=1000):
    return _fetch_sra_xml(
        search_term=search_term,
        retmax=retmax,
        verbose=True,
        timestamp_logs=True,
        progress_label='Retrieving SRA XML',
    )

def metadata_main(args):
    Entrez.email = args.entrez_email
    metadata_dir = os.path.join(args.out_dir, 'metadata')
    metadata_outfile_path = os.path.join(metadata_dir, "metadata.tsv")
    if os.path.exists(metadata_outfile_path):
        if not os.path.isfile(metadata_outfile_path):
            raise NotADirectoryError(
                'Output path exists but is not a file: {}'.format(metadata_outfile_path)
            )
        if not args.redo:
            raise AmalgkitExit(
                'Output file already exists (set --redo yes to overwrite): {}'.format(metadata_outfile_path),
                exit_code=0,
                use_stderr=False,
            )
    for path_dir in [args.out_dir, metadata_dir]:
        if os.path.exists(path_dir) and (not os.path.isdir(path_dir)):
            raise NotADirectoryError('Output path exists but is not a directory: {}'.format(path_dir))
        if not os.path.exists(path_dir):
            print('Creating directory: {}'.format(path_dir))
        os.makedirs(path_dir, exist_ok=True)
    search_term_raw = getattr(args, 'search_string', '')
    if search_term_raw is None:
        search_term = ''
    else:
        search_term = str(search_term_raw).strip()
    if search_term == '':
        raise ValueError('--search_string is required.')
    print('Entrez search term:', search_term)
    record_ids = search_sra_record_ids(search_term)
    metadata = Metadata.from_xml_roots(iter_sra_xml_chunks(record_ids=record_ids))
    metadata.df['tissue'] = metadata.df['tissue'].astype(str)
    metadata.df.loc[(metadata.df['tissue']=='nan'), 'tissue'] = ''
    metadata.df.loc[:, 'sample_group'] = metadata.df.loc[:, 'tissue'].str.lower()
    metadata.df['taxid'] = pandas.to_numeric(metadata.df['taxid'], errors='coerce').astype('Int64')
    metadata.add_standard_rank_taxids(args=args)
    if args.resolve_names:
        metadata.resolve_scientific_names(args=args)
    metadata.reorder(omit_misc=False)
    atomic_write_dataframe(metadata.df, metadata_outfile_path, sep="\t", index=False)
    if metadata.df.shape[0]==0:
        txt = 'No entry was found/survived in the metadata processing. Please reconsider the --search_string specification.\n'
        sys.stderr.write(txt)
