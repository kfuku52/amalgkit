from Bio import Entrez
import datetime
from http.client import IncompleteRead
import re
import time
import xml.etree.ElementTree as ET

from urllib.error import HTTPError, URLError

SRA_ACCESSION_PATTERN = re.compile(r'\b(?:[SED](?:RR|RP|RS|RX)\d+)\b', re.IGNORECASE)


def merge_xml_chunk(root, chunk):
    if (chunk.tag == root.tag) and root.tag.endswith('_SET'):
        root.extend(list(chunk))
        return root
    if root.tag.endswith('_SET'):
        root.append(chunk)
        return root
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


def esearch_sra_with_retry(search_term):
    try:
        sra_handle = Entrez.esearch(db='sra', term=search_term, retmax=10000000)
    except (HTTPError, URLError) as exc:
        print(exc, '- Trying Entrez.esearch() again...')
        sra_handle = Entrez.esearch(db='sra', term=search_term, retmax=10000000)
    return Entrez.read(sra_handle)


def fetch_sra_xml_chunk(record_ids, start, end, retmax, max_retry=10, verbose=True, retry_sleep_second=60):
    for _ in range(max_retry):
        try:
            handle = Entrez.efetch(db='sra', id=record_ids[start:end], rettype='full', retmode='xml', retmax=retmax)
        except (HTTPError, URLError) as exc:
            if verbose:
                print(
                    '{} - Trying Entrez.efetch() again after {:,} seconds...'.format(exc, int(retry_sleep_second)),
                    flush=True,
                )
            time.sleep(retry_sleep_second)
            continue
        try:
            return ET.parse(handle).getroot()
        except (ET.ParseError, IncompleteRead):
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


def raise_if_xml_has_error(root, search_term=None):
    error_node = root.find('.//Error')
    if error_node is None:
        return
    error_text = ''.join(error_node.itertext()).strip()
    if error_text != '':
        print(error_text)
    suffix = ''
    if search_term is not None:
        suffix = ' Search term: {}'.format(search_term)
    raise RuntimeError('Error found in Entrez XML response.{}'.format(suffix))


def search_sra_record_ids(search_term, verbose=True):
    sra_record = esearch_sra_with_retry(search_term)
    record_ids = sra_record['IdList']
    if verbose:
        print('Number of SRA records: {:,}'.format(len(record_ids)))
    return record_ids


def extract_sra_accessions(search_term, max_count=None):
    if search_term in [None, '']:
        return []
    accessions = []
    seen = set()
    for match in SRA_ACCESSION_PATTERN.finditer(str(search_term)):
        accession = match.group(0).upper()
        if accession in seen:
            continue
        seen.add(accession)
        accessions.append(accession)
        if (max_count is not None) and (len(accessions) >= max_count):
            break
    return accessions


def _extract_sra_summary_fields(root):
    scientific_name_node = root.find('.//SAMPLE_NAME/SCIENTIFIC_NAME')
    scientific_name = ''
    if (scientific_name_node is not None) and (scientific_name_node.text is not None):
        scientific_name = scientific_name_node.text.strip()
    platform = ''
    platform_node = root.find('.//PLATFORM')
    if platform_node is not None:
        platform_children = list(platform_node)
        if len(platform_children) > 0:
            platform_child = platform_children[0]
            platform = platform_child.tag.strip()
            instrument_model = platform_child.attrib.get('instrument_model', '').strip()
            if (instrument_model != '') and (instrument_model != platform):
                platform = '{} ({})'.format(platform, instrument_model)
    def get_text(path):
        node = root.find(path)
        if (node is None) or (node.text is None):
            return ''
        return node.text.strip()
    return {
        'scientific_name': scientific_name,
        'platform': platform,
        'library_strategy': get_text('.//LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY'),
        'library_source': get_text('.//LIBRARY_DESCRIPTOR/LIBRARY_SOURCE'),
        'library_selection': get_text('.//LIBRARY_DESCRIPTOR/LIBRARY_SELECTION'),
    }


def summarize_sra_record(record_id):
    root = fetch_sra_xml_chunk(
        record_ids=[record_id],
        start=0,
        end=1,
        retmax=1,
        max_retry=2,
        verbose=False,
    )
    raise_if_xml_has_error(root)
    return _extract_sra_summary_fields(root)


def inspect_accession_search_mismatches(search_term, max_accessions=3):
    diagnostics = []
    for accession in extract_sra_accessions(search_term, max_count=max_accessions):
        try:
            record_ids = search_sra_record_ids(accession, verbose=False)
        except Exception:
            continue
        if len(record_ids) == 0:
            continue
        summary = {}
        try:
            summary = summarize_sra_record(record_ids[0])
        except Exception:
            summary = {}
        summary.update({
            'accession': accession,
            'record_id': record_ids[0],
            'matched_record_count': len(record_ids),
        })
        diagnostics.append(summary)
    return diagnostics


def iter_sra_xml_chunks(record_ids, retmax=1000, verbose=True, timestamp_logs=True, progress_label='Retrieving SRA XML'):
    num_record = len(record_ids)
    if num_record == 0:
        return
    start_time = time.time()
    if verbose and timestamp_logs:
        print('{}: SRA XML retrieval started.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    for start in range(0, num_record, retmax):
        end = min(start + retmax, num_record)
        if verbose:
            if timestamp_logs:
                now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                print(
                    '{}: {}: {:,}-{:,} of {:,} records'.format(now, progress_label, start, end - 1, num_record),
                    flush=True,
                )
            else:
                print('{}: {} - {}'.format(progress_label, start, end - 1), flush=True)
        chunk = fetch_sra_xml_chunk(record_ids, start, end, retmax, max_retry=10, verbose=verbose)
        yield chunk
    if verbose and timestamp_logs:
        elapsed_time = int(time.time() - start_time)
        print('{}: SRA XML retrieval ended.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        print('SRA XML retrieval time: {:,.1f} sec'.format(elapsed_time), flush=True)


def fetch_sra_xml(search_term, retmax=1000, verbose=True, timestamp_logs=True, progress_label='Retrieving SRA XML'):
    record_ids = search_sra_record_ids(search_term, verbose=verbose)
    if len(record_ids) == 0:
        return ET.Element('EXPERIMENT_PACKAGE_SET')
    root = None
    for chunk in iter_sra_xml_chunks(
        record_ids=record_ids,
        retmax=retmax,
        verbose=verbose,
        timestamp_logs=timestamp_logs,
        progress_label=progress_label,
    ):
        raise_if_xml_has_error(chunk, search_term=search_term)
        if root is None:
            root = chunk
        else:
            root = merge_xml_chunk(root, chunk)
    return root
