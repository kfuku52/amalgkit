from Bio import Entrez
from contextlib import contextmanager
from contextvars import ContextVar
import datetime
import email.utils
from defusedxml.ElementTree import parse as parse_untrusted_xml
from http.client import IncompleteRead
import math
import random
import re
import time
import xml.etree.ElementTree as ET

from urllib.error import HTTPError, URLError

from amalgkit.download_utils import maybe_acquire_download_semaphore
from amalgkit.subprocess_utils import resolve_timeout_seconds

SRA_ACCESSION_PATTERN = re.compile(r'\b(?:[SED](?:RR|RP|RS|RX)\d+)\b', re.IGNORECASE)
NCBI_METADATA_TIMEOUT_SECONDS = 300
NCBI_RETRY_DELAY_CAP_SECONDS = 900
SRA_ESEARCH_PAGE_SIZE = 100_000
_ENTREZ_TIMEOUT_CONTEXT = ContextVar('amalgkit_entrez_timeout_seconds', default=None)
_ENTREZ_ORIGINAL_URLOPEN = getattr(
    Entrez.urlopen,
    '_amalgkit_original_urlopen',
    Entrez.urlopen,
)


def _entrez_urlopen_with_timeout(*args, **kwargs):
    timeout_seconds = _ENTREZ_TIMEOUT_CONTEXT.get()
    if timeout_seconds is not None:
        kwargs.setdefault('timeout', float(timeout_seconds))
    return _ENTREZ_ORIGINAL_URLOPEN(*args, **kwargs)


_entrez_urlopen_with_timeout._amalgkit_original_urlopen = _ENTREZ_ORIGINAL_URLOPEN
if not hasattr(Entrez.urlopen, '_amalgkit_original_urlopen'):
    # Bio.Entrez does not expose a request-timeout argument. Its private request
    # helper resolves this module-level urlopen symbol, so a context-local
    # wrapper lets concurrent amalgkit workers use their own timeout without
    # changing the process-wide socket default.
    Entrez.urlopen = _entrez_urlopen_with_timeout


def resolve_ncbi_metadata_timeout_seconds(args, default_seconds=NCBI_METADATA_TIMEOUT_SECONDS):
    return resolve_timeout_seconds(
        args=args,
        attribute_name='ncbi_metadata_timeout_seconds',
        default_seconds=default_seconds,
    )


@contextmanager
def entrez_request_timeout(args=None):
    timeout_seconds = resolve_ncbi_metadata_timeout_seconds(args)
    token = _ENTREZ_TIMEOUT_CONTEXT.set(timeout_seconds)
    try:
        yield timeout_seconds
    finally:
        _ENTREZ_TIMEOUT_CONTEXT.reset(token)


def _close_entrez_handle(handle):
    close = getattr(handle, 'close', None)
    if callable(close):
        close()


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


def _parse_retry_after_seconds(exc):
    headers = getattr(exc, 'headers', None)
    if headers is None:
        return None
    retry_after = headers.get('Retry-After')
    if retry_after is None:
        return None
    retry_after = str(retry_after).strip()
    try:
        seconds = float(retry_after)
    except ValueError:
        try:
            retry_at = email.utils.parsedate_to_datetime(retry_after)
        except (TypeError, ValueError, OverflowError):
            return None
        if retry_at.tzinfo is None:
            retry_at = retry_at.replace(tzinfo=datetime.timezone.utc)
        seconds = (
            retry_at - datetime.datetime.now(datetime.timezone.utc)
        ).total_seconds()
    if (not math.isfinite(seconds)) or seconds < 0:
        return None
    return min(float(seconds), float(NCBI_RETRY_DELAY_CAP_SECONDS))


def _calculate_entrez_retry_delay(exc, retry_sleep_second, failure_index):
    retry_after = _parse_retry_after_seconds(exc)
    if retry_after is not None:
        return retry_after
    base_delay = float(retry_sleep_second)
    if (not math.isfinite(base_delay)) or base_delay < 0:
        raise ValueError('retry_sleep_second must be a finite value >= 0.')
    exponential_delay = min(
        base_delay * (2 ** min(30, max(0, int(failure_index)))),
        float(NCBI_RETRY_DELAY_CAP_SECONDS),
    )
    if exponential_delay <= 0:
        return 0.0
    jitter_cap = min(1.0, exponential_delay * 0.1)
    return min(
        exponential_delay + random.uniform(0.0, jitter_cap),  # noqa: S311 - retry jitter is not security-sensitive
        float(NCBI_RETRY_DELAY_CAP_SECONDS),
    )


def esearch_sra_with_retry(
    search_term,
    args=None,
    max_retry=2,
    retry_sleep_second=1,
    retstart=0,
    retmax=SRA_ESEARCH_PAGE_SIZE,
):
    max_retry = int(max_retry)
    if max_retry <= 0:
        raise ValueError('max_retry must be > 0.')
    for failure_index in range(max_retry):
        try:
            with maybe_acquire_download_semaphore(
                args=args,
                limit_attr='ncbi_metadata_max_concurrency',
                semaphore_name='ncbi_metadata',
                lock_label='NCBI metadata download',
            ), entrez_request_timeout(args):
                sra_handle = Entrez.esearch(
                    db='sra',
                    term=search_term,
                    retstart=int(retstart),
                    retmax=int(retmax),
                )
                try:
                    return Entrez.read(sra_handle)
                finally:
                    _close_entrez_handle(sra_handle)
        except (HTTPError, URLError, TimeoutError) as exc:
            if failure_index + 1 >= max_retry:
                raise
            delay_seconds = _calculate_entrez_retry_delay(
                exc=exc,
                retry_sleep_second=retry_sleep_second,
                failure_index=failure_index,
            )
            print(
                '{} - Trying Entrez.esearch() again after {:.1f} seconds...'.format(
                    exc,
                    delay_seconds,
                ),
                flush=True,
            )
            time.sleep(delay_seconds)


def fetch_sra_xml_chunk(record_ids, start, end, retmax, max_retry=10, verbose=True, retry_sleep_second=60, args=None):
    max_retry = int(max_retry)
    if max_retry <= 0:
        raise ValueError('max_retry must be > 0.')
    last_exception = None
    for failure_index in range(max_retry):
        try:
            with maybe_acquire_download_semaphore(
                args=args,
                limit_attr='ncbi_metadata_max_concurrency',
                semaphore_name='ncbi_metadata',
                lock_label='NCBI metadata download',
            ), entrez_request_timeout(args):
                handle = Entrez.efetch(db='sra', id=record_ids[start:end], rettype='full', retmode='xml', retmax=retmax)
                try:
                    return parse_untrusted_xml(handle).getroot()
                finally:
                    _close_entrez_handle(handle)
        except (HTTPError, URLError, TimeoutError) as exc:
            last_exception = exc
            if failure_index + 1 >= max_retry:
                break
            delay_seconds = _calculate_entrez_retry_delay(
                exc=exc,
                retry_sleep_second=retry_sleep_second,
                failure_index=failure_index,
            )
            if verbose:
                print(
                    '{} - Trying Entrez.efetch() again after {:.1f} seconds...'.format(
                        exc,
                        delay_seconds,
                    ),
                    flush=True,
                )
            time.sleep(delay_seconds)
            continue
        except (ET.ParseError, IncompleteRead) as exc:
            last_exception = exc
            if failure_index + 1 >= max_retry:
                break
            delay_seconds = _calculate_entrez_retry_delay(
                exc=exc,
                retry_sleep_second=retry_sleep_second,
                failure_index=failure_index,
            )
            if verbose:
                print(
                    'XML may be truncated. Retrying after {:.1f} seconds...'.format(
                        delay_seconds,
                    ),
                    flush=True,
                )
            time.sleep(delay_seconds)
            continue
    error = RuntimeError(
        'Failed to parse Entrez XML chunk after {} retries (records {}-{}).'.format(
            max_retry,
            start,
            end - 1,
        )
    )
    raise error from last_exception


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


def search_sra_record_ids(search_term, verbose=True, args=None):
    sra_record = esearch_sra_with_retry(
        search_term,
        args=args,
        retstart=0,
        retmax=SRA_ESEARCH_PAGE_SIZE,
    )
    record_ids = [str(record_id) for record_id in sra_record.get('IdList', [])]
    raw_count = sra_record.get('Count', len(record_ids))
    try:
        total_count = int(raw_count)
    except (TypeError, ValueError) as exc:
        raise RuntimeError('Entrez.esearch returned an invalid Count: {}'.format(raw_count)) from exc
    if total_count < 0:
        raise RuntimeError('Entrez.esearch returned a negative Count: {}'.format(total_count))
    retstart = len(record_ids)
    while retstart < total_count:
        page = esearch_sra_with_retry(
            search_term,
            args=args,
            retstart=retstart,
            retmax=min(SRA_ESEARCH_PAGE_SIZE, total_count - retstart),
        )
        page_ids = [str(record_id) for record_id in page.get('IdList', [])]
        if len(page_ids) == 0:
            raise RuntimeError(
                'Entrez.esearch returned no IDs at offset {:,} of {:,}; refusing a truncated result.'.format(
                    retstart,
                    total_count,
                )
            )
        record_ids.extend(page_ids)
        retstart += len(page_ids)
    if len(record_ids) != total_count:
        raise RuntimeError(
            'Entrez.esearch Count/IdList mismatch: expected {:,}, received {:,}.'.format(
                total_count,
                len(record_ids),
            )
        )
    if len(set(record_ids)) != len(record_ids):
        raise RuntimeError('Entrez.esearch returned duplicate IDs across result pages.')
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


def _call_fetch_sra_xml_chunk(
    record_ids,
    start,
    end,
    retmax,
    max_retry=10,
    verbose=True,
    retry_sleep_second=60,
    args=None,
):
    kwargs = {
        'record_ids': record_ids,
        'start': start,
        'end': end,
        'retmax': retmax,
        'max_retry': max_retry,
        'verbose': verbose,
        'retry_sleep_second': retry_sleep_second,
    }
    if args is not None:
        kwargs['args'] = args
    try:
        return fetch_sra_xml_chunk(**kwargs)
    except TypeError as exc:
        if ('unexpected keyword argument' not in str(exc)) or ('args' not in str(exc)):
            raise
    kwargs.pop('args', None)
    return fetch_sra_xml_chunk(**kwargs)


def summarize_sra_record(record_id):
    root = _call_fetch_sra_xml_chunk(
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
        except Exception:  # noqa: S112 - optional mismatch diagnostics are best-effort
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


def iter_sra_xml_chunks(record_ids, retmax=1000, verbose=True, timestamp_logs=True, progress_label='Retrieving SRA XML', args=None):
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
        chunk = _call_fetch_sra_xml_chunk(
            record_ids=record_ids,
            start=start,
            end=end,
            retmax=retmax,
            max_retry=10,
            verbose=verbose,
            args=args,
        )
        yield chunk
    if verbose and timestamp_logs:
        elapsed_time = int(time.time() - start_time)
        print('{}: SRA XML retrieval ended.'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        print('SRA XML retrieval time: {:,.1f} sec'.format(elapsed_time), flush=True)


def fetch_sra_xml(search_term, retmax=1000, verbose=True, timestamp_logs=True, progress_label='Retrieving SRA XML', args=None):
    record_ids = search_sra_record_ids(search_term, verbose=verbose, args=args)
    if len(record_ids) == 0:
        return ET.Element('EXPERIMENT_PACKAGE_SET')
    root = None
    for chunk in iter_sra_xml_chunks(
        record_ids=record_ids,
        retmax=retmax,
        verbose=verbose,
        timestamp_logs=timestamp_logs,
        progress_label=progress_label,
        args=args,
    ):
        raise_if_xml_has_error(chunk, search_term=search_term)
        if root is None:
            root = chunk
        else:
            root = merge_xml_chunk(root, chunk)
    return root
