import xml.etree.ElementTree as ET
from io import BytesIO
from types import SimpleNamespace
from urllib.error import URLError

import pytest
from defusedxml.common import EntitiesForbidden

from amalgkit.sra import (
    _entrez_urlopen_with_timeout,
    entrez_request_timeout,
    esearch_sra_with_retry,
    extract_sra_accessions,
    fetch_sra_xml_chunk,
    inspect_accession_search_mismatches,
)
import amalgkit.sra as sra


class RecordingBytesIO(BytesIO):
    def __init__(self, value):
        super().__init__(value)
        self.was_closed = False

    def close(self):
        self.was_closed = True
        super().close()


class TestExtractSraAccessions:
    def test_deduplicates_and_preserves_order(self):
        search_term = '(srr29779594) AND PRJNA1 AND SRR29779594 AND ERX123'
        assert extract_sra_accessions(search_term) == ['SRR29779594', 'ERX123']

    def test_returns_empty_for_blank_query(self):
        assert extract_sra_accessions('') == []
        assert extract_sra_accessions(None) == []


def test_fetch_sra_xml_chunk_rejects_xml_entities(monkeypatch):
    malicious_xml = b'<!DOCTYPE root [<!ENTITY xxe SYSTEM "file:///etc/passwd">]><root>&xxe;</root>'
    handle = RecordingBytesIO(malicious_xml)
    monkeypatch.setattr(
        'amalgkit.sra.Entrez.efetch',
        lambda **_kwargs: handle,
    )

    with pytest.raises(EntitiesForbidden):
        fetch_sra_xml_chunk(
            record_ids=['ID1'],
            start=0,
            end=1,
            retmax=1,
            max_retry=1,
            retry_sleep_second=0,
        )
    assert handle.was_closed


def test_esearch_closes_entrez_handle(monkeypatch):
    handle = RecordingBytesIO(b'ignored')
    monkeypatch.setattr('amalgkit.sra.Entrez.esearch', lambda **_kwargs: handle)
    monkeypatch.setattr('amalgkit.sra.Entrez.read', lambda _handle: {'IdList': ['ID1']})

    assert esearch_sra_with_retry('SRR000001') == {'IdList': ['ID1']}
    assert handle.was_closed


def test_esearch_retry_uses_backoff_and_closes_successful_handle(monkeypatch):
    handle = RecordingBytesIO(b'ignored')
    calls = []
    sleeps = []

    def fake_esearch(**_kwargs):
        calls.append(True)
        if len(calls) == 1:
            raise URLError('temporary failure')
        return handle

    monkeypatch.setattr('amalgkit.sra.Entrez.esearch', fake_esearch)
    monkeypatch.setattr('amalgkit.sra.Entrez.read', lambda _handle: {'IdList': ['ID1']})
    monkeypatch.setattr('amalgkit.sra._calculate_entrez_retry_delay', lambda **_kwargs: 2.5)
    monkeypatch.setattr('amalgkit.sra.time.sleep', sleeps.append)

    assert esearch_sra_with_retry('SRR000001') == {'IdList': ['ID1']}
    assert len(calls) == 2
    assert sleeps == [2.5]
    assert handle.was_closed


def test_esearch_final_failure_does_not_sleep(monkeypatch):
    sleeps = []
    monkeypatch.setattr(
        'amalgkit.sra.Entrez.esearch',
        lambda **_kwargs: (_ for _ in ()).throw(URLError('unavailable')),
    )
    monkeypatch.setattr('amalgkit.sra.time.sleep', sleeps.append)

    with pytest.raises(URLError):
        esearch_sra_with_retry(
            'SRR000001',
            max_retry=1,
            retry_sleep_second=60,
        )

    assert sleeps == []


def test_search_sra_record_ids_fetches_all_counted_pages(monkeypatch):
    calls = []

    def fake_esearch(search_term, args=None, max_retry=2, retry_sleep_second=1, retstart=0, retmax=100000):
        _ = (search_term, args, max_retry, retry_sleep_second, retmax)
        calls.append(retstart)
        pages = {
            0: {'Count': '5', 'IdList': ['ID1', 'ID2']},
            2: {'Count': '5', 'IdList': ['ID3', 'ID4']},
            4: {'Count': '5', 'IdList': ['ID5']},
        }
        return pages[retstart]

    monkeypatch.setattr(sra, 'esearch_sra_with_retry', fake_esearch)

    assert sra.search_sra_record_ids('query', verbose=False) == [
        'ID1', 'ID2', 'ID3', 'ID4', 'ID5',
    ]
    assert calls == [0, 2, 4]


def test_search_sra_record_ids_rejects_truncated_page(monkeypatch):
    monkeypatch.setattr(
        sra,
        'esearch_sra_with_retry',
        lambda *_args, **_kwargs: {'Count': '3', 'IdList': ['ID1', 'ID2']}
        if _kwargs.get('retstart', 0) == 0
        else {'Count': '3', 'IdList': []},
    )

    with pytest.raises(RuntimeError, match='refusing a truncated result'):
        sra.search_sra_record_ids('query', verbose=False)


def test_fetch_final_failure_does_not_sleep(monkeypatch):
    sleeps = []
    monkeypatch.setattr(
        'amalgkit.sra.Entrez.efetch',
        lambda **_kwargs: (_ for _ in ()).throw(URLError('unavailable')),
    )
    monkeypatch.setattr('amalgkit.sra.time.sleep', sleeps.append)

    with pytest.raises(RuntimeError, match='after 1 retries'):
        fetch_sra_xml_chunk(
            record_ids=['ID1'],
            start=0,
            end=1,
            retmax=1,
            max_retry=1,
            retry_sleep_second=60,
            verbose=False,
        )

    assert sleeps == []


def test_retry_after_header_takes_precedence_over_exponential_delay():
    exc = SimpleNamespace(headers={'Retry-After': '7'})

    assert sra._calculate_entrez_retry_delay(
        exc=exc,
        retry_sleep_second=60,
        failure_index=3,
    ) == 7.0


def test_entrez_timeout_context_passes_cli_value_to_urlopen(monkeypatch):
    observed = {}

    def fake_urlopen(*args, **kwargs):
        observed['args'] = args
        observed['kwargs'] = kwargs
        return object()

    monkeypatch.setattr(sra, '_ENTREZ_ORIGINAL_URLOPEN', fake_urlopen)
    with entrez_request_timeout(SimpleNamespace(ncbi_metadata_timeout_seconds=17)):
        _entrez_urlopen_with_timeout('request')

    assert observed['args'] == ('request',)
    assert observed['kwargs']['timeout'] == 17.0


def test_entrez_timeout_can_be_disabled(monkeypatch):
    observed = {}

    def fake_urlopen(*_args, **kwargs):
        observed['kwargs'] = kwargs
        return object()

    monkeypatch.setattr(sra, '_ENTREZ_ORIGINAL_URLOPEN', fake_urlopen)
    with entrez_request_timeout(SimpleNamespace(ncbi_metadata_timeout_seconds=0)):
        _entrez_urlopen_with_timeout('request')

    assert 'timeout' not in observed['kwargs']


class TestInspectAccessionSearchMismatches:
    def test_collects_platform_and_library_fields(self, monkeypatch):
        def fake_search(search_term, verbose=True):
            if search_term == 'SRR29779594':
                return ['UID1']
            return []

        def fake_fetch(record_ids, start, end, retmax, max_retry=10, verbose=True, retry_sleep_second=60):
            root = ET.Element('EXPERIMENT_PACKAGE_SET')
            pkg = ET.SubElement(root, 'EXPERIMENT_PACKAGE')
            sample = ET.SubElement(pkg, 'SAMPLE')
            sample_name = ET.SubElement(sample, 'SAMPLE_NAME')
            ET.SubElement(sample_name, 'SCIENTIFIC_NAME').text = 'Acanthamoeba astronyxis'
            experiment = ET.SubElement(pkg, 'EXPERIMENT')
            platform = ET.SubElement(experiment, 'PLATFORM')
            ET.SubElement(platform, 'CAPILLARY', instrument_model='AB 3730xL Genetic Analyzer')
            lib = ET.SubElement(experiment, 'LIBRARY_DESCRIPTOR')
            ET.SubElement(lib, 'LIBRARY_STRATEGY').text = 'CLONE'
            ET.SubElement(lib, 'LIBRARY_SOURCE').text = 'TRANSCRIPTOMIC'
            ET.SubElement(lib, 'LIBRARY_SELECTION').text = 'cDNA'
            return root

        monkeypatch.setattr('amalgkit.sra.search_sra_record_ids', fake_search)
        monkeypatch.setattr('amalgkit.sra.fetch_sra_xml_chunk', fake_fetch)

        diagnostics = inspect_accession_search_mismatches(
            '(SRR29779594) AND "Illumina"[Platform] AND ("RNA-seq"[Strategy] OR "EST"[Strategy])'
        )

        assert diagnostics == [
            {
                'accession': 'SRR29779594',
                'record_id': 'UID1',
                'matched_record_count': 1,
                'scientific_name': 'Acanthamoeba astronyxis',
                'platform': 'CAPILLARY (AB 3730xL Genetic Analyzer)',
                'library_strategy': 'CLONE',
                'library_source': 'TRANSCRIPTOMIC',
                'library_selection': 'cDNA',
            }
        ]
