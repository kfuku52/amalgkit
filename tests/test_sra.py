import xml.etree.ElementTree as ET
from io import BytesIO

import pytest
from defusedxml.common import EntitiesForbidden

from amalgkit.sra import (
    esearch_sra_with_retry,
    extract_sra_accessions,
    fetch_sra_xml_chunk,
    inspect_accession_search_mismatches,
)


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
