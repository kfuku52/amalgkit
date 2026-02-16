import xml.etree.ElementTree as ET

import pytest

from amalgkit.metadata import fetch_sra_xml


class TestFetchSraXml:
    class _DummyTree:
        def __init__(self, root):
            self._root = root

        def getroot(self):
            return self._root

    def test_returns_empty_root_when_no_records(self, monkeypatch):
        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': []})

        efetch_calls = {'n': 0}

        def fake_efetch(**kwargs):
            efetch_calls['n'] += 1
            return object()

        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', fake_efetch)

        root = fetch_sra_xml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE_SET'
        assert efetch_calls['n'] == 0

    def test_batches_without_extra_request_on_exact_multiple(self, monkeypatch):
        id_list = ['ID{}'.format(i) for i in range(2000)]
        efetch_calls = []

        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': id_list})

        def fake_efetch(**kwargs):
            efetch_calls.append(list(kwargs['id']))
            return object()

        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', fake_efetch)
        monkeypatch.setattr(
            'amalgkit.metadata.ET.parse',
            lambda handle: self._DummyTree(ET.Element('EXPERIMENT_PACKAGE'))
        )

        root = fetch_sra_xml(search_term='SRR000000', retmax=1000)

        assert root is not None
        assert len(efetch_calls) == 2
        assert [len(c) for c in efetch_calls] == [1000, 1000]

    def test_merges_package_set_chunks_without_nested_container(self, monkeypatch):
        id_list = ['ID{}'.format(i) for i in range(2000)]
        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': id_list})
        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', lambda **kwargs: object())

        def fake_parse(_handle):
            root = ET.Element('EXPERIMENT_PACKAGE_SET')
            ET.SubElement(root, 'EXPERIMENT_PACKAGE')
            return self._DummyTree(root)

        monkeypatch.setattr('amalgkit.metadata.ET.parse', fake_parse)

        root = fetch_sra_xml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE_SET'
        assert len(root.findall('./EXPERIMENT_PACKAGE')) == 2
        assert len(root.findall('./EXPERIMENT_PACKAGE_SET')) == 0

    def test_raises_when_xml_chunk_parsing_never_succeeds(self, monkeypatch):
        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.ET.parse', lambda handle: (_ for _ in ()).throw(Exception('broken xml')))

        with pytest.raises(RuntimeError, match='Failed to parse Entrez XML chunk'):
            fetch_sra_xml(search_term='SRR000000', retmax=1000)

    def test_raises_when_error_tag_present(self, monkeypatch):
        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', lambda **kwargs: object())
        err_root = ET.Element('EXPERIMENT_PACKAGE')
        err = ET.SubElement(err_root, 'Error')
        err.text = 'Entrez error'
        monkeypatch.setattr(
            'amalgkit.metadata.ET.parse',
            lambda handle: self._DummyTree(err_root)
        )

        with pytest.raises(Exception, match='Error\\. found in the xml'):
            fetch_sra_xml(search_term='SRR000000', retmax=1000)
