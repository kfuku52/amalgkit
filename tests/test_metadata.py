import xml.etree.ElementTree as ET

import pytest
import pandas
import urllib.error

from types import SimpleNamespace

from amalgkit.metadata import fetch_sra_xml, metadata_main
from amalgkit.util import Metadata


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

    def test_retries_esearch_once_on_urlerror(self, monkeypatch):
        esearch_calls = {'n': 0}

        def flaky_esearch(**_kwargs):
            esearch_calls['n'] += 1
            if esearch_calls['n'] == 1:
                raise urllib.error.URLError('temporary network failure')
            return object()

        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', flaky_esearch)
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': []})
        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', lambda **kwargs: object())

        root = fetch_sra_xml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE_SET'
        assert esearch_calls['n'] == 2

    def test_retries_efetch_once_on_urlerror(self, monkeypatch):
        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': ['ID1']})
        efetch_calls = {'n': 0}

        def flaky_efetch(**_kwargs):
            efetch_calls['n'] += 1
            if efetch_calls['n'] == 1:
                raise urllib.error.URLError('temporary network failure')
            return object()

        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', flaky_efetch)
        monkeypatch.setattr(
            'amalgkit.metadata.ET.parse',
            lambda handle: self._DummyTree(ET.Element('EXPERIMENT_PACKAGE'))
        )
        monkeypatch.setattr('amalgkit.metadata.time.sleep', lambda *_args, **_kwargs: None)

        root = fetch_sra_xml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE'
        assert efetch_calls['n'] == 2

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

    def test_wraps_non_set_chunks_to_preserve_multiple_records(self, monkeypatch):
        id_list = ['ID{}'.format(i) for i in range(2000)]
        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': id_list})
        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', lambda **kwargs: object())

        def fake_parse(_handle):
            root = ET.Element('EXPERIMENT_PACKAGE')
            ET.SubElement(root, 'RUN_SET')
            return self._DummyTree(root)

        monkeypatch.setattr('amalgkit.metadata.ET.parse', fake_parse)

        root = fetch_sra_xml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE_SET'
        assert len(root.findall('./EXPERIMENT_PACKAGE')) == 2

    def test_raises_when_xml_chunk_parsing_never_succeeds(self, monkeypatch):
        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.ET.parse', lambda handle: (_ for _ in ()).throw(ET.ParseError('broken xml')))

        with pytest.raises(RuntimeError, match='Failed to parse Entrez XML chunk'):
            fetch_sra_xml(search_term='SRR000000', retmax=1000)

    def test_propagates_unexpected_xml_parse_error(self, monkeypatch):
        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.ET.parse', lambda handle: (_ for _ in ()).throw(ValueError('unexpected')))

        with pytest.raises(ValueError, match='unexpected'):
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

        with pytest.raises(RuntimeError, match='Error found in Entrez XML response'):
            fetch_sra_xml(search_term='SRR000000', retmax=1000)


class TestMetadataMain:
    @staticmethod
    def _args(out_dir):
        return SimpleNamespace(
            entrez_email='dummy@example.com',
            out_dir=str(out_dir),
            redo=False,
            search_string='dummy',
            resolve_names=False,
        )

    def test_rejects_out_dir_file_path(self, tmp_path):
        out_file = tmp_path / 'out'
        out_file.write_text('not a directory')
        args = self._args(out_file)

        with pytest.raises(NotADirectoryError, match='not a directory'):
            metadata_main(args)

    def test_rejects_metadata_dir_file_path(self, tmp_path):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        metadata_path = out_dir / 'metadata'
        metadata_path.write_text('not a directory')
        args = self._args(out_dir)

        with pytest.raises(NotADirectoryError, match='not a directory'):
            metadata_main(args)

    def test_rejects_metadata_tsv_directory_when_redo_disabled(self, tmp_path):
        out_dir = tmp_path / 'out'
        metadata_dir = out_dir / 'metadata'
        metadata_dir.mkdir(parents=True)
        (metadata_dir / 'metadata.tsv').mkdir()
        args = self._args(out_dir)
        args.redo = False

        with pytest.raises(NotADirectoryError, match='not a file'):
            metadata_main(args)

    def test_rejects_metadata_tsv_directory_when_redo_enabled(self, tmp_path):
        out_dir = tmp_path / 'out'
        metadata_dir = out_dir / 'metadata'
        metadata_dir.mkdir(parents=True)
        (metadata_dir / 'metadata.tsv').mkdir()
        args = self._args(out_dir)
        args.redo = True

        with pytest.raises(NotADirectoryError, match='not a file'):
            metadata_main(args)

    def test_empty_result_guidance_mentions_search_string(self, tmp_path, monkeypatch, capsys):
        args = self._args(tmp_path / 'out')
        empty_metadata = Metadata.from_DataFrame(pandas.DataFrame(columns=Metadata.column_names))

        monkeypatch.setattr('amalgkit.metadata.fetch_sra_xml', lambda search_term: ET.Element('EXPERIMENT_PACKAGE_SET'))
        monkeypatch.setattr('amalgkit.metadata.Metadata.from_xml', lambda xml_root: empty_metadata)
        monkeypatch.setattr('amalgkit.metadata.Metadata.add_standard_rank_taxids', lambda self: None)

        metadata_main(args)

        captured = capsys.readouterr()
        assert '--search_string' in captured.err

    def test_rejects_missing_search_string(self, tmp_path, monkeypatch):
        args = self._args(tmp_path / 'out')
        args.search_string = None
        monkeypatch.setattr(
            'amalgkit.metadata.fetch_sra_xml',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('fetch_sra_xml should not be called')),
        )

        with pytest.raises(ValueError, match='--search_string is required'):
            metadata_main(args)

    def test_rejects_blank_search_string(self, tmp_path, monkeypatch):
        args = self._args(tmp_path / 'out')
        args.search_string = '   '
        monkeypatch.setattr(
            'amalgkit.metadata.fetch_sra_xml',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('fetch_sra_xml should not be called')),
        )

        with pytest.raises(ValueError, match='--search_string is required'):
            metadata_main(args)
