import xml.etree.ElementTree as ET
import json
from contextlib import contextmanager

import pytest
import pandas
import urllib.error
from http.client import IncompleteRead

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

    def test_retries_chunk_parse_once_on_incomplete_read(self, monkeypatch):
        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', lambda **kwargs: object())
        parse_calls = {'n': 0}

        def flaky_parse(_handle):
            parse_calls['n'] += 1
            if parse_calls['n'] == 1:
                raise IncompleteRead(b'partial-xml', len(b'partial-xml') + 10)
            return self._DummyTree(ET.Element('EXPERIMENT_PACKAGE'))

        monkeypatch.setattr('amalgkit.metadata.ET.parse', flaky_parse)
        monkeypatch.setattr('amalgkit.metadata.time.sleep', lambda *_args, **_kwargs: None)

        root = fetch_sra_xml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE'
        assert parse_calls['n'] == 2

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

    def test_acquires_ncbi_metadata_semaphore_for_esearch_and_efetch(self, monkeypatch, tmp_path):
        observed = []

        @contextmanager
        def fake_maybe_acquire(args, limit_attr, semaphore_name, lock_label, resolve_download_dir_fn=None):
            observed.append((limit_attr, semaphore_name, lock_label, args.ncbi_metadata_max_concurrency))
            yield 'slot'

        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir='inferred',
            download_lock_dir='inferred',
            ncbi_metadata_max_concurrency=2,
        )

        monkeypatch.setattr('amalgkit.sra.maybe_acquire_download_semaphore', fake_maybe_acquire)
        monkeypatch.setattr('amalgkit.metadata.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.metadata.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.metadata.Entrez.efetch', lambda **kwargs: object())
        monkeypatch.setattr(
            'amalgkit.metadata.ET.parse',
            lambda handle: self._DummyTree(ET.Element('EXPERIMENT_PACKAGE'))
        )

        root = fetch_sra_xml(search_term='SRR000000', retmax=1000, args=args)

        assert root.tag == 'EXPERIMENT_PACKAGE'
        assert observed == [
            ('ncbi_metadata_max_concurrency', 'ncbi_metadata', 'NCBI metadata download', 2),
            ('ncbi_metadata_max_concurrency', 'ncbi_metadata', 'NCBI metadata download', 2),
        ]


class TestMetadataMain:
    @staticmethod
    def _args(out_dir):
        return SimpleNamespace(
            entrez_email='dummy@example.com',
            out_dir=str(out_dir),
            redo=False,
            search_string='dummy',
            species_tsv=None,
            mode='base',
            title_terms='flower,leaf,root',
            organ_terms_tsv=None,
            species_limit=None,
            merge=True,
            resolve_names=False,
            threads='auto',
            internal_jobs='auto',
            internal_cpu_budget='auto',
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

        monkeypatch.setattr('amalgkit.metadata.search_sra_record_ids', lambda search_term: [])
        monkeypatch.setattr('amalgkit.metadata.Metadata.from_xml_roots', lambda xml_roots: empty_metadata)
        monkeypatch.setattr('amalgkit.metadata.Metadata.add_standard_rank_taxids', lambda self, args=None: None)

        metadata_main(args)

        captured = capsys.readouterr()
        assert '--search_string' in captured.err

    def test_empty_result_guidance_reports_accession_filter_mismatch(self, tmp_path, monkeypatch, capsys):
        args = self._args(tmp_path / 'out')
        args.search_string = '(SRR29779594) AND "Illumina"[Platform] AND ("RNA-seq"[Strategy] OR "EST"[Strategy])'
        empty_metadata = Metadata.from_DataFrame(pandas.DataFrame(columns=Metadata.column_names))

        monkeypatch.setattr('amalgkit.metadata.search_sra_record_ids', lambda search_term: [])
        monkeypatch.setattr('amalgkit.metadata.Metadata.from_xml_roots', lambda xml_roots: empty_metadata)
        monkeypatch.setattr('amalgkit.metadata.Metadata.add_standard_rank_taxids', lambda self, args=None: None)
        monkeypatch.setattr(
            'amalgkit.metadata.inspect_accession_search_mismatches',
            lambda search_term: [
                {
                    'accession': 'SRR29779594',
                    'scientific_name': 'Acanthamoeba astronyxis',
                    'platform': 'CAPILLARY (AB 3730xL Genetic Analyzer)',
                    'library_strategy': 'CLONE',
                    'library_source': 'TRANSCRIPTOMIC',
                    'library_selection': 'cDNA',
                }
            ],
        )

        metadata_main(args)

        captured = capsys.readouterr()
        assert 'SRR29779594' in captured.err
        assert 'platform=CAPILLARY (AB 3730xL Genetic Analyzer)' in captured.err
        assert 'library_strategy=CLONE' in captured.err

    def test_redo_keeps_existing_metadata_when_generation_fails(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        metadata_dir = out_dir / 'metadata'
        metadata_dir.mkdir(parents=True)
        metadata_path = metadata_dir / 'metadata.tsv'
        metadata_path.write_text('old\tvalue\n', encoding='utf-8')
        args = self._args(out_dir)
        args.redo = True

        monkeypatch.setattr('amalgkit.metadata.search_sra_record_ids', lambda _search_term: ['ID1'])
        monkeypatch.setattr(
            'amalgkit.metadata.Metadata.from_xml_roots',
            lambda _xml_roots: (_ for _ in ()).throw(RuntimeError('generation failed')),
        )

        with pytest.raises(RuntimeError, match='generation failed'):
            metadata_main(args)

        assert metadata_path.read_text(encoding='utf-8') == 'old\tvalue\n'

    def test_rejects_missing_search_string(self, tmp_path, monkeypatch):
        args = self._args(tmp_path / 'out')
        args.search_string = None
        monkeypatch.setattr(
            'amalgkit.metadata.search_sra_record_ids',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('search_sra_record_ids should not be called')),
        )

        with pytest.raises(ValueError, match='--search_string is required'):
            metadata_main(args)

    def test_rejects_blank_search_string(self, tmp_path, monkeypatch):
        args = self._args(tmp_path / 'out')
        args.search_string = '   '
        monkeypatch.setattr(
            'amalgkit.metadata.search_sra_record_ids',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('search_sra_record_ids should not be called')),
        )

        with pytest.raises(ValueError, match='--search_string is required'):
            metadata_main(args)

    def test_writes_tsv_safe_metadata_when_text_contains_carriage_returns(self, tmp_path, monkeypatch):
        args = self._args(tmp_path / 'out')
        metadata = Metadata.from_DataFrame(
            pandas.DataFrame(
                [
                    {
                        'scientific_name': 'Musa acuminata AAA Group',
                        'sample_group': 'root',
                        'tissue': 'root',
                        'biosample': 'SAMN000001',
                        'experiment': 'SRX000001',
                        'run': 'SRR000001',
                        'sample_title': 'banana root sample',
                        'sample_description': 'Line 1\rLine 2\nLine 3\tTabbed',
                        'taxid': '4641',
                    }
                ],
                columns=Metadata.column_names,
            )
        )
        monkeypatch.setattr('amalgkit.metadata.search_sra_record_ids', lambda _search_term: ['ID1'])
        monkeypatch.setattr('amalgkit.metadata.Metadata.from_xml_roots', lambda _xml_roots: metadata)
        monkeypatch.setattr('amalgkit.metadata.Metadata.add_standard_rank_taxids', lambda self, args=None: None)

        metadata_main(args)

        metadata_path = tmp_path / 'out' / 'metadata' / 'metadata.tsv'
        raw = metadata_path.read_bytes()
        assert b'\r' not in raw

        df = pandas.read_csv(metadata_path, sep='\t')
        assert df.shape[0] == 1
        assert df.loc[0, 'run'] == 'SRR000001'
        assert df.loc[0, 'sample_description'] == 'Line 1 Line 2 Line 3 Tabbed'

    def test_writes_query_info_and_summary_sidecars(self, tmp_path, monkeypatch):
        args = self._args(tmp_path / 'out')
        metadata = Metadata.from_DataFrame(
            pandas.DataFrame(
                [
                    {
                        'scientific_name': 'Arabidopsis thaliana',
                        'sample_group': 'leaf',
                        'tissue': 'leaf',
                        'bioproject': 'PRJNA1',
                        'biosample': 'SAMN1',
                        'experiment': 'SRX1',
                        'run': 'SRR1',
                        'sample_title': 'leaf sample',
                        'source_name': 'rosette leaf',
                        'taxid': '3702',
                    }
                ],
                columns=Metadata.column_names,
            )
        )
        monkeypatch.setattr('amalgkit.metadata.search_sra_record_ids', lambda _search_term: ['ID1'])
        monkeypatch.setattr('amalgkit.metadata.Metadata.from_xml_roots', lambda _xml_roots: metadata)
        monkeypatch.setattr('amalgkit.metadata.Metadata.add_standard_rank_taxids', lambda self, args=None: None)

        metadata_main(args)

        query_info_path = tmp_path / 'out' / 'metadata' / 'metadata.query_info.json'
        summary_path = tmp_path / 'out' / 'metadata' / 'metadata.summary.tsv'
        assert query_info_path.exists()
        assert summary_path.exists()
        query_info = json.loads(query_info_path.read_text(encoding='utf-8'))
        assert query_info['search_string'] == 'dummy'
        assert query_info['record_id_count'] == 1
        assert query_info['metadata_row_count'] == 1

        summary = pandas.read_csv(summary_path, sep='\t')
        row_count = summary.loc[summary['name'] == 'row_count', 'value'].iloc[0]
        assert row_count == '1'

    def test_drops_missing_run_id_rows_in_generated_metadata(self, tmp_path, monkeypatch):
        args = self._args(tmp_path / 'out')
        metadata = Metadata.from_DataFrame(
            pandas.DataFrame(
                [
                    {
                        'scientific_name': 'Arabidopsis thaliana',
                        'sample_group': 'leaf',
                        'tissue': 'leaf',
                        'experiment': 'SRX1',
                        'run': 'SRR1',
                        'sample_title': 'leaf sample',
                        'taxid': '3702',
                    },
                    {
                        'scientific_name': 'Arabidopsis thaliana',
                        'sample_group': 'leaf',
                        'tissue': 'leaf',
                        'experiment': 'SRX2',
                        'run': '',
                        'sample_title': 'broken sample',
                        'taxid': '3702',
                    }
                ],
                columns=Metadata.column_names,
            )
        )
        monkeypatch.setattr('amalgkit.metadata.search_sra_record_ids', lambda _search_term: ['ID1'])
        monkeypatch.setattr('amalgkit.metadata.Metadata.from_xml_roots', lambda _xml_roots: metadata)
        monkeypatch.setattr('amalgkit.metadata.Metadata.add_standard_rank_taxids', lambda self, args=None: None)

        metadata_main(args)

        metadata_path = tmp_path / 'out' / 'metadata' / 'metadata.tsv'
        query_info_path = tmp_path / 'out' / 'metadata' / 'metadata.query_info.json'

        df = pandas.read_csv(metadata_path, sep='\t')
        assert df.shape[0] == 1
        assert df.loc[0, 'run'] == 'SRR1'

        query_info = json.loads(query_info_path.read_text(encoding='utf-8'))
        assert query_info['missing_run_drop_count'] == 1
        assert query_info['missing_run_drop_examples'][0]['sample_title'] == 'broken sample'

    def test_rejects_duplicate_run_id_in_generated_metadata(self, tmp_path, monkeypatch):
        args = self._args(tmp_path / 'out')
        metadata = Metadata.from_DataFrame(
            pandas.DataFrame(
                [
                    {
                        'scientific_name': 'Arabidopsis thaliana',
                        'sample_group': 'leaf',
                        'tissue': 'leaf',
                        'experiment': 'SRX1',
                        'run': 'SRR1',
                        'taxid': '3702',
                    },
                    {
                        'scientific_name': 'Arabidopsis thaliana',
                        'sample_group': 'root',
                        'tissue': 'root',
                        'experiment': 'SRX2',
                        'run': 'SRR1',
                        'taxid': '3702',
                    },
                ],
                columns=Metadata.column_names,
            )
        )
        monkeypatch.setattr('amalgkit.metadata.search_sra_record_ids', lambda _search_term: ['ID1'])
        monkeypatch.setattr('amalgkit.metadata.Metadata.from_xml_roots', lambda _xml_roots: metadata)
        monkeypatch.setattr('amalgkit.metadata.Metadata.add_standard_rank_taxids', lambda self, args=None: None)

        with pytest.raises(ValueError, match='Duplicate run ID'):
            metadata_main(args)

    def test_species_batch_mode_writes_query_and_merged_outputs(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        args = self._args(out_dir)
        args.search_string = None
        args.species_tsv = str(tmp_path / 'species.tsv')
        pandas.DataFrame(
            [
                {'scientific_name': 'Species alpha'},
                {'scientific_name': 'Species beta'},
            ]
        ).to_csv(args.species_tsv, sep='\t', index=False)

        monkeypatch.setattr('amalgkit.metadata.search_sra_record_ids', lambda search_term: [search_term])
        monkeypatch.setattr(
            'amalgkit.metadata.iter_sra_xml_chunks',
            lambda record_ids, retmax=1000: iter([{'search_term': record_ids[0]}]),
        )

        def fake_from_xml_roots(xml_roots):
            payload = list(xml_roots)[0]
            search_term = payload['search_term']
            species_name = 'Species alpha' if 'Species alpha' in search_term else 'Species beta'
            run_id = 'SRR_ALPHA' if species_name == 'Species alpha' else 'SRR_BETA'
            taxid = '1001' if species_name == 'Species alpha' else '1002'
            return Metadata.from_DataFrame(
                pandas.DataFrame(
                    [
                        {
                            'scientific_name': species_name,
                            'sample_group': 'leaf',
                            'tissue': 'leaf',
                            'bioproject': 'PRJNA_' + species_name.split()[-1].upper(),
                            'biosample': 'SAMN_' + species_name.split()[-1].upper(),
                            'experiment': 'SRX_' + species_name.split()[-1].upper(),
                            'run': run_id,
                            'sample_title': species_name + ' leaf sample',
                            'taxid': taxid,
                        }
                    ],
                    columns=Metadata.column_names,
                )
            )

        monkeypatch.setattr('amalgkit.metadata.Metadata.from_xml_roots', fake_from_xml_roots)
        monkeypatch.setattr('amalgkit.metadata.Metadata.add_standard_rank_taxids', lambda self, args=None: None)

        metadata_main(args)

        alpha_query = out_dir / 'metadata_specieswise' / 'Species_alpha' / 'base' / 'metadata' / 'metadata.tsv'
        alpha_merged = out_dir / 'metadata_specieswise' / 'Species_alpha' / 'Species_alpha.metadata.tsv'
        alpha_merged_info = out_dir / 'metadata_specieswise' / 'Species_alpha' / 'Species_alpha.query_info.json'
        beta_query = out_dir / 'metadata_specieswise' / 'Species_beta' / 'base' / 'metadata' / 'metadata.tsv'
        assert alpha_query.exists()
        assert alpha_merged.exists()
        assert alpha_merged_info.exists()
        assert beta_query.exists()
        merged_info = json.loads(alpha_merged_info.read_text(encoding='utf-8'))
        assert merged_info['kind'] == 'merged'
        assert merged_info['source_query_labels'] == ['base']

    def test_species_batch_mode_uses_parallel_workers(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        args = self._args(out_dir)
        args.search_string = None
        args.species_tsv = str(tmp_path / 'species.tsv')
        args.internal_jobs = 2
        args.threads = 2
        pandas.DataFrame(
            [
                {'scientific_name': 'Species alpha'},
                {'scientific_name': 'Species beta'},
            ]
        ).to_csv(args.species_tsv, sep='\t', index=False)

        observed = {
            'max_workers': None,
            'species': [],
            'merged': [],
        }

        def fake_run_single_query(args, out_dir=None, search_string=None, species_name=None, query_label=None, mode='single', allow_cached=False):
            _ = (args, out_dir, search_string, query_label, mode, allow_cached)
            observed['species'].append(species_name)
            token = species_name.split()[-1].upper()
            metadata = Metadata.from_DataFrame(
                pandas.DataFrame(
                    [
                        {
                            'scientific_name': species_name,
                            'sample_group': 'leaf',
                            'tissue': 'leaf',
                            'bioproject': 'PRJNA_' + token,
                            'biosample': 'SAMN_' + token,
                            'experiment': 'SRX_' + token,
                            'run': 'SRR_' + token,
                            'sample_title': species_name + ' leaf sample',
                            'taxid': '1001',
                        }
                    ],
                    columns=Metadata.column_names,
                )
            )
            return {
                'metadata': metadata,
                'paths': {
                    'query_info_path': str(tmp_path / (species_name.replace(' ', '_') + '.query_info.json')),
                    'metadata_path': str(tmp_path / (species_name.replace(' ', '_') + '.metadata.tsv')),
                },
                'query_info': {
                    'record_id_count': 1,
                    'missing_run_drop_count': 0,
                },
                'query_label': 'base',
            }

        def fake_write_species_merged_metadata(args, species_name, species_dir, species_token, query_results):
            _ = (args, species_dir, species_token, query_results)
            observed['merged'].append(species_name)

        def fake_run_tasks(task_items, task_fn, max_workers=1):
            observed['max_workers'] = max_workers
            results = {}
            failures = []
            for task in task_items:
                results[task] = task_fn(task)
            return results, failures

        monkeypatch.setattr('amalgkit.metadata._run_single_query', fake_run_single_query)
        monkeypatch.setattr('amalgkit.metadata._write_species_merged_metadata', fake_write_species_merged_metadata)
        monkeypatch.setattr('amalgkit.metadata.run_tasks_with_optional_threads', fake_run_tasks)

        metadata_main(args)

        assert observed['max_workers'] == 2
        assert set(observed['species']) == {'Species alpha', 'Species beta'}
        assert set(observed['merged']) == {'Species alpha', 'Species beta'}

    def test_species_batch_mode_cpu_budget_caps_parallel_workers_to_serial(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        args = self._args(out_dir)
        args.search_string = None
        args.species_tsv = str(tmp_path / 'species.tsv')
        args.internal_jobs = 4
        args.threads = 4
        args.internal_cpu_budget = 1
        pandas.DataFrame(
            [
                {'scientific_name': 'Species alpha'},
                {'scientific_name': 'Species beta'},
            ]
        ).to_csv(args.species_tsv, sep='\t', index=False)

        observed = {
            'species': [],
            'merged': [],
        }

        def fake_run_single_query(args, out_dir=None, search_string=None, species_name=None, query_label=None, mode='single', allow_cached=False):
            _ = (args, out_dir, search_string, query_label, mode, allow_cached)
            observed['species'].append(species_name)
            metadata = Metadata.from_DataFrame(
                pandas.DataFrame(
                    [
                        {
                            'scientific_name': species_name,
                            'sample_group': 'leaf',
                            'tissue': 'leaf',
                            'experiment': 'SRX1',
                            'run': 'SRR1_' + species_name.split()[-1].upper(),
                            'taxid': '1001',
                        }
                    ],
                    columns=Metadata.column_names,
                )
            )
            return {
                'metadata': metadata,
                'paths': {
                    'query_info_path': str(tmp_path / (species_name.replace(' ', '_') + '.query_info.json')),
                    'metadata_path': str(tmp_path / (species_name.replace(' ', '_') + '.metadata.tsv')),
                },
                'query_info': {
                    'record_id_count': 1,
                    'missing_run_drop_count': 0,
                },
                'query_label': 'base',
            }

        def fake_write_species_merged_metadata(args, species_name, species_dir, species_token, query_results):
            _ = (args, species_dir, species_token, query_results)
            observed['merged'].append(species_name)

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('run_tasks_with_optional_threads should not be used when --internal_cpu_budget caps internal_jobs to 1.')

        monkeypatch.setattr('amalgkit.metadata._run_single_query', fake_run_single_query)
        monkeypatch.setattr('amalgkit.metadata._write_species_merged_metadata', fake_write_species_merged_metadata)
        monkeypatch.setattr('amalgkit.metadata.run_tasks_with_optional_threads', fail_if_called)

        metadata_main(args)

        assert observed['species'] == ['Species alpha', 'Species beta']
        assert observed['merged'] == ['Species alpha', 'Species beta']


class TestMetadataFromXmlRoots:
    def test_preserves_colliding_sample_attributes(self):
        xml_root = ET.fromstring(
            '''
            <EXPERIMENT_PACKAGE_SET>
              <EXPERIMENT_PACKAGE>
                <SUBMISSION center_name="center">
                  <IDENTIFIERS>
                    <PRIMARY_ID>SRA000001</PRIMARY_ID>
                    <SUBMITTER_ID>submitter</SUBMITTER_ID>
                  </IDENTIFIERS>
                </SUBMISSION>
                <STUDY>
                  <DESCRIPTOR>
                    <STUDY_TITLE>study title</STUDY_TITLE>
                  </DESCRIPTOR>
                </STUDY>
                <SAMPLE>
                  <IDENTIFIERS>
                    <PRIMARY_ID>SRS000001</PRIMARY_ID>
                  </IDENTIFIERS>
                  <TITLE>sample title</TITLE>
                  <DESCRIPTION>sample description</DESCRIPTION>
                  <SAMPLE_NAME>
                    <SCIENTIFIC_NAME>Arabidopsis thaliana</SCIENTIFIC_NAME>
                    <TAXON_ID>3702</TAXON_ID>
                  </SAMPLE_NAME>
                  <SAMPLE_ATTRIBUTES>
                    <SAMPLE_ATTRIBUTE>
                      <TAG>tissue</TAG>
                      <VALUE>leaf</VALUE>
                    </SAMPLE_ATTRIBUTE>
                    <SAMPLE_ATTRIBUTE>
                      <TAG>leaf</TAG>
                      <VALUE>foliage</VALUE>
                    </SAMPLE_ATTRIBUTE>
                    <SAMPLE_ATTRIBUTE>
                      <TAG>leaf</TAG>
                      <VALUE>leaf blade</VALUE>
                    </SAMPLE_ATTRIBUTE>
                  </SAMPLE_ATTRIBUTES>
                </SAMPLE>
                <EXPERIMENT>
                  <IDENTIFIERS>
                    <PRIMARY_ID>SRX000001</PRIMARY_ID>
                  </IDENTIFIERS>
                  <TITLE>experiment title</TITLE>
                  <STUDY_REF>
                    <IDENTIFIERS>
                      <PRIMARY_ID>SRP000001</PRIMARY_ID>
                    </IDENTIFIERS>
                  </STUDY_REF>
                  <DESIGN>
                    <DESIGN_DESCRIPTION>design</DESIGN_DESCRIPTION>
                    <LIBRARY_DESCRIPTOR>
                      <LIBRARY_NAME>lib</LIBRARY_NAME>
                      <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
                      <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
                      <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
                      <LIBRARY_LAYOUT>
                        <SINGLE />
                      </LIBRARY_LAYOUT>
                    </LIBRARY_DESCRIPTOR>
                  </DESIGN>
                  <PLATFORM>
                    <ILLUMINA>
                      <INSTRUMENT_MODEL>NovaSeq</INSTRUMENT_MODEL>
                    </ILLUMINA>
                  </PLATFORM>
                </EXPERIMENT>
                <RUN_SET>
                  <RUN published="2024-01-01" total_spots="100" total_bases="10000" size="1000">
                    <IDENTIFIERS>
                      <PRIMARY_ID>SRR000001</PRIMARY_ID>
                    </IDENTIFIERS>
                  </RUN>
                </RUN_SET>
              </EXPERIMENT_PACKAGE>
            </EXPERIMENT_PACKAGE_SET>
            '''
        )

        metadata = Metadata.from_xml_roots([xml_root])

        assert metadata.df.loc[0, 'sample_attribute_tissue'] == 'leaf'
        assert metadata.df.loc[0, 'leaf'] == 'foliage | leaf blade'
        assert metadata.sample_attribute_collision_count >= 2

    def test_from_xml_roots_extracts_generic_platform_and_instrument(self):
        xml_root = ET.fromstring(
            '''
            <EXPERIMENT_PACKAGE_SET>
              <EXPERIMENT_PACKAGE>
                <SAMPLE>
                  <SAMPLE_NAME>
                    <SCIENTIFIC_NAME>Arabidopsis thaliana</SCIENTIFIC_NAME>
                    <TAXON_ID>3702</TAXON_ID>
                  </SAMPLE_NAME>
                </SAMPLE>
                <EXPERIMENT>
                  <IDENTIFIERS>
                    <PRIMARY_ID>SRX000002</PRIMARY_ID>
                  </IDENTIFIERS>
                  <DESIGN>
                    <LIBRARY_DESCRIPTOR>
                      <LIBRARY_LAYOUT>
                        <SINGLE />
                      </LIBRARY_LAYOUT>
                    </LIBRARY_DESCRIPTOR>
                  </DESIGN>
                  <PLATFORM>
                    <PACBIO_SMRT>
                      <INSTRUMENT_MODEL>Sequel II</INSTRUMENT_MODEL>
                    </PACBIO_SMRT>
                  </PLATFORM>
                </EXPERIMENT>
                <RUN_SET>
                  <RUN total_spots="10" total_bases="10000" size="1000">
                    <IDENTIFIERS>
                      <PRIMARY_ID>SRR000002</PRIMARY_ID>
                    </IDENTIFIERS>
                  </RUN>
                </RUN_SET>
              </EXPERIMENT_PACKAGE>
            </EXPERIMENT_PACKAGE_SET>
            '''
        )

        metadata = Metadata.from_xml_roots([xml_root])

        assert metadata.df.loc[0, 'platform'] == 'PACBIO_SMRT'
        assert metadata.df.loc[0, 'instrument'] == 'Sequel II'
