import os
import threading
import time

import pytest
import pandas
from types import SimpleNamespace

from amalgkit.util import (
    Metadata,
    get_ete_ncbitaxa,
    run_tasks_with_optional_threads,
)


# ---------------------------------------------------------------------------
# strtobool
# ---------------------------------------------------------------------------


class TestMetadataTaxidValidation:
    def test_add_standard_rank_taxids_requires_nullable_int64_taxid(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
            'taxid': ['9606'],
        }))

        def fail_if_called():
            raise AssertionError('NCBITaxa should not be initialized for invalid taxid dtype.')

        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', fail_if_called)

        with pytest.raises(TypeError, match='taxid column must be Int64 dtype'):
            metadata.add_standard_rank_taxids()

    def test_resolve_scientific_names_requires_nullable_int64_taxid(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['H. sapiens'],
            'exclusion': ['no'],
            'taxid': ['9606'],
        }))

        def fail_if_called():
            raise AssertionError('NCBITaxa should not be initialized for invalid taxid dtype.')

        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', fail_if_called)

        with pytest.raises(TypeError, match='taxid column must be Int64 dtype'):
            metadata.resolve_scientific_names()

    def test_add_standard_rank_taxids_does_not_swallow_keyboard_interrupt(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
            'taxid': [9606],
        }))
        metadata.df['taxid'] = metadata.df['taxid'].astype('Int64')

        class InterruptingNcbi:
            def get_lineage(self, _taxid):
                raise KeyboardInterrupt()

            def get_rank(self, _lineage):
                return {}

        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', lambda: InterruptingNcbi())

        with pytest.raises(KeyboardInterrupt):
            metadata.add_standard_rank_taxids()

    def test_add_standard_rank_taxids_warns_on_lineage_lookup_error(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
            'taxid': [9606],
        }))
        metadata.df['taxid'] = metadata.df['taxid'].astype('Int64')

        class FailingNcbi:
            def get_lineage(self, _taxid):
                raise RuntimeError('boom')

            def get_rank(self, _lineage):
                return {}

        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', lambda: FailingNcbi())

        with pytest.warns(UserWarning, match='Failed to resolve NCBI lineage'):
            metadata.add_standard_rank_taxids()
        assert 'taxid_species' in metadata.df.columns
        assert metadata.df['taxid_species'].isna().all()

    def test_add_standard_rank_taxids_batches_lineage_lookup_when_supported(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Homo sapiens', 'Mus musculus'],
            'exclusion': ['no', 'no'],
            'taxid': [9606, 10090],
        }))
        metadata.df['taxid'] = metadata.df['taxid'].astype('Int64')
        captured = {'lineage_translator_calls': 0, 'rank_calls': 0}

        class BatchNcbi:
            def get_lineage_translator(self, taxids):
                captured['lineage_translator_calls'] += 1
                assert sorted(taxids) == [9606, 10090]
                return {
                    9606: [1, 2759, 9605, 9606],
                    10090: [1, 2759, 10088, 10090],
                }

            def get_lineage(self, _taxid):
                raise AssertionError('Per-taxid lineage lookup should not be used when batch API is available.')

            def get_rank(self, taxids):
                captured['rank_calls'] += 1
                assert set(taxids) == {1, 2759, 9605, 9606, 10088, 10090}
                return {
                    1: 'domain',
                    2759: 'kingdom',
                    9605: 'genus',
                    9606: 'species',
                    10088: 'genus',
                    10090: 'species',
                }

        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', lambda: BatchNcbi())

        metadata.add_standard_rank_taxids()

        assert captured['lineage_translator_calls'] == 1
        assert captured['rank_calls'] == 1
        assert metadata.df.loc[0, 'taxid_species'] == 9606
        assert metadata.df.loc[1, 'taxid_species'] == 10090
        assert metadata.df.loc[0, 'taxid_genus'] == 9605
        assert metadata.df.loc[1, 'taxid_genus'] == 10088

    def test_add_standard_rank_taxids_reuses_cached_lineage_and_rank_results(self, monkeypatch):
        metadata1 = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
            'taxid': [9606],
        }))
        metadata2 = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR002'],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
            'taxid': [9606],
        }))
        metadata1.df['taxid'] = metadata1.df['taxid'].astype('Int64')
        metadata2.df['taxid'] = metadata2.df['taxid'].astype('Int64')
        captured = {'lineage_translator_calls': 0, 'rank_calls': 0}

        class BatchNcbi:
            def get_lineage_translator(self, taxids):
                captured['lineage_translator_calls'] += 1
                assert list(taxids) == [9606]
                return {9606: [1, 2759, 9605, 9606]}

            def get_lineage(self, _taxid):
                raise AssertionError('Per-taxid lineage lookup should not be needed.')

            def get_rank(self, taxids):
                captured['rank_calls'] += 1
                assert set(taxids) == {1, 2759, 9605, 9606}
                return {
                    1: 'domain',
                    2759: 'kingdom',
                    9605: 'genus',
                    9606: 'species',
                }

        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', lambda: BatchNcbi())

        metadata1.add_standard_rank_taxids()
        metadata2.add_standard_rank_taxids()

        assert captured['lineage_translator_calls'] == 1
        assert captured['rank_calls'] == 1

    def test_add_standard_rank_taxids_uses_download_dir_for_shared_ete_taxonomy_cache(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
            'taxid': [9606],
        }))
        metadata.df['taxid'] = metadata.df['taxid'].astype('Int64')
        captured = {'lock_path': None, 'downloaded': None}

        class DummyLock:
            def __init__(self, lock_path, lock_label='Lock', poll_seconds=5, timeout_seconds=3600):
                _ = (lock_label, poll_seconds, timeout_seconds)
                captured['lock_path'] = lock_path

            def __enter__(self):
                return None

            def __exit__(self, exc_type, exc, tb):
                return False

        class RecordingNcbi:
            def __init__(self, **kwargs):
                captured['kwargs'] = kwargs

            def get_lineage(self, _taxid):
                return [1, 9606]

            def get_rank(self, _lineage):
                return {1: 'domain', 9606: 'species'}

        def fake_download(url, output_path, timeout_seconds, urlopen_fn=None):
            _ = (timeout_seconds, urlopen_fn)
            captured['downloaded'] = (url, output_path)
            with open(output_path, 'wb') as fout:
                fout.write(b'taxdump')

        monkeypatch.setattr('amalgkit.download_utils.acquire_exclusive_lock', DummyLock)
        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', RecordingNcbi)
        monkeypatch.setattr('amalgkit.download_utils.download_url_to_regular_file', fake_download)
        monkeypatch.setattr('amalgkit.download_utils.validate_taxonomy_dump', lambda path: path)
        monkeypatch.setattr('amalgkit.download_utils.read_published_md5', lambda *args, **kwargs: 'd41d8cd98f00b204e9800998ecf8427e')
        monkeypatch.setattr('amalgkit.download_utils.calculate_file_md5', lambda path: 'd41d8cd98f00b204e9800998ecf8427e')
        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(tmp_path / 'shared_downloads'),
        )
        metadata.add_standard_rank_taxids(args=args)

        expected_download_dir = os.path.realpath(args.download_dir)
        expected_ete_dir = os.path.join(expected_download_dir, 'ete_taxonomy')
        expected_lock_path = os.path.join(expected_download_dir, 'locks', 'ete_taxonomy.lock')
        assert captured['kwargs']['dbfile'] == os.path.join(expected_ete_dir, 'taxa.sqlite')
        assert captured['kwargs']['taxdump_file'] == os.path.join(expected_ete_dir, 'taxdump.tar.gz')
        assert os.path.isdir(expected_ete_dir)
        assert os.path.isfile(os.path.join(expected_ete_dir, 'taxdump.tar.gz'))
        assert captured['lock_path'] == expected_lock_path
        # The hardened download path (timeout-aware, checksum-verified) is
        # used in production, not the untimed urllib.request.urlretrieve.
        assert captured['downloaded'][0].endswith('/taxdump.tar.gz')
        assert os.path.dirname(captured['downloaded'][1]) == expected_ete_dir
        assert os.path.basename(captured['downloaded'][1]).startswith('taxdump.tar.gz.')
        assert captured['downloaded'][1].endswith('.tmp')

    def test_add_standard_rank_taxids_replaces_existing_lineage_columns(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
            'taxid': [9606],
            'taxid_domain': [999],
            'taxid_kingdom': [999],
            'taxid_phylum': [pandas.NA],
            'taxid_class': [pandas.NA],
            'taxid_order': [pandas.NA],
            'taxid_family': [pandas.NA],
            'taxid_genus': [999],
            'taxid_species': [999],
        }))
        metadata.df['taxid'] = metadata.df['taxid'].astype('Int64')
        for col in [col for col in metadata.df.columns if col.startswith('taxid_') and col != 'taxid']:
            metadata.df[col] = pandas.to_numeric(metadata.df[col], errors='coerce').astype('Int64')

        class RecordingNcbi:
            def get_lineage(self, _taxid):
                return [1, 2759, 9605, 9606]

            def get_rank(self, taxids):
                assert set(taxids) == {1, 2759, 9605, 9606}
                return {
                    1: 'domain',
                    2759: 'kingdom',
                    9605: 'genus',
                    9606: 'species',
                }

        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', lambda: RecordingNcbi())

        metadata.add_standard_rank_taxids()

        assert 'taxid_domain_x' not in metadata.df.columns
        assert 'taxid_domain_y' not in metadata.df.columns
        assert metadata.df.loc[0, 'taxid_domain'] == 1
        assert metadata.df.loc[0, 'taxid_kingdom'] == 2759
        assert metadata.df.loc[0, 'taxid_genus'] == 9605
        assert metadata.df.loc[0, 'taxid_species'] == 9606

    def test_resolve_scientific_names_reuses_cached_translator_results(self, monkeypatch):
        metadata1 = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['H. sapiens'],
            'exclusion': ['no'],
            'taxid': [9606],
        }))
        metadata2 = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR002'],
            'scientific_name': ['human'],
            'exclusion': ['no'],
            'taxid': [9606],
        }))
        metadata1.df['taxid'] = metadata1.df['taxid'].astype('Int64')
        metadata2.df['taxid'] = metadata2.df['taxid'].astype('Int64')
        captured = {'translator_calls': 0}

        class RecordingNcbi:
            def get_taxid_translator(self, taxids):
                captured['translator_calls'] += 1
                assert list(taxids) == [9606]
                return {9606: 'Homo sapiens'}

        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', lambda: RecordingNcbi())

        metadata1.resolve_scientific_names()
        metadata2.resolve_scientific_names()

        assert captured['translator_calls'] == 1
        assert metadata1.df.loc[0, 'scientific_name'] == 'Homo sapiens'
        assert metadata2.df.loc[0, 'scientific_name'] == 'Homo sapiens'

    def test_get_ete_ncbitaxa_bootstraps_fresh_custom_download_dir(self, tmp_path, monkeypatch):
        captured = {'lock_path': None, 'downloaded': None}

        class DummyLock:
            def __init__(self, lock_path, lock_label='Lock', poll_seconds=5, timeout_seconds=3600):
                _ = (lock_label, poll_seconds, timeout_seconds)
                captured['lock_path'] = lock_path

            def __enter__(self):
                return None

            def __exit__(self, exc_type, exc, tb):
                return False

        class RecordingNcbi:
            def __init__(self, **kwargs):
                captured['kwargs'] = kwargs

        def fake_download(url, output_path, timeout_seconds, urlopen_fn=None):
            _ = (timeout_seconds, urlopen_fn)
            captured['downloaded'] = (url, output_path)
            with open(output_path, 'wb') as fout:
                fout.write(b'taxdump')

        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(tmp_path / 'fresh_downloads'),
        )
        expected_download_dir = os.path.realpath(args.download_dir)
        expected_ete_dir = os.path.join(expected_download_dir, 'ete_taxonomy')
        expected_lock_path = os.path.join(expected_download_dir, 'locks', 'ete_taxonomy.lock')

        assert not os.path.exists(args.download_dir)
        monkeypatch.setattr('amalgkit.util.acquire_exclusive_lock', DummyLock)
        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', RecordingNcbi)
        monkeypatch.setattr('amalgkit.download_utils.download_url_to_regular_file', fake_download)
        monkeypatch.setattr('amalgkit.download_utils.validate_taxonomy_dump', lambda path: path)
        monkeypatch.setattr('amalgkit.download_utils.read_published_md5', lambda *args, **kwargs: 'd41d8cd98f00b204e9800998ecf8427e')
        monkeypatch.setattr('amalgkit.download_utils.calculate_file_md5', lambda path: 'd41d8cd98f00b204e9800998ecf8427e')

        result = get_ete_ncbitaxa(args=args)

        assert isinstance(result, RecordingNcbi)
        assert os.path.isdir(expected_ete_dir)
        assert captured['kwargs']['dbfile'] == os.path.join(expected_ete_dir, 'taxa.sqlite')
        assert captured['kwargs']['taxdump_file'] == os.path.join(expected_ete_dir, 'taxdump.tar.gz')
        assert os.path.isfile(os.path.join(expected_ete_dir, 'taxdump.tar.gz'))
        assert captured['lock_path'] == expected_lock_path
        # The hardened download path (timeout-aware, checksum-verified) is
        # used in production, not the untimed urllib.request.urlretrieve.
        assert captured['downloaded'][0].endswith('/taxdump.tar.gz')
        assert os.path.dirname(captured['downloaded'][1]) == expected_ete_dir
        assert os.path.basename(captured['downloaded'][1]).startswith('taxdump.tar.gz.')
        assert captured['downloaded'][1].endswith('.tmp')

    def test_get_ete_ncbitaxa_reuses_existing_custom_db_without_taxdump_refresh(self, tmp_path, monkeypatch):
        captured = {'lock_path': None}

        class DummyLock:
            def __init__(self, lock_path, lock_label='Lock', poll_seconds=5, timeout_seconds=3600):
                _ = (lock_label, poll_seconds, timeout_seconds)
                captured['lock_path'] = lock_path

            def __enter__(self):
                return None

            def __exit__(self, exc_type, exc, tb):
                return False

        class RecordingNcbi:
            def __init__(self, **kwargs):
                captured['kwargs'] = kwargs

        def fail_urlretrieve(*_args, **_kwargs):
            raise AssertionError('urlretrieve should not be called when a usable custom ETE DB already exists.')

        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(tmp_path / 'shared_downloads'),
        )
        expected_download_dir = os.path.realpath(args.download_dir)
        expected_ete_dir = os.path.join(expected_download_dir, 'ete_taxonomy')
        expected_lock_path = os.path.join(expected_download_dir, 'locks', 'ete_taxonomy.lock')
        os.makedirs(expected_ete_dir, exist_ok=True)
        dbfile = os.path.join(expected_ete_dir, 'taxa.sqlite')
        with open(dbfile, 'wb') as fout:
            fout.write(b'sqlite')

        monkeypatch.setattr('amalgkit.util.acquire_exclusive_lock', DummyLock)
        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', RecordingNcbi)
        monkeypatch.setattr('amalgkit.download_utils.is_taxonomy_database_compatible', lambda path: path == dbfile)
        monkeypatch.setattr('amalgkit.download_utils.urllib.request.urlretrieve', fail_urlretrieve)

        result = get_ete_ncbitaxa(args=args)

        assert isinstance(result, RecordingNcbi)
        assert captured['kwargs']['dbfile'] == dbfile
        assert captured['kwargs']['update'] is False
        assert 'taxdump_file' not in captured['kwargs']
        assert captured['lock_path'] == expected_lock_path

    def test_get_ete_ncbitaxa_caches_instances_per_custom_db(self, tmp_path, monkeypatch):
        captured = {'init_calls': 0, 'lock_calls': 0}

        class DummyLock:
            def __init__(self, lock_path, lock_label='Lock', poll_seconds=5, timeout_seconds=3600):
                _ = (lock_path, lock_label, poll_seconds, timeout_seconds)

            def __enter__(self):
                captured['lock_calls'] += 1
                return None

            def __exit__(self, exc_type, exc, tb):
                return False

        class RecordingNcbi:
            def __init__(self, **kwargs):
                captured['init_calls'] += 1
                captured.setdefault('kwargs', kwargs)

        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(tmp_path / 'shared_downloads'),
        )
        expected_download_dir = os.path.realpath(args.download_dir)
        expected_ete_dir = os.path.join(expected_download_dir, 'ete_taxonomy')
        dbfile = os.path.join(expected_ete_dir, 'taxa.sqlite')
        os.makedirs(expected_ete_dir, exist_ok=True)
        with open(dbfile, 'wb') as fout:
            fout.write(b'sqlite')

        monkeypatch.setattr('amalgkit.util.acquire_exclusive_lock', DummyLock)
        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', RecordingNcbi)
        monkeypatch.setattr('amalgkit.download_utils.is_taxonomy_database_compatible', lambda path: path == dbfile)

        first = get_ete_ncbitaxa(args=args)
        second = get_ete_ncbitaxa(args=args)

        assert first is second
        assert captured['init_calls'] == 1
        assert captured['lock_calls'] == 1

    def test_get_ete_ncbitaxa_default_cache_ignores_constructor_id_collisions(self, monkeypatch):
        import builtins
        import amalgkit.download_utils as download_utils

        class FirstNcbi:
            pass

        class SecondNcbi:
            pass

        cache = download_utils._get_thread_local_ete_ncbitaxa_cache()
        cache.clear()
        original_id = builtins.id

        def fake_id(value):
            if value in {FirstNcbi, SecondNcbi}:
                return 424242
            return original_id(value)

        monkeypatch.setattr(builtins, 'id', fake_id)
        try:
            first = download_utils.get_ete_ncbitaxa(ncbitaxa_cls=FirstNcbi)
            second = download_utils.get_ete_ncbitaxa(ncbitaxa_cls=SecondNcbi)
        finally:
            cache.clear()

        assert isinstance(first, FirstNcbi)
        assert isinstance(second, SecondNcbi)
        assert first is not second

    def test_get_ete_ncbitaxa_serializes_custom_constructor_across_threads(self, tmp_path, monkeypatch):
        captured = {'active': 0, 'max_active': 0, 'init_calls': 0}
        state_lock = threading.Lock()

        class DummyLock:
            def __init__(self, lock_path, lock_label='Lock', poll_seconds=5, timeout_seconds=3600):
                _ = (lock_path, lock_label, poll_seconds, timeout_seconds)

            def __enter__(self):
                return None

            def __exit__(self, exc_type, exc, tb):
                return False

        class RecordingNcbi:
            def __init__(self, **kwargs):
                _ = kwargs
                with state_lock:
                    captured['active'] += 1
                    captured['max_active'] = max(captured['max_active'], captured['active'])
                    captured['init_calls'] += 1
                    current_active = captured['active']
                if current_active > 1:
                    raise RuntimeError('constructor entered concurrently')
                time.sleep(0.05)
                with state_lock:
                    captured['active'] -= 1

        monkeypatch.setattr('amalgkit.util.acquire_exclusive_lock', DummyLock)
        monkeypatch.setattr('amalgkit.download_utils.NcbiTaxonomy', RecordingNcbi)
        monkeypatch.setattr('amalgkit.download_utils.should_build_ncbi_taxonomy_db', lambda *args, **kwargs: False)

        args_by_label = {
            'first': SimpleNamespace(out_dir=str(tmp_path / 'out1'), download_dir=str(tmp_path / 'shared1')),
            'second': SimpleNamespace(out_dir=str(tmp_path / 'out2'), download_dir=str(tmp_path / 'shared2')),
        }
        results, failures = run_tasks_with_optional_threads(
            task_items=list(args_by_label),
            task_fn=lambda label: get_ete_ncbitaxa(args=args_by_label[label]),
            max_workers=2,
        )

        assert failures == []
        assert len(results) == 2
        assert captured['init_calls'] == 2
        assert captured['max_active'] == 1
