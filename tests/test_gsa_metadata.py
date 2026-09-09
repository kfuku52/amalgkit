import io
import json
import urllib.error
import urllib.parse

import pandas
import pytest

from amalgkit import gsa
from amalgkit.main import build_main_parser
from amalgkit.metadata import metadata_main, _build_metadata_cache_fingerprint
from tests.support.gsa import fake_pages, search_payload, EXP_URL, RUN_URL


@pytest.mark.parametrize('language', ['en', 'cn'])
@pytest.mark.parametrize('paired', [True, False])
def test_native_crr_metadata_normalizes_pages_without_local_fastq(monkeypatch, language, paired):
    pages = fake_pages(language, paired)
    monkeypatch.setattr(gsa.GsaClient, 'read', lambda self, url: pages[url])
    metadata = gsa.fetch_gsa_metadata('CRR0001')
    row = metadata.df.iloc[0]
    assert row['run'] == 'CRR0001'
    assert row['tissue'] == 'leaf'
    assert row['scientific_name'] == 'Arabidopsis thaliana'
    assert row['taxid'] == '3702'
    assert row['data_source'] == 'gsa'
    assert row['private_file'] == 'no'
    assert row['total_spots'] == ''
    assert row['read_count_status'] == 'unknown'
    files = json.loads(row['gsa_fastq_files'])
    assert [entry['mate'] for entry in files] == ([1, 2] if paired else [0])
    assert all(entry['sources'][0]['url'].startswith('https://download.cncb.ac.cn/') for entry in files)


def test_search_uses_all_pages_and_checks_returned_total(monkeypatch):
    calls = []
    def read(self, url):
        query = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
        calls.append(int(query['start'][0]))
        if calls[-1] == 0:
            return search_payload([{'id': 'CRA0001'}, {'id': 'CRX0001'}], total=3)
        return search_payload([{'id': 'CRX0002'}], total=3)
    monkeypatch.setattr(gsa.GsaClient, 'read', read)
    assert [row['id'] for row in gsa.GsaClient().search('leaf')] == ['CRA0001', 'CRX0001', 'CRX0002']
    assert calls == [0, 2]


@pytest.mark.parametrize('second', [search_payload([], 2), search_payload([{'id': 'CRX0001'}], 2),
                                   search_payload([{'id': 'CRX0002'}], 3), '{}', '<html>outage</html>'])
def test_incomplete_search_never_becomes_zero_results(monkeypatch, second):
    pages = iter([search_payload([{'id': 'CRX0001'}], 2), second])
    monkeypatch.setattr(gsa.GsaClient, 'read', lambda *args: next(pages))
    with pytest.raises(RuntimeError, match='Invalid/incomplete'):
        gsa.GsaClient().search('leaf')


@pytest.mark.parametrize('accession', ['CRA0001', 'CRX0001', 'PRJCA0001'])
def test_archive_accessions_resolve_experiments(monkeypatch, accession):
    pages = fake_pages()
    monkeypatch.setattr(gsa.GsaClient, 'read', lambda self, url: pages[url] if url in pages else search_payload([
        {'id': 'CRX0001', 'url': EXP_URL, 'attrs': {'BioProject': 'PRJCA0001'}},
        {'id': 'CRA0001', 'url': gsa.GSA_ROOT + 'browse/CRA0001'},
    ]))
    assert gsa.fetch_gsa_metadata(accession).df['run'].tolist() == ['CRR0001']


@pytest.mark.parametrize('query', ['HRA001', 'leaf[Title]', ''])
def test_unsupported_search_is_explicit(query):
    with pytest.raises(ValueError):
        gsa.fetch_gsa_metadata(query)


def test_changed_html_or_missing_public_files_is_error(monkeypatch):
    pages = fake_pages()
    pages[RUN_URL] = pages[RUN_URL].replace('download.cncb.ac.cn', 'not-download.example')
    monkeypatch.setattr(gsa.GsaClient, 'read', lambda self, url: pages[url])
    with pytest.raises(ValueError, match='No public supported FASTQ'):
        gsa.fetch_gsa_metadata('CRR0001')


def test_missing_mate_and_ambiguous_lane_names_are_errors():
    with pytest.raises(ValueError, match='Missing GSA FASTQ mate'):
        gsa.assign_file_mates([{'filename': 'run_f1.fq.gz'}], 'paired')
    with pytest.raises(ValueError, match='Cannot establish'):
        gsa.assign_file_mates([{'filename': 'run_a.fq.gz'}, {'filename': 'run_b.fq.gz'}], 'paired')
    entries = [{'filename': name} for name in ['run_L2_r2.fq.gz', 'run_L1_f1.fq.gz', 'run_L1_r2.fq.gz', 'run_L2_f1.fq.gz']]
    assert [(entry['group'], entry['mate']) for entry in gsa.assign_file_mates(entries, 'paired')] == [(0, 1), (0, 2), (1, 1), (1, 2)]


def test_request_retries_transient_errors_and_sets_user_agent(monkeypatch):
    calls = []
    def open_request(request, timeout):
        calls.append(request)
        if len(calls) < 3:
            raise urllib.error.URLError('temporary outage')
        return io.BytesIO(b'public metadata')
    monkeypatch.setattr(gsa, '_open_metadata', open_request)
    monkeypatch.setattr(gsa.time, 'sleep', lambda _: None)
    client = gsa.GsaClient()
    assert client.read(EXP_URL) == 'public metadata'
    assert client.read(EXP_URL) == 'public metadata'
    assert len(calls) == 3
    assert calls[0].get_header('User-agent').startswith('amalgkit/')
    with pytest.raises(ValueError, match='Unexpected'):
        client.read('file:///etc/passwd')


def test_metadata_cli_writes_portable_native_manifest_and_source_cache(monkeypatch, tmp_path):
    pages = fake_pages()
    monkeypatch.setattr(gsa.GsaClient, 'read', lambda self, url: pages[url])
    args = build_main_parser().parse_args(['metadata', '--resolve_names', 'no', '--source', 'gsa', '--accession', 'CRR0001', '--out_dir', str(tmp_path)])
    metadata_main(args)
    df = pandas.read_csv(tmp_path / 'metadata' / 'metadata.tsv', sep='\t')
    assert df.loc[0, 'run'] == 'CRR0001'
    assert df.loc[0, 'sample_group'] == 'leaf'
    assert json.loads(df.loc[0, 'gsa_fastq_files'])[0]['mate'] == 1
    info = json.loads((tmp_path / 'metadata' / 'metadata.query_info.json').read_text())
    assert info['source'] == 'gsa'
    fingerprint = _build_metadata_cache_fingerprint(args, 'leaf', None, None, 'single')
    args.source = 'ncbi'
    assert fingerprint != _build_metadata_cache_fingerprint(args, 'leaf', None, None, 'single')


def test_species_and_title_matching_are_postfiltered(monkeypatch):
    pages = fake_pages()
    monkeypatch.setattr(gsa.GsaClient, 'read', lambda self, url: pages[url])
    monkeypatch.setattr(gsa.GsaClient, 'experiment_urls', lambda *args: {'CRX0001': EXP_URL})
    assert gsa.fetch_gsa_metadata('leaf', species_name='Arabidopsis thaliana', title_terms=['leaf']).df.shape[0] == 1
    assert gsa.fetch_gsa_metadata('leaf', species_name='Arabidopsis lyrata').df.empty
    assert gsa.fetch_gsa_metadata('leaf', title_terms=['root']).df.empty


def test_species_batch_uses_big_search_and_preserves_outputs(monkeypatch, tmp_path):
    pages = fake_pages()
    queries = []
    def experiments(self, query):
        queries.append(query)
        return {'CRX0001': EXP_URL}
    monkeypatch.setattr(gsa.GsaClient, 'read', lambda self, url: pages[url])
    monkeypatch.setattr(gsa.GsaClient, 'experiment_urls', experiments)
    species = tmp_path / 'species.tsv'
    species.write_text('scientific_name\nArabidopsis thaliana\n')
    args = build_main_parser().parse_args(['metadata', '--resolve_names', 'no', '--source', 'gsa', '--species_tsv', str(species),
                                         '--mode', 'title_split', '--title_terms', 'leaf,root', '--out_dir', str(tmp_path)])
    metadata_main(args)
    assert all(query == '"Arabidopsis thaliana" AND "RNA-Seq"' for query in queries)
    merged = pandas.read_csv(tmp_path / 'metadata_specieswise/Arabidopsis_thaliana/Arabidopsis_thaliana.metadata.tsv', sep='\t')
    assert merged['run'].tolist() == ['CRR0001']


def test_metadata_redirects_cannot_leave_provider():
    handler = gsa._GsaRedirectHandler()
    with pytest.raises(ValueError, match='Unexpected GSA metadata URL'):
        handler.redirect_request(None, None, 302, '', {}, 'https://127.0.0.1/private')
    with pytest.raises(ValueError, match='Unexpected GSA metadata URL'):
        gsa._validate_metadata_url('https://user:pass@ngdc.cncb.ac.cn/gsa/')


def test_provider_download_hosts_are_https_only():
    from amalgkit.sra_sources import is_allowed_download_url
    assert is_allowed_download_url('https://download.cncb.ac.cn/gsa2/CRA0001/CRR0001/file.fastq.gz')
    assert not is_allowed_download_url('ftp://download.cncb.ac.cn/file.fastq.gz')
    assert not is_allowed_download_url('https://download.cncb.ac.cn.evil.example/file.fastq.gz')
    assert not is_allowed_download_url('https://127.0.0.1/file.fastq.gz')


@pytest.mark.parametrize('original,replacement,reason', [
    ('<td>fastq</td>', '<td>bam</td>', 'unsupported_gsa_format'),
    ('Illumina NovaSeq 6000', 'Oxford Nanopore', 'unsupported_gsa_long_read'),
    ('<td>PAIRED</td>', '<td>OTHER</td>', 'unsupported_gsa_layout'),
])
def test_known_unsupported_runs_remain_visible_but_excluded(monkeypatch, original, replacement, reason):
    pages = fake_pages()
    pages[RUN_URL] = pages[RUN_URL].replace(original, replacement)
    monkeypatch.setattr(gsa.GsaClient, 'read', lambda self, url: pages[url])
    row = gsa.fetch_gsa_metadata('CRR0001').df.iloc[0]
    assert row['exclusion'] == reason
    assert row['run'] == 'CRR0001'
    assert json.loads(row['gsa_fastq_files']) == []
