import gzip
import json
from types import SimpleNamespace

import pandas
import pytest

from amalgkit import getfastq, gsa
from amalgkit.gsa_select import defer_unknown_count_filter, defer_unknown_count_dedup, DEFERRED_COLUMN
from amalgkit.gsa_snapshot import publish_gsa_snapshot, preferred_gsa_snapshot
from amalgkit.main import build_main_parser
from amalgkit.metadata_utils import Metadata, load_metadata
from amalgkit.quant import build_quant_tasks
from tests.support.gsa import manifest_row, payloads_for_row, fake_pages, fastq_bytes
from tests.test_gsa_fastq import downloader

pytestmark = pytest.mark.integration


def write_metadata(tmp_path, rows):
    directory = tmp_path / 'metadata'
    directory.mkdir(exist_ok=True)
    table = Metadata.from_DataFrame(pandas.DataFrame(rows))
    table.df.to_csv(directory / 'metadata.tsv', sep='\t', index=False)
    return table


def native_args(tmp_path, *extra):
    return build_main_parser().parse_args(['getfastq', '--out_dir', str(tmp_path), '--fastp', 'no',
                                          '--min_read_length', '1', '--threads', '1', '--internal_jobs', '1', *extra])


def snapshot_args(tmp_path):
    args = native_args(tmp_path)
    args._capture_gsa_source = True
    load_metadata(args)
    return args


def forbid_sra(*args, **kwargs):
    raise AssertionError('Native GSA must not use Entrez, SRA download, fasterq-dump, or private FASTQ')


def install_native_inputs(monkeypatch, row, payloads=None):
    calls = []
    monkeypatch.setattr(getfastq, 'check_getfastq_dependency', lambda args: None)
    monkeypatch.setattr(getfastq, 'download_sra', forbid_sra)
    monkeypatch.setattr(getfastq, 'run_fasterq_dump', forbid_sra)
    monkeypatch.setattr(getfastq, 'sequence_extraction_private', forbid_sra)
    monkeypatch.setattr(getfastq, '_call_getfastq_getxml', forbid_sra)
    monkeypatch.setattr(getfastq, 'download_file_from_candidate_sources', downloader(payloads or payloads_for_row(row), calls))
    return calls


@pytest.mark.parametrize('paired', [True, False])
def test_metadata_to_getfastq_produces_standard_outputs_and_resumes(tmp_path, monkeypatch, paired):
    row = manifest_row(paired)
    write_metadata(tmp_path, [row])
    original = (tmp_path / 'metadata/metadata.tsv').read_bytes()
    calls = install_native_inputs(monkeypatch, row)
    args = native_args(tmp_path, '--max_bp', '1000')
    getfastq.getfastq_main(args)
    suffix = '_1' if paired else ''
    output = tmp_path / f'getfastq/CRR0001/CRR0001{suffix}.amalgkit.fastq.gz'
    assert gzip.decompress(output.read_bytes()).count(b'@read') == 4
    assert (tmp_path / 'metadata/metadata.tsv').read_bytes() == original
    stats = pandas.read_csv(tmp_path / 'getfastq/CRR0001/getfastq_stats.tsv', sep='\t')
    assert stats.loc[0, 'data_source'] == 'gsa'
    assert stats.loc[0, 'total_spots'] == 4
    assert stats.loc[0, 'num_written'] == 4
    state = json.loads((tmp_path / 'getfastq/CRR0001/getfastq_run_state.json').read_text())
    assert state['phase'] == 'complete'
    timestamp = output.stat().st_mtime_ns
    getfastq.getfastq_main(args)
    assert output.stat().st_mtime_ns == timestamp
    assert len(calls) == (2 if paired else 1)
    quant_args = SimpleNamespace(out_dir=str(tmp_path), metadata='inferred', _prefer_gsa_snapshot=True)
    metadata = load_metadata(quant_args)
    assert metadata.df.loc[0, 'total_spots'] == '4.0' or metadata.df.loc[0, 'total_spots'] == '4'
    assert build_quant_tasks(metadata) == [('CRR0001', 'Arabidopsis thaliana')]


def test_gsa_getfastq_id_resolves_remote_metadata_without_entrez(tmp_path, monkeypatch):
    pages = fake_pages()
    monkeypatch.setattr(gsa.GsaClient, 'read', lambda self, url: pages[url])
    resolved = gsa.fetch_gsa_metadata('CRR0001').df.iloc[0].to_dict()
    calls = install_native_inputs(monkeypatch, resolved)
    getfastq.getfastq_main(native_args(tmp_path, '--id', 'CRR0001'))
    assert len(calls) == 2
    assert (tmp_path / 'getfastq/CRR0001/CRR0001_1.amalgkit.fastq.gz').exists()


def test_native_minimum_read_length_preserves_pairs_and_second_round_reuses_inputs(tmp_path, monkeypatch):
    row = manifest_row()
    payloads = {}
    for entry in json.loads(row['gsa_fastq_files']):
        mate = entry['mate']
        raw = fastq_bytes(9999, mate=mate)
        raw += fastq_bytes(2, mate=mate, offset=9999, length=5)
        raw += fastq_bytes(19, mate=mate, offset=10001)
        payloads[entry['filename']] = gzip.compress(raw, mtime=0)
    write_metadata(tmp_path, [row])
    calls = install_native_inputs(monkeypatch, row, payloads)
    args = native_args(tmp_path, '--min_read_length', '8', '--max_bp', '80')
    getfastq.getfastq_main(args)
    stats = pandas.read_csv(tmp_path / 'getfastq/CRR0001/getfastq_stats.tsv', sep='\t').iloc[0]
    assert stats['num_rejected'] == 2
    assert stats['num_written'] >= 4
    assert stats['spot_start_2nd'] > stats['spot_end_1st']
    assert len(calls) == 2
    outputs = [gzip.decompress((tmp_path / f'getfastq/CRR0001/CRR0001_{mate}.amalgkit.fastq.gz').read_bytes()) for mate in (1, 2)]
    for output in outputs:
        headers = [line for line in output.splitlines() if line.startswith(b'@')]
        assert len(headers) == len(set(headers))
        assert b'@read9999/' not in output and b'@read10000/' not in output
    assert outputs[0].count(b'@') == outputs[1].count(b'@')


def test_native_dependency_probe_does_not_require_fasterq(monkeypatch):
    probes = []
    def probe(**kwargs):
        probes.append(kwargs['label'])
        return SimpleNamespace(stdout=b'', stderr=b'', returncode=0), '', ''
    monkeypatch.setattr(getfastq, 'probe_dependency_command', probe)
    args = SimpleNamespace(_requires_sra_toolkit=False, fastp=False, rrna_filter=False, contam_filter=False)
    getfastq.check_getfastq_dependency(args)
    assert 'seqkit' in probes
    assert 'fasterq-dump' not in probes


def deferred_table(rows, minimum=3, dedup=False):
    metadata = Metadata.from_DataFrame(pandas.DataFrame(rows))
    protected = pandas.Series(False, index=metadata.df.index)
    mask = defer_unknown_count_filter(metadata.df, 'total_spots', minimum, 'exclusion', 'low_nspots', protected)
    assert mask.any()
    if dedup:
        defer_unknown_count_dedup(metadata.df, {'target_column': 'total_spots', 'columns': ['bioproject', 'biosample'],
                                              'outcome': 'redundant_biosample'}, ~protected)
    return metadata


def test_unknown_counts_defer_only_gsa_and_keep_ncbi_unchanged():
    rows = [manifest_row(), dict(manifest_row(run='CRR0002'), data_source='ncbi')]
    table = deferred_table(rows)
    assert json.loads(table.df.loc[0, DEFERRED_COLUMN])[0]['threshold'] == 3
    assert table.df.loc[1, DEFERRED_COLUMN] == ''


def test_array_snapshots_aggregate_counts_without_changing_source_or_batch_ids(tmp_path):
    source = deferred_table([manifest_row(run='CRR0001'), manifest_row(run='CRR0002')], dedup=True)
    write_metadata(tmp_path, source.df.to_dict('records'))
    source_path = tmp_path / 'metadata/metadata.tsv'
    original = source_path.read_bytes()
    args = snapshot_args(tmp_path)
    first = Metadata.from_DataFrame(source.df.iloc[[0]])
    first.df['total_spots'] = 10
    first.df['total_bases'] = 200
    first.df['spot_length'] = 20
    first.df['read_count_status'] = 'measured'
    publish_gsa_snapshot(args, first)
    with pytest.raises(ValueError, match='pending input counts'):
        preferred_gsa_snapshot(args, str(source_path))
    second = Metadata.from_DataFrame(source.df.iloc[[1]])
    second.df['total_spots'] = 20
    second.df['total_bases'] = 400
    second.df['spot_length'] = 20
    second.df['read_count_status'] = 'measured'
    publish_gsa_snapshot(args, second)
    snapshot = pandas.read_csv(tmp_path / 'getfastq/metadata.tsv', sep='\t')
    assert snapshot['run'].tolist() == ['CRR0001', 'CRR0002']
    assert snapshot['exclusion'].tolist() == ['redundant_biosample', 'no']
    assert snapshot['is_sampled'].tolist() == ['yes', 'yes']
    # Repeating an earlier array job must not resurrect a redundant sample.
    publish_gsa_snapshot(args, first)
    repeated = pandas.read_csv(tmp_path / 'getfastq/metadata.tsv', sep='\t')
    assert repeated['exclusion'].tolist() == ['redundant_biosample', 'no']
    assert source_path.read_bytes() == original
    quant_args = SimpleNamespace(out_dir=str(tmp_path), metadata='inferred', _prefer_gsa_snapshot=True, batch=2)
    loaded = load_metadata(quant_args)
    assert loaded.df['run'].tolist() == ['CRR0002']
    source_path.write_text(source_path.read_text().replace('leaf', 'root'))
    assert preferred_gsa_snapshot(args, str(source_path)) == str(source_path)


def test_pending_count_threshold_is_applied_after_real_measurement(tmp_path, monkeypatch):
    source = deferred_table([manifest_row()], minimum=100)
    write_metadata(tmp_path, source.df.to_dict('records'))
    install_native_inputs(monkeypatch, manifest_row())
    with pytest.raises(ValueError, match='No eligible getfastq'):
        getfastq.getfastq_main(native_args(tmp_path))
    snapshot = pandas.read_csv(tmp_path / 'getfastq/metadata.tsv', sep='\t')
    assert snapshot.loc[0, 'total_spots'] == 4
    assert snapshot.loc[0, 'exclusion'] == 'low_nspots'
    assert not list((tmp_path / 'getfastq').glob('CRR*/*.amalgkit.fastq.gz'))


def test_native_workflow_uses_real_http_transfer(tmp_path, monkeypatch):
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    import threading
    import urllib.parse
    import urllib.request
    row = manifest_row()
    payloads = payloads_for_row(row)
    requests = []
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            name = self.path.rsplit('/', 1)[-1]
            requests.append(name)
            payload = payloads[name]
            self.send_response(200)
            self.send_header('Content-Length', str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
        def log_message(self, *args):
            pass
    server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    worker = threading.Thread(target=server.serve_forever, kwargs={'poll_interval': 0.01})
    worker.start()
    real_transfer = getfastq.download_file_from_candidate_sources
    def open_fixture(url, timeout):
        assert urllib.parse.urlparse(url).hostname == 'download.cncb.ac.cn'
        path = urllib.parse.urlparse(url).path
        return urllib.request.urlopen(f'http://127.0.0.1:{server.server_port}{path}', timeout=timeout)
    try:
        write_metadata(tmp_path, [row])
        install_native_inputs(monkeypatch, row)
        monkeypatch.setattr(getfastq, 'download_file_from_candidate_sources', real_transfer)
        monkeypatch.setattr(getfastq, 'build_allowed_host_opener', lambda: SimpleNamespace(open=open_fixture))
        getfastq.getfastq_main(native_args(tmp_path, '--sra_download_method', 'urllib'))
    finally:
        server.shutdown()
        worker.join()
        server.server_close()
    assert sorted(requests) == sorted(payloads)
    output = tmp_path / 'getfastq/CRR0001/CRR0001_1.amalgkit.fastq.gz'
    assert gzip.decompress(output.read_bytes()).count(b'@read') == 4


def test_ncbi_second_round_preserves_allocated_range_after_checkpoint_restore(tmp_path, monkeypatch):
    from contextlib import nullcontext
    row = {'run': 'SRR001', 'scientific_name': 'Arabidopsis thaliana', 'lib_layout': 'paired',
           'total_spots': 20000, 'total_bases': 400000, 'spot_length': 20,
           'spot_start_2nd': 10010, 'spot_end_2nd': 10020}
    metadata = Metadata.from_DataFrame(pandas.DataFrame([row]))
    def restore(**kwargs):
        current = kwargs['metadata']
        current.df['spot_start_2nd'] = 0
        current.df['spot_end_2nd'] = 0
        return current, {'sra_id': 'SRR001', 'getfastq_sra_dir': str(tmp_path)}, {'phase': 'first_round'}, False
    ranges = []
    def extract(args, sra_stat, metadata, g, runtime_context=None):
        ranges.append(tuple(metadata.df.loc[0, ['spot_start_2nd', 'spot_end_2nd']]))
        return metadata
    monkeypatch.setattr(getfastq, '_inspect_getfastq_run_after_lock', restore)
    monkeypatch.setattr(getfastq, 'acquire_exclusive_lock', lambda **kwargs: nullcontext())
    monkeypatch.setattr(getfastq, 'write_getfastq_run_state', lambda *args, **kwargs: None)
    monkeypatch.setattr(getfastq, 'write_getfastq_stats', lambda *args, **kwargs: None)
    monkeypatch.setattr(getfastq, 'sequence_extraction_2nd_round', extract)
    getfastq._process_getfastq_second_round_run(native_args(tmp_path), 0, 'SRR001', metadata.df, {})
    assert ranges == [(10010, 10020)]


def test_corrupt_measured_snapshot_blocks_consumers_and_can_be_rebuilt(tmp_path):
    row = manifest_row()
    write_metadata(tmp_path, [row])
    measured = dict(row, total_spots=4, total_bases=80, spot_length=20, read_count_status='measured')
    metadata = Metadata.from_DataFrame(pandas.DataFrame([measured]))
    args = snapshot_args(tmp_path)
    publish_gsa_snapshot(args, metadata)
    snapshot = tmp_path / 'getfastq/metadata.tsv'
    snapshot.write_text('corrupt')
    with pytest.raises(ValueError, match='snapshot is invalid'):
        preferred_gsa_snapshot(args, str(tmp_path / 'metadata/metadata.tsv'))
    publish_gsa_snapshot(args, metadata)
    assert preferred_gsa_snapshot(args, str(tmp_path / 'metadata/metadata.tsv'), read_table=True).loc[0, 'run'] == 'CRR0001'
