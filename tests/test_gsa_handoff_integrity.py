import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import threading
import urllib.request
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pandas
import pytest

from amalgkit import getfastq, gsa_fastq
from amalgkit.gsa_snapshot import publish_gsa_snapshot
from amalgkit.gsa_select import resolve_deferred_selection
from amalgkit.metadata_utils import Metadata, load_metadata
from amalgkit.quant import build_quant_tasks
from amalgkit.select import apply_select_filter_rules
from tests.support.gsa import manifest_row, payloads_for_row, fastq_bytes
from tests.test_gsa_workflow import native_args, snapshot_args, write_metadata, deferred_table, install_native_inputs
from tests.test_gsa_regressions import measured

pytestmark = pytest.mark.integration


def test_short_http_body_cannot_become_complete_gsa_input(tmp_path, monkeypatch):
    row = manifest_row()
    parts = payloads_for_row(row)
    complete = {entry['filename']: parts[entry['filename']] + gzip.compress(fastq_bytes(4, entry['mate'], offset=4))
                for entry in json.loads(row['gsa_fastq_files'])}
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            name = self.path.rsplit('/', 1)[-1]
            self.send_response(200)
            self.send_header('Content-Length', str(len(complete[name])))
            self.end_headers()
            self.wfile.write(parts[name])
        def log_message(self, *args):
            pass
    server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    worker = threading.Thread(target=server.serve_forever, kwargs={'poll_interval': 0.01})
    worker.start()
    def opener(url, timeout):
        assert urllib.parse.urlparse(url).hostname == 'download.cncb.ac.cn'
        return urllib.request.urlopen(f'http://127.0.0.1:{server.server_port}' + urllib.parse.urlparse(url).path, timeout=timeout)
    try:
        monkeypatch.setattr(getfastq, 'build_allowed_host_opener', lambda: SimpleNamespace(open=opener))
        args = native_args(tmp_path, '--sra_download_method', 'urllib')
        try:
            stats = gsa_fastq.prepare_run(args, row, getfastq.download_file_from_candidate_sources)
        except (OSError, ValueError):
            return
        assert stats['total_spots'] == 8, 'Truncated HTTP response certified as a complete 4-spot input'
    finally:
        server.shutdown()
        worker.join()
        server.server_close()


def test_partial_content_range_cannot_become_complete_gsa_input(tmp_path, monkeypatch):
    import subprocess
    row = manifest_row(False)
    payload = next(iter(payloads_for_row(row).values()))
    def run(command, **kwargs):
        Path(command[command.index('--dump-header') + 1]).write_bytes(
            f'HTTP/1.1 206 Partial Content\r\nContent-Range: bytes 0-{len(payload)-1}/{len(payload)*2}\r\n\r\n'.encode())
        Path(command[command.index('-o') + 1]).write_bytes(payload)
        return subprocess.CompletedProcess(command, 0, stdout=b'AMALGKIT_CURL_STATUS:206\nAMALGKIT_CURL_REDIRECT:\n', stderr=b'')
    monkeypatch.setattr(getfastq.subprocess, 'run', run)
    output = tmp_path / 'download.part'
    assert not getfastq.download_with_curl(
        source_url=json.loads(row['gsa_fastq_files'])[0]['sources'][0]['url'], output_path=str(output),
        args=native_args(tmp_path), sra_source_name='GSA', artifact_label='GSA original FASTQ', resume_existing=True)


def test_explicit_snapshot_blocks_pending_selection_like_inferred(tmp_path):
    rows = [manifest_row(), manifest_row(run='CRR0002')]
    source = deferred_table(rows, dedup=True)
    write_metadata(tmp_path, source.df.to_dict('records'))
    args = snapshot_args(tmp_path)
    first = Metadata.from_DataFrame(pandas.DataFrame([measured(source.df.iloc[0].to_dict())]))
    publish_gsa_snapshot(args, first)
    args.metadata = str(tmp_path / 'getfastq/metadata.tsv')
    args._prefer_gsa_snapshot = True
    with pytest.raises(ValueError, match='pending'):
        build_quant_tasks(load_metadata(args))


def test_reselecting_pending_input_replaces_old_threshold():
    source = deferred_table([manifest_row()], minimum=100)
    rules = [{'stage': 'filter', 'target_column': 'exclusion', 'action': 'exclude_if_lt_parameter',
              'parameter_name': 'min_nspots', 'rule_id': 'low_depth', 'columns': ['total_spots'],
              'outcome': 'low_nspots', 'stop_on_match': True}]
    apply_select_filter_rules(source, SimpleNamespace(min_nspots=3), rules)
    source.df['read_count_status'] = 'measured'
    source.df['total_spots'] = 4
    resolve_deferred_selection(source)
    assert source.df.loc[0, 'exclusion'] == 'no', 'Old threshold 100 survives a new selection threshold 3'


def test_whitespace_source_marker_keeps_supported_workflow(tmp_path, monkeypatch):
    row = dict(manifest_row(), data_source=' GSA ')
    write_metadata(tmp_path, [row])
    install_native_inputs(monkeypatch, row)
    getfastq.getfastq_main(native_args(tmp_path))


def test_snapshot_as_input_allows_parallel_array_publications(tmp_path):
    rows = [measured(manifest_row(run=run)) for run in ['CRR0001', 'CRR0002']]
    for row, mean in zip(rows, [150, 350]):
        row.update(fragment_length_mean=mean, fragment_length_sd=7,
                   fragment_length_source='measured', fragment_length_source_detail='insert assay L1')
    write_metadata(tmp_path, rows)
    publish_gsa_snapshot(snapshot_args(tmp_path), Metadata.from_DataFrame(pandas.DataFrame(rows)))
    jobs = []
    for batch in (1, 2):
        args = native_args(tmp_path, '--metadata', str(tmp_path / 'getfastq/metadata.tsv'), '--batch', str(batch))
        args._capture_gsa_source = True
        table = load_metadata(args)
        table.df['gsa_input_seconds'] = float(batch)
        jobs.append((args, table))
    for args, table in jobs:
        publish_gsa_snapshot(args, table)
    result = pandas.read_csv(tmp_path / 'getfastq/metadata.tsv', sep='\t')
    assert result['gsa_input_seconds'].tolist() == [1.0, 2.0]
    assert result['fragment_length_mean'].tolist() == [150, 350]
    assert result['fragment_length_sd'].tolist() == [7, 7]
    assert result['fragment_length_source_detail'].tolist() == ['insert assay L1', 'insert assay L1']


def test_full_select_keeps_gsa_manifest_usable_by_getfastq(tmp_path, monkeypatch):
    from amalgkit.main import build_main_parser
    from amalgkit.select import select_main
    row = dict(manifest_row(), taxid_species=3702)
    write_metadata(tmp_path, [row])
    install_native_inputs(monkeypatch, row)
    base = Path(__file__).resolve().parents[1] / 'amalgkit/select_rule_sets/base/select_rules.tsv'
    config = pandas.read_csv(base, sep='\t', dtype=str, keep_default_na=False)
    config.loc[config['parameter_name'].eq('min_nspots'), 'parameter_value'] = '3'
    rules = tmp_path / 'select_rules.tsv'
    config.to_csv(rules, sep='\t', index=False)
    args = build_main_parser().parse_args(['select', '--out_dir', str(tmp_path), '--select_rules_tsv', str(rules)])
    select_main(args)
    selected = pandas.read_csv(tmp_path / 'metadata/metadata.tsv', sep='\t')
    assert selected.loc[0, 'exclusion'] == 'no'
    assert selected.loc[0, 'is_sampled'] == 'yes'
    assert json.loads(selected.loc[0, 'gsa_fastq_files']) == json.loads(row['gsa_fastq_files'])
    getfastq.getfastq_main(native_args(tmp_path))
    assert (tmp_path / 'getfastq/CRR0001/CRR0001_1.amalgkit.fastq.gz').exists()


def test_self_snapshot_late_job_shares_immutable_generation(tmp_path):
    rows = [measured(manifest_row(run=run)) for run in ['CRR0001', 'CRR0002', 'CRR0003']]
    write_metadata(tmp_path, rows)
    publish_gsa_snapshot(snapshot_args(tmp_path), Metadata.from_DataFrame(pandas.DataFrame(rows)))
    def capture(batch):
        args = native_args(tmp_path, '--metadata', str(tmp_path / 'getfastq/metadata.tsv'), '--batch', str(batch))
        args._capture_gsa_source = True
        table = load_metadata(args)
        table.df['gsa_input_seconds'] = float(batch)
        return args, table
    early, delayed = capture(1), capture(2)
    publish_gsa_snapshot(*early)
    late = capture(3)
    assert early[0]._gsa_snapshot_source['digest'] == late[0]._gsa_snapshot_source['digest']
    publish_gsa_snapshot(*late)
    publish_gsa_snapshot(*delayed)
    result = pandas.read_csv(tmp_path / 'getfastq/metadata.tsv', sep='\t')
    assert result['gsa_input_seconds'].tolist() == [1.0, 2.0, 3.0]


def test_self_snapshot_still_rejects_user_edits(tmp_path):
    rows = [measured(manifest_row())]
    write_metadata(tmp_path, rows)
    publish_gsa_snapshot(snapshot_args(tmp_path), Metadata.from_DataFrame(pandas.DataFrame(rows)))
    path = tmp_path / 'getfastq/metadata.tsv'
    args = native_args(tmp_path, '--metadata', str(path))
    args._capture_gsa_source = True
    table = load_metadata(args)
    path.write_text(path.read_text().replace('leaf', 'root'))
    with pytest.raises(ValueError, match='source changed'):
        publish_gsa_snapshot(args, table)


@pytest.mark.parametrize('status', [200, 206])
def test_complete_http_response_passes_length_validation(tmp_path, status):
    import io
    response = io.BytesIO(b'ACGT')
    response.headers = {'Content-Length': '4', 'Content-Range': 'bytes 0-3/4'}
    response.status = status
    path = tmp_path / 'payload'
    getfastq.download_with_urllib('https://download.cncb.ac.cn/file', str(path), 30,
                                 urlopen_fn=lambda *args, **kwargs: response)
    assert path.read_bytes() == b'ACGT'


def test_reselection_can_remove_all_deferred_rules():
    from amalgkit.select import apply_select_filters
    table = deferred_table([manifest_row()], minimum=100, dedup=True)
    args = SimpleNamespace(mark_redundant_biosamples=False, max_sample=99999, random_seed=0,
                           sampling_strategy='maximize_bioproject_diversity')
    apply_select_filters(table, args, [])
    assert json.loads(table.df.loc[0, 'gsa_deferred_select']) == []


def test_explicit_pending_merge_input_is_rejected(tmp_path):
    table = deferred_table([manifest_row()])
    path = tmp_path / 'pending.tsv'
    table.df.to_csv(path, sep='\t', index=False)
    args = SimpleNamespace(out_dir=str(tmp_path), metadata=str(path), _prefer_gsa_snapshot=True)
    with pytest.raises(ValueError, match='pending'):
        load_metadata(args)
