import gzip
from pathlib import Path
import json

import pandas
import pytest

from amalgkit import gsa, gsa_fastq, getfastq, merge
from amalgkit.gsa_snapshot import publish_gsa_snapshot, preferred_gsa_snapshot, capture_gsa_accession_source
from amalgkit.metadata_utils import Metadata
from tests.support.gsa import manifest_row, payloads_for_row
from tests.test_gsa_fastq import downloader
from tests.test_gsa_workflow import write_metadata, native_args, snapshot_args, install_native_inputs, deferred_table

pytestmark = pytest.mark.integration


def measured(row, count=4, fingerprint='a' * 64):
    return dict(row, total_spots=count, total_bases=count * 20, spot_length=20,
                read_count_status='measured', gsa_input_fingerprint=fingerprint)


def test_snapshot_rejects_source_changed_since_loading(tmp_path):
    old = manifest_row()
    write_metadata(tmp_path, [old])
    args = snapshot_args(tmp_path)
    loaded = Metadata.from_DataFrame(pandas.DataFrame([measured(old)]))
    new = dict(old, scientific_name='Oryza sativa')
    write_metadata(tmp_path, [new])
    with pytest.raises(ValueError, match='source changed since loading'):
        publish_gsa_snapshot(args, loaded)
    assert not (tmp_path / 'getfastq/metadata.tsv').exists()


def test_merge_rejects_previous_content_fingerprint(tmp_path):
    current = measured(manifest_row(), count=8, fingerprint='b' * 64)
    metadata = Metadata.from_DataFrame(pandas.DataFrame([current]))
    directory = tmp_path / 'getfastq/CRR0001'
    directory.mkdir(parents=True)
    pandas.DataFrame([measured(manifest_row(), count=4)]).to_csv(directory / 'getfastq_stats.tsv', sep='\t', index=False)
    with pytest.raises(ValueError, match='current input fingerprint'):
        merge.merge_fastp_stats_into_metadata(metadata, str(tmp_path), max_workers=1)
    assert metadata.df.loc[0, 'total_spots'] == 8, 'Stale stats replaced current measured counts despite different input digest'


@pytest.mark.parametrize('prefix', ['sample_1', 'sample_2'])
def test_explicit_mate_marker_takes_priority_over_sample_number(prefix):
    files = [{'filename': f'{prefix}_R{mate}.fastq.gz'} for mate in (1, 2)]
    assigned = gsa.assign_file_mates(files, 'paired')
    assert [entry['mate'] for entry in assigned] == [1, 2]


def test_each_corrupt_mate_gets_one_retry(tmp_path):
    row = manifest_row()
    args = native_args(tmp_path)
    directory = Path(gsa_fastq.cache_directory(args, row))
    payloads = payloads_for_row(row)
    for name in payloads:
        (directory / name).write_bytes(b'broken gzip')
    calls = []
    stats = gsa_fastq.prepare_run(args, row, downloader(payloads, calls))
    assert stats['total_spots'] == 4
    assert sorted(calls) == sorted(payloads)


def test_no_terminal_newline_does_not_corrupt_multilane_output(tmp_path):
    row = manifest_row(False, groups=2)
    args = native_args(tmp_path)
    payloads = {name: gzip.compress(gzip.decompress(data).rstrip(b'\n'))
                for name, data in payloads_for_row(row).items()}
    row.update(gsa_fastq.prepare_run(args, row, downloader(payloads, [])))
    output = tmp_path / 'output'
    output.mkdir()
    gsa_fastq.extract_run(args, row, str(output), 1, 8)
    with gzip.open(output / 'CRR0001.fastq.gz', 'rb') as handle:
        assert len(list(gsa_fastq._records(handle, 'output'))) == 8


def test_identical_paired_conversion_with_second_round(tmp_path, monkeypatch):
    row = manifest_row()
    write_metadata(tmp_path, [row])
    install_native_inputs(monkeypatch, row, payloads_for_row(row, count=10040))
    args = native_args(tmp_path, '--max_bp', '200', '--treat_identical_paired_as_single', 'yes')
    getfastq.getfastq_main(args)
    output = tmp_path / 'getfastq/CRR0001/CRR0001.amalgkit.fastq.gz'
    assert output.exists()
    stats = pandas.read_csv(output.parent / 'getfastq_stats.tsv', sep='\t').iloc[0]
    assert stats['bp_written'] >= 200


def test_completed_snapshot_input_can_be_reused_explicitly(tmp_path, monkeypatch):
    row = manifest_row()
    write_metadata(tmp_path, [row])
    install_native_inputs(monkeypatch, row)
    getfastq.getfastq_main(native_args(tmp_path))
    args = native_args(tmp_path, '--metadata', str(tmp_path / 'getfastq/metadata.tsv'))
    getfastq.getfastq_main(args)
    assert (tmp_path / 'getfastq/CRR0001/CRR0001_1.amalgkit.fastq.gz').exists()


def test_extraction_rejects_different_validated_content(tmp_path):
    row = manifest_row()
    args = native_args(tmp_path)
    first = gsa_fastq.prepare_run(args, row, downloader(payloads_for_row(row), []))
    old_row = dict(row, **first)
    directory = Path(gsa_fastq.cache_directory(args, row))
    for name, payload in payloads_for_row(row, length=11).items():
        (directory / name).write_bytes(payload)
    second = gsa_fastq.prepare_run(args, row, downloader({}, []))
    assert first['gsa_input_fingerprint'] != second['gsa_input_fingerprint']
    output = tmp_path / 'output'
    output.mkdir()
    with pytest.raises(ValueError):
        gsa_fastq.extract_run(args, old_row, str(output), 1, 4)


def test_parallel_snapshot_publications_preserve_all_measurements(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    import threading
    rows = [dict(manifest_row(run=f'CRR{i:04d}'), biosample=f'SAMC{i:04d}') for i in range(1, 7)]
    source = deferred_table(rows, minimum=3)
    write_metadata(tmp_path, source.df.to_dict('records'))
    args = snapshot_args(tmp_path)
    barrier = threading.Barrier(len(rows))
    def publish(index):
        table = Metadata.from_DataFrame(pandas.DataFrame([measured(source.df.iloc[index].to_dict(), index + 4)]))
        barrier.wait()
        return publish_gsa_snapshot(args, table)
    with ThreadPoolExecutor(max_workers=len(rows)) as executor:
        list(executor.map(publish, range(len(rows))))
    table = preferred_gsa_snapshot(native_args(tmp_path), str(tmp_path / 'metadata/metadata.tsv'), read_table=True)
    assert table['total_spots'].astype(float).tolist() == list(range(4, 10))
    assert table['gsa_selection_status'].eq('resolved').all()


def test_multiple_direct_id_jobs_preserve_snapshot_union(tmp_path):
    args = native_args(tmp_path, '--id_list', str(tmp_path / 'ids.txt'), '--batch', '1')
    (tmp_path / 'ids.txt').write_text('CRR0001\nCRR0002\n')
    source = Metadata.from_DataFrame(pandas.DataFrame([manifest_row(run=run) for run in ['CRR0001', 'CRR0002']]))
    capture_gsa_accession_source(args, source, ['CRR0001', 'CRR0002'])
    for batch, run in enumerate(['CRR0001', 'CRR0002'], 1):
        args.batch = batch
        table = Metadata.from_DataFrame(pandas.DataFrame([measured(manifest_row(run=run))]))
        publish_gsa_snapshot(args, table)
    snapshot = pandas.read_csv(tmp_path / 'getfastq/metadata.tsv', sep='\t')
    assert snapshot['run'].tolist() == ['CRR0001', 'CRR0002'], 'Later --id_list --batch replaces earlier job metadata'


@pytest.mark.parametrize('order', [(1, 2), (2, 1)])
def test_accession_array_workflow_keeps_all_runs_for_quant(tmp_path, monkeypatch, order):
    rows = [manifest_row(run=run) for run in ['CRR0001', 'CRR0002']]
    ids = tmp_path / 'ids.txt'
    ids.write_text('CRR0001\nCRR0002\n')
    payloads = {}
    for row in rows:
        payloads.update(payloads_for_row(row))
    calls = install_native_inputs(monkeypatch, rows[0], payloads)
    def fetch(run, args):
        row = dict(next(row for row in rows if row['run'] == run), gsa_retrieved_at=str(len(calls)))
        return Metadata.from_DataFrame(pandas.DataFrame([row]))
    monkeypatch.setattr(getfastq, 'fetch_gsa_metadata', fetch)
    for batch in order:
        getfastq.getfastq_main(native_args(tmp_path, '--id_list', str(ids), '--batch', str(batch)))
    path = tmp_path / 'getfastq/metadata.tsv'
    table = pandas.read_csv(path, sep='\t')
    assert table['run'].tolist() == ['CRR0001', 'CRR0002']
    assert table['read_count_status'].eq('measured').all()
    from amalgkit.quant import build_quant_tasks
    assert [run for run, _ in build_quant_tasks(Metadata.from_DataFrame(table))] == ['CRR0001', 'CRR0002']
    for row in rows:
        assert (tmp_path / f'getfastq/{row["run"]}/{row["run"]}_1.amalgkit.fastq.gz').exists()
    assert len(calls) == 4


def test_source_edit_during_download_blocks_workflow_publication(tmp_path, monkeypatch):
    row = manifest_row()
    write_metadata(tmp_path, [row])
    install_native_inputs(monkeypatch, row)
    original_transfer = getfastq.download_file_from_candidate_sources
    def transfer(**kwargs):
        result = original_transfer(**kwargs)
        write_metadata(tmp_path, [dict(row, scientific_name='Oryza sativa')])
        return result
    monkeypatch.setattr(getfastq, 'download_file_from_candidate_sources', transfer)
    with pytest.raises(ValueError, match='source changed since loading'):
        getfastq.getfastq_main(native_args(tmp_path))
    assert not (tmp_path / 'getfastq/metadata.tsv').exists()
    assert not list((tmp_path / 'getfastq').glob('CRR*/*.amalgkit.fastq.gz'))


def test_accession_list_edit_during_processing_is_rejected(tmp_path):
    ids = tmp_path / 'ids.txt'
    ids.write_text('CRR0001\n')
    args = native_args(tmp_path, '--id_list', str(ids))
    source = Metadata.from_DataFrame(pandas.DataFrame([manifest_row()]))
    capture_gsa_accession_source(args, source, ['CRR0001'])
    ids.write_text('CRR0002\n')
    with pytest.raises(ValueError, match='accession list changed'):
        publish_gsa_snapshot(args, Metadata.from_DataFrame(pandas.DataFrame([measured(manifest_row())])))


def test_measured_overlay_preserves_captured_annotations(tmp_path):
    row = manifest_row()
    write_metadata(tmp_path, [row])
    args = snapshot_args(tmp_path)
    table = Metadata.from_DataFrame(pandas.DataFrame([dict(measured(row), scientific_name='changed annotation')]))
    publish_gsa_snapshot(args, table)
    result = pandas.read_csv(tmp_path / 'getfastq/metadata.tsv', sep='\t')
    assert result.loc[0, 'scientific_name'] == row['scientific_name']


@pytest.mark.parametrize('measured_source', [False, True])
def test_matching_gsa_stats_preserve_content_provenance(tmp_path, measured_source):
    row = manifest_row()
    metadata = Metadata.from_DataFrame(pandas.DataFrame([measured(row) if measured_source else row]))
    directory = tmp_path / 'getfastq/CRR0001'
    directory.mkdir(parents=True)
    pandas.DataFrame([measured(row)]).to_csv(directory / 'getfastq_stats.tsv', sep='\t', index=False)
    merge.merge_fastp_stats_into_metadata(metadata, str(tmp_path), max_workers=1)
    assert metadata.df.loc[0, 'gsa_input_fingerprint'] == 'a' * 64
    assert metadata.df.loc[0, 'read_count_status'] == 'measured'
    assert metadata.df.loc[0, 'total_spots'] == 4


def test_identical_mate_stats_keep_original_input_spot_length(tmp_path, monkeypatch):
    row = manifest_row()
    write_metadata(tmp_path, [row])
    install_native_inputs(monkeypatch, row)
    getfastq.getfastq_main(native_args(tmp_path, '--treat_identical_paired_as_single', 'yes'))
    table = Metadata.from_DataFrame(pandas.read_csv(tmp_path / 'getfastq/metadata.tsv', sep='\t'))
    merge.merge_fastp_stats_into_metadata(table, str(tmp_path), max_workers=1)
    assert table.df.loc[0, 'spot_length'] == 20


def test_persistent_corruption_has_bounded_file_retries(tmp_path):
    row = manifest_row()
    payloads = {name: b'broken gzip' for name in payloads_for_row(row)}
    calls = []
    with pytest.raises(gsa_fastq.GsaCorruptInputError):
        gsa_fastq.prepare_run(native_args(tmp_path), row, downloader(payloads, calls))
    first_name = json.loads(row['gsa_fastq_files'])[0]['filename']
    assert calls.count(first_name) == 2
    assert len(calls) <= 4


@pytest.mark.parametrize('names', [
    ['sample_1.fastq.gz', 'sample_2.fastq.gz'],
    ['sample_1_read1.fq.gz', 'sample_1_read2.fq.gz'],
    ['sample-2-f1.fq.gz', 'sample-2-r2.fq.gz'],
])
def test_explicit_and_unambiguous_bare_mate_markers(names):
    result = gsa.assign_file_mates([{'filename': name} for name in names], 'paired')
    assert [entry['mate'] for entry in result] == [1, 2]


def test_conflicting_explicit_mate_markers_are_rejected():
    with pytest.raises(ValueError, match='Ambiguous GSA FASTQ mate'):
        gsa.assign_file_mates([{'filename': 'sample_R1_R2.fastq.gz'}], 'paired')


@pytest.mark.parametrize('bad_values', [
    {'gsa_input_fingerprint': ''},
    {'read_count_status': 'unknown'},
    {'total_spots': 4.5},
    {'total_bases': float('inf')},
    {'spot_length': 21},
    {'data_source': 'ncbi'},
])
def test_merge_rejects_invalid_gsa_input_provenance_without_overwriting_counts(tmp_path, bad_values):
    row = manifest_row()
    metadata = Metadata.from_DataFrame(pandas.DataFrame([measured(row)]))
    directory = tmp_path / 'getfastq/CRR0001'
    directory.mkdir(parents=True)
    pandas.DataFrame([dict(measured(row), **bad_values)]).to_csv(directory / 'getfastq_stats.tsv', sep='\t', index=False)
    with pytest.raises(ValueError):
        merge.merge_fastp_stats_into_metadata(metadata, str(tmp_path), max_workers=1)
    assert metadata.df.loc[0, 'total_spots'] == 4
    assert metadata.df.loc[0, 'total_bases'] == 80
