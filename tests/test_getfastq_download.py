import pytest
import pandas
import numpy
import gzip
import http.client
import ssl
import os
import stat
import subprocess
from contextlib import contextmanager

import urllib.error
from types import SimpleNamespace


from amalgkit.getfastq import (
    initialize_columns,
    sequence_extraction_2nd_round,
    build_sra_source_candidates,
    resolve_sra_download_sources,
    download_file_from_candidate_sources,
    download_file_from_source,
    download_public_original_fastq_files,
    download_with_curl,
    download_verified_original_fastq,
    _public_fastq_resume_directory,
    download_sra,
    run_fasterq_dump,
    count_fastq_records,
    apply_first_round_getfastq_results,
)
from amalgkit.util import Metadata



class TestRunFasterqDumpOutputValidation:
    @staticmethod
    def _metadata_for_extraction(sra_id):
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'num_dumped': [0],
            'num_rejected': [0],
            'num_written': [0],
            'bp_dumped': [0],
            'bp_rejected': [0],
            'bp_written': [0],
            'layout_amalgkit': ['single'],
        }))

    @staticmethod
    def _args_for_fasterq_dump():
        class Args:
            threads = 1
            min_read_length = 25
            dump_print = False
            fasterq_dump_exe = 'fasterq-dump'
            fasterq_size_check = True
            fasterq_disk_limit = None
            fasterq_disk_limit_tmp = None
        return Args()

    def test_retries_then_fails_when_exit_zero_generates_no_fastq_files(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        (tmp_path / '{}.sra'.format(sra_id)).write_text('initial')
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()

        run_calls = {'count': 0}

        def fake_run(cmd, stdout=None, stderr=None):
            run_calls['count'] += 1
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        def fake_download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
            _ = (metadata, args)
            assert overwrite is True
            (tmp_path / '{}.sra'.format(sra_stat['sra_id'])).write_text('fresh')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.download_sra', fake_download_sra)
        monkeypatch.setattr(
            'amalgkit.getfastq.download_public_original_fastq_files',
            lambda **_kwargs: None,
        )

        with pytest.raises(RuntimeError, match='after re-download'):
            run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 2

    def test_retries_then_fails_when_exit_zero_generates_empty_fastq(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        (tmp_path / '{}.sra'.format(sra_id)).write_text('initial')
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()
        run_calls = {'count': 0}
        fallback_calls = {'count': 0}

        def fake_run(cmd, stdout=None, stderr=None):
            run_calls['count'] += 1
            (tmp_path / '{}.fastq'.format(sra_id)).write_bytes(b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        def fake_download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
            _ = (metadata, args, work_dir)
            assert overwrite is True
            (tmp_path / '{}.sra'.format(sra_stat['sra_id'])).write_text('fresh')

        def fake_original_fastq_fallback(**_kwargs):
            fallback_calls['count'] += 1
            return None

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.download_sra', fake_download_sra)
        monkeypatch.setattr(
            'amalgkit.getfastq.download_public_original_fastq_files',
            fake_original_fastq_fallback,
        )

        with pytest.raises(RuntimeError, match='after re-download'):
            run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert run_calls['count'] == 2
        assert fallback_calls['count'] == 1
        assert not (tmp_path / '{}.fastq'.format(sra_id)).exists()

    def test_raises_when_fasterq_dump_output_path_is_directory(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()
        (tmp_path / '{}.fastq'.format(sra_id)).mkdir()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(IsADirectoryError, match='FASTQ artifact exists but is not a file'):
            run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)


class TestSequenceExtractionSecondRound:
    def test_extracts_single_inclusive_spot_when_start_equals_end(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'

        def write_fastq_gz(path, read_id):
            with gzip.open(path, 'wt') as handle:
                handle.write('@{}\nACGT\n+\nIIII\n'.format(read_id))

        write_fastq_gz(tmp_path / '{}.amalgkit.fastq.gz'.format(sra_id), 'first')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'spot_start_2nd': [2],
            'spot_end_2nd': [2],
            'bp_written': [100],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        observed = {}

        def fake_sequence_extraction(args, sra_stat, metadata, g, start, end, runtime_context=None):
            observed['range'] = (start, end)
            write_fastq_gz(
                tmp_path / '{}.amalgkit.fastq.gz'.format(sra_id),
                'second',
            )
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.sequence_extraction', fake_sequence_extraction)

        sequence_extraction_2nd_round(
            args=SimpleNamespace(threads=1),
            sra_stat=sra_stat,
            metadata=metadata,
            g={},
        )

        assert observed['range'] == (2, 2)
        assert count_fastq_records(str(tmp_path / '{}.amalgkit.fastq.gz'.format(sra_id))) == 2

    def test_raises_when_second_round_fastq_is_missing(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        (tmp_path / '{}.amalgkit.fastq.gz'.format(sra_id)).write_text('AAAA\n')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
            'spot_start_2nd': [2],
            'spot_end_2nd': [5],
            'bp_written': [100],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        args = SimpleNamespace()
        g = {}

        monkeypatch.setattr('amalgkit.getfastq.sequence_extraction', lambda *a, **k: metadata)

        with pytest.raises(FileNotFoundError, match='Dumped fastq not found'):
            sequence_extraction_2nd_round(args=args, sra_stat=sra_stat, metadata=metadata, g=g)

    def test_raises_when_second_round_fastq_path_is_directory(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        (tmp_path / '{}.amalgkit.fastq.gz'.format(sra_id)).write_text('AAAA\n')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
            'spot_start_2nd': [2],
            'spot_end_2nd': [5],
            'bp_written': [100],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        args = SimpleNamespace()
        g = {}

        def fake_sequence_extraction(*_args, **_kwargs):
            tmp_first = tmp_path / '{}.amalgkit_1st.fastq.gz'.format(sra_id)
            if tmp_first.exists():
                tmp_first.unlink()
            tmp_first.mkdir()
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.sequence_extraction', fake_sequence_extraction)

        with pytest.raises(IsADirectoryError, match='Dumped fastq path exists but is not a file'):
            sequence_extraction_2nd_round(args=args, sra_stat=sra_stat, metadata=metadata, g=g)

    def test_paired_2nd_round_merges_subexts_with_parallel_workers(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'

        def write_fastq_gz(path, read_name, seq):
            with gzip.open(path, 'wt') as handle:
                handle.write('@{}\n'.format(read_name))
                handle.write(seq + '\n')
                handle.write('+\n')
                handle.write('I' * len(seq) + '\n')

        def read_fastq_gz(path):
            with gzip.open(path, 'rt') as handle:
                return handle.read()

        write_fastq_gz(
            tmp_path / '{}_1.amalgkit.fastq.gz'.format(sra_id),
            'first 1:N:0:ACGT',
            'FIRST1',
        )
        write_fastq_gz(
            tmp_path / '{}_2.amalgkit.fastq.gz'.format(sra_id),
            'first 2:N:0:ACGT',
            'FIRST2',
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['paired'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
            'spot_start_2nd': [2],
            'spot_end_2nd': [5],
            'bp_written': [100],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        args = SimpleNamespace(threads=2)
        g = {}
        observed = {'max_workers': None, 'task_items': None}

        def fake_sequence_extraction(*_args, **_kwargs):
            write_fastq_gz(
                tmp_path / '{}_1.amalgkit.fastq.gz'.format(sra_id),
                'second 1:N:0:ACGT',
                'SECOND1',
            )
            write_fastq_gz(
                tmp_path / '{}_2.amalgkit.fastq.gz'.format(sra_id),
                'second 2:N:0:ACGT',
                'SECOND2',
            )
            return metadata

        def fake_run_tasks(task_items, task_fn, max_workers):
            observed['max_workers'] = max_workers
            observed['task_items'] = list(task_items)
            results = {}
            failures = []
            for task_item in task_items:
                try:
                    results[task_item] = task_fn(task_item)
                except Exception as exc:
                    failures.append((task_item, exc))
            return results, failures

        monkeypatch.setattr('amalgkit.getfastq.sequence_extraction', fake_sequence_extraction)
        monkeypatch.setattr('amalgkit.getfastq.run_tasks_with_optional_threads', fake_run_tasks)

        sequence_extraction_2nd_round(args=args, sra_stat=sra_stat, metadata=metadata, g=g)

        assert observed['max_workers'] == 2
        assert set(observed['task_items']) == {'_1', '_2'}
        assert read_fastq_gz(tmp_path / '{}_1.amalgkit.fastq.gz'.format(sra_id)) == (
            '@first 1:N:0:ACGT\nFIRST1\n+\nIIIIII\n'
            '@second 1:N:0:ACGT\nSECOND1\n+\nIIIIIII\n'
        )
        assert read_fastq_gz(tmp_path / '{}_2.amalgkit.fastq.gz'.format(sra_id)) == (
            '@first 2:N:0:ACGT\nFIRST2\n+\nIIIIII\n'
            '@second 2:N:0:ACGT\nSECOND2\n+\nIIIIIII\n'
        )


class TestDownloadSraUrlSchemes:
    @staticmethod
    def _make_metadata(
        sra_id,
        aws_link,
        gcp_link,
        ncbi_link,
        experiment='',
        ena_sra_link='',
        ddbj_sra_link='',
    ):
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'experiment': [experiment],
            'AWS_Link': [aws_link],
            'GCP_Link': [gcp_link],
            'NCBI_Link': [ncbi_link],
            'ENA_SRA_Link': [ena_sra_link],
            'DDBJ_SRA_Link': [ddbj_sra_link],
        }))

    @staticmethod
    def _make_args(gcp_project='', sra_download_method='urllib', ena=False, ddbj=False):
        args = type('Args', (), {})()
        args.aws = True
        args.gcp = True
        args.ncbi = True
        args.ena = ena
        args.ddbj = ddbj
        args.gcp_project = gcp_project
        args.sra_download_method = sra_download_method
        args.ena_download_max_concurrency = 'auto'
        args.ddbj_download_max_concurrency = 'auto'
        return args

    def test_build_sra_source_candidates_uses_explicit_priority(self):
        source_candidates = build_sra_source_candidates({
            'DDBJ': 'https://example.invalid/ddbj.sra',
            'NCBI': 'https://example.invalid/ncbi.sra',
            'AWS': 'https://example.invalid/aws.sra',
            'ENA': 'https://example.invalid/ena.sra',
            'GCP': 'https://example.invalid/gcp.sra',
        })

        assert [source['source_name'] for source in source_candidates] == [
            'AWS',
            'GCP',
            'NCBI',
            'ENA',
            'DDBJ',
        ]

    def test_malformed_candidate_url_does_not_block_valid_fallback(self, tmp_path, monkeypatch):
        observed_urls = []

        def fake_download(source_url, output_path, timeout_seconds):
            observed_urls.append((source_url, timeout_seconds))
            with open(output_path, 'wb') as handle:
                handle.write(b'ok')

        monkeypatch.setattr('amalgkit.getfastq.download_with_urllib', fake_download)
        args = SimpleNamespace(
            sra_download_method='urllib',
            gcp_project='',
            sra_download_transfer_timeout_seconds=7,
        )

        downloaded = download_file_from_candidate_sources(
            sra_id='SRR001',
            source_candidates=[
                {'source_name': 'AWS', 'url': 'https://[invalid'},
                {'source_name': 'NCBI', 'url': 'https://example.org/valid.sra'},
            ],
            output_path=str(tmp_path / 'SRR001.sra'),
            args=args,
            artifact_label='SRA file',
        )

        assert downloaded is True
        assert observed_urls == [('https://example.org/valid.sra', 7.0)]
        assert (tmp_path / 'SRR001.sra').read_bytes() == b'ok'

    def test_urllib_empty_artifact_is_rejected(self, tmp_path, monkeypatch):
        def fake_empty_download(source_url, output_path, timeout_seconds):
            _ = (source_url, timeout_seconds)
            open(output_path, 'wb').close()

        monkeypatch.setattr('amalgkit.getfastq.download_with_urllib', fake_empty_download)
        downloaded = download_file_from_candidate_sources(
            sra_id='SRR001',
            source_candidates=[
                {'source_name': 'NCBI', 'url': 'https://example.org/empty.sra'},
            ],
            output_path=str(tmp_path / 'SRR001.sra'),
            args=SimpleNamespace(sra_download_method='urllib', gcp_project=''),
            artifact_label='SRA file',
        )

        assert downloaded is False
        assert not (tmp_path / 'SRR001.sra').exists()

    @pytest.mark.parametrize(
        'recoverable_error',
        [
            ConnectionResetError('connection reset'),
            ssl.SSLError('TLS failure'),
            http.client.IncompleteRead(b'partial', 100),
        ],
        ids=['connection-reset', 'ssl-error', 'incomplete-read'],
    )
    def test_recoverable_transport_error_tries_next_candidate(
        self,
        recoverable_error,
        tmp_path,
        monkeypatch,
    ):
        attempted_urls = []

        def fake_download(source_url, output_path, timeout_seconds):
            _ = timeout_seconds
            attempted_urls.append(source_url)
            if len(attempted_urls) == 1:
                with open(output_path, 'wb') as handle:
                    handle.write(b'partial')
                raise recoverable_error
            with open(output_path, 'wb') as handle:
                handle.write(b'complete')

        monkeypatch.setattr('amalgkit.getfastq.download_with_urllib', fake_download)
        output_path = tmp_path / 'SRR001.sra'

        downloaded = download_file_from_candidate_sources(
            sra_id='SRR001',
            source_candidates=[
                {'source_name': 'AWS', 'url': 'https://example.org/aws.sra'},
                {'source_name': 'NCBI', 'url': 'https://example.org/ncbi.sra'},
            ],
            output_path=str(output_path),
            args=SimpleNamespace(sra_download_method='urllib', gcp_project=''),
            artifact_label='SRA file',
        )

        assert downloaded is True
        assert attempted_urls == [
            'https://example.org/aws.sra',
            'https://example.org/ncbi.sra',
        ]
        assert output_path.read_bytes() == b'complete'
        assert list(tmp_path.glob('SRR001.sra.urllibtmp.*')) == []

    def test_empty_cached_sra_is_redownloaded(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        output_path = tmp_path / '{}.sra'.format(sra_id)
        output_path.write_bytes(b'')
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='https://example.org/run.sra',
        )
        args = self._make_args()

        def fake_download(source_url, output_path, timeout_seconds):
            _ = (source_url, timeout_seconds)
            with open(output_path, 'wb') as handle:
                handle.write(b'valid')

        monkeypatch.setattr('amalgkit.getfastq.download_with_urllib', fake_download)
        download_sra(
            metadata=metadata,
            sra_stat={'sra_id': sra_id},
            args=args,
            work_dir=str(tmp_path),
            overwrite=False,
        )

        assert output_path.read_bytes() == b'valid'

    def test_curl_timeout_and_empty_artifact_guard(self, tmp_path, monkeypatch):
        observed = {}
        monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda name: '/usr/bin/curl')

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            output_path = cmd[cmd.index('-o') + 1]
            open(output_path, 'wb').close()
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        output_path = tmp_path / 'SRR001.sra'
        downloaded = download_with_curl(
            source_url='https://ftp.sra.ebi.ac.uk/empty.sra',
            output_path=str(output_path),
            args=SimpleNamespace(
                dump_print=False,
                sra_download_transfer_timeout_seconds=123,
            ),
            sra_source_name='NCBI',
            artifact_label='SRA file',
        )

        assert downloaded is False
        assert observed['cmd'][observed['cmd'].index('--max-time') + 1] == '123'
        assert not output_path.exists()

    def test_candidate_slot_wait_has_a_bounded_timeout(self, tmp_path, monkeypatch):
        @contextmanager
        def always_busy(*_args, **_kwargs):
            yield True, None

        monotonic_values = iter([10.0, 12.0])
        monkeypatch.setattr('amalgkit.getfastq.maybe_acquire_source_download_slot', always_busy)
        monkeypatch.setattr('amalgkit.getfastq.time.monotonic', lambda: next(monotonic_values))
        monkeypatch.setattr(
            'amalgkit.getfastq.time.sleep',
            lambda _seconds: pytest.fail('timeout should be detected before sleeping'),
        )

        downloaded = download_file_from_candidate_sources(
            sra_id='SRR001',
            source_candidates=[
                {'source_name': 'NCBI', 'url': 'https://example.org/run.sra'},
            ],
            output_path=str(tmp_path / 'SRR001.sra'),
            args=SimpleNamespace(sra_download_wait_timeout_seconds=1),
            artifact_label='SRA file',
        )

        assert downloaded is False

    def test_single_source_api_slot_wait_has_a_bounded_timeout(self, tmp_path, monkeypatch):
        acquire_wait_values = []

        @contextmanager
        def always_busy(*_args, **kwargs):
            acquire_wait_values.append(kwargs.get('wait'))
            yield True, None

        monotonic_values = iter([20.0, 22.0])
        monkeypatch.setattr('amalgkit.getfastq.maybe_acquire_source_download_slot', always_busy)
        monkeypatch.setattr('amalgkit.getfastq.time.monotonic', lambda: next(monotonic_values))
        monkeypatch.setattr(
            'amalgkit.getfastq.time.sleep',
            lambda _seconds: pytest.fail('timeout should be detected before sleeping'),
        )

        downloaded = download_file_from_source(
            sra_id='SRR001',
            sra_source_name='NCBI',
            source_url_original='https://example.org/run.sra',
            output_path=str(tmp_path / 'SRR001.sra'),
            args=SimpleNamespace(sra_download_wait_timeout_seconds=1),
            artifact_label='SRA file',
        )

        assert downloaded is False
        assert acquire_wait_values == [False]

    def test_raises_when_existing_sra_path_is_directory(self, tmp_path):
        sra_id = 'SRR_DIR'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='https://example.invalid/path/to.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args()
        (tmp_path / '{}.sra'.format(sra_id)).mkdir()

        with pytest.raises(IsADirectoryError, match='SRA path exists but is not a file'):
            download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

    def test_skips_non_http_schemes_before_urllib(self, tmp_path, monkeypatch, capsys):
        sra_id = 'SRR_SCHEME'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='s3://bucket/path/to.sra',
            gcp_link='gs://bucket/path/to.sra',
            ncbi_link='https://example.invalid/path/to.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args()
        called_urls = []

        def fake_urlretrieve(url, _path):
            called_urls.append(url)
            raise urllib.error.URLError('network down')

        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )

        with pytest.raises(FileNotFoundError):
            download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert called_urls == [
            'https://storage.googleapis.com/bucket/path/to.sra',
            'https://example.invalid/path/to.sra',
        ]
        stderr = capsys.readouterr().err
        assert 'unsupported URL scheme for urllib: s3' in stderr

    def test_downloads_when_http_source_succeeds(self, tmp_path, monkeypatch):
        sra_id = 'SRR_OK'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='https://example.invalid/path/to.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args()
        out_path = tmp_path / '{}.sra'.format(sra_id)

        def fake_urlretrieve(_url, path):
            with open(path, 'w') as f:
                f.write('ok')
            return (path, None)

        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )
        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)
        assert out_path.exists()

    def test_appends_user_project_for_gcp_requester_pays(self, tmp_path, monkeypatch):
        sra_id = 'SRR_GCP_PROJECT'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='gs://bucket/path/to.sra',
            ncbi_link='',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args(gcp_project='test-project')
        called_urls = []

        def fake_urlretrieve(url, _path):
            called_urls.append(url)
            raise urllib.error.URLError('network down')

        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )

        with pytest.raises(FileNotFoundError):
            download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert called_urls == ['https://storage.googleapis.com/bucket/path/to.sra?userProject=test-project']

    def test_uses_curl_when_requested(self, tmp_path, monkeypatch):
        sra_id = 'SRR_CURL'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='https://sra-downloadb.be-md.ncbi.nlm.nih.gov/path/to.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args(sra_download_method='curl')
        observed = {'cmd': None}

        monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda exe: '/usr/bin/curl' if exe == 'curl' else None)

        def fail_urlretrieve(*_args, **_kwargs):
            raise AssertionError('urllib.request.urlretrieve should not be called when curl succeeds.')

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            out_path = cmd[cmd.index('-o') + 1]
            with open(out_path, 'w') as fh:
                fh.write('ok')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fail_urlretrieve(source_url, output_path),
        )
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert observed['cmd'][0] == '/usr/bin/curl'
        assert (tmp_path / '{}.sra'.format(sra_id)).exists()

    def test_acquires_ncbi_download_semaphore_for_ncbi_source(self, tmp_path, monkeypatch):
        sra_id = 'SRR_LIMITED'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='https://example.invalid/path/to.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args()
        args.download_lock_dir = str(tmp_path / 'locks')
        args.ncbi_download_max_concurrency = 4
        args.aws_download_max_concurrency = 'auto'
        args.gcp_download_max_concurrency = 'auto'
        observed = {}

        @contextmanager
        def fake_maybe_acquire(args, limit_attr, semaphore_name, lock_label, resolve_download_dir_fn=None, wait=True):
            observed['limit_attr'] = limit_attr
            observed['semaphore_name'] = semaphore_name
            observed['lock_label'] = lock_label
            observed['max_concurrency'] = args.ncbi_download_max_concurrency
            observed['wait'] = wait
            yield 'slot'

        def fake_urlretrieve(_url, path):
            with open(path, 'w') as f:
                f.write('ok')
            return (path, None)

        monkeypatch.setattr('amalgkit.getfastq.maybe_acquire_download_semaphore', fake_maybe_acquire)
        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )

        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert observed == {
            'limit_attr': 'ncbi_download_max_concurrency',
            'semaphore_name': 'ncbi_download',
            'lock_label': 'NCBI download',
            'max_concurrency': 4,
            'wait': False,
        }

    def test_tries_next_source_when_first_source_slot_is_unavailable(self, tmp_path, monkeypatch):
        sra_id = 'SRR_FALLTHROUGH'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='https://example.invalid/aws.sra',
            gcp_link='',
            ncbi_link='https://example.invalid/ncbi.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args()
        args.download_lock_dir = str(tmp_path / 'locks')
        args.aws_download_max_concurrency = 1
        args.ncbi_download_max_concurrency = 1
        args.gcp_download_max_concurrency = 'auto'
        observed = {
            'acquire_calls': [],
            'download_urls': [],
        }

        @contextmanager
        def fake_maybe_acquire(args, limit_attr, semaphore_name, lock_label, resolve_download_dir_fn=None, wait=True):
            observed['acquire_calls'].append((limit_attr, semaphore_name, wait))
            if limit_attr == 'aws_download_max_concurrency':
                yield None
                return
            yield 'slot'

        def fake_urlretrieve(url, path):
            observed['download_urls'].append(url)
            with open(path, 'w') as f:
                f.write('ok')
            return (path, None)

        monkeypatch.setattr('amalgkit.getfastq.maybe_acquire_download_semaphore', fake_maybe_acquire)
        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )

        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert observed['acquire_calls'] == [
            ('aws_download_max_concurrency', 'aws_download', False),
            ('ncbi_download_max_concurrency', 'ncbi_download', False),
        ]
        assert observed['download_urls'] == ['https://example.invalid/ncbi.sra']
        assert (tmp_path / '{}.sra'.format(sra_id)).exists()

    def test_tries_ena_when_higher_priority_sources_are_at_concurrency_limit(self, tmp_path, monkeypatch):
        sra_id = 'SRR123456'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='https://example.invalid/aws.sra',
            gcp_link='https://example.invalid/gcp.sra',
            ncbi_link='https://example.invalid/ncbi.sra',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args(ena=True)
        args.download_lock_dir = str(tmp_path / 'locks')
        args.aws_download_max_concurrency = 1
        args.gcp_download_max_concurrency = 1
        args.ncbi_download_max_concurrency = 1
        args.ena_download_max_concurrency = 1
        observed = {
            'acquire_calls': [],
            'download_urls': [],
        }

        @contextmanager
        def fake_maybe_acquire(args, limit_attr, semaphore_name, lock_label, resolve_download_dir_fn=None, wait=True):
            observed['acquire_calls'].append((limit_attr, semaphore_name, wait))
            if limit_attr in {
                'aws_download_max_concurrency',
                'gcp_download_max_concurrency',
                'ncbi_download_max_concurrency',
            }:
                yield None
                return
            yield 'slot'

        def fake_urlretrieve(url, path):
            observed['download_urls'].append(url)
            with open(path, 'w') as f:
                f.write('ok')
            return (path, None)

        monkeypatch.setattr('amalgkit.getfastq.maybe_acquire_download_semaphore', fake_maybe_acquire)
        monkeypatch.setattr(
            'amalgkit.getfastq.fetch_ena_run_file_report',
            lambda run_accession, timeout: {
                'sra_urls': [
                    'https://ftp.sra.ebi.ac.uk/vol1/srr/SRR123/SRR123456/SRR123456.sra'
                ],
                'fastq_urls': [],
            },
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )

        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert observed['acquire_calls'] == [
            ('aws_download_max_concurrency', 'aws_download', False),
            ('gcp_download_max_concurrency', 'gcp_download', False),
            ('ncbi_download_max_concurrency', 'ncbi_download', False),
            ('ena_download_max_concurrency', 'ena_download', False),
        ]
        assert observed['download_urls'] == [
            'https://ftp.sra.ebi.ac.uk/vol1/srr/SRR123/SRR123456/SRR123456.sra'
        ]
        assert (tmp_path / '{}.sra'.format(sra_id)).exists()

    def test_tries_ddbj_when_higher_priority_sources_are_at_concurrency_limit(self, tmp_path, monkeypatch):
        sra_id = 'DRR000001'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='https://example.invalid/aws.sra',
            gcp_link='https://example.invalid/gcp.sra',
            ncbi_link='https://example.invalid/ncbi.sra',
            experiment='DRX000001',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args(ddbj=True)
        args.download_lock_dir = str(tmp_path / 'locks')
        args.aws_download_max_concurrency = 1
        args.gcp_download_max_concurrency = 1
        args.ncbi_download_max_concurrency = 1
        args.ddbj_download_max_concurrency = 1
        observed = {
            'acquire_calls': [],
            'download_urls': [],
        }

        @contextmanager
        def fake_maybe_acquire(args, limit_attr, semaphore_name, lock_label, resolve_download_dir_fn=None, wait=True):
            observed['acquire_calls'].append((limit_attr, semaphore_name, wait))
            if limit_attr in {
                'aws_download_max_concurrency',
                'gcp_download_max_concurrency',
                'ncbi_download_max_concurrency',
            }:
                yield None
                return
            yield 'slot'

        def fake_urlretrieve(url, path):
            observed['download_urls'].append(url)
            with open(path, 'w') as f:
                f.write('ok')
            return (path, None)

        monkeypatch.setattr('amalgkit.getfastq.maybe_acquire_download_semaphore', fake_maybe_acquire)
        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )

        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert observed['acquire_calls'] == [
            ('aws_download_max_concurrency', 'aws_download', False),
            ('gcp_download_max_concurrency', 'gcp_download', False),
            ('ncbi_download_max_concurrency', 'ncbi_download', False),
            ('ddbj_download_max_concurrency', 'ddbj_download', False),
        ]
        assert observed['download_urls'] == [
            'https://ddbj.nig.ac.jp/public/ddbj_database/dra/sra/ByExp/sra/DRX/'
            'DRX000/DRX000001/DRR000001/DRR000001.sra'
        ]
        assert (tmp_path / '{}.sra'.format(sra_id)).exists()

    def test_waits_until_any_deferred_source_slot_opens(self, tmp_path, monkeypatch):
        args = SimpleNamespace(
            aws_download_max_concurrency=1,
            ena_download_max_concurrency=1,
        )
        observed = {
            'acquire_calls': [],
            'download_sources': [],
            'sleep_calls': 0,
        }
        retry_state = {'allow_ena': False}

        @contextmanager
        def fake_maybe_acquire(args, limit_attr, semaphore_name, lock_label, resolve_download_dir_fn=None, wait=True):
            observed['acquire_calls'].append((limit_attr, wait, retry_state['allow_ena']))
            if limit_attr == 'ena_download_max_concurrency' and retry_state['allow_ena']:
                yield 'slot'
                return
            yield None

        def fake_download(*, sra_id, sra_source_name, source_url_original, output_path, args, artifact_label):
            observed['download_sources'].append(sra_source_name)
            return True

        def fake_sleep(_seconds):
            observed['sleep_calls'] += 1
            retry_state['allow_ena'] = True

        monkeypatch.setattr('amalgkit.getfastq.maybe_acquire_download_semaphore', fake_maybe_acquire)
        monkeypatch.setattr('amalgkit.getfastq._download_file_from_source_without_semaphore', fake_download)
        monkeypatch.setattr('amalgkit.getfastq.time.sleep', fake_sleep)

        is_downloaded = download_file_from_candidate_sources(
            sra_id='SRR_WAIT',
            source_candidates=[
                {'source_name': 'AWS', 'url': 'https://example.invalid/aws.sra'},
                {'source_name': 'ENA', 'url': 'https://ftp.sra.ebi.ac.uk/vol1/srr/SRR000/SRR_WAIT/SRR_WAIT.sra'},
            ],
            output_path=str(tmp_path / 'SRR_WAIT.sra'),
            args=args,
            artifact_label='SRA file',
        )

        assert is_downloaded is True
        assert observed['sleep_calls'] == 1
        assert observed['download_sources'] == ['ENA']

    def test_uses_ena_filereport_sra_url_when_enabled(self, tmp_path, monkeypatch):
        sra_id = 'SRR000001'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args(ena=True)
        args.aws = False
        args.gcp = False
        args.ncbi = False
        called_urls = []

        def fake_urlretrieve(url, path):
            called_urls.append(url)
            with open(path, 'w') as f:
                f.write('ok')
            return (path, None)

        monkeypatch.setattr(
            'amalgkit.getfastq.fetch_ena_run_file_report',
            lambda run_accession, timeout: {
                'sra_urls': [
                    'https://ftp.sra.ebi.ac.uk/vol1/srr/SRR000/SRR000001/SRR000001.sra'
                ],
                'fastq_urls': [],
            },
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )

        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert called_urls == ['https://ftp.sra.ebi.ac.uk/vol1/srr/SRR000/SRR000001/SRR000001.sra']
        assert metadata.df.loc[0, 'ENA_SRA_Link'] == called_urls[0]
        assert (tmp_path / '{}.sra'.format(sra_id)).exists()

    def test_uses_ena_filereport_when_link_column_is_float_nan(self, tmp_path, monkeypatch):
        sra_id = 'SRR000001'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link=numpy.nan,
            gcp_link=numpy.nan,
            ncbi_link=numpy.nan,
            ena_sra_link=numpy.nan,
        )
        metadata.df['AWS_Link'] = pandas.Series([numpy.nan], dtype='float64')
        metadata.df['GCP_Link'] = pandas.Series([numpy.nan], dtype='float64')
        metadata.df['NCBI_Link'] = pandas.Series([numpy.nan], dtype='float64')
        metadata.df['ENA_SRA_Link'] = pandas.Series([numpy.nan], dtype='float64')
        sra_stat = {'sra_id': sra_id}
        args = self._make_args(ena=True)
        args.aws = True
        args.gcp = True
        args.ncbi = True
        called_urls = []

        def fake_urlretrieve(url, path):
            called_urls.append(url)
            with open(path, 'w') as f:
                f.write('ok')
            return (path, None)

        monkeypatch.setattr(
            'amalgkit.getfastq.fetch_ena_run_file_report',
            lambda run_accession, timeout: {
                'sra_urls': [
                    'https://ftp.sra.ebi.ac.uk/vol1/srr/SRR000/SRR000001/SRR000001.sra'
                ],
                'fastq_urls': [],
            },
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )

        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert called_urls == ['https://ftp.sra.ebi.ac.uk/vol1/srr/SRR000/SRR000001/SRR000001.sra']
        assert metadata.df.loc[0, 'ENA_SRA_Link'] == called_urls[0]
        assert metadata.df['ENA_SRA_Link'].dtype == object
        assert (tmp_path / '{}.sra'.format(sra_id)).exists()

    def test_does_not_reuse_cached_ena_link_absent_from_filereport(self, monkeypatch):
        sra_id = 'SRR000001'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='',
            ena_sra_link=(
                'https://ftp.sra.ebi.ac.uk/vol1/srr/SRR000/'
                'SRR000001/SRR000001.sra'
            ),
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args(ena=True)
        args.aws = False
        args.gcp = False
        args.ncbi = False
        args.ddbj = False
        monkeypatch.setattr(
            'amalgkit.getfastq.fetch_ena_run_file_report',
            lambda run_accession, timeout: {
                'sra_urls': [],
                'fastq_urls': [],
            },
        )

        sources = resolve_sra_download_sources(
            metadata=metadata,
            sra_stat=sra_stat,
            args=args,
        )

        assert sources == {}
        assert metadata.df.loc[0, 'ENA_SRA_Link'] == ''

    def test_ena_fastq_only_report_reaches_original_fastq_fallback(self, tmp_path, monkeypatch):
        sra_id = 'SRR000001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'experiment': ['SRX000001'],
            'AWS_Link': [''],
            'GCP_Link': [''],
            'NCBI_Link': [''],
            'ENA_SRA_Link': [''],
            'DDBJ_SRA_Link': [''],
            'lib_layout': ['paired'],
            'total_spots': [2],
            'total_bases': [16],
            'size': [16],
            'nominal_length': [4],
            'nominal_sdev': [0],
            'spot_length': [8],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        metadata = initialize_columns(metadata, {'num_bp_per_sra': 16})
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'spot_length': 8,
            'total_spot': 2,
            'getfastq_sra_dir': str(tmp_path),
            'metadata_idx': 0,
        }
        args = self._make_args(ena=True)
        args.aws = False
        args.gcp = False
        args.ncbi = False
        args.ddbj = False
        args.min_read_length = 1
        args.threads = 1
        args.fasterq_size_check = True
        args.fasterq_disk_limit = None
        args.fasterq_disk_limit_tmp = None
        args.fastp = False
        args.rrna_filter = False
        args.contam_filter = False
        args.dump_print = False

        fastq_urls = [
            'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/SRR000001/SRR000001_1.fastq.gz',
            'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/SRR000001/SRR000001_2.fastq.gz',
        ]
        monkeypatch.setattr(
            'amalgkit.getfastq.fetch_ena_run_file_report',
            lambda run_accession, timeout: {
                'sra_urls': [],
                'fastq_urls': fastq_urls,
            },
        )

        def fake_download(source_url, output_path, timeout_seconds):
            _ = (source_url, timeout_seconds)
            with gzip.open(output_path, 'wt') as handle:
                for spot in [1, 2]:
                    handle.write('@{}.{}\n'.format(sra_id, spot))
                    handle.write('ACGT\n+\nIIII\n')

        monkeypatch.setattr('amalgkit.getfastq.download_with_urllib', fake_download)
        monkeypatch.setattr(
            'amalgkit.getfastq.run_fasterq_dump_with_retry',
            lambda **_kwargs: pytest.fail('fasterq-dump must not run without an SRA file'),
        )

        result = download_sra(
            metadata=metadata,
            sra_stat=sra_stat,
            args=args,
            work_dir=str(tmp_path),
            overwrite=False,
        )
        metadata, sra_stat_out = run_fasterq_dump(
            sra_stat=sra_stat,
            args=args,
            metadata=metadata,
            start=1,
            end=2,
        )

        assert result == 'original-fastq'
        assert not (tmp_path / '{}.sra'.format(sra_id)).exists()
        assert count_fastq_records(str(tmp_path / '{}_1.fastq.gz'.format(sra_id))) == 2
        assert count_fastq_records(str(tmp_path / '{}_2.fastq.gz'.format(sra_id))) == 2
        assert metadata.df.loc[0, 'num_dumped'] == 2
        assert metadata.df.loc[0, 'num_written'] == 2
        assert sra_stat_out['layout'] == 'paired'

    def test_derives_ddbj_sra_url_when_enabled_for_drr(self, tmp_path, monkeypatch):
        sra_id = 'DRR000001'
        metadata = self._make_metadata(
            sra_id=sra_id,
            aws_link='',
            gcp_link='',
            ncbi_link='',
            experiment='DRX000001',
        )
        sra_stat = {'sra_id': sra_id}
        args = self._make_args(ddbj=True)
        args.aws = False
        args.gcp = False
        args.ncbi = False
        called_urls = []

        def fake_urlretrieve(url, path):
            called_urls.append(url)
            with open(path, 'w') as f:
                f.write('ok')
            return (path, None)

        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )

        download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert called_urls == [
            'https://ddbj.nig.ac.jp/public/ddbj_database/dra/sra/ByExp/sra/DRX/'
            'DRX000/DRX000001/DRR000001/DRR000001.sra'
        ]
        assert metadata.df.loc[0, 'DDBJ_SRA_Link'] == called_urls[0]
        assert (tmp_path / '{}.sra'.format(sra_id)).exists()

    def test_apply_first_round_results_writes_dynamic_sra_links_to_float_columns(self):
        sra_id = 'SRR000001'
        ena_url = 'https://ftp.sra.ebi.ac.uk/vol1/srr/SRR000/SRR000001/SRR000001.sra'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['paired'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
            'ENA_SRA_Link': pandas.Series([numpy.nan], dtype='float64'),
        }))
        row_series = metadata.df.iloc[0].copy()
        row_series['ENA_SRA_Link'] = ena_url
        run_results_by_id = {
            sra_id: {
                'row': row_series,
                'flag_any_output_file_present': False,
                'flag_private_file': False,
                'getfastq_sra_dir': '/tmp/getfastq/{}'.format(sra_id),
            },
        }

        metadata, flag_private_file, flag_any_output_file_present, last_getfastq_sra_dir = (
            apply_first_round_getfastq_results(metadata, [(0, sra_id)], run_results_by_id)
        )

        assert metadata.df.loc[0, 'ENA_SRA_Link'] == ena_url
        assert metadata.df['ENA_SRA_Link'].dtype == object
        assert flag_private_file is False
        assert flag_any_output_file_present is False
        assert last_getfastq_sra_dir == '/tmp/getfastq/{}'.format(sra_id)

def test_download_with_curl_rejects_disallowed_host(tmp_path, monkeypatch):
    monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda name: '/usr/bin/curl')
    with pytest.raises(ValueError, match='not an allowed SRA/ENA/cloud download endpoint'):
        download_with_curl(
            source_url='http://169.254.169.254/latest/meta-data/',
            output_path=str(tmp_path / 'x.sra'),
            args=SimpleNamespace(sra_download_transfer_timeout_seconds=60),
            sra_source_name='NCBI',
            artifact_label='SRA file',
        )

def test_download_with_curl_rejects_redirect_to_disallowed_host(tmp_path, monkeypatch):
    monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda name: '/usr/bin/curl')
    requested_urls = []

    def fake_run(cmd, stdout=None, stderr=None):
        requested_urls.append(cmd[-1])
        out_path = cmd[cmd.index('-o') + 1]
        with open(out_path, 'wb') as fh:
            fh.write(b'data')
        return subprocess.CompletedProcess(
            cmd,
            0,
            stdout=(
                b'AMALGKIT_CURL_STATUS:302\n'
                b'AMALGKIT_CURL_REDIRECT:http://169.254.169.254/latest/meta-data/\n'
            ),
            stderr=b'',
        )

    monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
    output_path = tmp_path / 'SRR001.sra'
    downloaded = download_with_curl(
        source_url='https://ftp.sra.ebi.ac.uk/run.sra',
        output_path=str(output_path),
        args=SimpleNamespace(sra_download_transfer_timeout_seconds=60),
        sra_source_name='ENA',
        artifact_label='SRA file',
    )
    assert downloaded is False
    assert requested_urls == ['https://ftp.sra.ebi.ac.uk/run.sra']
    assert not output_path.exists()


def test_download_with_curl_validates_each_allowed_redirect_before_connecting(tmp_path, monkeypatch):
    monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda name: '/usr/bin/curl')
    requested_urls = []

    def fake_run(cmd, stdout=None, stderr=None):
        requested_urls.append(cmd[-1])
        out_path = cmd[cmd.index('-o') + 1]
        with open(out_path, 'wb') as handle:
            handle.write(b'final' if len(requested_urls) == 2 else b'redirect')
        if len(requested_urls) == 1:
            response = (
                b'AMALGKIT_CURL_STATUS:302\n'
                b'AMALGKIT_CURL_REDIRECT:https://sra-download.be-md.ncbi.nlm.nih.gov/final.sra\n'
            )
        else:
            response = b'AMALGKIT_CURL_STATUS:200\nAMALGKIT_CURL_REDIRECT:\n'
        return subprocess.CompletedProcess(cmd, 0, stdout=response, stderr=b'')

    monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
    output_path = tmp_path / 'SRR001.sra'
    assert download_with_curl(
        source_url='https://ftp.sra.ebi.ac.uk/run.sra',
        output_path=str(output_path),
        args=SimpleNamespace(sra_download_transfer_timeout_seconds=60),
        sra_source_name='ENA',
        artifact_label='SRA file',
    )
    assert requested_urls == [
        'https://ftp.sra.ebi.ac.uk/run.sra',
        'https://sra-download.be-md.ncbi.nlm.nih.gov/final.sra',
    ]
    assert output_path.read_bytes() == b'final'


def test_download_with_curl_resumes_only_a_matching_http_byte_range(tmp_path, monkeypatch):
    monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda name: '/usr/bin/curl')
    output_path = tmp_path / 'SRR001.fastq.gz'
    output_path.write_bytes(b'ab')

    def fake_run(cmd, stdout=None, stderr=None):
        assert cmd[cmd.index('--continue-at') + 1] == '-'
        header_path = cmd[cmd.index('--dump-header') + 1]
        with open(header_path, 'wb') as handle:
            handle.write(b'HTTP/1.1 206 Partial Content\r\nContent-Range: bytes 2-5/6\r\n\r\n')
        with open(cmd[cmd.index('-o') + 1], 'ab') as handle:
            handle.write(b'cdef')
        return subprocess.CompletedProcess(
            cmd,
            0,
            stdout=b'AMALGKIT_CURL_STATUS:206\nAMALGKIT_CURL_REDIRECT:\n',
            stderr=b'',
        )

    monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
    assert download_with_curl(
        source_url='https://ftp.sra.ebi.ac.uk/run.fastq.gz',
        output_path=str(output_path),
        args=SimpleNamespace(sra_download_transfer_timeout_seconds=60),
        sra_source_name='ENA',
        artifact_label='Original FASTQ file',
        resume_existing=True,
    )
    assert output_path.read_bytes() == b'abcdef'


def test_download_with_curl_restarts_when_server_ignores_range(tmp_path, monkeypatch):
    monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda name: '/usr/bin/curl')
    output_path = tmp_path / 'SRR001.fastq.gz'
    output_path.write_bytes(b'ab')
    commands = []

    def fake_run(cmd, stdout=None, stderr=None):
        commands.append(list(cmd))
        header_path = cmd[cmd.index('--dump-header') + 1]
        with open(header_path, 'wb') as handle:
            handle.write(b'HTTP/1.1 200 OK\r\n\r\n')
        if len(commands) == 1:
            return subprocess.CompletedProcess(
                cmd,
                33,
                stdout=b'AMALGKIT_CURL_STATUS:200\nAMALGKIT_CURL_REDIRECT:\n',
                stderr=b'',
            )
        with open(cmd[cmd.index('-o') + 1], 'wb') as handle:
            handle.write(b'fresh')
        return subprocess.CompletedProcess(
            cmd,
            0,
            stdout=b'AMALGKIT_CURL_STATUS:200\nAMALGKIT_CURL_REDIRECT:\n',
            stderr=b'',
        )

    monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
    assert download_with_curl(
        source_url='https://ftp.sra.ebi.ac.uk/run.fastq.gz',
        output_path=str(output_path),
        args=SimpleNamespace(sra_download_transfer_timeout_seconds=60),
        sra_source_name='ENA',
        artifact_label='Original FASTQ file',
        resume_existing=True,
    )
    assert '--continue-at' in commands[0]
    assert '--continue-at' not in commands[1]
    assert output_path.read_bytes() == b'fresh'


def test_download_with_curl_retains_partial_after_transfer_failure(tmp_path, monkeypatch):
    monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda name: '/usr/bin/curl')
    output_path = tmp_path / 'SRR001.fastq.gz'
    output_path.write_bytes(b'ab')

    def fake_run(cmd, stdout=None, stderr=None):
        header_path = cmd[cmd.index('--dump-header') + 1]
        with open(header_path, 'wb') as handle:
            handle.write(b'HTTP/1.1 206 Partial Content\r\nContent-Range: bytes 2-2/6\r\n\r\n')
        with open(cmd[cmd.index('-o') + 1], 'ab') as handle:
            handle.write(b'c')
        return subprocess.CompletedProcess(
            cmd,
            28,
            stdout=b'AMALGKIT_CURL_STATUS:206\nAMALGKIT_CURL_REDIRECT:\n',
            stderr=b'',
        )

    monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
    assert not download_with_curl(
        source_url='https://ftp.sra.ebi.ac.uk/run.fastq.gz',
        output_path=str(output_path),
        args=SimpleNamespace(sra_download_transfer_timeout_seconds=60),
        sra_source_name='ENA',
        artifact_label='Original FASTQ file',
        resume_existing=True,
    )
    assert output_path.read_bytes() == b'abc'


def test_original_fastq_resume_directory_survives_interruption(tmp_path):
    args = SimpleNamespace(sra_download_wait_timeout_seconds=5)
    with pytest.raises(RuntimeError, match='interrupt'):
        with _public_fastq_resume_directory(str(tmp_path), args) as resume_dir:
            payload = os.path.join(resume_dir, 'SRR001_1.fastq.gz')
            with open(payload, 'wb') as handle:
                handle.write(b'partial')
            raise RuntimeError('interrupt')
    with _public_fastq_resume_directory(str(tmp_path), args) as resume_dir:
        with open(os.path.join(resume_dir, 'SRR001_1.fastq.gz'), 'rb') as handle:
            assert handle.read() == b'partial'


def test_paired_original_fastq_retry_reuses_the_verified_first_mate(tmp_path, monkeypatch):
    sra_id = 'SRR001'
    download_calls = []
    mate2_attempts = {'count': 0}

    def fake_download(**kwargs):
        output_path = kwargs['output_path']
        download_calls.append(os.path.basename(output_path))
        if output_path.endswith('_2.fastq.gz'):
            mate2_attempts['count'] += 1
            if mate2_attempts['count'] == 1:
                return False
            payload = b'pair'
        else:
            payload = b'mate'
        with open(output_path, 'wb') as handle:
            handle.write(payload)
        return True

    monkeypatch.setattr('amalgkit.getfastq.download_file_from_candidate_sources', fake_download)
    monkeypatch.setattr('amalgkit.getfastq.validate_original_fastq_download', lambda path, entry: None)
    entries = [
        {
            'filename': '{}_1.fastq.gz'.format(sra_id),
            'expected_bytes': 4,
            'expected_md5': None,
            'sources': [{'source_name': 'ENA', 'url': 'https://ftp.sra.ebi.ac.uk/read1'}],
        },
        {
            'filename': '{}_2.fastq.gz'.format(sra_id),
            'expected_bytes': 4,
            'expected_md5': None,
            'sources': [{'source_name': 'ENA', 'url': 'https://ftp.sra.ebi.ac.uk/read2'}],
        },
    ]
    sra_stat = {
        'sra_id': sra_id,
        'getfastq_sra_dir': str(tmp_path),
        'total_spot': 1,
    }
    args = SimpleNamespace(sra_download_wait_timeout_seconds=5)

    assert download_public_original_fastq_files(
        sra_stat=sra_stat,
        args=args,
        start=1,
        end=1,
        public_original_fastqs=entries,
        source_description='test',
    ) is None
    assert (tmp_path / '.amalgkit-original-fastq-resume' / '{}_1.fastq.gz'.format(sra_id)).read_bytes() == b'mate'

    assert download_public_original_fastq_files(
        sra_stat=sra_stat,
        args=args,
        start=1,
        end=1,
        public_original_fastqs=entries,
        source_description='test',
    ) is not None
    assert download_calls.count('{}_1.fastq.gz'.format(sra_id)) == 1
    assert (tmp_path / '{}_1.fastq.gz'.format(sra_id)).read_bytes() == b'mate'
    assert (tmp_path / '{}_2.fastq.gz'.format(sra_id)).read_bytes() == b'pair'
    assert not (tmp_path / '.amalgkit-original-fastq-resume').exists()


def test_verified_original_fastq_continues_short_partial(tmp_path, monkeypatch):
    output_path = tmp_path / 'SRR001_1.fastq.gz'
    output_path.write_bytes(b'ab')
    observed = {}

    def fake_download(**kwargs):
        observed.update(kwargs)
        with open(kwargs['output_path'], 'ab') as handle:
            handle.write(b'cdef')
        return True

    monkeypatch.setattr('amalgkit.getfastq.download_file_from_candidate_sources', fake_download)
    monkeypatch.setattr('amalgkit.getfastq.validate_original_fastq_download', lambda path, entry: None)
    download_verified_original_fastq(
        'SRR001',
        {
            'expected_bytes': 6,
            'expected_md5': None,
            'sources': [{'source_name': 'ENA', 'url': 'https://ftp.sra.ebi.ac.uk/run.fastq.gz'}],
        },
        str(output_path),
        SimpleNamespace(),
    )
    assert observed['resume_existing'] is True
    assert output_path.read_bytes() == b'abcdef'


def test_download_with_urllib_refuses_preexisting_symlink(tmp_path):
    from io import BytesIO

    from amalgkit.getfastq import download_with_urllib

    victim_path = tmp_path / 'victim'
    victim_path.write_bytes(b'original')
    output_path = tmp_path / 'output'
    output_path.symlink_to(victim_path)

    class _Response(BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            self.close()

    with pytest.raises(FileExistsError):
        download_with_urllib(
            source_url='https://ftp.sra.ebi.ac.uk/run.sra',
            output_path=str(output_path),
            timeout_seconds=30,
            urlopen_fn=lambda _url, timeout: _Response(b'attacker-controlled'),
        )

    assert victim_path.read_bytes() == b'original'


def test_download_temp_path_is_inside_private_unpredictable_directory(tmp_path):
    from amalgkit.getfastq import _secure_download_temp_path

    output = str(tmp_path / 'SRR001.sra')
    with _secure_download_temp_path(output, '.urllibtmp') as first_path:
        first_parent = os.path.dirname(first_path)
        assert os.path.basename(first_path) == 'payload'
        assert not os.path.exists(first_path)
        assert stat.S_IMODE(os.stat(first_parent).st_mode) == 0o700
    assert not os.path.exists(first_parent)

    with _secure_download_temp_path(output, '.urllibtmp') as second_path:
        assert os.path.dirname(second_path) != first_parent
