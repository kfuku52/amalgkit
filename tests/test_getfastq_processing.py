import pytest
import pandas
import numpy
import gzip
import subprocess
import xml.etree.ElementTree as ET

import os
from types import SimpleNamespace


from amalgkit.getfastq import (
    initialize_columns,
    get_identical_paired_ratio,
    maybe_treat_paired_as_single,
    parse_fastp_metrics,
    parse_fastp_summary_counts,
    update_fastp_metrics,
    write_fastp_stats,
    write_getfastq_stats,
    print_read_stats,
    run_fastp,
    guard_fasterq_full_dump_disk_space,
    sequence_extraction,
    estimate_num_written_spots_from_fastq,
    process_getfastq_run,
    run_fasterq_dump,
    count_fastq_records,
    remove_sra_files,
    remove_sra_path,
    resolve_public_original_fastq_sources_from_xml_root,
    assign_public_original_fastq_suffixes,
    GETFASTQ_PHASE_COMPLETE,
)
from amalgkit.util import Metadata



class TestFastpMetrics:
    def test_parse_fastp_metrics_extracts_duplication_and_insert_size(self):
        stderr_txt = '\n'.join([
            'Filtering result:',
            'reads passed filter: 274442354',
            'Duplication rate: 22.9565%',
            'Insert size peak (evaluated by paired-end reads): 138',
            'JSON report: /dev/null',
        ])
        duplication_rate, insert_size_peak = parse_fastp_metrics(stderr_txt)
        assert duplication_rate == pytest.approx(22.9565)
        assert insert_size_peak == pytest.approx(138.0)

    def test_parse_fastp_metrics_returns_nan_when_not_present(self):
        duplication_rate, insert_size_peak = parse_fastp_metrics('no matching lines')
        assert numpy.isnan(duplication_rate)
        assert numpy.isnan(insert_size_peak)

    def test_parse_fastp_summary_counts_raises_on_truncated_section(self):
        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            # missing total bases line
        ])
        with pytest.raises(RuntimeError, match='Unexpected fastp stderr format'):
            parse_fastp_summary_counts(stderr_txt)

    def test_update_fastp_metrics_weighted_average(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_fastp_in': [100],
            'fastp_duplication_rate': [10.0],
            'fastp_insert_size_peak': [200.0],
        }))
        update_fastp_metrics(
            metadata=metadata,
            ind_sra=0,
            current_num_in=300,
            duplication_rate=20.0,
            insert_size_peak=100.0,
        )
        assert metadata.df.loc[0, 'fastp_duplication_rate'] == pytest.approx(17.5)
        assert metadata.df.loc[0, 'fastp_insert_size_peak'] == pytest.approx(125.0)

    def test_update_fastp_metrics_sets_initial_values(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_fastp_in': [0],
            'fastp_duplication_rate': [numpy.nan],
            'fastp_insert_size_peak': [numpy.nan],
        }))
        update_fastp_metrics(
            metadata=metadata,
            ind_sra=0,
            current_num_in=50,
            duplication_rate=33.0,
            insert_size_peak=222.0,
        )
        assert metadata.df.loc[0, 'fastp_duplication_rate'] == pytest.approx(33.0)
        assert metadata.df.loc[0, 'fastp_insert_size_peak'] == pytest.approx(222.0)

    def test_write_fastp_stats_writes_getfastq_stats_only(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'fastp_duplication_rate': [22.9565],
            'fastp_insert_size_peak': [138.0],
            'num_fastp_in': [1000],
            'num_fastp_out': [900],
            'bp_fastp_in': [100000],
            'bp_fastp_out': [90000],
        }))
        sra_stat = {'sra_id': 'SRR001'}
        write_fastp_stats(sra_stat=sra_stat, metadata=metadata, output_dir=str(tmp_path))
        out_path = tmp_path / 'getfastq_stats.tsv'
        assert out_path.exists()
        assert not (tmp_path / 'fastp_stats.tsv').exists()
        out = pandas.read_csv(out_path, sep='\t')
        assert out.loc[0, 'run'] == 'SRR001'
        assert out.loc[0, 'fastp_duplication_rate'] == pytest.approx(22.9565)
        assert out.loc[0, 'fastp_insert_size_peak'] == pytest.approx(138.0)
        assert out.loc[0, 'num_fastp_in'] == 1000
        assert out.loc[0, 'num_fastp_out'] == 900
        assert out.loc[0, 'bp_fastp_in'] == 100000
        assert out.loc[0, 'bp_fastp_out'] == 90000

    def test_write_getfastq_stats_writes_tsv(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_dumped': [10],
            'num_rejected': [3],
            'num_written': [7],
            'num_fastp_in': [7],
            'num_fastp_out': [6],
            'num_rrna_in': [6],
            'num_rrna_out': [5],
            'num_contam_in': [5],
            'num_contam_out': [4],
            'bp_dumped': [1000],
            'bp_rejected': [300],
            'bp_written': [700],
            'bp_fastp_in': [700],
            'bp_fastp_out': [600],
            'bp_rrna_in': [600],
            'bp_rrna_out': [500],
            'bp_contam_in': [500],
            'bp_contam_out': [450],
            'bp_discarded': [500],
            'sec_sra_download': [0.75],
            'sec_fasterq_dump': [1.25],
            'sec_fastp': [2.5],
            'sec_rrna_filter': [3.75],
            'sec_rrna_search': [3.0],
            'sec_rrna_rewrite': [0.75],
            'sec_contam_filter': [4.5],
            'sec_ete_taxonomy': [1.5],
            'fastp_duplication_rate': [12.0],
            'fastp_insert_size_peak': [250.0],
        }))
        sra_stat = {'sra_id': 'SRR001'}
        write_getfastq_stats(sra_stat=sra_stat, metadata=metadata, output_dir=str(tmp_path))
        out_path = tmp_path / 'getfastq_stats.tsv'
        assert out_path.exists()
        out = pandas.read_csv(out_path, sep='\t')
        assert out.loc[0, 'run'] == 'SRR001'
        assert out.loc[0, 'num_dumped'] == 10
        assert out.loc[0, 'num_rejected'] == 3
        assert out.loc[0, 'bp_rejected'] == 300
        assert out.loc[0, 'bp_discarded'] == 500
        assert out.loc[0, 'sec_sra_download'] == pytest.approx(0.75)
        assert out.loc[0, 'sec_fasterq_dump'] == pytest.approx(1.25)
        assert out.loc[0, 'sec_fastp'] == pytest.approx(2.5)
        assert out.loc[0, 'sec_rrna_filter'] == pytest.approx(3.75)
        assert out.loc[0, 'sec_rrna_search'] == pytest.approx(3.0)
        assert out.loc[0, 'sec_rrna_rewrite'] == pytest.approx(0.75)
        assert out.loc[0, 'sec_contam_filter'] == pytest.approx(4.5)
        assert out.loc[0, 'sec_ete_taxonomy'] == pytest.approx(1.5)
        assert out.loc[0, 'fastp_duplication_rate'] == pytest.approx(12.0)
        assert out.loc[0, 'fastp_insert_size_peak'] == pytest.approx(250.0)
        assert out.loc[0, 'percent_fastp_filtered'] == pytest.approx((100.0 / 700.0) * 100.0)
        assert out.loc[0, 'percent_rrna_filtered'] == pytest.approx((100.0 / 600.0) * 100.0)
        assert out.loc[0, 'percent_contam_filtered'] == pytest.approx(10.0)

    def test_write_getfastq_stats_keeps_existing_file_when_atomic_write_fails(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_dumped': [10],
            'num_rejected': [3],
            'num_written': [7],
            'num_fastp_in': [7],
            'num_fastp_out': [6],
            'num_rrna_in': [6],
            'num_rrna_out': [5],
            'bp_dumped': [1000],
            'bp_rejected': [300],
            'bp_written': [700],
            'bp_fastp_in': [700],
            'bp_fastp_out': [600],
            'bp_rrna_in': [600],
            'bp_rrna_out': [500],
            'bp_discarded': [500],
            'sec_sra_download': [0.75],
            'sec_fasterq_dump': [1.25],
            'sec_fastp': [2.5],
            'sec_rrna_filter': [3.75],
            'sec_rrna_search': [3.0],
            'sec_rrna_rewrite': [0.75],
            'sec_contam_filter': [4.5],
            'sec_ete_taxonomy': [1.5],
            'fastp_duplication_rate': [12.0],
            'fastp_insert_size_peak': [250.0],
        }))
        sra_stat = {'sra_id': 'SRR001'}
        out_path = tmp_path / 'getfastq_stats.tsv'
        out_path.write_text('old\n')
        original_to_csv = pandas.DataFrame.to_csv

        def fake_to_csv(self, path_or_buf=None, *args, **kwargs):
            with open(path_or_buf, 'w') as handle:
                handle.write('partial\n')
            raise RuntimeError('boom')

        monkeypatch.setattr(pandas.DataFrame, 'to_csv', fake_to_csv)

        with pytest.raises(RuntimeError, match='boom'):
            write_getfastq_stats(sra_stat=sra_stat, metadata=metadata, output_dir=str(tmp_path))

        assert out_path.read_text() == 'old\n'
        assert list(tmp_path.glob('amalgkit_atomic_*')) == []
        monkeypatch.setattr(pandas.DataFrame, 'to_csv', original_to_csv)


class TestPrintReadStats:
    def test_includes_stage_duration_lines(self, capsys):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'bp_dumped': [1000, 500],
            'bp_rejected': [100, 50],
            'bp_written': [900, 450],
            'bp_fastp_in': [900, 450],
            'bp_fastp_out': [800, 400],
            'bp_rrna_in': [800, 400],
            'bp_rrna_out': [700, 350],
            'bp_contam_in': [700, 350],
            'bp_contam_out': [650, 300],
            'sec_sra_download': [0.5, 1.0],
            'sec_fasterq_dump': [1.5, 2.0],
            'sec_fastp': [3.0, 4.0],
            'sec_rrna_filter': [5.0, 6.0],
            'sec_rrna_search': [4.0, 5.0],
            'sec_rrna_rewrite': [1.0, 1.0],
            'sec_contam_filter': [7.0, 8.0],
            'sec_ete_taxonomy': [0.75, 1.25],
        }))
        args = SimpleNamespace(fastp=True, rrna_filter=True, contam_filter=True)

        print_read_stats(args=args, metadata=metadata, g={'max_bp': 10000}, sra_stat=None, individual=True)

        out = capsys.readouterr().out
        assert 'Sum of SRA download wall time: 1.5 sec' in out
        assert 'Sum of fasterq-dump wall time: 3.5 sec' in out
        assert 'Sum of fastp wall time: 7.0 sec' in out
        assert 'Sum of MMseqs2 rRNA filter wall time: 11.0 sec' in out
        assert 'Sum of MMseqs rRNA search wall time: 9.0 sec' in out
        assert 'Sum of rRNA FASTQ rewrite wall time: 2.0 sec' in out
        assert 'Sum of contaminant-filter wall time: 15.0 sec' in out
        assert 'Sum of NCBI taxonomy wall time: 2.0 sec' in out
        assert 'Individual SRA download wall time (sec): 0.5 1.0' in out
        assert 'Individual fasterq-dump wall time (sec): 1.5 2.0' in out
        assert 'Individual fastp wall time (sec): 3.0 4.0' in out
        assert 'Individual MMseqs rRNA search wall time (sec): 4.0 5.0' in out
        assert 'Individual rRNA FASTQ rewrite wall time (sec): 1.0 1.0' in out


class TestProcessGetfastqRun:
    def test_rechecks_resume_state_only_after_acquiring_run_lock(self, tmp_path, monkeypatch):
        args = SimpleNamespace(
            redo=False,
            out_dir=str(tmp_path),
            download_lock_dir=str(tmp_path / 'locks'),
        )
        run_row_df = pandas.DataFrame({
            'run': ['SRR001'],
            'total_bases': [1000],
        })
        run_dir = tmp_path / 'getfastq' / 'SRR001'
        run_dir.mkdir(parents=True)
        lock_state = {'held': False, 'path': None}

        class DummyLock:
            def __init__(
                self,
                lock_path,
                lock_label='Lock',
                poll_seconds=5,
                timeout_seconds=3600,
            ):
                _ = (lock_label, poll_seconds, timeout_seconds)
                lock_state['path'] = lock_path

            def __enter__(self):
                lock_state['held'] = True
                return None

            def __exit__(self, exc_type, exc, tb):
                lock_state['held'] = False
                return False

        monkeypatch.setattr('amalgkit.getfastq.acquire_exclusive_lock', DummyLock)
        monkeypatch.setattr(
            'amalgkit.getfastq.get_sra_stat',
            lambda sra_id, _metadata, _num_bp_per_sra: {
                'sra_id': sra_id,
                'layout': 'single',
                'total_spot': 10,
                'spot_length': 100,
                'metadata_idx': 0,
            },
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.get_getfastq_run_dir',
            lambda _args, _sra_id: str(run_dir),
        )

        def fake_inspect(_args, _sra_stat, _g, run_metadata):
            assert lock_state['held']
            return {
                'phase': GETFASTQ_PHASE_COMPLETE,
                'stats_row': run_metadata.df.iloc[0],
            }

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('completed resume state must skip fresh extraction')

        monkeypatch.setattr('amalgkit.getfastq.inspect_getfastq_resume_output', fake_inspect)
        monkeypatch.setattr('amalgkit.getfastq.download_sra', fail_if_called)
        monkeypatch.setattr('amalgkit.getfastq.sequence_extraction_1st_round', fail_if_called)

        result = process_getfastq_run(
            args=args,
            row_index=0,
            sra_id='SRR001',
            run_row_df=run_row_df,
            g={'num_bp_per_sra': 500},
        )

        assert result['completion_phase'] == GETFASTQ_PHASE_COMPLETE
        assert not lock_state['held']
        assert lock_state['path'] == str(tmp_path / 'locks' / 'getfastq_runs' / 'SRR001.lock')

    def test_tracks_sra_download_wall_time(self, tmp_path, monkeypatch, capsys):
        args = SimpleNamespace(redo=False, out_dir=str(tmp_path))
        run_row_df = pandas.DataFrame({
            'run': ['SRR001'],
            'total_bases': [1000],
        })
        g = {'num_bp_per_sra': 500}
        run_dir = tmp_path / 'SRR001'
        run_dir.mkdir()

        monkeypatch.setattr(
            'amalgkit.getfastq.get_sra_stat',
            lambda sra_id, _metadata, _num_bp_per_sra: {
                'sra_id': sra_id,
                'layout': 'single',
                'total_spot': 10,
                'spot_length': 100,
            },
        )
        monkeypatch.setattr('amalgkit.getfastq.get_getfastq_run_dir', lambda _args, _sra_id: str(run_dir))
        monkeypatch.setattr('amalgkit.getfastq.remove_old_intermediate_files', lambda **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.download_sra', lambda *_args, **_kwargs: None)
        monkeypatch.setattr(
            'amalgkit.getfastq.sequence_extraction_1st_round',
            lambda _args, _sra_stat, metadata, _g, runtime_context=None: metadata,
        )
        monkeypatch.setattr('amalgkit.getfastq.write_getfastq_run_state', lambda *_args, **_kwargs: None)
        perf_counter_values = iter([10.0, 12.5])
        monkeypatch.setattr('amalgkit.getfastq.time.perf_counter', lambda: next(perf_counter_values))

        out = process_getfastq_run(
            args=args,
            row_index=0,
            sra_id='SRR001',
            run_row_df=run_row_df,
            g=g,
        )

        assert out['flag_private_file'] is False
        assert out['flag_any_output_file_present'] is False
        assert out['row']['sec_sra_download'] == pytest.approx(2.5)
        stdout = capsys.readouterr().out
        assert 'Time elapsed for SRA download (SRR001):' in stdout


class TestRunFastp:
    @staticmethod
    def _build_metadata():
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'num_fastp_in': [0],
            'num_fastp_out': [0],
            'bp_fastp_in': [0],
            'bp_fastp_out': [0],
            'fastp_duplication_rate': [numpy.nan],
            'fastp_insert_size_peak': [numpy.nan],
        }))

    @staticmethod
    def _materialize_fastp_outputs(cmd):
        for flag in ['--out1', '--out2']:
            if flag not in cmd:
                continue
            out_path = cmd[cmd.index(flag) + 1]
            with open(out_path, 'wb') as fout:
                fout.write(b'@r0\nAAAA\n+\nIIII\n')

    def test_uses_configured_fastp_exe_and_shlex_parsing(self, tmp_path, monkeypatch, capsys):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        captured = {}

        class Args:
            threads = 4
            min_read_length = 25
            fastp_option = '--adapter_sequence "A B"'
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            captured['cmd'] = cmd
            self._materialize_fastp_outputs(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        perf_counter_values = iter([30.0, 31.75])
        monkeypatch.setattr('amalgkit.getfastq.time.perf_counter', lambda: next(perf_counter_values))

        run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

        assert captured['cmd'][0] == 'fastp-custom'
        assert '--adapter_sequence' in captured['cmd']
        assert 'A B' in captured['cmd']
        assert metadata.df.loc[0, 'num_fastp_in'] == 10
        assert metadata.df.loc[0, 'num_fastp_out'] == 8
        assert metadata.df.loc[0, 'bp_fastp_in'] == 100
        assert metadata.df.loc[0, 'bp_fastp_out'] == 80
        assert metadata.df.loc[0, 'sec_fastp'] == pytest.approx(1.75)
        assert not (tmp_path / 'fastp_stats.tsv').exists()
        assert not (tmp_path / 'getfastq_stats.tsv').exists()
        out = capsys.readouterr().out
        assert 'Time elapsed for fastp (SRR001):' in out

    def test_raises_when_fastp_exits_nonzero(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 1, stdout=b'', stderr=b'fastp failed')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(RuntimeError, match='fastp did not finish safely'):
            run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

    def test_raises_when_fastp_stderr_is_truncated(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        def fake_run(cmd, stdout=None, stderr=None):
            self._materialize_fastp_outputs(cmd)
            bad_stderr = '\n'.join([
                ' before filtering:',
                'total reads: 10',
            ]).encode('utf8')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=bad_stderr)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(RuntimeError, match='Unexpected fastp stderr format'):
            run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

    def test_removes_empty_outputs_when_fastp_filters_everything(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 0',
            'total bases: 0',
            'Duplication rate: 0.0%',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            self._materialize_fastp_outputs(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

        assert metadata.df.loc[0, 'num_fastp_out'] == 0
        assert metadata.df.loc[0, 'bp_fastp_out'] == 0
        assert not (tmp_path / 'SRR001.fastp.fastq.gz').exists()

    def test_uses_cached_extension_without_redetection(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single', 'current_ext': '.fastq.gz'}
        captured = {}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('get_newest_intermediate_file_extension should not be called with cached current_ext.')

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            fail_if_called,
        )

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            captured['cmd'] = cmd
            self._materialize_fastp_outputs(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

        assert '--in1' in captured['cmd']
        in1_arg = captured['cmd'][captured['cmd'].index('--in1') + 1]
        assert in1_arg.endswith('SRR001.fastq.gz')
        assert sra_stat['current_ext'] == '.fastp.fastq.gz'

    def test_uses_prefetched_files_for_extension_detection_without_rescan(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        captured = {'seen_files': None}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        def fake_get_newest_intermediate_file_extension(sra_stat, work_dir, files=None):
            captured['seen_files'] = set(files)
            return '.fastq.gz'

        def fail_if_listed(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when files are provided.')

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            fake_get_newest_intermediate_file_extension,
        )
        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_if_listed)

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            self._materialize_fastp_outputs(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        metadata, files_out = run_fastp(
            sra_stat=sra_stat,
            args=Args(),
            output_dir=str(tmp_path),
            metadata=metadata,
            files={'SRR001.fastq.gz'},
            return_files=True,
        )

        assert captured['seen_files'] == {'SRR001.fastq.gz'}
        assert 'SRR001.fastp.fastq.gz' in files_out
        assert metadata.df.loc[0, 'num_fastp_out'] == 8

    def test_tolerates_non_utf8_fastp_stderr(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = True
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        stderr_bytes = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8') + b'\xff'

        def fake_run(cmd, stdout=None, stderr=None):
            self._materialize_fastp_outputs(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'\xff', stderr=stderr_bytes)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

        assert metadata.df.loc[0, 'num_fastp_in'] == 10
        assert metadata.df.loc[0, 'num_fastp_out'] == 8

    def test_raises_when_fastp_output_file_is_missing(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            # Intentionally do not materialize --out1.
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(FileNotFoundError, match='fastp output file was not generated'):
            run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

    def test_raises_when_fastp_output_path_is_directory(self, tmp_path, monkeypatch):
        metadata = self._build_metadata()
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        class Args:
            threads = 1
            min_read_length = 25
            fastp_option = ''
            fastp_print = False
            remove_tmp = False
            fastp_exe = 'fastp-custom'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            lambda sra_stat, work_dir, files=None: '.fastq.gz'
        )

        stderr_txt = '\n'.join([
            ' before filtering:',
            'total reads: 10',
            'total bases: 100',
            ' after filtering:',
            'total reads: 8',
            'total bases: 80',
            'Duplication rate: 20.0%',
            'Insert size peak (evaluated by paired-end reads): 150',
        ]).encode('utf8')

        def fake_run(cmd, stdout=None, stderr=None):
            out1 = cmd[cmd.index('--out1') + 1]
            os.makedirs(out1, exist_ok=True)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(IsADirectoryError, match='fastp output path exists but is not a file'):
            run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)


class TestIdenticalPairedReads:
    @staticmethod
    def _write_fastq_gz(path, seqs):
        with gzip.open(path, 'wt') as out:
            for i, seq in enumerate(seqs, start=1):
                out.write('@read{}\n{}\n+\n{}\n'.format(i, seq, 'I' * len(seq)))

    @staticmethod
    def _write_fastq_plain(path, seqs):
        with open(path, 'wt') as out:
            for i, seq in enumerate(seqs, start=1):
                out.write('@read{}\n{}\n+\n{}\n'.format(i, seq, 'I' * len(seq)))

    def test_get_identical_paired_ratio(self, tmp_path):
        read1 = tmp_path / 'read1.fastq.gz'
        read2 = tmp_path / 'read2.fastq.gz'
        self._write_fastq_gz(read1, ['AAAA', 'CCCC', 'GGGG'])
        self._write_fastq_gz(read2, ['AAAA', 'TTTT', 'GGGG'])
        ratio, num_checked, read_length = get_identical_paired_ratio(str(read1), str(read2), num_checked_reads=3)
        assert ratio == pytest.approx(2 / 3)
        assert num_checked == 3
        assert read_length == 4

    def test_get_identical_paired_ratio_plain_fastq(self, tmp_path):
        read1 = tmp_path / 'read1.fastq'
        read2 = tmp_path / 'read2.fastq'
        self._write_fastq_plain(read1, ['AAAA', 'CCCC', 'GGGG'])
        self._write_fastq_plain(read2, ['AAAA', 'TTTT', 'GGGG'])
        ratio, num_checked, read_length = get_identical_paired_ratio(str(read1), str(read2), num_checked_reads=3)
        assert ratio == pytest.approx(2 / 3)
        assert num_checked == 3
        assert read_length == 4

    def test_maybe_treat_paired_as_single_converts_files(self, tmp_path):
        sra_id = 'SRR001'
        read1 = tmp_path / '{}_1.fastq.gz'.format(sra_id)
        read2 = tmp_path / '{}_2.fastq.gz'.format(sra_id)
        self._write_fastq_gz(read1, ['AAAA', 'CCCC', 'GGGG'])
        self._write_fastq_gz(read2, ['AAAA', 'CCCC', 'GGGG'])
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['paired'],
            'layout_amalgkit': ['paired'],
            'spot_length': [8],
            'total_spots': [3],
            'total_bases': [24],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 8,
        }
        metadata, sra_stat = maybe_treat_paired_as_single(
            sra_stat=sra_stat,
            metadata=metadata,
            work_dir=str(tmp_path),
            threshold=0.99,
            num_checked_reads=3,
        )
        assert sra_stat['layout'] == 'single'
        assert os.path.exists(str(tmp_path / '{}.fastq.gz'.format(sra_id)))
        assert not os.path.exists(str(read1))
        assert not os.path.exists(str(read2))
        assert metadata.df.loc[0, 'layout_amalgkit'] == 'single'
        assert metadata.df.loc[0, 'spot_length'] == 4

    def test_maybe_treat_paired_as_single_keeps_paired_when_ratio_low(self, tmp_path):
        sra_id = 'SRR001'
        read1 = tmp_path / '{}_1.fastq.gz'.format(sra_id)
        read2 = tmp_path / '{}_2.fastq.gz'.format(sra_id)
        self._write_fastq_gz(read1, ['AAAA', 'CCCC', 'GGGG'])
        self._write_fastq_gz(read2, ['TTTT', 'CCCC', 'AAAA'])
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['paired'],
            'layout_amalgkit': ['paired'],
            'spot_length': [8],
            'total_spots': [3],
            'total_bases': [24],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 8,
        }
        metadata, sra_stat = maybe_treat_paired_as_single(
            sra_stat=sra_stat,
            metadata=metadata,
            work_dir=str(tmp_path),
            threshold=0.99,
            num_checked_reads=3,
        )
        assert sra_stat['layout'] == 'paired'
        assert os.path.exists(str(read1))
        assert os.path.exists(str(read2))
        assert metadata.df.loc[0, 'layout_amalgkit'] == 'paired'

    def test_maybe_treat_paired_as_single_uses_prefetched_files_without_rescan(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        read1 = tmp_path / '{}_1.fastq.gz'.format(sra_id)
        read2 = tmp_path / '{}_2.fastq.gz'.format(sra_id)
        self._write_fastq_gz(read1, ['AAAA', 'CCCC'])
        self._write_fastq_gz(read2, ['AAAA', 'CCCC'])
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['paired'],
            'layout_amalgkit': ['paired'],
            'spot_length': [8],
            'total_spots': [2],
            'total_bases': [16],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 8,
        }

        def fail_if_called(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when files are provided.')

        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_if_called)
        metadata, sra_stat = maybe_treat_paired_as_single(
            sra_stat=sra_stat,
            metadata=metadata,
            work_dir=str(tmp_path),
            threshold=0.99,
            num_checked_reads=2,
            files={'{}_1.fastq.gz'.format(sra_id), '{}_2.fastq.gz'.format(sra_id)},
        )
        assert sra_stat['layout'] == 'single'


class TestWrittenSpotEstimation:
    @staticmethod
    def _write_fastq(path, seqs):
        with open(path, 'wt') as out:
            for i, seq in enumerate(seqs):
                out.write('@r{}\n{}\n+\n{}\n'.format(i, seq, 'I' * len(seq)))

    def test_uses_prefetched_files_without_rescan(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        fastq_path = tmp_path / '{}.fastq'.format(sra_id)
        self._write_fastq(fastq_path, ['AAAA', 'CCCC', 'GGGG'])
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
        }

        def fail_if_called(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when files are provided.')

        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_if_called)
        observed = estimate_num_written_spots_from_fastq(sra_stat, files={'{}.fastq'.format(sra_id)})
        assert observed == 3


class TestSraRecovery:
    def test_full_dump_fallback_is_guarded_by_disk_estimate(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            'amalgkit.getfastq.shutil.disk_usage',
            lambda _path: SimpleNamespace(total=1000, used=900, free=100),
        )

        with pytest.raises(OSError, match='Insufficient disk space'):
            guard_fasterq_full_dump_disk_space({
                'sra_id': 'SRR001',
                'getfastq_sra_dir': str(tmp_path),
                'total_spot': 100,
                'spot_length': 100,
            })

    def test_full_dump_fails_closed_without_estimate_when_size_check_is_off(self, tmp_path):
        with pytest.raises(RuntimeError, match='could not be estimated'):
            guard_fasterq_full_dump_disk_space(
                {
                    'sra_id': 'SRR001',
                    'getfastq_sra_dir': str(tmp_path),
                    'total_spot': 0,
                    'spot_length': 0,
                },
                size_check='off',
            )

    def test_non_range_fasterq_stops_before_full_dump_when_disk_is_insufficient(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        (tmp_path / '{}.sra'.format(sra_id)).write_bytes(b'sra')
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'total_spot': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()
        args._fasterq_supports_spot_range = False
        monkeypatch.setattr(
            'amalgkit.getfastq.shutil.disk_usage',
            lambda _path: SimpleNamespace(total=1000, used=900, free=100),
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.subprocess.run',
            lambda *_args, **_kwargs: pytest.fail('fasterq-dump must not start after disk guard failure'),
        )

        with pytest.raises(OSError, match='Insufficient disk space'):
            run_fasterq_dump(sra_stat, args, metadata, start=1, end=2)

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
            'layout_amalgkit': ['paired'],
        }))

    @staticmethod
    def _args_for_fasterq_dump():
        class Args:
            threads = 2
            min_read_length = 25
            dump_print = False
            fasterq_dump_exe = 'fasterq-dump'
            fasterq_size_check = True
            fasterq_disk_limit = None
            fasterq_disk_limit_tmp = None
        return Args()

    @pytest.fixture(autouse=True)
    def _disable_fasterq_output_validation(self, monkeypatch):
        monkeypatch.setattr('amalgkit.getfastq.ensure_fasterq_output_files_exist', lambda **_kwargs: None)

    def test_resolve_public_original_fastq_sources_from_xml_root_extracts_original_fastqs(self):
        xml_root = ET.fromstring(
            """
            <RunBundle>
              <RUN accession="SRR001">
                <SRAFiles>
                  <SRAFile filename="forward_reads.fastq.gz" semantic_name="fastq" supertype="Original"
                           url="https://trace.example/forward.fastq.gz">
                    <Alternatives url="https://aws.example/forward.fastq.gz" org="AWS"/>
                  </SRAFile>
                  <SRAFile filename="reverse_reads.fastq.gz" semantic_name="fastq" supertype="Original"
                           url="https://trace.example/reverse.fastq.gz">
                    <Alternatives url="https://aws.example/reverse.fastq.gz" org="AWS"/>
                  </SRAFile>
                  <SRAFile filename="SRR001" semantic_name="SRA Normalized" supertype="Primary ETL"
                           url="https://aws.example/SRR001"/>
                </SRAFiles>
              </RUN>
            </RunBundle>
            """
        )

        observed = resolve_public_original_fastq_sources_from_xml_root(xml_root)
        assigned = assign_public_original_fastq_suffixes(observed, 'SRR001')

        assert [entry['filename'] for entry in observed] == [
            'forward_reads.fastq.gz',
            'reverse_reads.fastq.gz',
        ]
        assert observed[0]['sources'][0]['source_name'] == 'AWS'
        assert observed[0]['sources'][0]['url'] == 'https://aws.example/forward.fastq.gz'
        assert [(entry['filename'], suffix) for entry, suffix in assigned] == [
            ('forward_reads.fastq.gz', '_1'),
            ('reverse_reads.fastq.gz', '_2'),
        ]

    def test_assigns_split3_original_fastq_files(self):
        entries = [
            {'filename': 'SRR001.fastq.gz', 'sources': []},
            {'filename': 'SRR001_1.fastq.gz', 'sources': []},
            {'filename': 'SRR001_2.fastq.gz', 'sources': []},
        ]

        assigned = assign_public_original_fastq_suffixes(entries, 'SRR001')

        assert [(entry['filename'], suffix) for entry, suffix in assigned] == [
            ('SRR001.fastq.gz', ''),
            ('SRR001_1.fastq.gz', '_1'),
            ('SRR001_2.fastq.gz', '_2'),
        ]

    def test_remove_sra_files_deletes_matching_sra_files(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001.sra').write_text('a')
        (sra_dir / 'SRR001.sra.vdbcache').write_text('b')
        (sra_dir / 'other.txt').write_text('keep')

        remove_sra_files(metadata, str(tmp_path))

        assert not (sra_dir / 'SRR001.sra').exists()
        assert not (sra_dir / 'SRR001.sra.vdbcache').exists()
        assert (sra_dir / 'other.txt').exists()

    def test_remove_sra_files_avoids_root_listdir_scan(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001.sra').write_text('a')

        def fail_if_listdir_called(_path):
            raise AssertionError('remove_sra_files should not call os.listdir on getfastq root.')

        monkeypatch.setattr('amalgkit.getfastq.os.listdir', fail_if_listdir_called)
        remove_sra_files(metadata, str(tmp_path))

        assert not (sra_dir / 'SRR001.sra').exists()

    def test_remove_sra_files_ignores_non_directory_entries(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        getfastq_root = tmp_path / 'getfastq'
        getfastq_root.mkdir(parents=True)
        (getfastq_root / 'SRR001').write_text('not a directory')

        remove_sra_files(metadata, str(tmp_path))

        assert (getfastq_root / 'SRR001').exists()

    def test_remove_sra_files_handles_missing_getfastq_root(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))

        remove_sra_files(metadata, str(tmp_path))

    def test_remove_sra_files_skips_missing_run_ids(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', numpy.nan, ''],
            'scientific_name': ['Sp1', 'Sp1', 'Sp1'],
            'exclusion': ['no', 'no', 'no'],
        }))
        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001.sra').write_text('a')

        remove_sra_files(metadata, str(tmp_path))

        assert not (sra_dir / 'SRR001.sra').exists()

    def test_remove_sra_files_keeps_non_artifact_suffix(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001.sra').write_text('a')
        (sra_dir / 'SRR001.sra.vdbcache').write_text('b')
        (sra_dir / 'SRR001.sra2').write_text('keep')

        remove_sra_files(metadata, str(tmp_path))

        assert not (sra_dir / 'SRR001.sra').exists()
        assert not (sra_dir / 'SRR001.sra.vdbcache').exists()
        assert (sra_dir / 'SRR001.sra2').exists()

    def test_remove_sra_files_rejects_unsafe_run_before_deleting(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', '../escape'],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
        }))
        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        sra_path = sra_dir / 'SRR001.sra'
        sra_path.write_text('keep')

        with pytest.raises(ValueError, match='Run ID'):
            remove_sra_files(metadata, str(tmp_path))

        assert sra_path.exists()

    def test_remove_sra_files_rejects_symlink_run_directory(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        external_dir = tmp_path / 'external'
        external_dir.mkdir()
        external_sra = external_dir / 'SRR001.sra'
        external_sra.write_text('keep')
        getfastq_root = tmp_path / 'getfastq'
        getfastq_root.mkdir()
        os.symlink(external_dir, getfastq_root / 'SRR001')

        with pytest.raises(ValueError, match='symbolic-link'):
            remove_sra_files(metadata, str(tmp_path))

        assert external_sra.exists()

    def test_remove_sra_files_rejects_symlink_artifact(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        external_sra = tmp_path / 'external.sra'
        external_sra.write_text('keep')
        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        os.symlink(external_sra, sra_dir / 'SRR001.sra')

        with pytest.raises(ValueError, match='symbolic-link SRA artifact'):
            remove_sra_files(metadata, str(tmp_path))

        assert external_sra.exists()

    def test_remove_sra_path_file_and_directory(self, tmp_path):
        file_path = tmp_path / 'SRR001.sra'
        file_path.write_text('dummy')
        remove_sra_path(str(file_path))
        assert not file_path.exists()

        dir_path = tmp_path / 'SRR002.sra'
        (dir_path / 'tbl').mkdir(parents=True)
        (dir_path / 'tbl' / 'x').write_text('dummy')
        remove_sra_path(str(dir_path))
        assert not dir_path.exists()

    def test_run_fasterq_dump_retries_once_after_redownload(self, tmp_path, monkeypatch, capsys):
        sra_id = 'SRR001'
        sra_path = tmp_path / '{}.sra'.format(sra_id)
        sra_path.mkdir()
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        run_calls = {'count': 0}

        def fake_run(cmd, stdout=None, stderr=None):
            run_calls['count'] += 1
            if run_calls['count'] == 1:
                (tmp_path / '{}_1.fastq'.format(sra_id)).write_text('partial')
                return subprocess.CompletedProcess(cmd, 1, stdout=b'', stderr=b'fasterq-dump failed')
            assert not (tmp_path / '{}_1.fastq'.format(sra_id)).exists()
            (tmp_path / '{}_1.fastq'.format(sra_id)).write_text('retry-r1')
            (tmp_path / '{}_2.fastq'.format(sra_id)).write_text('retry-r2')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        redownload_calls = []

        def fake_download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
            redownload_calls.append(overwrite)
            with open(os.path.join(work_dir, sra_stat['sra_id'] + '.sra'), 'w') as fh:
                fh.write('fresh')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.download_sra', fake_download_sra)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 4)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, sra_stat_out = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 2
        assert redownload_calls == [True]
        assert sra_path.is_file()
        assert metadata.df.loc[0, 'num_written'] == 4
        assert metadata.df.loc[0, 'num_dumped'] == 10
        assert metadata.df.loc[0, 'num_rejected'] == 6
        assert metadata.df.loc[0, 'bp_written'] == 400
        assert metadata.df.loc[0, 'bp_dumped'] == 1000
        assert metadata.df.loc[0, 'bp_rejected'] == 600
        assert sra_stat_out['layout'] == 'paired'
        out = capsys.readouterr().out
        assert 'Time elapsed for SRA re-download (SRR001):' in out

    @pytest.mark.slow
    def test_run_fasterq_dump_exits_when_retry_fails(self, tmp_path, monkeypatch, capsys):
        sra_id = 'SRR001'
        sra_path = tmp_path / '{}.sra'.format(sra_id)
        sra_path.write_text('broken')
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        run_calls = {'count': 0}

        def fake_run(cmd, stdout=None, stderr=None):
            run_calls['count'] += 1
            return subprocess.CompletedProcess(cmd, 1, stdout=b'', stderr=b'fasterq-dump failed')

        redownload_calls = []

        def fake_download_sra(metadata, sra_stat, args, work_dir, overwrite=False):
            redownload_calls.append(overwrite)
            with open(os.path.join(work_dir, sra_stat['sra_id'] + '.sra'), 'w') as fh:
                fh.write('fresh')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.download_sra', fake_download_sra)

        with pytest.raises(RuntimeError, match='fasterq-dump did not finish safely after re-download'):
            run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 2
        assert redownload_calls == [True]
        captured = capsys.readouterr()
        assert 'Command failed with exit code 1:' in captured.err
        assert 'Retry command failed with exit code 1:' in captured.err
        assert 'fasterq-dump stderr:' in captured.err
        assert 'fasterq-dump failed' in captured.err

    def test_run_fasterq_dump_falls_back_to_public_original_fastq_when_retry_fails(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
            'total_spot': 10,
        }
        args = self._args_for_fasterq_dump()
        args.sra_download_method = 'urllib'

        def fail_retry(**_kwargs):
            raise RuntimeError('fasterq-dump did not finish safely after re-download.')

        def fake_fetch_public_original_fastq_sources(sra_id):
            assert sra_id == 'SRR001'
            return [
                {
                    'filename': 'forward_reads.fastq.gz',
                    'sources': [{'source_name': 'AWS', 'url': 'https://example.org/forward.fastq.gz'}],
                },
                {
                    'filename': 'reverse_reads.fastq.gz',
                    'sources': [{'source_name': 'AWS', 'url': 'https://example.org/reverse.fastq.gz'}],
                },
            ]

        def _write_fastq_gz(path, suffix):
            with gzip.open(path, 'wt') as handle:
                for i in range(10):
                    handle.write('@read{}_{}\n'.format(i + 1, suffix))
                    handle.write('ACGT\n')
                    handle.write('+\n')
                    handle.write('!!!!\n')

        def fake_urlretrieve(url, out_path):
            if 'forward' in url:
                _write_fastq_gz(out_path, '1')
            elif 'reverse' in url:
                _write_fastq_gz(out_path, '2')
            else:
                raise AssertionError('Unexpected URL {}'.format(url))

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump_with_retry', fail_retry)
        monkeypatch.setattr('amalgkit.getfastq.fetch_public_original_fastq_sources', fake_fetch_public_original_fastq_sources)
        monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda _name: None)
        monkeypatch.setattr(
            'amalgkit.getfastq.download_with_urllib',
            lambda source_url, output_path, timeout_seconds: fake_urlretrieve(source_url, output_path),
        )

        metadata, sra_stat_out = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert count_fastq_records(str(tmp_path / 'SRR001_1.fastq.gz')) == 10
        assert count_fastq_records(str(tmp_path / 'SRR001_2.fastq.gz')) == 10
        assert metadata.df.loc[0, 'num_written'] == 10
        assert metadata.df.loc[0, 'num_dumped'] == 10
        assert metadata.df.loc[0, 'num_rejected'] == 0
        assert metadata.df.loc[0, 'bp_written'] == 1000
        assert metadata.df.loc[0, 'layout_amalgkit'] == 'paired'
        assert sra_stat_out['layout'] == 'paired'

    def test_run_fasterq_dump_raises_when_public_original_fastq_fallback_is_unavailable(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
            'total_spot': 10,
        }
        args = self._args_for_fasterq_dump()

        def fail_retry(**_kwargs):
            raise RuntimeError('fasterq-dump did not finish safely after re-download.')

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump_with_retry', fail_retry)
        monkeypatch.setattr('amalgkit.getfastq.fetch_public_original_fastq_sources', lambda _sra_id: [])

        with pytest.raises(RuntimeError, match='fasterq-dump did not finish safely after re-download'):
            run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

    def test_run_fasterq_dump_no_redownload_when_first_attempt_succeeds(self, tmp_path, monkeypatch, capsys):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        run_calls = {'count': 0}

        def fake_run(cmd, stdout=None, stderr=None):
            run_calls['count'] += 1
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        def fail_download(*_args, **_kwargs):
            raise AssertionError('download_sra should not be called when first extraction succeeds.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.download_sra', fail_download)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 5)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)
        perf_counter_values = iter([100.0, 103.25])
        monkeypatch.setattr('amalgkit.getfastq.time.perf_counter', lambda: next(perf_counter_values))

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 1
        assert metadata.df.loc[0, 'num_written'] == 5
        assert metadata.df.loc[0, 'num_dumped'] == 10
        assert metadata.df.loc[0, 'num_rejected'] == 5
        assert metadata.df.loc[0, 'bp_written'] == 500
        assert metadata.df.loc[0, 'bp_dumped'] == 1000
        assert metadata.df.loc[0, 'sec_fasterq_dump'] == pytest.approx(3.25)
        out = capsys.readouterr().out
        assert 'Time elapsed for fasterq-dump (SRR001):' in out

    def test_run_fasterq_dump_uses_reported_spots_without_recount_for_partial_range(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'spots written      : 7\n',
                stderr=b'',
            )

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when spots written is reported.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=2, end=8)

        assert metadata.df.loc[0, 'num_written'] == 7
        assert metadata.df.loc[0, 'num_dumped'] == 7
        assert metadata.df.loc[0, 'num_rejected'] == 0
        assert metadata.df.loc[0, 'bp_written'] == 700
        assert metadata.df.loc[0, 'bp_dumped'] == 700
        assert metadata.df.loc[0, 'bp_rejected'] == 0

    def test_run_fasterq_dump_ignores_written_total_line(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 0, stdout=b'written markers should be ignored', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 6667)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert metadata.df.loc[0, 'num_written'] == 6667
        assert metadata.df.loc[0, 'num_dumped'] == 6667
        assert metadata.df.loc[0, 'num_rejected'] == 0

    def test_run_fasterq_dump_infers_spots_from_reported_written_reads_for_paired_layout(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'',
                stderr=b'reads written   : 20\n',
            )

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when spots are inferred from reads written.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert metadata.df.loc[0, 'num_written'] == 10
        assert metadata.df.loc[0, 'num_dumped'] == 10
        assert metadata.df.loc[0, 'num_rejected'] == 0

    def test_run_fasterq_dump_falls_back_to_fastq_count_when_singletons_exist(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        (tmp_path / '{}.fastq'.format(sra_id)).write_text('@r0\nA\n+\nI\n')
        observed = {'estimate_called': False}

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'',
                stderr=b'reads written   : 21\n',
            )

        def fake_estimate(*_args, **_kwargs):
            observed['estimate_called'] = True
            return 3

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fake_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert observed['estimate_called']
        assert metadata.df.loc[0, 'num_written'] == 3

    def test_run_fasterq_dump_skips_pre_fastp_compression(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        args.fastp = True
        args.rrna_filter = False

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'',
                stderr=b'reads written   : 20\n',
            )

        def fail_compress(*_args, **_kwargs):
            raise AssertionError('compress_fasterq_output_files should be skipped when fastp is the next filter.')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when reads written is usable.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', fail_compress)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, sra_stat_out = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert sra_stat_out['current_ext'] == '.fastq'
        assert metadata.df.loc[0, 'num_written'] == 10

    def test_run_fasterq_dump_skips_pre_rrna_compression(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
        }
        args = self._args_for_fasterq_dump()
        args.fastp = False
        args.rrna_filter = True

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'',
                stderr=b'reads written   : 20\n',
            )

        def fail_compress(*_args, **_kwargs):
            raise AssertionError('compress_fasterq_output_files should be skipped when a downstream filter exists.')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when reads written is usable.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', fail_compress)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, sra_stat_out = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert sra_stat_out['current_ext'] == '.fastq'
        assert metadata.df.loc[0, 'num_written'] == 10

    def test_run_fasterq_dump_skips_trim_for_full_range(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
            'total_spot': 10,
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        def fail_trim(*_args, **_kwargs):
            raise AssertionError('trim_fasterq_output_files should be skipped for full-range extraction.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', fail_trim)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 10)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert metadata.df.loc[0, 'num_written'] == 10

    def test_run_fasterq_dump_skips_trim_for_partial_range(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
            'total_spot': 20,
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'spots written      : 7\n',
                stderr=b'',
            )

        def fail_trim(*_args, **_kwargs):
            raise AssertionError('trim_fasterq_output_files should not be called when fasterq-dump is spot-range limited.')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when spots written is reported.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', fail_trim)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=2, end=8)

        assert metadata.df.loc[0, 'num_written'] == 7

    def test_run_fasterq_dump_trims_when_spot_range_flags_are_not_supported(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'paired',
            'total_spot': 20,
        }
        args = self._args_for_fasterq_dump()
        args.fastp = True
        args._fasterq_supports_spot_range = False
        observed = {'cmd': None, 'trim': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'spots written : 100\n', stderr=b'')

        def fake_trim(*_args, **kwargs):
            observed['trim'] = (kwargs.get('start'), kwargs.get('end'), kwargs.get('compress_to_gz'))
            return {'': 0, '_1': 3, '_2': 3}, kwargs.get('file_state')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when trim counts are available.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', fake_trim)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=2, end=8)

        cmd = observed['cmd']
        assert '-N' not in cmd
        assert '-X' not in cmd
        assert observed['trim'] == (2, 8, False)
        assert metadata.df.loc[0, 'num_written'] == 3
        assert metadata.df.loc[0, 'num_dumped'] == 7
        assert metadata.df.loc[0, 'num_rejected'] == 4

    def test_run_fasterq_dump_full_range_uses_reported_spots_without_recount(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
            'total_spot': 10,
        }
        args = self._args_for_fasterq_dump()

        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(
                cmd,
                0,
                stdout=b'spots written      : 7\n',
                stderr=b'',
            )

        def fail_trim(*_args, **_kwargs):
            raise AssertionError('trim_fasterq_output_files should be skipped for full-range extraction.')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called when spots written is reported.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', fail_trim)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', fail_estimate)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        assert metadata.df.loc[0, 'num_written'] == 7
        assert metadata.df.loc[0, 'num_dumped'] == 10
        assert metadata.df.loc[0, 'num_rejected'] == 3

    def test_run_fasterq_dump_passes_size_check_and_disk_limits(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()
        args.fasterq_size_check = False
        args.fasterq_disk_limit = '10G'
        args.fasterq_disk_limit_tmp = '20G'
        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 1)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        run_fasterq_dump(sra_stat, args, metadata, start=1, end=1)

        cmd = observed['cmd']
        assert '--size-check' in cmd
        assert cmd[cmd.index('--size-check') + 1] == 'off'
        assert '--disk-limit' in cmd
        assert cmd[cmd.index('--disk-limit') + 1] == '10G'
        assert '--disk-limit-tmp' in cmd
        assert cmd[cmd.index('--disk-limit-tmp') + 1] == '20G'
        assert '-N' in cmd
        assert cmd[cmd.index('-N') + 1] == '1'
        assert '-X' in cmd
        assert cmd[cmd.index('-X') + 1] == '1'

    def test_run_fasterq_dump_passes_requested_spot_range(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
        }
        args = self._args_for_fasterq_dump()
        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 8)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        run_fasterq_dump(sra_stat, args, metadata, start=5, end=12)

        cmd = observed['cmd']
        assert '-N' in cmd
        assert cmd[cmd.index('-N') + 1] == '5'
        assert '-X' in cmd
        assert cmd[cmd.index('-X') + 1] == '12'

    def test_run_fasterq_dump_skips_spot_range_flags_for_full_range(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = self._metadata_for_extraction(sra_id)
        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
            'spot_length': 100,
            'layout': 'single',
            'total_spot': 10,
        }
        args = self._args_for_fasterq_dump()
        observed = {'cmd': None}

        def fake_run(cmd, stdout=None, stderr=None):
            observed['cmd'] = cmd
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'spots written : 10')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.compress_fasterq_output_files', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 10)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)

        cmd = observed['cmd']
        assert '-N' not in cmd
        assert '-X' not in cmd

    def test_sequence_extraction_reuses_run_files_without_rescan_when_fastp_disabled(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = False
            read_name = 'trinity'

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            if return_files:
                return metadata, sra_stat, {'{}.fastq.gz'.format(sra_id)}
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if file_state is not None:
                files = file_state.to_set()
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=files)
            if return_files:
                return metadata, sra_stat, set(files)
            return metadata, sra_stat

        captured = {'files': None}

        def fake_rename_reads(sra_stat, args, output_dir, files=None, file_state=None, return_file_state=False):
            if file_state is not None:
                files = file_state.to_set()
            captured['files'] = set(files)
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return RunFileState(work_dir=output_dir, files=files)
            return set(files)

        def fail_list_run_dir_files(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when run file snapshot is already available.')

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.rename_reads', fake_rename_reads)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_list_run_dir_files)

        metadata = sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert captured['files'] == {'{}.fastq.gz'.format(sra_id)}
        assert metadata.df.loc[0, 'num_written'] == 3
        assert metadata.df.loc[0, 'num_dumped'] == 3
        assert metadata.df.loc[0, 'bp_specified_for_extraction'] == 400

    def test_sequence_extraction_fastp_trinity_reuses_run_files_without_rescan(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            read_name = 'trinity'

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            if return_files:
                return metadata, sra_stat, {'{}.fastq.gz'.format(sra_id)}
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if file_state is not None:
                files = file_state.to_set()
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=files)
            if return_files:
                return metadata, sra_stat, set(files)
            return metadata, sra_stat

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if file_state is not None:
                files = file_state.to_set()
            assert set(files) == {'{}.fastq.gz'.format(sra_id)}
            idx = 0
            metadata.df.at[idx, 'bp_fastp_out'] += 250
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.fastp.fastq.gz'.format(sra_id)})
            if return_files:
                return metadata, {'{}.fastp.fastq.gz'.format(sra_id)}
            return metadata

        captured = {'files': None}

        def fake_rename_reads(sra_stat, args, output_dir, files=None, file_state=None, return_file_state=False):
            if file_state is not None:
                files = file_state.to_set()
            captured['files'] = set(files)
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return RunFileState(work_dir=output_dir, files=files)
            return set(files)

        def fail_list_run_dir_files(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when run file snapshot is already available.')

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.rename_reads', fake_rename_reads)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_list_run_dir_files)

        metadata = sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert captured['files'] == {'{}.fastp.fastq.gz'.format(sra_id)}
        assert metadata.df.loc[0, 'bp_fastp_out'] == 250
        assert metadata.df.loc[0, 'bp_amalgkit'] == 250

    def test_sequence_extraction_uses_rrna_filtered_bp_for_target_size(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            rrna_filter = True
            read_name = 'default'

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=file_state.to_set())
            return metadata, sra_stat

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            idx = 0
            metadata.df.at[idx, 'bp_fastp_out'] += 250
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.fastp.fastq.gz'.format(sra_id)})
            return metadata

        def fake_run_rrna_filter(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            known_input_counts=None,
            return_files=False,
            return_file_state=False,
        ):
            idx = 0
            metadata.df.at[idx, 'bp_rrna_in'] += 250
            metadata.df.at[idx, 'bp_rrna_out'] += 180
            sra_stat['current_ext'] = '.rrna-filtered.fastq.gz'
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.rrna-filtered.fastq.gz'.format(sra_id)})
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.run_rrna_filter', fake_run_rrna_filter)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)

        metadata = sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert metadata.df.loc[0, 'bp_fastp_out'] == 250
        assert metadata.df.loc[0, 'bp_rrna_out'] == 180
        assert metadata.df.loc[0, 'bp_amalgkit'] == 180

    def test_sequence_extraction_passes_fastp_output_counts_to_rrna(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            rrna_filter = True
            read_name = 'default'

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=file_state.to_set())
            return metadata, sra_stat

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            idx = 0
            metadata.df.at[idx, 'num_fastp_out'] += 2
            metadata.df.at[idx, 'bp_fastp_out'] += 250
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.fastp.fastq.gz'.format(sra_id)})
            return metadata

        observed = {'known_input_counts': None}

        def fake_run_rrna_filter(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            known_input_counts=None,
            return_files=False,
            return_file_state=False,
        ):
            observed['known_input_counts'] = known_input_counts
            idx = 0
            metadata.df.at[idx, 'num_rrna_in'] += 2
            metadata.df.at[idx, 'num_rrna_out'] += 1
            metadata.df.at[idx, 'bp_rrna_in'] += 250
            metadata.df.at[idx, 'bp_rrna_out'] += 180
            sra_stat['current_ext'] = '.rrna-filtered.fastq.gz'
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.rrna-filtered.fastq.gz'.format(sra_id)})
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.run_rrna_filter', fake_run_rrna_filter)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)

        sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert observed['known_input_counts'] == {'num_spots': 2, 'bp_total': 250}

    def test_sequence_extraction_stops_when_fastp_retains_no_reads(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            rrna_filter = True
            read_name = 'default'

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=file_state.to_set())
            return metadata, sra_stat

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files=set())
            return metadata

        def fail_rrna(*_args, **_kwargs):
            raise AssertionError('rrna filter should not run when fastp retains no reads.')

        def fail_rename(*_args, **_kwargs):
            raise AssertionError('rename_fastq should not run when fastp retains no reads.')

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.run_rrna_filter', fail_rrna)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', fail_rename)

        metadata = sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        stats_path = tmp_path / 'getfastq_stats.tsv'
        assert stats_path.exists()
        stats_df = pandas.read_csv(stats_path, sep='\t')
        assert stats_df.loc[0, 'bp_fastp_out'] == 0
        assert metadata.df.loc[0, 'bp_fastp_out'] == 0

    def test_sequence_extraction_paired_normalizes_fastp_read_counts_to_spots_for_rrna(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['paired'],
            'total_spots': [10],
            'total_bases': [2000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            rrna_filter = True
            read_name = 'default'

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 600
            metadata.df.at[idx, 'bp_written'] += 600
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(
                    work_dir=sra_stat['getfastq_sra_dir'],
                    files={'{}_1.fastq.gz'.format(sra_id), '{}_2.fastq.gz'.format(sra_id)},
                )
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=file_state.to_set())
            return metadata, sra_stat

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            idx = 0
            metadata.df.at[idx, 'num_fastp_out'] += 6
            metadata.df.at[idx, 'bp_fastp_out'] += 580
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(
                    work_dir=output_dir,
                    files={'{}_1.fastp.fastq.gz'.format(sra_id), '{}_2.fastp.fastq.gz'.format(sra_id)},
                )
            return metadata

        observed = {'known_input_counts': None}

        def fake_run_rrna_filter(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            known_input_counts=None,
            return_files=False,
            return_file_state=False,
        ):
            observed['known_input_counts'] = known_input_counts
            idx = 0
            metadata.df.at[idx, 'num_rrna_in'] += 3
            metadata.df.at[idx, 'num_rrna_out'] += 2
            metadata.df.at[idx, 'bp_rrna_in'] += 580
            metadata.df.at[idx, 'bp_rrna_out'] += 400
            sra_stat['current_ext'] = '.rrna-filtered.fastq.gz'
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(
                    work_dir=output_dir,
                    files={'{}_1.rrna-filtered.fastq.gz'.format(sra_id), '{}_2.rrna-filtered.fastq.gz'.format(sra_id)},
                )
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.run_rrna_filter', fake_run_rrna_filter)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)

        sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert observed['known_input_counts'] == {'num_spots': 3, 'bp_total': 580}

    def test_sequence_extraction_respects_rrna_first_order_for_bp_target(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            rrna_filter = True
            filter_order = 'rrna,fastp'
            read_name = 'default'

        observed_order = []

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=file_state.to_set())
            return metadata, sra_stat

        def fake_run_rrna_filter(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            known_input_counts=None,
            return_files=False,
            return_file_state=False,
        ):
            observed_order.append('rrna')
            idx = 0
            metadata.df.at[idx, 'bp_rrna_in'] += 300
            metadata.df.at[idx, 'bp_rrna_out'] += 210
            sra_stat['current_ext'] = '.rrna-filtered.fastq.gz'
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.rrna-filtered.fastq.gz'.format(sra_id)})
            return metadata

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            observed_order.append('fastp')
            idx = 0
            metadata.df.at[idx, 'bp_fastp_in'] += 210
            metadata.df.at[idx, 'bp_fastp_out'] += 180
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.fastp.fastq.gz'.format(sra_id)})
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_rrna_filter', fake_run_rrna_filter)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)

        metadata = sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert observed_order == ['rrna', 'fastp']
        assert metadata.df.loc[0, 'bp_rrna_out'] == 210
        assert metadata.df.loc[0, 'bp_fastp_out'] == 180
        assert metadata.df.loc[0, 'bp_amalgkit'] == 180

    def test_sequence_extraction_respects_custom_three_filter_order(self, tmp_path, monkeypatch):
        sra_id = 'SRR001'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        g = {'num_bp_per_sra': 1000}
        metadata = initialize_columns(metadata, g)
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'spot_length': 100,
            'total_spot': 10,
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            fastp = True
            rrna_filter = True
            contam_filter = True
            filter_order = 'contam,fastp,rrna'
            read_name = 'default'

        observed_order = []

        def fake_run_fasterq_dump(sra_stat, args, metadata, start, end, return_files=False, return_file_state=False):
            idx = 0
            metadata.df.at[idx, 'num_dumped'] += 3
            metadata.df.at[idx, 'num_written'] += 3
            metadata.df.at[idx, 'bp_dumped'] += 300
            metadata.df.at[idx, 'bp_written'] += 300
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=sra_stat['getfastq_sra_dir'], files={'{}.fastq.gz'.format(sra_id)})
            return metadata, sra_stat

        def fake_maybe_treat_paired_as_single(
            sra_stat,
            metadata,
            work_dir,
            threshold=0.99,
            num_checked_reads=2000,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, sra_stat, RunFileState(work_dir=work_dir, files=file_state.to_set())
            return metadata, sra_stat

        def fake_run_contam_filter(
            sra_stat,
            args,
            output_dir,
            metadata,
            runtime_context=None,
            files=None,
            file_state=None,
            known_input_counts=None,
            return_files=False,
            return_file_state=False,
        ):
            _ = (runtime_context, files, known_input_counts, return_files)
            observed_order.append('contam')
            idx = 0
            metadata.df.at[idx, 'bp_contam_in'] += 300
            metadata.df.at[idx, 'bp_contam_out'] += 240
            metadata.df.at[idx, 'num_contam_out'] += 3
            sra_stat['current_ext'] = '.contam-filtered.fastq.gz'
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.contam-filtered.fastq.gz'.format(sra_id)})
            return metadata

        def fake_run_fastp(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            return_files=False,
            return_file_state=False,
        ):
            _ = (files, file_state, return_files)
            observed_order.append('fastp')
            idx = 0
            metadata.df.at[idx, 'bp_fastp_in'] += 240
            metadata.df.at[idx, 'bp_fastp_out'] += 180
            metadata.df.at[idx, 'num_fastp_out'] += 3
            sra_stat['current_ext'] = '.fastp.fastq.gz'
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.fastp.fastq.gz'.format(sra_id)})
            return metadata

        def fake_run_rrna_filter(
            sra_stat,
            args,
            output_dir,
            metadata,
            files=None,
            file_state=None,
            known_input_counts=None,
            return_files=False,
            return_file_state=False,
        ):
            _ = (files, file_state, known_input_counts, return_files)
            observed_order.append('rrna')
            idx = 0
            metadata.df.at[idx, 'bp_rrna_in'] += 180
            metadata.df.at[idx, 'bp_rrna_out'] += 120
            metadata.df.at[idx, 'num_rrna_out'] += 2
            sra_stat['current_ext'] = '.rrna-filtered.fastq.gz'
            if return_file_state:
                from amalgkit.getfastq import RunFileState
                return metadata, RunFileState(work_dir=output_dir, files={'{}.rrna-filtered.fastq.gz'.format(sra_id)})
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.run_fasterq_dump', fake_run_fasterq_dump)
        monkeypatch.setattr('amalgkit.getfastq.maybe_treat_paired_as_single', fake_maybe_treat_paired_as_single)
        monkeypatch.setattr('amalgkit.getfastq.run_mmseqs_contam_filter', fake_run_contam_filter)
        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.run_rrna_filter', fake_run_rrna_filter)
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *args, **kwargs: None)

        metadata = sequence_extraction(args=Args(), sra_stat=sra_stat, metadata=metadata, g=g, start=1, end=4)

        assert observed_order == ['contam', 'fastp', 'rrna']
        assert metadata.df.loc[0, 'bp_contam_out'] == 240
        assert metadata.df.loc[0, 'bp_fastp_out'] == 180
        assert metadata.df.loc[0, 'bp_rrna_out'] == 120
        assert metadata.df.loc[0, 'bp_amalgkit'] == 120
