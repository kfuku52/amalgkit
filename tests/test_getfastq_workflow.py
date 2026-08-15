import pytest
import pandas
import numpy
import gzip
import hashlib
import subprocess
from io import BytesIO

import os
import re
from types import SimpleNamespace


from amalgkit.command_context import GetfastqRuntimeContext
from amalgkit.getfastq import (
    rename_reads,
    sequence_extraction_private,
    compress_fasterq_output_files,
    getfastq_main,
    check_getfastq_dependency,
    ensure_rrna_reference_files_exist,
    ensure_mmseqs_rrna_reference_db_exists,
    ensure_mmseqs_rrna_search_index_exists,
    download_rrna_reference_gz,
    GETFASTQ_PHASE_COMPLETE,
    GETFASTQ_PHASE_FIRST_ROUND,
)
from amalgkit.util import Metadata



class TestSequenceExtractionPrivate:
    def test_uses_run_id_for_private_symlink_names(self, tmp_path, monkeypatch):
        read1_path = tmp_path / 'input_R1.fastq.gz'
        read1_path.write_text('dummy')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['Homo_sapiens_sample1'],
            'read1_path': [str(read1_path)],
            'read2_path': [numpy.nan],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'work'
        sra_dir.mkdir()
        sra_stat = {
            'sra_id': 'Homo_sapiens_sample1',
            'layout': 'single',
            'getfastq_sra_dir': str(sra_dir),
        }
        args = SimpleNamespace(fastp=False)

        monkeypatch.setattr('amalgkit.getfastq.set_current_intermediate_extension', lambda *_args, **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.get_or_detect_intermediate_extension', lambda *_args, **_kwargs: '.fastq.gz')
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *_args, **_kwargs: None)

        sequence_extraction_private(metadata=metadata, sra_stat=sra_stat, args=args)

        assert (sra_dir / 'Homo_sapiens_sample1.fastq.gz').exists()

    def test_compresses_plain_private_fastq_before_processing(self, tmp_path, monkeypatch):
        read1_path = tmp_path / 'input_R1.fastq'
        fastq = '@read1\nACGT\n+\nIIII\n'
        read1_path.write_text(fastq)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'read1_path': [str(read1_path)],
            'read2_path': [numpy.nan],
            'lib_layout': ['single'],
            'total_spots': [1],
            'total_bases': [4],
            'spot_length': [4],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'work'
        sra_dir.mkdir()
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(sra_dir),
        }
        args = SimpleNamespace(fastp=False)

        monkeypatch.setattr('amalgkit.getfastq.set_current_intermediate_extension', lambda *_args, **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.get_or_detect_intermediate_extension', lambda *_args, **_kwargs: '.fastq.gz')
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *_args, **_kwargs: None)

        sequence_extraction_private(metadata=metadata, sra_stat=sra_stat, args=args)

        output_path = sra_dir / 'SRR001.fastq.gz'
        assert output_path.is_file()
        assert not output_path.is_symlink()
        with gzip.open(output_path, 'rt') as handle:
            assert handle.read() == fastq

    def test_handles_missing_read2_path_without_type_error(self, tmp_path, monkeypatch):
        read1_path = tmp_path / 'input_R1.fastq.gz'
        read1_path.write_text('dummy')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'read1_path': [str(read1_path)],
            'read2_path': [numpy.nan],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'work'
        sra_dir.mkdir()
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(sra_dir),
        }
        args = SimpleNamespace(fastp=False)

        monkeypatch.setattr('amalgkit.getfastq.set_current_intermediate_extension', lambda *_args, **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.get_or_detect_intermediate_extension', lambda *_args, **_kwargs: '.fastq.gz')
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *_args, **_kwargs: None)

        sequence_extraction_private(metadata=metadata, sra_stat=sra_stat, args=args)

        assert (sra_dir / 'SRR001.fastq.gz').exists()

    def test_warns_when_private_path_is_directory(self, tmp_path, monkeypatch, capsys):
        read1_dir = tmp_path / 'read1_dir'
        read1_dir.mkdir()
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'read1_path': [str(read1_dir)],
            'read2_path': ['unavailable'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'work'
        sra_dir.mkdir()
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(sra_dir),
        }
        args = SimpleNamespace(fastp=False)

        monkeypatch.setattr('amalgkit.getfastq.set_current_intermediate_extension', lambda *_args, **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.get_or_detect_intermediate_extension', lambda *_args, **_kwargs: '.fastq.gz')
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *_args, **_kwargs: None)

        sequence_extraction_private(metadata=metadata, sra_stat=sra_stat, args=args)
        err = capsys.readouterr().err

        assert 'exists but is not a file' in err

    def test_raises_when_private_output_path_is_directory(self, tmp_path, monkeypatch):
        read1_path = tmp_path / 'input_R1.fastq.gz'
        read1_path.write_text('dummy')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'read1_path': [str(read1_path)],
            'read2_path': ['unavailable'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp'],
            'exclusion': ['no'],
        }))
        sra_dir = tmp_path / 'work'
        sra_dir.mkdir()
        (sra_dir / 'SRR001.fastq.gz').mkdir()
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(sra_dir),
        }
        args = SimpleNamespace(fastp=False)

        monkeypatch.setattr('amalgkit.getfastq.set_current_intermediate_extension', lambda *_args, **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.get_or_detect_intermediate_extension', lambda *_args, **_kwargs: '.fastq.gz')
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *_args, **_kwargs: None)

        with pytest.raises(IsADirectoryError, match='Private output path exists but is not a file/symlink'):
            sequence_extraction_private(metadata=metadata, sra_stat=sra_stat, args=args)

    def test_paired_fastp_counts_are_normalized_before_rrna(self, tmp_path, monkeypatch):
        read1_path = tmp_path / 'input_R1.fastq.gz'
        read2_path = tmp_path / 'input_R2.fastq.gz'
        read1_path.write_text('dummy')
        read2_path.write_text('dummy')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'read1_path': [str(read1_path)],
            'read2_path': [str(read2_path)],
            'lib_layout': ['paired'],
            'num_fastp_out': [0],
            'bp_fastp_out': [0],
            'num_rrna_in': [0],
            'num_rrna_out': [0],
            'bp_rrna_in': [0],
            'bp_rrna_out': [0],
        }))
        sra_dir = tmp_path / 'work'
        sra_dir.mkdir()
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(sra_dir),
        }
        args = SimpleNamespace(fastp=True, rrna_filter=True)
        observed = {'known_input_counts': None}

        def fake_run_fastp(sra_stat, args, output_dir, metadata, **kwargs):
            metadata.df.at[0, 'num_fastp_out'] += 6
            metadata.df.at[0, 'bp_fastp_out'] += 580
            return metadata

        def fake_run_rrna_filter(sra_stat, args, output_dir, metadata, known_input_counts=None, **kwargs):
            observed['known_input_counts'] = known_input_counts
            metadata.df.at[0, 'num_rrna_in'] += 3
            metadata.df.at[0, 'num_rrna_out'] += 2
            metadata.df.at[0, 'bp_rrna_in'] += 580
            metadata.df.at[0, 'bp_rrna_out'] += 400
            sra_stat['current_ext'] = '.rrna-filtered.fastq.gz'
            return metadata

        monkeypatch.setattr('amalgkit.getfastq.run_fastp', fake_run_fastp)
        monkeypatch.setattr('amalgkit.getfastq.run_rrna_filter', fake_run_rrna_filter)
        monkeypatch.setattr('amalgkit.getfastq.get_or_detect_intermediate_extension', lambda *_args, **_kwargs: '.rrna-filtered.fastq.gz')
        monkeypatch.setattr('amalgkit.getfastq.rename_fastq', lambda *_args, **_kwargs: None)

        sequence_extraction_private(metadata=metadata, sra_stat=sra_stat, args=args)

        assert observed['known_input_counts'] == {'num_spots': 3, 'bp_total': 580}


class TestGetfastqMainJobs:
    def test_rejects_nonpositive_jobs(self):
        args = SimpleNamespace(internal_jobs=0)
        with pytest.raises(ValueError, match='--internal_jobs must be > 0'):
            getfastq_main(args)

    @pytest.mark.parametrize(
        ('option_name', 'option_value', 'message'),
        [
            ('rrna_filter_jobs', 3, '--rrna_filter_jobs must be 1, 2'),
            ('rrna_filter_chunk_spots', 0, '--rrna_filter_chunk_spots must be an integer > 0'),
            ('rrna_filter_memory_limit', 'lots', '--rrna_filter_memory_limit must be "auto"'),
        ],
    )
    def test_rejects_invalid_rrna_resource_options_before_metadata(self, option_name, option_value, message):
        args = SimpleNamespace(
            internal_jobs='auto',
            rrna_filter=True,
            **{option_name: option_value},
        )
        with pytest.raises(ValueError, match=re.escape(message)):
            getfastq_main(args)

    def test_rrna_memory_limit_is_validated_before_metadata_or_dependencies(self, monkeypatch):
        reached = {'metadata': False, 'dependency': False}

        def fail_metadata(_args):
            reached['metadata'] = True
            raise AssertionError('metadata must not be loaded before memory validation')

        def fail_dependency(_args):
            reached['dependency'] = True
            raise AssertionError('dependencies must not be probed before memory validation')

        monkeypatch.setattr('amalgkit.getfastq.getfastq_metadata', fail_metadata)
        monkeypatch.setattr('amalgkit.getfastq.check_getfastq_dependency', fail_dependency)

        with pytest.raises(ValueError, match='uppercase'):
            getfastq_main(SimpleNamespace(
                threads=1,
                internal_jobs=1,
                rrna_filter=True,
                rrna_filter_memory_limit='32g',
            ))

        assert reached == {'metadata': False, 'dependency': False}

    def test_auto_jobs_cap_to_single_run_and_preserve_threads(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        args = SimpleNamespace(
            threads=2,
            internal_jobs='auto',
            redo=False,
            out_dir=str(tmp_path),
            remove_sra=False,
            tol=1,
            fastp=False,
            read_name='default',
        )
        observed = {}

        def fake_check_dependency(runtime_args):
            observed['checked_threads'] = runtime_args.threads
            observed['checked_jobs'] = runtime_args.internal_jobs

        def fake_process_getfastq_run(args, row_index, sra_id, run_row_df, g, runtime_context=None):
            _ = (row_index, g, runtime_context)
            observed['runtime_threads'] = args.threads
            observed['runtime_jobs'] = args.internal_jobs
            return {
                'row_index': 0,
                'sra_id': sra_id,
                'row': run_row_df.iloc[0].copy(),
                'flag_any_output_file_present': True,
                'flag_private_file': False,
                'getfastq_sra_dir': os.path.join(args.out_dir, 'getfastq', sra_id),
            }

        monkeypatch.setattr('amalgkit.getfastq.check_getfastq_dependency', fake_check_dependency)
        monkeypatch.setattr('amalgkit.getfastq.getfastq_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.getfastq.remove_experiment_without_run', lambda m: m)
        monkeypatch.setattr('amalgkit.getfastq.check_metadata_validity', lambda m: m)
        monkeypatch.setattr(
            'amalgkit.getfastq.initialize_global_params',
            lambda _args, _metadata: {
                'start_time': 0.0,
                'max_bp': 1000,
                'num_sra': 1,
                'num_bp_per_sra': 1000,
                'total_sra_bp': 1000,
            },
        )
        monkeypatch.setattr('amalgkit.getfastq.initialize_columns', lambda m, _g: m)
        monkeypatch.setattr('amalgkit.getfastq.process_getfastq_run', fake_process_getfastq_run)
        monkeypatch.setattr('amalgkit.getfastq.write_getfastq_completion_manifest', lambda **_kwargs: None)

        getfastq_main(args)

        assert observed['checked_threads'] == 2
        assert observed['checked_jobs'] == 1
        assert observed['runtime_threads'] == 2
        assert observed['runtime_jobs'] == 1

    def test_parallel_jobs_process_all_runs(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'lib_layout': ['single', 'single'],
            'total_spots': [10, 10],
            'total_bases': [1000, 1000],
            'spot_length': [100, 100],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        args = SimpleNamespace(
            internal_jobs=2,
            redo=False,
            out_dir=str(tmp_path),
            remove_sra=False,
            tol=1,
            fastp=False,
            read_name='default',
        )
        processed = []
        observed = {}
        runtime_context_ids = []

        def fake_process_getfastq_run(args, row_index, sra_id, run_row_df, g, runtime_context=None):
            assert isinstance(runtime_context, GetfastqRuntimeContext)
            observed.setdefault('runtime_threads', args.threads)
            observed.setdefault('runtime_jobs', args.internal_jobs)
            processed.append(sra_id)
            runtime_context_ids.append(id(runtime_context))
            return {
                'row_index': row_index,
                'sra_id': sra_id,
                'row': run_row_df.iloc[0].copy(),
                'flag_any_output_file_present': True,
                'flag_private_file': False,
                'getfastq_sra_dir': os.path.join(args.out_dir, 'getfastq', sra_id),
            }

        monkeypatch.setattr('amalgkit.getfastq.check_getfastq_dependency', lambda _args: None)
        monkeypatch.setattr('amalgkit.getfastq.getfastq_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.getfastq.remove_experiment_without_run', lambda m: m)
        monkeypatch.setattr('amalgkit.getfastq.check_metadata_validity', lambda m: m)
        monkeypatch.setattr(
            'amalgkit.getfastq.initialize_global_params',
            lambda _args, _metadata: {
                'start_time': 0.0,
                'max_bp': 2000,
                'num_sra': 2,
                'num_bp_per_sra': 1000,
                'total_sra_bp': 2000,
            },
        )
        monkeypatch.setattr('amalgkit.getfastq.initialize_columns', lambda m, _g: m)
        monkeypatch.setattr('amalgkit.getfastq.process_getfastq_run', fake_process_getfastq_run)
        monkeypatch.setattr('amalgkit.getfastq.write_getfastq_completion_manifest', lambda **_kwargs: None)

        getfastq_main(args)

        assert set(processed) == {'SRR001', 'SRR002'}
        assert len(set(runtime_context_ids)) == 2

    def test_cpu_budget_caps_jobs_to_serial(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'lib_layout': ['single', 'single'],
            'total_spots': [10, 10],
            'total_bases': [1000, 1000],
            'spot_length': [100, 100],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        args = SimpleNamespace(
            internal_jobs=4,
            threads=4,
            internal_cpu_budget=1,
            redo=False,
            out_dir=str(tmp_path),
            remove_sra=False,
            tol=1,
            fastp=False,
            read_name='default',
        )
        processed = []
        observed = {}

        def fake_process_getfastq_run(args, row_index, sra_id, run_row_df, g, runtime_context=None):
            _ = runtime_context
            observed.setdefault('runtime_threads', args.threads)
            observed.setdefault('runtime_jobs', args.internal_jobs)
            processed.append(sra_id)
            return {
                'row_index': row_index,
                'sra_id': sra_id,
                'row': run_row_df.iloc[0].copy(),
                'flag_any_output_file_present': True,
                'flag_private_file': False,
                'getfastq_sra_dir': os.path.join(args.out_dir, 'getfastq', sra_id),
            }

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('run_tasks_with_optional_threads should not be used when --internal_cpu_budget caps internal_jobs to 1.')

        monkeypatch.setattr('amalgkit.getfastq.check_getfastq_dependency', lambda _args: None)
        monkeypatch.setattr('amalgkit.getfastq.getfastq_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.getfastq.remove_experiment_without_run', lambda m: m)
        monkeypatch.setattr('amalgkit.getfastq.check_metadata_validity', lambda m: m)
        monkeypatch.setattr(
            'amalgkit.getfastq.initialize_global_params',
            lambda _args, _metadata: {
                'start_time': 0.0,
                'max_bp': 2000,
                'num_sra': 2,
                'num_bp_per_sra': 1000,
                'total_sra_bp': 2000,
            },
        )
        monkeypatch.setattr('amalgkit.getfastq.initialize_columns', lambda m, _g: m)
        monkeypatch.setattr('amalgkit.getfastq.process_getfastq_run', fake_process_getfastq_run)
        monkeypatch.setattr('amalgkit.getfastq.write_getfastq_completion_manifest', lambda **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.run_tasks_with_optional_threads', fail_if_called)

        getfastq_main(args)

        assert set(processed) == {'SRR001', 'SRR002'}
        assert observed['runtime_threads'] == 1
        assert observed['runtime_jobs'] == 1
        assert args.threads == 4
        assert args.internal_jobs == 4

    def test_complete_private_run_does_not_suppress_pending_public_run(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'lib_layout': ['single', 'single'],
            'total_spots': [10, 10],
            'total_bases': [1000, 1000],
            'spot_length': [100, 100],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
            'bp_amalgkit': [500, 500],
        }))
        args = SimpleNamespace(
            internal_jobs=2,
            redo=False,
            out_dir=str(tmp_path),
            remove_sra=False,
            tol=1,
            fastp=False,
            read_name='default',
        )

        def fake_process_getfastq_run(args, row_index, sra_id, run_row_df, g, runtime_context=None):
            _ = runtime_context
            return {
                'row_index': row_index,
                'sra_id': sra_id,
                'row': run_row_df.iloc[0].copy(),
                'flag_any_output_file_present': False,
                'flag_private_file': (sra_id == 'SRR001'),
                'getfastq_sra_dir': os.path.join(args.out_dir, 'getfastq', sra_id),
                'completion_phase': (
                    GETFASTQ_PHASE_COMPLETE
                    if sra_id == 'SRR001'
                    else GETFASTQ_PHASE_FIRST_ROUND
                ),
            }

        monkeypatch.setattr('amalgkit.getfastq.check_getfastq_dependency', lambda _args: None)
        monkeypatch.setattr('amalgkit.getfastq.getfastq_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.getfastq.remove_experiment_without_run', lambda m: m)
        monkeypatch.setattr('amalgkit.getfastq.check_metadata_validity', lambda m: m)
        monkeypatch.setattr(
            'amalgkit.getfastq.initialize_global_params',
            lambda _args, _metadata: {
                'start_time': 0.0,
                'max_bp': 2000,
                'num_sra': 2,
                'num_bp_per_sra': 1000,
                'total_sra_bp': 2000,
            },
        )
        monkeypatch.setattr('amalgkit.getfastq.initialize_columns', lambda m, _g: m)
        monkeypatch.setattr('amalgkit.getfastq.process_getfastq_run', fake_process_getfastq_run)
        monkeypatch.setattr('amalgkit.getfastq.print_read_stats', lambda *_args, **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.write_getfastq_run_state', lambda *_args, **_kwargs: None)
        monkeypatch.setattr('amalgkit.getfastq.write_getfastq_completion_manifest', lambda **_kwargs: None)

        def fake_inspect_after_lock(
            args,
            metadata,
            row_index,
            sra_id,
            g,
            runtime_context=None,
        ):
            _ = (args, g, runtime_context)
            phase = (
                GETFASTQ_PHASE_COMPLETE
                if sra_id == 'SRR001'
                else GETFASTQ_PHASE_FIRST_ROUND
            )
            return (
                metadata,
                {
                    'sra_id': sra_id,
                    'layout': 'single',
                    'metadata_idx': row_index,
                    'getfastq_sra_dir': str(tmp_path / 'getfastq' / sra_id),
                },
                {
                    'phase': phase,
                    'stats_row': metadata.df.loc[row_index, :],
                },
                False,
            )

        monkeypatch.setattr(
            'amalgkit.getfastq._inspect_getfastq_run_after_lock',
            fake_inspect_after_lock,
        )

        calls = {'second_round_decision': 0}

        def no_second_round(*_args, **_kwargs):
            calls['second_round_decision'] += 1
            return False

        monkeypatch.setattr('amalgkit.getfastq.is_2nd_round_needed', no_second_round)

        getfastq_main(args)
        assert calls['second_round_decision'] == 1


class TestGetfastqDependencyChecks:
    def test_handles_missing_fastp_attributes(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'3.0.10\n', stderr=b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'] == [['fasterq-dump', '--version'], ['fasterq-dump', '-h'], ['seqkit', '--help']]

    def test_requires_fasterq_dump_version_3_or_newer(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'

        def fake_run(cmd, stdout=None, stderr=None):
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'2.9.6\n', stderr=b'')
            if cmd == ['fasterq-dump', '-h']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'Usage:\n', stderr=b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        with pytest.raises(RuntimeError, match='sra-tools >= 3 is required'):
            check_getfastq_dependency(Args())

    def test_detects_missing_fasterq_spot_range_flags(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'

        def fake_run(cmd, stdout=None, stderr=None):
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'3.0.10\n', stderr=b'')
            if cmd == ['fasterq-dump', '-h']:
                return subprocess.CompletedProcess(
                    cmd,
                    0,
                    stdout=b'Usage:\\n  -x|--details\\n  -s|--split-spot\\n',
                    stderr=b'',
                )
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        args = Args()
        check_getfastq_dependency(args)
        assert args._fasterq_supports_spot_range is False

    def test_detects_fasterq_spot_range_flags_when_present(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'

        def fake_run(cmd, stdout=None, stderr=None):
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'3.0.10\n', stderr=b'')
            if cmd == ['fasterq-dump', '-h']:
                return subprocess.CompletedProcess(
                    cmd,
                    0,
                    stdout=b'Usage:\\n  -N|--minSpotId\\n  -X|--maxSpotId\\n',
                    stderr=b'',
                )
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        args = Args()
        check_getfastq_dependency(args)
        assert args._fasterq_supports_spot_range is True

    def test_uses_default_fastp_exe_when_fastp_enabled_without_override(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            fastp = True

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'3.0.10\n', stderr=b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'][0] == ['fasterq-dump', '--version']
        assert called['cmds'][1] == ['fasterq-dump', '-h']
        assert called['cmds'][2] == ['seqkit', '--help']
        assert called['cmds'][3] == ['fastp', '--help']

    def test_uses_fasterq_dump_dependency(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'default'

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'3.0.10\n', stderr=b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'][0][0] == 'fasterq-dump'

    def test_obsolete_flags_are_ignored_and_still_uses_fasterq_dump(self, monkeypatch):
        class Args:
            obsolete_pfd = True
            obsolete_pfd_exe = '/tmp/legacy_pfd_exe'
            obsolete_fastq_dump_exe = '/tmp/legacy_fastq_dump_exe'
            fasterq_dump_exe = 'fasterq-dump'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'default'

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'3.0.10\n', stderr=b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'][0][0] == 'fasterq-dump'

    def test_trinity_mode_uses_same_dependencies(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'trinity'

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'3.0.10\n', stderr=b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'][0][0] == 'fasterq-dump'

    def test_raises_clear_error_when_fasterq_dump_missing(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'missing-fasterq-dump'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'default'

        def fake_run(cmd, stdout=None, stderr=None):
            raise FileNotFoundError(cmd[0])

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        with pytest.raises(FileNotFoundError, match='fasterq-dump executable not found'):
            check_getfastq_dependency(Args())

    def test_raises_clear_error_when_seqkit_missing(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            seqkit_exe = 'missing-seqkit'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'default'

        def fake_run(cmd, stdout=None, stderr=None):
            if cmd[0] == 'missing-seqkit':
                raise FileNotFoundError(cmd[0])
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'3.0.10\n', stderr=b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        with pytest.raises(FileNotFoundError, match='seqkit executable not found'):
            check_getfastq_dependency(Args())

    def test_raises_clear_error_when_fastp_probe_fails(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            fastp = True
            fastp_exe = 'fastp'
            read_name = 'default'

        def fake_run(cmd, stdout=None, stderr=None):
            if cmd[0] == 'fastp':
                return subprocess.CompletedProcess(cmd, 127, stdout=b'', stderr=b'')
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'3.0.10\n', stderr=b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        with pytest.raises(RuntimeError, match='fastp dependency probe failed'):
            check_getfastq_dependency(Args())

    def test_rrna_filter_probes_mmseqs_and_prepares_rrna_db(self, tmp_path, monkeypatch):
        download_dir = tmp_path / 'downloads'
        args = SimpleNamespace(
            fasterq_dump_exe='fasterq-dump',
            fastp=False,
            rrna_filter=True,
            download_dir=str(download_dir),
            mmseqs_exe='mmseqs',
        )
        called = {'cmds': [], 'db_prepare': 0}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            if cmd == ['fasterq-dump', '--version']:
                return subprocess.CompletedProcess(cmd, 0, stdout=b'3.0.10\n', stderr=b'')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        def fake_prepare_db(_args):
            called['db_prepare'] += 1
            return '/tmp/mmseqs_rrna_db'

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.ensure_mmseqs_rrna_reference_db_exists', fake_prepare_db)
        check_getfastq_dependency(args)
        assert called['cmds'][0] == ['fasterq-dump', '--version']
        assert called['cmds'][1] == ['fasterq-dump', '-h']
        assert called['cmds'][2] == ['seqkit', '--help']
        assert called['cmds'][3] == ['mmseqs', '--help']
        assert called['db_prepare'] == 1
        assert download_dir.exists()
        assert download_dir.is_dir()

    def test_rrna_filter_raises_when_download_path_is_file(self, tmp_path, monkeypatch):
        bad_path = tmp_path / 'downloads_as_file'
        bad_path.write_text('not a directory')

        args = SimpleNamespace(
            fasterq_dump_exe='fasterq-dump',
            fastp=False,
            rrna_filter=True,
            download_dir=str(bad_path),
        )

        monkeypatch.setattr(
            'amalgkit.getfastq.subprocess.run',
            lambda cmd, stdout=None, stderr=None: subprocess.CompletedProcess(
                cmd,
                0,
                stdout=(b'3.0.10\n' if cmd == ['fasterq-dump', '--version'] else b''),
                stderr=b'',
            ),
        )
        with pytest.raises(NotADirectoryError, match='Download path exists but is not a directory'):
            check_getfastq_dependency(args)


class TestRrnaReferenceDownload:
    def test_downloads_silva_refs_to_default_out_dir_downloads(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        calls = {'n': 0}
        checksums = {}

        def fake_urlretrieve(url, out_path):
            calls['n'] += 1
            with gzip.open(out_path, 'wt') as fout:
                fout.write('>rRNA\nACGT\n')
            with open(out_path, 'rb') as fin:
                downloaded_bytes = fin.read()
            checksums[url + '.md5'] = hashlib.md5(  # noqa: S324 - mirrors the upstream integrity checksum
                downloaded_bytes,
                usedforsecurity=False,
            ).hexdigest()
            return (out_path, None)

        def fake_urlopen(url, timeout):
            assert timeout == 30
            filename = os.path.basename(url.removesuffix('.md5'))
            return BytesIO('{}  {}\n'.format(checksums[url], filename).encode('ascii'))

        def fake_stream_download(url, output_path, timeout_seconds, urlopen_fn=None):
            _ = (timeout_seconds, urlopen_fn)
            return fake_urlretrieve(url, output_path)[0]

        monkeypatch.setattr('amalgkit.getfastq.download_url_to_regular_file', fake_stream_download)
        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlopen', fake_urlopen)

        args = SimpleNamespace(
            out_dir=str(out_dir),
            download_dir='inferred',
        )

        refs_first = ensure_rrna_reference_files_exist(args)
        assert len(refs_first) == 2
        for ref_path in refs_first:
            assert ref_path.startswith(os.path.join(str(out_dir), 'downloads', 'silva'))
            assert os.path.exists(ref_path)
        assert os.path.isfile(out_dir / 'downloads' / 'silva' / 'silva.ready')
        assert calls['n'] == 2

        refs_second = ensure_rrna_reference_files_exist(args)
        assert refs_second == refs_first
        assert calls['n'] == 2

    def test_downloads_silva_refs_to_custom_download_dir(self, tmp_path, monkeypatch):
        custom_dir = tmp_path / 'custom_downloads'
        captured = {'lock_path': None}
        checksums = {}

        class DummyLock:
            def __init__(self, lock_path, lock_label='Lock', poll_seconds=5, timeout_seconds=3600):
                _ = (lock_label, poll_seconds, timeout_seconds)
                captured['lock_path'] = lock_path

            def __enter__(self):
                return None

            def __exit__(self, exc_type, exc, tb):
                return False

        def fake_urlretrieve(url, out_path):
            with gzip.open(out_path, 'wt') as fout:
                fout.write('>rRNA\nACGT\n')
            with open(out_path, 'rb') as fin:
                downloaded_bytes = fin.read()
            checksums[url + '.md5'] = hashlib.md5(  # noqa: S324 - mirrors the upstream integrity checksum
                downloaded_bytes,
                usedforsecurity=False,
            ).hexdigest()
            return (out_path, None)

        def fake_urlopen(url, timeout):
            assert timeout == 30
            filename = os.path.basename(url.removesuffix('.md5'))
            return BytesIO('{}  {}\n'.format(checksums[url], filename).encode('ascii'))

        monkeypatch.setattr('amalgkit.getfastq.acquire_exclusive_lock', DummyLock)
        def fake_stream_download(url, output_path, timeout_seconds, urlopen_fn=None):
            _ = (timeout_seconds, urlopen_fn)
            return fake_urlretrieve(url, output_path)[0]

        monkeypatch.setattr('amalgkit.getfastq.download_url_to_regular_file', fake_stream_download)
        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlopen', fake_urlopen)

        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(custom_dir),
        )

        refs = ensure_rrna_reference_files_exist(args)
        assert len(refs) == 2
        for ref_path in refs:
            assert ref_path.startswith(str(custom_dir / 'silva'))
            assert os.path.exists(ref_path)
        assert os.path.isfile(custom_dir / 'silva' / 'silva.ready')
        assert captured['lock_path'] == os.path.join(str(custom_dir), 'locks', 'silva_refs.lock')

    def test_checksum_mismatch_is_not_cached(self, tmp_path):
        gz_path = tmp_path / 'silva_ssu.fasta.gz'
        url = 'https://example.test/silva_ssu.fasta.gz'

        def fake_urlretrieve(_url, out_path):
            with gzip.open(out_path, 'wt') as fout:
                fout.write('>rRNA\nACGT\n')

        def fake_urlopen(_url, timeout):
            assert timeout == 30
            return BytesIO(b'00000000000000000000000000000000  silva_ssu.fasta.gz\n')

        with pytest.raises(ValueError, match='checksum mismatch'):
            download_rrna_reference_gz(
                url=url,
                gz_path=str(gz_path),
                label='SILVA SSU',
                urlretrieve_fn=fake_urlretrieve,
                checksum_urlopen_fn=fake_urlopen,
            )

        assert not gz_path.exists()
        assert list(tmp_path.glob('silva_ssu.fasta.gz.tmp.*')) == []

    def test_default_streaming_download_uses_bounded_timeout_and_checksum(self, tmp_path):
        gz_path = tmp_path / 'silva_ssu.fasta.gz'
        url = 'https://example.test/silva_ssu.fasta.gz'
        payload = gzip.compress(b'>rRNA\nACGT\n')
        expected_md5 = hashlib.md5(  # noqa: S324 - mirrors the upstream integrity checksum
            payload,
            usedforsecurity=False,
        ).hexdigest()
        observed = []

        def fake_download_urlopen(requested_url, timeout):
            observed.append((requested_url, timeout))
            return BytesIO(payload)

        def fake_checksum_urlopen(requested_url, timeout):
            observed.append((requested_url, timeout))
            return BytesIO(
                '{}  silva_ssu.fasta.gz\n'.format(expected_md5).encode('ascii')
            )

        download_rrna_reference_gz(
            url=url,
            gz_path=str(gz_path),
            label='SILVA SSU',
            checksum_urlopen_fn=fake_checksum_urlopen,
            download_timeout_seconds=7,
            download_urlopen_fn=fake_download_urlopen,
        )

        assert gz_path.read_bytes() == payload
        assert observed == [
            (url, 7.0),
            (url + '.md5', 30),
        ]

    def test_legacy_download_callback_rejects_empty_reference(self, tmp_path):
        gz_path = tmp_path / 'silva_ssu.fasta.gz'

        def fake_urlretrieve(_url, out_path):
            open(out_path, 'wb').close()

        with pytest.raises(ValueError, match='reference is empty'):
            download_rrna_reference_gz(
                url='https://example.test/silva_ssu.fasta.gz',
                gz_path=str(gz_path),
                label='SILVA SSU',
                urlretrieve_fn=fake_urlretrieve,
            )

        assert not gz_path.exists()
        assert list(tmp_path.glob('silva_ssu.fasta.gz.tmp.*')) == []


class TestMmseqsRrnaDbPreparation:
    def test_builds_mmseqs_rrna_db_with_lock(self, tmp_path, monkeypatch):
        custom_dir = tmp_path / 'custom_downloads'
        captured = {'lock_paths': [], 'createdb_cmd': None, 'createindex_cmd': None}

        class DummyLock:
            def __init__(self, lock_path, lock_label='Lock', poll_seconds=5, timeout_seconds=3600):
                _ = (lock_label, poll_seconds, timeout_seconds)
                captured['lock_paths'].append(lock_path)

            def __enter__(self):
                return None

            def __exit__(self, exc_type, exc, tb):
                return False

        ref1 = custom_dir / 'silva' / 'silva_ssu.fasta'
        ref2 = custom_dir / 'silva' / 'silva_lsu.fasta'

        def fake_prepare_refs(_args):
            os.makedirs(str(custom_dir / 'silva'), exist_ok=True)
            ref1.write_text('>ssu\nACGT\n')
            ref2.write_text('>lsu\nTGCA\n')
            return [str(ref1), str(ref2)]

        def fake_run(cmd, stdout=None, stderr=None):
            _ = (stdout, stderr)
            if cmd[1] == 'createdb':
                captured['createdb_cmd'] = cmd
                with open(cmd[3], 'wt') as fout:
                    fout.write('mmseqsdb')
                with open(cmd[3] + '.dbtype', 'wt') as fout:
                    fout.write('nucleotide')
            elif cmd[1] == 'createindex':
                captured['createindex_cmd'] = cmd
                with open(cmd[2] + '.idx', 'wt') as fout:
                    fout.write('index data')
                with open(cmd[2] + '.idx.index', 'wt') as fout:
                    fout.write('offsets')
                with open(cmd[2] + '.idx.dbtype', 'wt') as fout:
                    fout.write('index type')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.acquire_exclusive_lock', DummyLock)
        monkeypatch.setattr('amalgkit.getfastq.ensure_rrna_reference_files_exist', fake_prepare_refs)
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(custom_dir),
            download_lock_dir=str(tmp_path / 'custom_locks'),
            mmseqs_exe='mmseqs',
            dump_print=False,
        )
        db_path = ensure_mmseqs_rrna_reference_db_exists(args)
        assert captured['lock_paths'] == [
            os.path.join(
                str(tmp_path / 'custom_locks'),
                'mmseqs_db_SILVA_DB.lock',
            ),
            os.path.join(
                str(tmp_path / 'custom_locks'),
                'mmseqs_index_SILVA_DB.idx.lock',
            ),
        ]
        assert captured['createdb_cmd'][0:2] == ['mmseqs', 'createdb']
        assert captured['createdb_cmd'][2] == os.path.join(str(custom_dir), 'silva', 'silva_refs.fasta')
        assert captured['createdb_cmd'][3] == db_path
        assert db_path == os.path.join(str(custom_dir), 'mmseqs2', 'SILVA_DB')
        assert os.path.isfile(db_path + '.dbtype')
        assert os.path.isfile(db_path)
        assert os.path.isfile(db_path + '.ready')
        assert captured['createindex_cmd'][0:3] == ['mmseqs', 'createindex', db_path]
        assert captured['createindex_cmd'][captured['createindex_cmd'].index('--split-memory-limit') + 1] == '32G'
        assert os.path.isfile(db_path + '.idx')
        assert os.path.isfile(db_path + '.idx.index')
        assert os.path.isfile(db_path + '.idx.dbtype')
        assert os.path.isfile(db_path + '.idx.ready')

    def test_markerless_rrna_db_is_rebuilt(self, tmp_path, monkeypatch):
        custom_dir = tmp_path / 'custom_downloads'
        db_path = custom_dir / 'mmseqs2' / 'SILVA_DB'
        db_path.parent.mkdir(parents=True)
        db_path.write_text('untrusted-db')
        (db_path.parent / 'SILVA_DB.dbtype').write_text('untrusted-type')
        calls = []

        def fake_prepare_refs(_args):
            silva_dir = custom_dir / 'silva'
            silva_dir.mkdir(parents=True, exist_ok=True)
            ref_path = silva_dir / 'silva_ssu.fasta'
            ref_path.write_text('>ssu\nACGT\n')
            return [str(ref_path)]

        def fake_run(cmd, stdout=None, stderr=None):
            _ = (stdout, stderr)
            calls.append(cmd)
            with open(cmd[3], 'wt') as fout:
                fout.write('rebuilt-db')
            with open(cmd[3] + '.dbtype', 'wt') as fout:
                fout.write('rebuilt-type')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr(
            'amalgkit.getfastq.ensure_rrna_reference_files_exist',
            fake_prepare_refs,
        )
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr(
            'amalgkit.getfastq.ensure_mmseqs_rrna_search_index_exists',
            lambda args, db_path: db_path,
        )
        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(custom_dir),
            download_lock_dir=str(tmp_path / 'locks'),
            mmseqs_exe='mmseqs',
            dump_print=False,
        )

        observed = ensure_mmseqs_rrna_reference_db_exists(args)

        assert observed == str(db_path)
        assert len(calls) == 1
        assert calls[0][0:2] == ['mmseqs', 'createdb']
        assert db_path.read_text() == 'rebuilt-db'
        assert (db_path.parent / 'SILVA_DB.ready').is_file()

    def test_complete_markerless_rrna_index_is_rebuilt(self, tmp_path, monkeypatch):
        db_path = str(tmp_path / 'mmseqs2' / 'SILVA_DB')
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        for suffix, content in [
            ('', 'db'),
            ('.dbtype', 'nucleotide'),
            ('.idx', 'untrusted-index'),
            ('.idx.index', 'untrusted-offsets'),
            ('.idx.dbtype', 'untrusted-type'),
        ]:
            with open(db_path + suffix, 'wt') as fout:
                fout.write(content)
        calls = []

        def fake_run(cmd, stdout=None, stderr=None):
            _ = (stdout, stderr)
            calls.append(cmd)
            with open(db_path + '.idx', 'wt') as fout:
                fout.write('rebuilt-index')
            with open(db_path + '.idx.index', 'wt') as fout:
                fout.write('rebuilt-offsets')
            with open(db_path + '.idx.dbtype', 'wt') as fout:
                fout.write('rebuilt-type')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        args = SimpleNamespace(
            mmseqs_exe='mmseqs',
            threads=2,
            dump_print=False,
            rrna_filter_sensitivity=1.0,
            rrna_filter_max_seqs=20,
            rrna_filter_memory_limit='16G',
        )

        ensure_mmseqs_rrna_search_index_exists(args=args, db_path=db_path)

        assert len(calls) == 1
        assert calls[0][0:3] == ['mmseqs', 'createindex', db_path]
        with open(db_path + '.idx', 'rt') as fin:
            assert fin.read() == 'rebuilt-index'
        assert os.path.isfile(db_path + '.idx.ready')

    def test_index_guard_reuses_complete_index_and_repairs_missing_table(self, tmp_path, monkeypatch):
        db_path = str(tmp_path / 'mmseqs2' / 'SILVA_DB')
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        with open(db_path, 'wt') as fout:
            fout.write('db')
        with open(db_path + '.dbtype', 'wt') as fout:
            fout.write('nucleotide')
        calls = []

        def fake_run(cmd, stdout=None, stderr=None):
            _ = (stdout, stderr)
            calls.append(cmd)
            with open(db_path + '.idx', 'wt') as fout:
                fout.write('index data')
            with open(db_path + '.idx.index', 'wt') as fout:
                fout.write('offset table')
            with open(db_path + '.idx.dbtype', 'wt') as fout:
                fout.write('index type')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        args = SimpleNamespace(
            mmseqs_exe='mmseqs',
            threads=2,
            dump_print=False,
            rrna_filter_sensitivity=1.0,
            rrna_filter_max_seqs=20,
            rrna_filter_memory_limit='16G',
        )

        ensure_mmseqs_rrna_search_index_exists(args=args, db_path=db_path)
        ensure_mmseqs_rrna_search_index_exists(args=args, db_path=db_path)
        assert len(calls) == 1
        assert calls[0][0:3] == ['mmseqs', 'createindex', db_path]
        assert calls[0][calls[0].index('--split-memory-limit') + 1] == '16G'

        os.remove(db_path + '.idx.index')
        ensure_mmseqs_rrna_search_index_exists(args=args, db_path=db_path)
        assert len(calls) == 2
        assert os.path.isfile(db_path + '.idx.index')

        args.rrna_filter_sensitivity = 2.0
        ensure_mmseqs_rrna_search_index_exists(args=args, db_path=db_path)
        assert len(calls) == 3
        assert calls[-1][calls[-1].index('-s') + 1] == '2'


class TestFasterqSeqkitCompression:
    def _write_fastq(self, path, reads):
        with open(path, 'wt') as out:
            for i, seq in enumerate(reads):
                out.write('@r{}\n'.format(i))
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_uses_seqkit_for_compression(self, tmp_path, monkeypatch):
        sra_id = 'SRR999'
        fastq_path = tmp_path / '{}.fastq'.format(sra_id)
        self._write_fastq(str(fastq_path), ['ACGT', 'TGCA'])

        class Args:
            threads = 2
            dump_print = False
            seqkit_exe = 'seqkit'

        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
        }

        def fake_seqkit_run(cmd, stdout=None, stderr=None):
            assert cmd[0] == 'seqkit'
            assert cmd[1] == 'seq'
            out_index = cmd.index('-o')
            out_path = cmd[out_index + 1]
            in_path = cmd[out_index + 2]
            with open(in_path, 'rb') as fin, gzip.open(out_path, 'wb') as fout:
                fout.write(fin.read())
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_seqkit_run)
        compress_fasterq_output_files(sra_stat=sra_stat, args=Args())

        gz_path = tmp_path / '{}.fastq.gz'.format(sra_id)
        assert not fastq_path.exists()
        assert gz_path.exists()
        with gzip.open(str(gz_path), 'rt') as f:
            content = f.read()
        assert '@r0' in content
        assert '@r1' in content


class TestRenameReadsSeqkitCompression:
    def _write_fastq_gz(self, path, headers_and_seqs):
        with gzip.open(path, 'wt') as out:
            for header, seq in headers_and_seqs:
                out.write(header + '\n')
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_single_end_trinity_rename(self, tmp_path, monkeypatch):
        sra_id = 'SRR777'
        in_path = tmp_path / '{}.fastq.gz'.format(sra_id)
        self._write_fastq_gz(str(in_path), [('@r0 comment', 'ACGT'), ('@r1 x', 'TGCA')])
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            threads = 1
            remove_tmp = False
            dump_print = False
            seqkit_exe = 'seqkit'

        def fake_seqkit_run(cmd, stdout=None, stderr=None):
            assert cmd[0] == 'seqkit'
            assert cmd[1] == 'replace'
            assert '-p' in cmd
            assert '-r' in cmd
            out_index = cmd.index('-o')
            out_path = cmd[out_index + 1]
            in_path = cmd[-1]
            replacement = cmd[cmd.index('-r') + 1]
            suffix = replacement[len('$1'):] if replacement.startswith('$1') else ''
            in_open = gzip.open if str(in_path).endswith('.gz') else open
            out_open = gzip.open if str(out_path).endswith('.gz') else open
            with in_open(in_path, 'rt') as fin, out_open(out_path, 'wt') as fout:
                while True:
                    line1 = fin.readline()
                    if line1 == '':
                        break
                    line2 = fin.readline()
                    line3 = fin.readline()
                    line4 = fin.readline()
                    header = line1.rstrip('\n').split()[0]
                    fout.write(header + suffix + '\n')
                    fout.write(line2)
                    fout.write(line3)
                    fout.write(line4)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_seqkit_run)

        rename_reads(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path))
        out_path = tmp_path / '{}.rename.fastq.gz'.format(sra_id)
        assert out_path.exists()
        with gzip.open(str(out_path), 'rt') as f:
            lines = [f.readline().rstrip('\n') for _ in range(4)]
        assert lines[0] == '@r0/1'

    def test_single_end_trinity_rename_falls_back_to_python_when_seqkit_fails(self, tmp_path, monkeypatch):
        sra_id = 'SRR779'
        in_path = tmp_path / '{}.fastq.gz'.format(sra_id)
        self._write_fastq_gz(str(in_path), [('@r0 comment', 'ACGT'), ('@r1 x', 'TGCA')])
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            threads = 1
            remove_tmp = False
            dump_print = False
            seqkit_exe = 'seqkit'

        def fake_seqkit_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 1, stdout=b'', stderr=b'seqkit failed')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_seqkit_run)

        rename_reads(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path))

        out_path = tmp_path / '{}.rename.fastq.gz'.format(sra_id)
        assert out_path.exists()
        with gzip.open(str(out_path), 'rt') as f:
            lines = [f.readline().rstrip('\n') for _ in range(8)]
        assert lines[0] == '@r0/1'
        assert lines[4] == '@r1/1'

    def test_paired_end_trinity_rename_uses_parallel_seqkit_workers(self, tmp_path, monkeypatch):
        sra_id = 'SRR778'
        in_path1 = tmp_path / '{}_1.fastq.gz'.format(sra_id)
        in_path2 = tmp_path / '{}_2.fastq.gz'.format(sra_id)
        self._write_fastq_gz(str(in_path1), [('@r0 comment1', 'ACGT')])
        self._write_fastq_gz(str(in_path2), [('@r0 comment2', 'TGCA')])
        sra_stat = {
            'sra_id': sra_id,
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }

        class Args:
            threads = 2
            remove_tmp = False
            dump_print = False
            seqkit_exe = 'seqkit'

        observed_threads = []

        def fake_seqkit_run(cmd, stdout=None, stderr=None):
            assert cmd[0] == 'seqkit'
            assert cmd[1] == 'replace'
            observed_threads.append(cmd[cmd.index('-j') + 1])
            out_index = cmd.index('-o')
            out_path = cmd[out_index + 1]
            in_path = cmd[-1]
            replacement = cmd[cmd.index('-r') + 1]
            suffix = replacement[len('$1'):] if replacement.startswith('$1') else ''
            with gzip.open(in_path, 'rt') as fin, gzip.open(out_path, 'wt') as fout:
                while True:
                    line1 = fin.readline()
                    if line1 == '':
                        break
                    line2 = fin.readline()
                    line3 = fin.readline()
                    line4 = fin.readline()
                    header = line1.rstrip('\n').split()[0]
                    fout.write(header + suffix + '\n')
                    fout.write(line2)
                    fout.write(line3)
                    fout.write(line4)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_seqkit_run)

        rename_reads(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path))
        out_path1 = tmp_path / '{}_1.rename.fastq.gz'.format(sra_id)
        out_path2 = tmp_path / '{}_2.rename.fastq.gz'.format(sra_id)
        assert out_path1.exists()
        assert out_path2.exists()
        with gzip.open(str(out_path1), 'rt') as f:
            first_header1 = f.readline().rstrip('\n')
        with gzip.open(str(out_path2), 'rt') as f:
            first_header2 = f.readline().rstrip('\n')
        assert first_header1 == '@r0/1'
        assert first_header2 == '@r0/2'
        assert sorted(observed_threads) == ['1', '1']
