import pytest
import pandas
import numpy
import gzip
import json
import subprocess

import os
from types import SimpleNamespace


from amalgkit.getfastq import (
    get_range,
    get_layout,
    detect_concat_input_files,
    concat_fastq,
    remove_experiment_without_run,
    filter_getfastq_eligible_metadata,
    check_metadata_validity,
    initialize_global_params,
    rename_fastq,
    remove_old_intermediate_files,
    remove_intermediate_files,
    is_getfastq_output_present,
    initialize_columns,
    calc_2nd_ranges,
    has_remaining_spots_after_first_round,
    is_2nd_round_needed,
    write_getfastq_stats,
    calculate_requested_spots,
    collect_valid_run_ids,
    maybe_run_getfastq_second_round,
    remove_stale_getfastq_completion_manifest,
    GETFASTQ_RESUME_SCHEMA_VERSION,
    GETFASTQ_TOOL_TIMEOUT_SECONDS,
    GETFASTQ_PHASE_COMPLETE,
    GETFASTQ_PHASE_FIRST_ROUND,
    GETFASTQ_PHASE_SECOND_ROUND_IN_PROGRESS,
    build_getfastq_run_fingerprint,
    inspect_getfastq_resume_output,
    read_getfastq_run_state,
    restore_getfastq_stats,
    validate_getfastq_resume_output,
    write_getfastq_completion_manifest,
    write_getfastq_run_state,
)
from amalgkit.util import Metadata



class TestGetRange:
    def test_total_within_max(self):
        sra_stat = {'total_spot': 1000, 'num_read_per_sra': 500}
        start, end = get_range(sra_stat, offset=0, total_sra_bp=100, max_bp=200)
        assert start == 1
        assert end == 1000

    def test_total_exceeds_max_with_offset(self):
        sra_stat = {'total_spot': 10000, 'num_read_per_sra': 5000}
        start, end = get_range(sra_stat, offset=100, total_sra_bp=2000, max_bp=1000)
        assert start == 100
        assert end == 5099
        assert (end - start + 1) == 5000

    def test_total_exceeds_max_offset_too_large(self):
        sra_stat = {'total_spot': 6000, 'num_read_per_sra': 5000}
        start, end = get_range(sra_stat, offset=2000, total_sra_bp=2000, max_bp=1000)
        # total_spot > num_read_per_sra but total_spot <= num_read_per_sra + offset
        assert start == 1001
        assert end == 6000
        assert (end - start + 1) == 5000

    def test_total_spot_less_than_num_reads(self):
        sra_stat = {'total_spot': 3000, 'num_read_per_sra': 5000}
        start, end = get_range(sra_stat, offset=0, total_sra_bp=2000, max_bp=1000)
        assert start == 1
        assert end == 3000

    def test_inclusive_range_count_distinguishes_empty_and_one_spot(self):
        assert calculate_requested_spots(start=5, end=5) == 1
        assert calculate_requested_spots(start=5, end=4) == 0


# ---------------------------------------------------------------------------
# get_layout (wiki: auto-detects paired/single from metadata)
# ---------------------------------------------------------------------------

class TestGetLayout:
    def test_auto_prefers_paired(self):
        """Wiki: auto layout prefers paired when multiple layouts exist."""
        class Args:
            layout = 'auto'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2', 'R3'],
            'lib_layout': ['paired', 'single', 'paired'],
            'exclusion': ['no', 'no', 'no'],
        }))
        result = get_layout(Args(), m)
        assert result == 'paired'

    def test_auto_single_only(self):
        """When all samples are single-end, auto returns single."""
        class Args:
            layout = 'auto'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2'],
            'lib_layout': ['single', 'single'],
            'exclusion': ['no', 'no'],
        }))
        result = get_layout(Args(), m)
        assert result == 'single'

    def test_explicit_override(self):
        """Explicit layout overrides metadata."""
        class Args:
            layout = 'single'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'lib_layout': ['paired'],
            'exclusion': ['no'],
        }))
        result = get_layout(Args(), m)
        assert result == 'single'

    def test_auto_normalizes_case_and_whitespace(self):
        class Args:
            layout = 'auto'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2'],
            'lib_layout': ['  PAIRED ', 'single'],
            'exclusion': ['no', 'no'],
        }))
        result = get_layout(Args(), m)
        assert result == 'paired'

    def test_auto_raises_when_no_valid_layout(self):
        class Args:
            layout = 'auto'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2'],
            'lib_layout': [None, ''],
            'exclusion': ['no', 'no'],
        }))
        with pytest.raises(ValueError, match='No valid lib_layout'):
            get_layout(Args(), m)

    def test_auto_raises_when_lib_layout_column_missing(self):
        class Args:
            layout = 'auto'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['no'],
        }))
        m.df = m.df.drop(columns=['lib_layout'])
        with pytest.raises(ValueError, match='Column \"lib_layout\" is required'):
            get_layout(Args(), m)

    def test_explicit_invalid_layout_raises(self):
        class Args:
            layout = 'paired_end'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'lib_layout': ['paired'],
            'exclusion': ['no'],
        }))
        with pytest.raises(ValueError, match='--layout must be one of'):
            get_layout(Args(), m)


# ---------------------------------------------------------------------------
# concat_fastq
# ---------------------------------------------------------------------------

class TestConcatFastq:
    def _metadata_single(self, runs):
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': runs,
            'lib_layout': ['single'] * len(runs),
            'total_spots': [1] * len(runs),
            'spot_length': [4] * len(runs),
            'total_bases': [4] * len(runs),
            'scientific_name': ['Sp'] * len(runs),
            'exclusion': ['no'] * len(runs),
        }))

    def _metadata_paired(self, runs):
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': runs,
            'lib_layout': ['paired'] * len(runs),
            'total_spots': [1] * len(runs),
            'spot_length': [4] * len(runs),
            'total_bases': [4] * len(runs),
            'scientific_name': ['Sp'] * len(runs),
            'exclusion': ['no'] * len(runs),
        }))

    def test_single_file_uses_single_directory_scan(self, tmp_path, monkeypatch):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('ACGT\n')
        metadata = self._metadata_single(['SRR001'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'NEWID_',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}
        calls = {'num': 0}

        def fake_list_run_dir_files(work_dir):
            calls['num'] += 1
            return set(os.listdir(work_dir))

        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fake_list_run_dir_files)

        concat_fastq(args, metadata, str(tmp_path), g)

        assert calls['num'] == 1
        assert (tmp_path / 'NEWID_SRR001.amalgkit.fastq.gz').exists()

    def test_single_file_with_run_id_does_not_duplicate_prefix(self, tmp_path):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('ACGT\n')
        metadata = self._metadata_single(['SRR001'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'SRR001',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        concat_fastq(args, metadata, str(tmp_path), g)

        assert (tmp_path / 'SRR001.amalgkit.fastq.gz').exists()
        assert not (tmp_path / 'SRR001SRR001.amalgkit.fastq.gz').exists()

    def test_remove_tmp_reuses_prefetched_file_set_for_extension_lookup(self, tmp_path, monkeypatch):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        (tmp_path / 'SRR002.amalgkit.fastq.gz').write_text('CCCC\n')
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': True,
        })()
        g = {'num_bp_per_sra': 4}
        seen = []

        def fake_get_newest_intermediate_file_extension(sra_stat, work_dir, files=None):
            assert files is not None
            seen.append((sra_stat['sra_id'], files))
            return '.amalgkit.fastq.gz'

        monkeypatch.setattr(
            'amalgkit.getfastq.get_newest_intermediate_file_extension',
            fake_get_newest_intermediate_file_extension,
        )
        monkeypatch.setattr('amalgkit.getfastq.remove_intermediate_files', lambda sra_stat, ext, work_dir: None)

        concat_fastq(args, metadata, str(tmp_path), g)

        assert (tmp_path / 'MERGED.amalgkit.fastq.gz').exists()
        assert [run_id for run_id, _ in seen] == ['SRR001', 'SRR002']
        assert seen[0][1] is seen[1][1]

    def test_concat_uses_system_cat_when_available(self, tmp_path, monkeypatch):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_bytes(b'AAAA\n')
        (tmp_path / 'SRR002.amalgkit.fastq.gz').write_bytes(b'CCCC\n')
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}
        captured = {}

        monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda name: '/bin/cat' if name == 'cat' else None)

        def fake_run(cmd, stdout=None, stderr=None):
            captured['cmd'] = cmd
            stdout.write(b'AAAA\nCCCC\n')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('append_file_binary should not be used when system cat succeeds.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.append_file_binary', fail_if_called)

        concat_fastq(args, metadata, str(tmp_path), g)

        assert captured['cmd'][0] == '/bin/cat'
        assert (tmp_path / 'MERGED.amalgkit.fastq.gz').exists()
        assert (tmp_path / 'MERGED.amalgkit.fastq.gz').read_bytes() == b'AAAA\nCCCC\n'

    def test_skips_missing_run_ids_while_concatenating(self, tmp_path):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        (tmp_path / 'SRR002.amalgkit.fastq.gz').write_text('CCCC\n')
        metadata = self._metadata_single(['SRR001', numpy.nan, 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        concat_fastq(args, metadata, str(tmp_path), g)

        merged = tmp_path / 'MERGED.amalgkit.fastq.gz'
        assert merged.exists()
        assert merged.read_text() == 'AAAA\nCCCC\n'

    def test_raises_when_expected_input_fastq_is_missing(self, tmp_path, monkeypatch):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        # Simulate a race where a second input disappears after an initial directory snapshot.
        monkeypatch.setattr(
            'amalgkit.getfastq.list_run_dir_files',
            lambda _work_dir: {'SRR001.amalgkit.fastq.gz', 'SRR002.amalgkit.fastq.gz'},
        )

        with pytest.raises(FileNotFoundError, match='Dumped fastq not found'):
            concat_fastq(args, metadata, str(tmp_path), g)

    def test_single_shortcut_does_not_hide_missing_runs(self, tmp_path):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        with pytest.raises(FileNotFoundError, match='Dumped fastq not found'):
            concat_fastq(args, metadata, str(tmp_path), g)

    def test_paired_shortcut_does_not_hide_missing_mates_across_runs(self, tmp_path):
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz').write_text('AAAA\n')
        (tmp_path / 'SRR002_2.amalgkit.fastq.gz').write_text('CCCC\n')
        metadata = self._metadata_paired(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        with pytest.raises(FileNotFoundError, match='Dumped fastq not found'):
            concat_fastq(args, metadata, str(tmp_path), g)

    def test_ignores_duplicate_run_ids_in_metadata(self, tmp_path):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        metadata = self._metadata_single(['SRR001', 'SRR001'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED_',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        concat_fastq(args, metadata, str(tmp_path), g)

        merged = tmp_path / 'MERGED_SRR001.amalgkit.fastq.gz'
        assert merged.exists()
        assert merged.read_text() == 'AAAA\n'

    def test_raises_when_concat_output_path_exists_as_directory(self, tmp_path):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        (tmp_path / 'SRR002.amalgkit.fastq.gz').write_text('CCCC\n')
        (tmp_path / 'MERGED.amalgkit.fastq.gz').mkdir()
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        with pytest.raises(IsADirectoryError, match='Concatenation output path exists but is not a file'):
            concat_fastq(args, metadata, str(tmp_path), g)

    def test_raises_when_concat_input_path_is_directory(self, tmp_path, monkeypatch):
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('AAAA\n')
        (tmp_path / 'SRR002.amalgkit.fastq.gz').mkdir()
        metadata = self._metadata_single(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
        })()
        g = {'num_bp_per_sra': 4}

        monkeypatch.setattr(
            'amalgkit.getfastq.list_run_dir_files',
            lambda _work_dir: {'SRR001.amalgkit.fastq.gz', 'SRR002.amalgkit.fastq.gz'},
        )

        with pytest.raises(IsADirectoryError, match='Concatenation input path exists but is not a file'):
            concat_fastq(args, metadata, str(tmp_path), g)

    def test_paired_concat_uses_parallel_workers_when_threads_gt_one(self, tmp_path, monkeypatch):
        metadata = self._metadata_paired(['SRR001', 'SRR002'])
        args = type('Args', (), {
            'layout': 'auto',
            'id': 'MERGED',
            'id_list': None,
            'remove_tmp': False,
            'threads': 2,
        })()
        g = {'num_bp_per_sra': 4}
        observed = {'max_workers': None, 'subexts': []}

        def fake_concat_fastq_files_for_subext(
            run_ids,
            subext,
            inext,
            output_dir,
            outfile_path,
            timeout_seconds,
        ):
            assert run_ids == ['SRR001', 'SRR002']
            assert inext == '.amalgkit.fastq.gz'
            assert timeout_seconds == float(GETFASTQ_TOOL_TIMEOUT_SECONDS)
            observed['subexts'].append(subext)
            with open(outfile_path, 'wt') as out:
                out.write('dummy-{}\n'.format(subext))

        def fake_run_tasks(task_items, task_fn, max_workers=1):
            observed['max_workers'] = max_workers
            results = {}
            for item in task_items:
                results[item] = task_fn(item)
            return results, []

        monkeypatch.setattr('amalgkit.getfastq.concat_fastq_files_for_subext', fake_concat_fastq_files_for_subext)
        monkeypatch.setattr('amalgkit.getfastq.run_tasks_with_optional_threads', fake_run_tasks)
        monkeypatch.setattr(
            'amalgkit.getfastq.list_run_dir_files',
            lambda _work_dir: {
                'SRR001_1.amalgkit.fastq.gz',
                'SRR001_2.amalgkit.fastq.gz',
                'SRR002_1.amalgkit.fastq.gz',
                'SRR002_2.amalgkit.fastq.gz',
            },
        )

        concat_fastq(args, metadata, str(tmp_path), g)

        assert observed['max_workers'] == 2
        assert set(observed['subexts']) == {'_1', '_2'}
        assert (tmp_path / 'MERGED_1.amalgkit.fastq.gz').exists()
        assert (tmp_path / 'MERGED_2.amalgkit.fastq.gz').exists()


class TestDetectConcatInputFiles:
    def test_does_not_match_prefix_only_run_ids(self):
        output_files = {
            'SRR1.amalgkit.fastq.gz',
            'SRR10.amalgkit.fastq.gz',
            'SRR100_1.amalgkit.fastq.gz',
            'SRR100_2.amalgkit.fastq.gz',
        }
        detected = detect_concat_input_files(
            output_files=output_files,
            run_ids=['SRR1', 'SRR100'],
            inext='.amalgkit.fastq.gz',
        )
        assert detected == [
            'SRR1.amalgkit.fastq.gz',
            'SRR100_1.amalgkit.fastq.gz',
            'SRR100_2.amalgkit.fastq.gz',
        ]


class TestCollectValidRunIds:
    def test_filters_blank_nan_and_duplicates(self):
        values = ['SRR001', '', '   ', numpy.nan, None, 'SRR001', 'SRR002  ']
        assert collect_valid_run_ids(values, unique=True) == ['SRR001', 'SRR002']


# ---------------------------------------------------------------------------
# getfastq eligibility
# ---------------------------------------------------------------------------

class TestFilterGetfastqEligibleMetadata:
    def test_uses_only_selected_nonexcluded_rows(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003', 'SRR004'],
            'exclusion': ['no', 'low_nspots', 'no', ''],
            'is_sampled': ['yes', 'yes', 'no', ''],
        }))

        filtered = filter_getfastq_eligible_metadata(metadata)

        assert filtered.df['run'].tolist() == ['SRR001']
        assert filtered.df.index.tolist() == [0]

    def test_preserves_legacy_rows_when_eligibility_columns_are_blank(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'exclusion': ['', numpy.nan],
            'is_sampled': [None, '  '],
        }))

        filtered = filter_getfastq_eligible_metadata(metadata)

        assert filtered.df['run'].tolist() == ['SRR001', 'SRR002']

    def test_rejects_invalid_sampling_flag(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'exclusion': ['no', 'excluded'],
            'is_sampled': ['yes', 'maybe'],
        }))

        with pytest.raises(ValueError, match='invalid flag'):
            filter_getfastq_eligible_metadata(metadata)

    def test_rejects_empty_eligible_result(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'exclusion': ['low_nspots', 'no'],
            'is_sampled': ['yes', 'no'],
        }))

        with pytest.raises(ValueError, match='No eligible getfastq rows'):
            filter_getfastq_eligible_metadata(metadata)


# ---------------------------------------------------------------------------
# remove_experiment_without_run
# ---------------------------------------------------------------------------

class TestRemoveExperimentWithoutRun:
    def test_removes_empty_run(self):
        """Experiments without run IDs should be filtered out."""
        data = {
            'run': ['SRR001', '', 'SRR003'],
            'scientific_name': ['Sp1', 'Sp1', 'Sp1'],
            'exclusion': ['no', 'no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = remove_experiment_without_run(m)
        assert m.df.shape[0] == 2
        assert '' not in m.df['run'].values

    def test_no_removal_needed(self):
        """All experiments have runs, nothing removed."""
        data = {
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = remove_experiment_without_run(m)
        assert m.df.shape[0] == 2

    def test_removes_nan_and_whitespace_runs(self):
        data = {
            'run': ['SRR001', numpy.nan, None, '   ', 'SRR002  '],
            'scientific_name': ['Sp1', 'Sp1', 'Sp1', 'Sp1', 'Sp1'],
            'exclusion': ['no', 'no', 'no', 'no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = remove_experiment_without_run(m)
        assert m.df['run'].tolist() == ['SRR001', 'SRR002']

    def test_handles_all_nan_run_column_without_typeerror(self):
        data = {
            'run': [numpy.nan, numpy.nan],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = remove_experiment_without_run(m)
        assert m.df.shape[0] == 0


# ---------------------------------------------------------------------------
# check_metadata_validity (issues #96, #110: empty total_bases/total_spots)
# ---------------------------------------------------------------------------

class TestCheckMetadataValidity:
    def test_rejects_empty_metadata_table(self):
        data = {
            'run': [],
            'scientific_name': [],
            'exclusion': [],
            'total_bases': [],
            'total_spots': [],
            'spot_length': [],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='No SRA entry found'):
            check_metadata_validity(m)

    def test_rejects_missing_run_ids(self):
        data = {
            'run': ['SRR001', '   '],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
            'total_bases': [1000000, 1000000],
            'total_spots': [5000, 5000],
            'spot_length': [200, 200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='Missing Run ID\\(s\\) were detected'):
            check_metadata_validity(m)

    @pytest.mark.parametrize('unsafe_run_id', [
        '../escaped',
        '/tmp/escaped',
        'nested/run',
        'nested\\run',
        '.',
        '..',
        'bad\nrun',
    ])
    def test_rejects_unsafe_run_path_components(self, unsafe_run_id):
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [unsafe_run_id],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000],
            'total_spots': [10],
            'spot_length': [100],
        }))

        with pytest.raises(ValueError, match='Unsafe Run ID'):
            check_metadata_validity(m)

    def test_fills_missing_total_bases(self):
        """Issue #96: Empty total_bases should be filled with placeholder 999999999999."""
        data = {
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
            'total_bases': [1000000, 0],
            'total_spots': [5000, 5000],
            'spot_length': [200, 200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        assert m.df.loc[m.df['run'] == 'SRR002', 'total_bases'].values[0] == 999999999999

    def test_fills_missing_total_spots(self):
        """Issue #110: Empty total_spots should be filled with placeholder."""
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000000],
            'total_spots': [0],
            'spot_length': [200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        # total_spots was 0, should be filled with a value
        assert m.df.loc[0, 'total_spots'] > 0

    def test_valid_metadata_unchanged(self):
        """Valid metadata should pass through unchanged."""
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [2000000000],
            'total_spots': [10000000],
            'spot_length': [200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        assert m.df.loc[0, 'total_bases'] == 2000000000
        assert m.df.loc[0, 'total_spots'] == 10000000

    def test_uses_placeholder_when_spot_length_is_zero(self):
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000000],
            'total_spots': [0],
            'spot_length': [0],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        assert int(m.df.loc[0, 'total_spots']) == 999999999999

    def test_fills_non_numeric_total_bases(self):
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': ['abc'],
            'total_spots': [100],
            'spot_length': [100],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        assert int(m.df.loc[0, 'total_bases']) == 999999999999

    def test_estimates_total_spots_from_non_numeric_value(self):
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000],
            'total_spots': ['abc'],
            'spot_length': [100],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m = check_metadata_validity(m)
        assert int(m.df.loc[0, 'total_spots']) == 10

    def test_rejects_duplicate_run_ids(self):
        data = {
            'run': ['SRR001', ' SRR001 '],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
            'total_bases': [1000000, 1000000],
            'total_spots': [5000, 5000],
            'spot_length': [200, 200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='Duplicate Run ID\\(s\\) were detected'):
            check_metadata_validity(m)

    def test_rejects_missing_required_columns(self):
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000000],
            'total_spots': [5000],
            'spot_length': [200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.df = m.df.drop(columns=['total_spots'])

        with pytest.raises(ValueError, match='Missing required metadata column\\(s\\) for getfastq: total_spots'):
            check_metadata_validity(m)


# ---------------------------------------------------------------------------
# initialize_global_params (wiki: calculates per-SRA bp targets)
# ---------------------------------------------------------------------------

class TestInitializeGlobalParams:
    def test_basic_calculation(self):
        """Calculates max_bp, num_sra, num_bp_per_sra, total_sra_bp."""
        class Args:
            max_bp = '1,000,000'
        data = {
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
            'total_bases': [500000, 500000],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        g = initialize_global_params(Args(), m)
        assert g['max_bp'] == 1000000
        assert g['num_sra'] == 2
        assert g['num_bp_per_sra'] == 500000
        assert g['total_sra_bp'] == 1000000

    def test_comma_removal(self):
        """Commas in max_bp string should be removed."""
        class Args:
            max_bp = '999,999,999,999,999'
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        g = initialize_global_params(Args(), m)
        assert g['max_bp'] == 999999999999999

    def test_rejects_non_positive_max_bp(self):
        class Args:
            max_bp = '0'

        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_bases': [1000],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='--max_bp must be > 0'):
            initialize_global_params(Args(), m)

    def test_rejects_max_bp_too_small_for_num_sra(self):
        class Args:
            max_bp = '1'

        data = {
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
            'total_bases': [500, 500],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='--max_bp \\(1\\) is too small for 2 SRA runs'):
            initialize_global_params(Args(), m)


# ---------------------------------------------------------------------------
# rename_fastq (renames fastq files by extension)
# ---------------------------------------------------------------------------

class TestRenameFastq:
    @staticmethod
    def _write_fastq_gz(path, seqs):
        with gzip.open(path, 'wt') as handle:
            for idx, seq in enumerate(seqs, start=1):
                handle.write('@read{}\n'.format(idx))
                handle.write(seq + '\n')
                handle.write('+\n')
                handle.write('I' * len(seq) + '\n')

    @staticmethod
    def _write_header_fastq_gz(path, headers):
        with gzip.open(path, 'wt') as handle:
            for header in headers:
                handle.write('{}\n'.format(header))
                handle.write('ACGT\n+\nIIII\n')

    def test_rename_single(self, tmp_path):
        """Renames single-end fastq file."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        rename_fastq(sra_stat, str(tmp_path), '.fastq.gz', '.amalgkit.fastq.gz')
        assert os.path.exists(str(tmp_path / 'SRR001.amalgkit.fastq.gz'))
        assert not os.path.exists(str(tmp_path / 'SRR001.fastq.gz'))
        assert sra_stat['current_ext'] == '.amalgkit.fastq.gz'

    def test_rename_paired(self, tmp_path):
        """Renames paired-end fastq files."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        (tmp_path / 'SRR001_1.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.fastq.gz').write_text('data')
        rename_fastq(sra_stat, str(tmp_path), '.fastq.gz', '.amalgkit.fastq.gz')
        assert os.path.exists(str(tmp_path / 'SRR001_1.amalgkit.fastq.gz'))
        assert os.path.exists(str(tmp_path / 'SRR001_2.amalgkit.fastq.gz'))
        assert sra_stat['current_ext'] == '.amalgkit.fastq.gz'

    def test_rejects_directory_as_source_path(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        (tmp_path / 'SRR001.fastq.gz').mkdir()

        with pytest.raises(IsADirectoryError, match='Intermediate path exists but is not a file'):
            rename_fastq(sra_stat, str(tmp_path), '.fastq.gz', '.amalgkit.fastq.gz')

    def test_rejects_invalid_fastq_when_validation_enabled(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        input_path = tmp_path / 'SRR001.fastq.gz'
        with gzip.open(input_path, 'wt') as handle:
            handle.write('not-a-fastq-header\n')
            handle.write('ACGT\n')
            handle.write('+\n')
            handle.write('IIII\n')

        with pytest.raises(ValueError, match='FASTQ validation failed'):
            rename_fastq(
                sra_stat,
                str(tmp_path),
                '.fastq.gz',
                '.amalgkit.fastq.gz',
                validate_fastq=True,
            )

        assert input_path.exists()
        assert not (tmp_path / 'SRR001.amalgkit.fastq.gz').exists()

    def test_rejects_empty_final_fastq(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        input_path = tmp_path / 'SRR001.fastq.gz'
        with gzip.open(input_path, 'wt'):
            pass

        with pytest.raises(ValueError, match='contains no records'):
            rename_fastq(
                sra_stat,
                str(tmp_path),
                '.fastq.gz',
                '.amalgkit.fastq.gz',
                validate_fastq=True,
            )

        assert input_path.exists()
        assert not (tmp_path / 'SRR001.amalgkit.fastq.gz').exists()

    def test_rejects_paired_record_count_mismatch_when_validation_enabled(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        path1 = tmp_path / 'SRR001_1.fastq.gz'
        path2 = tmp_path / 'SRR001_2.fastq.gz'
        self._write_fastq_gz(str(path1), ['AAAA', 'CCCC'])
        self._write_fastq_gz(str(path2), ['TTTT'])

        with pytest.raises(ValueError, match='Paired FASTQ read count mismatch'):
            rename_fastq(
                sra_stat,
                str(tmp_path),
                '.fastq.gz',
                '.amalgkit.fastq.gz',
                validate_fastq=True,
            )

        assert path1.exists()
        assert path2.exists()
        assert not (tmp_path / 'SRR001_1.amalgkit.fastq.gz').exists()
        assert not (tmp_path / 'SRR001_2.amalgkit.fastq.gz').exists()

    def test_rejects_paired_sequence_quality_length_mismatch_in_single_pass(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        path1 = tmp_path / 'SRR001_1.fastq.gz'
        path2 = tmp_path / 'SRR001_2.fastq.gz'
        with gzip.open(path1, 'wt') as handle:
            handle.write('@spotA/1\nACGT\n+\nIII\n')
        self._write_header_fastq_gz(path2, ['@spotA/2'])

        with pytest.raises(ValueError, match='sequence and quality lengths differ'):
            rename_fastq(
                sra_stat,
                str(tmp_path),
                '.fastq.gz',
                '.amalgkit.fastq.gz',
                validate_fastq=True,
            )

        assert path1.exists()
        assert path2.exists()

    def test_rejects_empty_paired_final_fastq(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        path1 = tmp_path / 'SRR001_1.fastq.gz'
        path2 = tmp_path / 'SRR001_2.fastq.gz'
        with gzip.open(path1, 'wt'), gzip.open(path2, 'wt'):
            pass

        with pytest.raises(ValueError, match='contains no records'):
            rename_fastq(
                sra_stat,
                str(tmp_path),
                '.fastq.gz',
                '.amalgkit.fastq.gz',
                validate_fastq=True,
            )

        assert path1.exists()
        assert path2.exists()

    def test_rejects_equal_count_paired_mate_order_mismatch(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        path1 = tmp_path / 'SRR001_1.fastq.gz'
        path2 = tmp_path / 'SRR001_2.fastq.gz'
        self._write_header_fastq_gz(path1, ['@spotA/1', '@spotB/1'])
        self._write_header_fastq_gz(path2, ['@spotB/2', '@spotA/2'])

        with pytest.raises(ValueError, match='mate ID/order mismatch'):
            rename_fastq(
                sra_stat,
                str(tmp_path),
                '.fastq.gz',
                '.amalgkit.fastq.gz',
                validate_fastq=True,
            )

    def test_preserves_sra_spot_suffix_when_comparing_mates(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        path1 = tmp_path / 'SRR001_1.fastq.gz'
        path2 = tmp_path / 'SRR001_2.fastq.gz'
        self._write_header_fastq_gz(path1, ['@SRR001.1'])
        self._write_header_fastq_gz(path2, ['@SRR001.2'])

        with pytest.raises(ValueError, match='mate ID/order mismatch'):
            rename_fastq(
                sra_stat,
                str(tmp_path),
                '.fastq.gz',
                '.amalgkit.fastq.gz',
                validate_fastq=True,
            )

    def test_accepts_whitespace_mate_fields(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        path1 = tmp_path / 'SRR001_1.fastq.gz'
        path2 = tmp_path / 'SRR001_2.fastq.gz'
        self._write_header_fastq_gz(path1, ['@SRR001.1 1:N:0:ACGT'])
        self._write_header_fastq_gz(path2, ['@SRR001.1 2:N:0:ACGT'])

        rename_fastq(
            sra_stat,
            str(tmp_path),
            '.fastq.gz',
            '.amalgkit.fastq.gz',
            validate_fastq=True,
        )

        assert (tmp_path / 'SRR001_1.amalgkit.fastq.gz').is_file()
        assert (tmp_path / 'SRR001_2.amalgkit.fastq.gz').is_file()

    def test_rejects_wrong_whitespace_mate_field(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        path1 = tmp_path / 'SRR001_1.fastq.gz'
        path2 = tmp_path / 'SRR001_2.fastq.gz'
        self._write_header_fastq_gz(path1, ['@SRR001.1 2:N:0:ACGT'])
        self._write_header_fastq_gz(path2, ['@SRR001.1 1:N:0:ACGT'])

        with pytest.raises(ValueError, match='mate field mismatch'):
            rename_fastq(
                sra_stat,
                str(tmp_path),
                '.fastq.gz',
                '.amalgkit.fastq.gz',
                validate_fastq=True,
            )


# ---------------------------------------------------------------------------
# remove_old_intermediate_files (removes old files but keeps .sra)
# ---------------------------------------------------------------------------

class TestRemoveOldIntermediateFiles:
    def test_removes_intermediate_keeps_sra(self, tmp_path):
        """Removes intermediate files but keeps .sra files."""
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_1.fastp.fastq.gz').write_text('data')
        (tmp_path / 'SRR001.sra').write_text('data')
        remove_old_intermediate_files('SRR001', str(tmp_path))
        assert not os.path.exists(str(tmp_path / 'SRR001.fastq.gz'))
        assert not os.path.exists(str(tmp_path / 'SRR001_1.fastp.fastq.gz'))
        assert os.path.exists(str(tmp_path / 'SRR001.sra'))

    def test_no_files_to_remove(self, tmp_path):
        """No matching intermediate files -> no error."""
        (tmp_path / 'other_file.txt').write_text('data')
        remove_old_intermediate_files('SRR001', str(tmp_path))

    def test_ignores_similar_prefix_files(self, tmp_path):
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        (tmp_path / 'SRR0010.fastq.gz').write_text('data')
        remove_old_intermediate_files('SRR001', str(tmp_path))
        assert not os.path.exists(str(tmp_path / 'SRR001.fastq.gz'))
        assert os.path.exists(str(tmp_path / 'SRR0010.fastq.gz'))


# ---------------------------------------------------------------------------
# remove_intermediate_files (removes single/paired intermediate files)
# ---------------------------------------------------------------------------

class TestRemoveIntermediateFiles:
    def test_removes_single(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        remove_intermediate_files(sra_stat, '.fastq.gz', str(tmp_path))
        assert not os.path.exists(str(tmp_path / 'SRR001.fastq.gz'))

    def test_removes_paired(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        (tmp_path / 'SRR001_1.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.fastq.gz').write_text('data')
        remove_intermediate_files(sra_stat, '.fastq.gz', str(tmp_path))
        assert not os.path.exists(str(tmp_path / 'SRR001_1.fastq.gz'))
        assert not os.path.exists(str(tmp_path / 'SRR001_2.fastq.gz'))

    def test_missing_files_no_error(self, tmp_path):
        """Missing files should not raise."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        remove_intermediate_files(sra_stat, '.fastq.gz', str(tmp_path))

    def test_raises_when_intermediate_path_is_directory(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        (tmp_path / 'SRR001.fastq.gz').mkdir()

        with pytest.raises(IsADirectoryError, match='Intermediate path exists but is not a file'):
            remove_intermediate_files(sra_stat, '.fastq.gz', str(tmp_path))


# ---------------------------------------------------------------------------
# is_getfastq_output_present (checks for output files)
# ---------------------------------------------------------------------------

class TestIsGetfastqOutputPresent:
    def test_output_present_single(self, tmp_path):
        """Single-end output detected."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('data')
        observed = is_getfastq_output_present(sra_stat)
        assert observed
        assert isinstance(observed, bool)

    def test_output_present_paired(self, tmp_path):
        """Paired-end output detected."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.amalgkit.fastq.gz').write_text('data')
        assert is_getfastq_output_present(sra_stat)

    def test_safely_removed_counts(self, tmp_path):
        """Safely removed output counts as present."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.amalgkit.fastq.gz.safely_removed').write_text('')
        assert is_getfastq_output_present(sra_stat)

    def test_output_missing(self, tmp_path):
        """No output files -> False."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        assert not is_getfastq_output_present(sra_stat)


class TestGetfastqResume:
    @staticmethod
    def _make_case(tmp_path, layout='single'):
        args = SimpleNamespace(
            out_dir=str(tmp_path),
            max_bp='1000',
            min_read_length=25,
            fastp=True,
            fastp_option='-j /dev/null -h /dev/null',
            rrna_filter=True,
            rrna_filter_sensitivity=1.0,
            rrna_filter_max_seqs=20,
            rrna_filter_chunk_spots=5000000,
            rrna_filter_memory_limit='32G',
            filter_order='fastp,rrna,contam',
            contam_filter=False,
            contam_filter_rank='superkingdom',
            contam_filter_db_name='UniRef90',
            contam_filter_db='inferred',
            contam_filter_sensitivity='auto',
            contam_filter_max_seqs='auto',
            read_name='trinity',
            tol=1,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'lib_layout': [layout],
            'total_spots': [10],
            'total_bases': [1000],
            'spot_length': [100],
            'size': [1000],
            'nominal_length': [100],
            'nominal_sdev': [0],
            'scientific_name': ['Species one'],
            'taxid': [123],
        }))
        g = {
            'start_time': 0.0,
            'max_bp': 1000,
            'num_sra': 1,
            'num_bp_per_sra': 1000,
            'total_sra_bp': 1000,
        }
        metadata = initialize_columns(metadata, g)
        metadata.df.at[0, 'num_written'] = 4
        metadata.df.at[0, 'bp_written'] = 400
        metadata.df.at[0, 'num_fastp_in'] = 4
        metadata.df.at[0, 'num_fastp_out'] = 4
        metadata.df.at[0, 'bp_fastp_in'] = 400
        metadata.df.at[0, 'bp_fastp_out'] = 400
        metadata.df.at[0, 'num_rrna_in'] = 4
        metadata.df.at[0, 'num_rrna_out'] = 4
        metadata.df.at[0, 'bp_rrna_in'] = 400
        metadata.df.at[0, 'bp_rrna_out'] = 400
        metadata.df.at[0, 'bp_amalgkit'] = 400
        metadata.df.at[0, 'spot_length_amalgkit'] = 100
        metadata.df.at[0, 'spot_start_1st'] = 1
        metadata.df.at[0, 'spot_end_1st'] = 4
        metadata.df.at[0, 'bp_still_available'] = 600
        metadata.df.at[0, 'bp_specified_for_extraction'] = 400
        metadata.df.at[0, 'bp_until_target_size'] = 600
        metadata.df.at[0, 'rate_obtained'] = 0.4
        metadata.df.at[0, 'layout_amalgkit'] = layout
        run_dir = tmp_path / 'getfastq' / 'SRR001'
        run_dir.mkdir(parents=True)
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': layout,
            'total_spot': 10,
            'spot_length': 100,
            'metadata_idx': 0,
            'getfastq_sra_dir': str(run_dir),
        }
        return args, metadata, g, sra_stat, run_dir

    @staticmethod
    def _write_fastq(path, num_records=1):
        with gzip.open(path, 'wt') as handle:
            for index in range(num_records):
                handle.write('@r{}\nACGT\n+\nIIII\n'.format(index))

    @staticmethod
    def _write_fastq_ids(path, read_ids):
        with gzip.open(path, 'wt') as handle:
            for read_id in read_ids:
                handle.write('@{}\nACGT\n+\nIIII\n'.format(read_id))

    def test_adopts_valid_unmarked_output_and_restores_stats(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        self._write_fastq(run_dir / 'SRR001.amalgkit.fastq.gz')
        write_getfastq_stats(sra_stat, metadata, str(run_dir))

        resume = inspect_getfastq_resume_output(args, sra_stat, g, metadata)

        assert resume['phase'] == GETFASTQ_PHASE_FIRST_ROUND
        state = read_getfastq_run_state(str(run_dir))
        assert state['fingerprint'] == build_getfastq_run_fingerprint(args, sra_stat, g, metadata)
        reset_metadata = initialize_columns(Metadata.from_DataFrame(metadata.df.iloc[:, :9].copy()), g)
        restore_getfastq_stats(reset_metadata, sra_stat, resume['stats_row'])
        assert reset_metadata.df.at[0, 'bp_amalgkit'] == 400
        assert reset_metadata.df.at[0, 'spot_end_1st'] == 4

    def test_corrupt_changed_output_is_discarded(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        output_path = run_dir / 'SRR001.amalgkit.fastq.gz'
        self._write_fastq(output_path)
        write_getfastq_stats(sra_stat, metadata, str(run_dir))
        write_getfastq_run_state(args, sra_stat, g, metadata, GETFASTQ_PHASE_FIRST_ROUND)
        with open(output_path, 'ab') as handle:
            handle.write(b'not-a-gzip-member')

        assert inspect_getfastq_resume_output(args, sra_stat, g, metadata) is None
        assert not output_path.exists()
        assert not (run_dir / 'getfastq_stats.tsv').exists()

    def test_option_change_invalidates_only_that_run(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        output_path = run_dir / 'SRR001.amalgkit.fastq.gz'
        self._write_fastq(output_path)
        write_getfastq_stats(sra_stat, metadata, str(run_dir))
        write_getfastq_run_state(args, sra_stat, g, metadata, GETFASTQ_PHASE_COMPLETE)
        args.min_read_length = 50

        assert inspect_getfastq_resume_output(args, sra_stat, g, metadata) is None
        assert not output_path.exists()

    def test_incomplete_second_round_is_restarted_from_first_round(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        output_path = run_dir / 'SRR001.amalgkit.fastq.gz'
        self._write_fastq(output_path)
        write_getfastq_stats(sra_stat, metadata, str(run_dir))
        write_getfastq_run_state(
            args,
            sra_stat,
            g,
            metadata,
            GETFASTQ_PHASE_SECOND_ROUND_IN_PROGRESS,
        )

        assert inspect_getfastq_resume_output(args, sra_stat, g, metadata) is None
        assert not output_path.exists()

    def test_rejects_mismatched_paired_outputs(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path, layout='paired')
        self._write_fastq(run_dir / 'SRR001_1.amalgkit.fastq.gz', num_records=1)
        self._write_fastq(run_dir / 'SRR001_2.amalgkit.fastq.gz', num_records=2)
        write_getfastq_stats(sra_stat, metadata, str(run_dir))

        with pytest.raises(ValueError, match='read count mismatch'):
            validate_getfastq_resume_output(sra_stat)

    def test_paired_full_validation_uses_only_synchronized_scan(self, tmp_path, monkeypatch):
        args, metadata, g, sra_stat, run_dir = self._make_case(
            tmp_path,
            layout='paired',
        )
        _ = (args, g)
        self._write_fastq(run_dir / 'SRR001_1.amalgkit.fastq.gz', num_records=2)
        self._write_fastq(run_dir / 'SRR001_2.amalgkit.fastq.gz', num_records=2)
        write_getfastq_stats(sra_stat, metadata, str(run_dir))

        monkeypatch.setattr(
            'amalgkit.getfastq.shared_validate_fastq_structure',
            lambda _path: (_ for _ in ()).throw(
                AssertionError('paired validation must not scan each file first')
            ),
        )

        validated = validate_getfastq_resume_output(sra_stat)

        assert len(validated['outputs']) == 2

    def test_rejects_equal_count_paired_outputs_with_mismatched_ids(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path, layout='paired')
        self._write_fastq_ids(
            run_dir / 'SRR001_1.amalgkit.fastq.gz',
            ['spot-a/1', 'spot-b/1'],
        )
        self._write_fastq_ids(
            run_dir / 'SRR001_2.amalgkit.fastq.gz',
            ['spot-a/2', 'spot-c/2'],
        )
        write_getfastq_stats(sra_stat, metadata, str(run_dir))

        with pytest.raises(ValueError, match='mate ID/order mismatch'):
            validate_getfastq_resume_output(sra_stat)

    def test_rejects_empty_resume_fastq(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        with gzip.open(run_dir / 'SRR001.amalgkit.fastq.gz', 'wt'):
            pass
        write_getfastq_stats(sra_stat, metadata, str(run_dir))

        with pytest.raises(ValueError, match='contains no reads'):
            validate_getfastq_resume_output(sra_stat)

    def test_schema_one_state_is_invalidated_by_schema_two(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        output_path = run_dir / 'SRR001.amalgkit.fastq.gz'
        state_path = run_dir / 'getfastq_run_state.json'
        self._write_fastq(output_path)
        write_getfastq_stats(sra_stat, metadata, str(run_dir))
        write_getfastq_run_state(args, sra_stat, g, metadata, GETFASTQ_PHASE_COMPLETE)
        with open(state_path, encoding='utf-8') as handle:
            state = json.load(handle)
        state['schema_version'] = 1
        with open(state_path, 'w', encoding='utf-8') as handle:
            json.dump(state, handle)

        assert GETFASTQ_RESUME_SCHEMA_VERSION == 2
        assert inspect_getfastq_resume_output(args, sra_stat, g, metadata) is None
        assert not output_path.exists()
        assert not state_path.exists()
        assert not (run_dir / 'getfastq_stats.tsv').exists()

    def test_state_snapshot_records_file_identity(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        self._write_fastq(run_dir / 'SRR001.amalgkit.fastq.gz')
        write_getfastq_stats(sra_stat, metadata, str(run_dir))

        state = write_getfastq_run_state(
            args,
            sra_stat,
            g,
            metadata,
            GETFASTQ_PHASE_COMPLETE,
        )

        snapshot = state['outputs'][0]
        for field_name in [
            'dev',
            'inode',
            'ctime_ns',
            'target_dev',
            'target_inode',
            'target_ctime_ns',
        ]:
            assert isinstance(snapshot[field_name], int)
            assert snapshot[field_name] >= 0

    def test_same_size_same_mtime_replacement_is_invalidated(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        output_path = run_dir / 'SRR001.amalgkit.fastq.gz'
        replacement_path = run_dir / 'replacement.fastq.gz'
        self._write_fastq(output_path)
        write_getfastq_stats(sra_stat, metadata, str(run_dir))
        state = write_getfastq_run_state(
            args,
            sra_stat,
            g,
            metadata,
            GETFASTQ_PHASE_COMPLETE,
        )
        original_stat = os.stat(output_path)
        replacement_path.write_bytes(output_path.read_bytes())
        os.utime(
            replacement_path,
            ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns),
        )
        os.replace(replacement_path, output_path)
        replacement_stat = os.stat(output_path)
        assert replacement_stat.st_size == state['outputs'][0]['size']
        assert replacement_stat.st_mtime_ns == state['outputs'][0]['mtime_ns']
        assert replacement_stat.st_ino != state['outputs'][0]['target_inode']

        assert inspect_getfastq_resume_output(args, sra_stat, g, metadata) is None
        assert not output_path.exists()

    def test_mutated_stats_are_invalidated(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        output_path = run_dir / 'SRR001.amalgkit.fastq.gz'
        stats_path = run_dir / 'getfastq_stats.tsv'
        self._write_fastq(output_path)
        write_getfastq_stats(sra_stat, metadata, str(run_dir))
        write_getfastq_run_state(args, sra_stat, g, metadata, GETFASTQ_PHASE_COMPLETE)
        stats = pandas.read_csv(stats_path, sep='\t')
        stats.at[0, 'bp_amalgkit'] += 1
        stats.to_csv(stats_path, sep='\t', index=False)

        assert inspect_getfastq_resume_output(args, sra_stat, g, metadata) is None
        assert not output_path.exists()
        assert not stats_path.exists()

    def test_private_symlink_output_can_be_snapshotted_and_resumed(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        metadata.df['private_file'] = 'yes'
        source_path = tmp_path / 'private-source.fastq.gz'
        output_path = run_dir / 'SRR001.amalgkit.fastq.gz'
        self._write_fastq(source_path)
        output_path.symlink_to(source_path)
        write_getfastq_stats(sra_stat, metadata, str(run_dir))

        state = write_getfastq_run_state(
            args,
            sra_stat,
            g,
            metadata,
            GETFASTQ_PHASE_COMPLETE,
        )
        resume = inspect_getfastq_resume_output(args, sra_stat, g, metadata)

        assert resume['phase'] == GETFASTQ_PHASE_COMPLETE
        assert state['outputs'][0]['inode'] == os.lstat(output_path).st_ino
        assert state['outputs'][0]['target_inode'] == os.stat(source_path).st_ino

    def test_writes_global_manifest_only_for_complete_runs(self, tmp_path):
        args, metadata, g, sra_stat, run_dir = self._make_case(tmp_path)
        self._write_fastq(run_dir / 'SRR001.amalgkit.fastq.gz')
        write_getfastq_stats(sra_stat, metadata, str(run_dir))
        write_getfastq_run_state(args, sra_stat, g, metadata, GETFASTQ_PHASE_COMPLETE)

        completion_path = write_getfastq_completion_manifest(args, metadata, [(0, 'SRR001')], g)

        with open(completion_path, encoding='utf-8') as handle:
            completion = json.load(handle)
        assert completion['status'] == GETFASTQ_PHASE_COMPLETE
        assert completion['run_count'] == 1
        assert completion['runs'][0]['run'] == 'SRR001'

    def test_completion_manifest_holds_all_run_locks_in_sorted_order(
        self,
        tmp_path,
        monkeypatch,
    ):
        args = SimpleNamespace(
            out_dir=str(tmp_path),
            download_lock_dir=str(tmp_path / 'locks'),
        )
        events = []
        held_paths = set()

        class DummyLock:
            def __init__(
                self,
                lock_path,
                lock_label='Lock',
                poll_seconds=5,
                timeout_seconds=3600,
            ):
                _ = (lock_label, poll_seconds, timeout_seconds)
                self.lock_path = lock_path

            def __enter__(self):
                events.append(('enter', os.path.basename(self.lock_path)))
                held_paths.add(self.lock_path)
                return None

            def __exit__(self, exc_type, exc, tb):
                held_paths.remove(self.lock_path)
                events.append(('exit', os.path.basename(self.lock_path)))
                return False

        def fake_write_locked(**_kwargs):
            assert {
                str(tmp_path / 'locks' / 'getfastq_runs' / 'SRR001.lock'),
                str(tmp_path / 'locks' / 'getfastq_runs' / 'SRR002.lock'),
            } == held_paths
            return str(tmp_path / 'getfastq' / 'getfastq_completion.json')

        monkeypatch.setattr('amalgkit.getfastq.acquire_exclusive_lock', DummyLock)
        monkeypatch.setattr(
            'amalgkit.getfastq._write_getfastq_completion_manifest_locked',
            fake_write_locked,
        )

        write_getfastq_completion_manifest(
            args=args,
            metadata=Metadata.from_DataFrame(pandas.DataFrame()),
            run_rows=[(0, 'SRR002'), (1, 'SRR001')],
            g={},
        )

        assert events == [
            ('enter', 'SRR001.lock'),
            ('enter', 'SRR002.lock'),
            ('exit', 'SRR002.lock'),
            ('exit', 'SRR001.lock'),
        ]

    def test_stale_manifest_cleanup_refuses_symlinked_parent(self, tmp_path):
        out_dir = tmp_path / 'out'
        external_dir = tmp_path / 'external-getfastq'
        external_manifest = external_dir / 'getfastq_completion.json'
        out_dir.mkdir()
        external_dir.mkdir()
        external_manifest.write_text('external\n')
        (out_dir / 'getfastq').symlink_to(external_dir, target_is_directory=True)
        args = SimpleNamespace(out_dir=str(out_dir))

        with pytest.raises(NotADirectoryError, match='not a regular directory'):
            remove_stale_getfastq_completion_manifest(args)

        assert external_manifest.read_text() == 'external\n'

    def test_uses_prefetched_file_set(self, tmp_path, monkeypatch):
        """When files set is provided, no directory re-scan is needed."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        files = {'SRR001.amalgkit.fastq.gz'}

        def fail_if_called(_work_dir):
            raise AssertionError('list_run_dir_files should not be called when files are provided.')

        monkeypatch.setattr('amalgkit.getfastq.list_run_dir_files', fail_if_called)
        assert is_getfastq_output_present(sra_stat, files=files)

    def test_ignores_present_output_when_stats_show_zero_fastp_output(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('data')
        pandas.DataFrame([{
            'run': 'SRR001',
            'bp_written': 1000,
            'num_written': 10,
            'bp_fastp_in': 1000,
            'bp_fastp_out': 0,
            'num_fastp_out': 0,
        }]).to_csv(tmp_path / 'getfastq_stats.tsv', sep='\t', index=False)

        assert not is_getfastq_output_present(sra_stat)


# ---------------------------------------------------------------------------
# initialize_columns (initializes tracking columns for getfastq metrics)
# ---------------------------------------------------------------------------

class TestInitializeColumns:
    def test_initializes_columns(self):
        """Adds all tracking columns to metadata DataFrame."""
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_spots': [10000000],
            'total_bases': [2000000000],
            'size': [500000000],
            'nominal_length': [200],
            'nominal_sdev': [0],
            'spot_length': [200],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        g = {'num_bp_per_sra': 1000000}
        m = initialize_columns(m, g)
        assert 'bp_amalgkit' in m.df.columns
        assert 'bp_dumped' in m.df.columns
        assert 'rate_obtained' in m.df.columns
        assert 'layout_amalgkit' in m.df.columns
        assert 'fastp_duplication_rate' in m.df.columns
        assert 'fastp_insert_size_peak' in m.df.columns
        assert 'num_rrna_in' in m.df.columns
        assert 'num_rrna_out' in m.df.columns
        assert 'bp_rrna_in' in m.df.columns
        assert 'bp_rrna_out' in m.df.columns
        assert 'spot_length_amalgkit' in m.df.columns
        assert 'sec_sra_download' in m.df.columns
        assert 'sec_fasterq_dump' in m.df.columns
        assert 'sec_fastp' in m.df.columns
        assert 'sec_rrna_filter' in m.df.columns
        assert 'sec_rrna_search' in m.df.columns
        assert 'sec_rrna_rewrite' in m.df.columns
        assert 'sec_contam_filter' in m.df.columns
        assert 'sec_ete_taxonomy' in m.df.columns
        assert m.df.loc[0, 'bp_until_target_size'] == 1000000
        assert m.df.loc[0, 'bp_dumped'] == 0
        assert m.df.loc[0, 'spot_length_amalgkit'] == 0
        assert m.df.loc[0, 'layout_amalgkit'] == ''
        assert m.df.loc[0, 'sec_sra_download'] == 0.0
        assert m.df.loc[0, 'sec_fasterq_dump'] == 0.0
        assert m.df.loc[0, 'sec_fastp'] == 0.0
        assert m.df.loc[0, 'sec_rrna_filter'] == 0.0
        assert m.df.loc[0, 'sec_rrna_search'] == 0.0
        assert m.df.loc[0, 'sec_rrna_rewrite'] == 0.0
        assert m.df.loc[0, 'sec_contam_filter'] == 0.0
        assert m.df.loc[0, 'sec_ete_taxonomy'] == 0.0
        assert numpy.isnan(m.df.loc[0, 'fastp_duplication_rate'])
        assert numpy.isnan(m.df.loc[0, 'fastp_insert_size_peak'])

    def test_layout_amalgkit_column_is_dtype_safe_when_preexisting_float(self):
        data = {
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
            'total_spots': [10000000],
            'total_bases': [2000000000],
            'size': [500000000],
            'nominal_length': [200],
            'nominal_sdev': [0],
            'spot_length': [200],
            'layout_amalgkit': [numpy.nan],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.df['layout_amalgkit'] = pandas.Series([numpy.nan], dtype='float64')

        m = initialize_columns(m, {'num_bp_per_sra': 1000000})

        assert m.df.loc[0, 'layout_amalgkit'] == ''
        assert m.df['layout_amalgkit'].dtype == object


# ---------------------------------------------------------------------------
# is_2nd_round_needed (triggers compensatory extraction based on tolerated loss)
# ---------------------------------------------------------------------------

class TestIs2ndRoundNeeded:
    def test_zero_obtained_requires_second_round(self):
        assert is_2nd_round_needed(rate_obtained_1st=0.0, tol=1.0)

    def test_half_obtained_requires_second_round_at_default_tol(self):
        assert is_2nd_round_needed(rate_obtained_1st=0.5, tol=1.0)

    def test_within_default_tolerance_skips_second_round(self):
        assert not is_2nd_round_needed(rate_obtained_1st=0.99, tol=1.0)

    def test_custom_tolerance_boundary(self):
        assert not is_2nd_round_needed(rate_obtained_1st=0.95, tol=5.0)
        assert is_2nd_round_needed(rate_obtained_1st=0.9499, tol=5.0)


class TestCalc2ndRanges:
    def test_non_zero_index_uses_label_based_access(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [300.0, 600.0, 900.0],
            'rate_obtained': [0.5, numpy.nan, 1.0],
            'spot_length_amalgkit': [100.0, 200.0, 300.0],
            'total_spots': [1000, 1000, 1000],
            'spot_end_1st': [100, 200, 300],
        }))
        metadata.df.index = [10, 20, 30]

        out = calc_2nd_ranges(metadata)

        assert list(out.df.index) == [10, 20, 30]
        assert out.df.loc[10, 'spot_start_2nd'] == 101
        assert out.df.loc[20, 'spot_start_2nd'] == 201
        assert out.df.loc[30, 'spot_start_2nd'] == 301
        assert out.df.loc[10, 'spot_end_2nd'] == 106
        assert out.df.loc[20, 'spot_end_2nd'] == 203
        assert out.df.loc[30, 'spot_end_2nd'] == 303

    def test_redistributes_missing_reads_when_other_sra_has_capacity(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [50.0, 0.0],
            'rate_obtained': [1.0, 1.0],
            'spot_length_amalgkit': [1.0, 1.0],
            'total_spots': [10, 100],
            'spot_end_1st': [0, 0],
        }))

        out = calc_2nd_ranges(metadata)

        assert out.df.loc[0, 'spot_start_2nd'] == 1
        assert out.df.loc[1, 'spot_start_2nd'] == 1
        assert out.df.loc[0, 'spot_end_2nd'] == 10
        assert out.df.loc[1, 'spot_end_2nd'] == 40

    def test_zero_rate_obtained_uses_base_target_without_inf(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [100.0],
            'rate_obtained': [0.0],
            'spot_length_amalgkit': [10.0],
            'total_spots': [1000],
            'spot_end_1st': [0],
        }))

        out = calc_2nd_ranges(metadata)

        assert out.df.loc[0, 'spot_start_2nd'] == 1
        assert out.df.loc[0, 'spot_end_2nd'] == 10

    def test_never_emits_end_before_start_when_total_spots_is_smaller(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [100.0],
            'rate_obtained': [1.0],
            'spot_length_amalgkit': [10.0],
            'total_spots': [5],
            'spot_end_1st': [10],
        }))

        out = calc_2nd_ranges(metadata)

        assert out.df.loc[0, 'spot_start_2nd'] == 11
        assert out.df.loc[0, 'spot_end_2nd'] == 10

    def test_uses_spot_length_when_spot_length_amalgkit_is_missing(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [100.0],
            'rate_obtained': [1.0],
            'spot_length': [10.0],
            'total_spots': [1000],
            'spot_end_1st': [0],
        }))

        out = calc_2nd_ranges(metadata)

        assert out.df.loc[0, 'spot_start_2nd'] == 1
        assert out.df.loc[0, 'spot_end_2nd'] == 10

    def test_zero_remaining_target_is_an_empty_range(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'bp_until_target_size': [0.0],
            'rate_obtained': [1.0],
            'spot_length_amalgkit': [10.0],
            'total_spots': [1000],
            'spot_end_1st': [25],
        }))

        out = calc_2nd_ranges(metadata)

        assert out.df.loc[0, 'spot_start_2nd'] == 26
        assert out.df.loc[0, 'spot_end_2nd'] == 25
        assert calculate_requested_spots(
            out.df.loc[0, 'spot_start_2nd'],
            out.df.loc[0, 'spot_end_2nd'],
        ) == 0


class TestHasRemainingSpotsAfterFirstRound:
    def test_returns_true_when_any_spot_remains(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'total_spots': [10, 20],
            'spot_end_1st': [10, 19],
        }))
        assert has_remaining_spots_after_first_round(metadata)

    def test_returns_false_when_all_spots_are_exhausted(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'total_spots': [10, 20],
            'spot_end_1st': [10, 20],
        }))
        assert not has_remaining_spots_after_first_round(metadata)

    def test_returns_true_when_required_columns_are_missing(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'total_spots': [10],
        }))
        assert has_remaining_spots_after_first_round(metadata)


class TestMaybeRunGetfastqSecondRound:
    def test_allocates_and_extracts_only_pending_run(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'bp_amalgkit': [1000, 100],
            'bp_until_target_size': [0, 900],
            'rate_obtained': [1.0, 0.1],
            'spot_length_amalgkit': [100, 100],
            'total_spots': [10, 100],
            'spot_end_1st': [10, 2],
            'spot_start_2nd': [777, 0],
            'spot_end_2nd': [888, 0],
        }))
        metadata.df.index = [10, 20]
        args = SimpleNamespace(
            out_dir=str(tmp_path),
            download_lock_dir=str(tmp_path / 'locks'),
            tol=1,
        )
        g = {
            'max_bp': 2000,
            'num_bp_per_sra': 1000,
        }
        inspected = []
        extracted = []
        written_phases = []

        def fake_inspect_after_lock(
            args,
            metadata,
            row_index,
            sra_id,
            g,
            runtime_context=None,
        ):
            _ = (args, g, runtime_context)
            inspected.append((sra_id, row_index, list(metadata.df.index)))
            return (
                metadata,
                {
                    'sra_id': sra_id,
                    'layout': 'single',
                    'metadata_idx': row_index,
                    'getfastq_sra_dir': str(tmp_path / 'getfastq' / sra_id),
                },
                {
                    'phase': GETFASTQ_PHASE_FIRST_ROUND,
                    'stats_row': metadata.df.loc[row_index, :],
                },
                False,
            )

        def fake_second_round(
            args,
            sra_stat,
            metadata,
            g,
            runtime_context=None,
        ):
            _ = (args, g, runtime_context)
            extracted.append(sra_stat['sra_id'])
            return metadata

        def fake_write_state(
            args,
            sra_stat,
            g,
            run_metadata,
            phase,
            full_validation=True,
        ):
            _ = (args, g, run_metadata, full_validation)
            written_phases.append((sra_stat['sra_id'], phase))

        monkeypatch.setattr(
            'amalgkit.getfastq._inspect_getfastq_run_after_lock',
            fake_inspect_after_lock,
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.sequence_extraction_2nd_round',
            fake_second_round,
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.write_getfastq_run_state',
            fake_write_state,
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.write_getfastq_stats',
            lambda **_kwargs: None,
        )

        observed = maybe_run_getfastq_second_round(
            args=args,
            metadata=metadata,
            run_rows=[(10, 'SRR001'), (20, 'SRR002')],
            g=g,
            flag_private_file=True,
            flag_any_output_file_present=False,
            completion_phase_by_run={
                'SRR001': GETFASTQ_PHASE_COMPLETE,
                'SRR002': GETFASTQ_PHASE_FIRST_ROUND,
            },
        )

        assert observed.df.loc[10, 'spot_start_2nd'] == 777
        assert observed.df.loc[10, 'spot_end_2nd'] == 888
        assert observed.df.loc[20, 'spot_start_2nd'] == 3
        assert observed.df.loc[20, 'spot_end_2nd'] == 92
        assert extracted == ['SRR002']
        assert [run_id for run_id, _row_index, _indices in inspected] == [
            'SRR002',
            'SRR002',
        ]
        assert all(indices == [10, 20] for _run_id, _row_index, indices in inspected)
        assert written_phases == [
            ('SRR002', GETFASTQ_PHASE_SECOND_ROUND_IN_PROGRESS),
            ('SRR002', GETFASTQ_PHASE_COMPLETE),
        ]

    def test_parallel_second_round_merges_isolated_run_metadata(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'bp_amalgkit': [100, 100],
            'bp_until_target_size': [900, 900],
            'rate_obtained': [0.1, 0.1],
            'spot_length_amalgkit': [100, 100],
            'total_spots': [100, 100],
            'spot_end_1st': [2, 2],
            'spot_start_2nd': [0, 0],
            'spot_end_2nd': [0, 0],
        }))
        metadata.df.index = [10, 20]
        args = SimpleNamespace(
            out_dir=str(tmp_path),
            download_lock_dir=str(tmp_path / 'locks'),
            tol=1,
            internal_jobs=2,
        )
        g = {'max_bp': 2000, 'num_bp_per_sra': 1000}
        observed_workers = []

        def fake_inspect_after_lock(args, metadata, row_index, sra_id, g, runtime_context=None):
            _ = (args, g, runtime_context)
            return (
                metadata,
                {
                    'sra_id': sra_id,
                    'layout': 'single',
                    'metadata_idx': row_index,
                    'getfastq_sra_dir': str(tmp_path / 'getfastq' / sra_id),
                },
                {
                    'phase': GETFASTQ_PHASE_FIRST_ROUND,
                    'stats_row': metadata.df.loc[row_index, :],
                },
                False,
            )

        def fake_second_round(args, sra_stat, metadata, g, runtime_context=None):
            _ = (args, g, runtime_context)
            row_index = sra_stat['metadata_idx']
            metadata.df.at[row_index, 'bp_amalgkit'] = (
                501 if sra_stat['sra_id'] == 'SRR001' else 502
            )
            return metadata

        def fake_run_tasks(
            task_items,
            task_fn,
            max_workers,
            fail_fast=False,
            stop_scheduling_on_failure=None,
        ):
            _ = (fail_fast, stop_scheduling_on_failure)
            observed_workers.append(max_workers)
            return {item: task_fn(item) for item in task_items}, []

        monkeypatch.setattr(
            'amalgkit.getfastq._inspect_getfastq_run_after_lock',
            fake_inspect_after_lock,
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.sequence_extraction_2nd_round',
            fake_second_round,
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.write_getfastq_run_state',
            lambda *_args, **_kwargs: None,
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.write_getfastq_stats',
            lambda **_kwargs: None,
        )
        monkeypatch.setattr(
            'amalgkit.getfastq.run_tasks_with_optional_threads',
            fake_run_tasks,
        )

        observed = maybe_run_getfastq_second_round(
            args=args,
            metadata=metadata,
            run_rows=[(10, 'SRR001'), (20, 'SRR002')],
            g=g,
            flag_private_file=False,
            flag_any_output_file_present=False,
            completion_phase_by_run={
                'SRR001': GETFASTQ_PHASE_FIRST_ROUND,
                'SRR002': GETFASTQ_PHASE_FIRST_ROUND,
            },
        )

        assert observed_workers == [2]
        assert observed.df.loc[10, 'bp_amalgkit'] == 501
        assert observed.df.loc[20, 'bp_amalgkit'] == 502
        assert observed.df.loc[10, 'spot_start_2nd'] == 3
        assert observed.df.loc[20, 'spot_start_2nd'] == 3
