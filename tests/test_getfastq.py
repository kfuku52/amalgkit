import pytest
import pandas
import numpy
import time

import os

from amalgkit.getfastq import (
    getfastq_search_term,
    get_range,
    get_layout,
    remove_experiment_without_run,
    check_metadata_validity,
    initialize_global_params,
    rename_fastq,
    remove_old_intermediate_files,
    remove_intermediate_files,
    is_getfastq_output_present,
    initialize_columns,
)
from amalgkit.util import Metadata


class TestGetfastqSearchTerm:
    def test_id_only(self):
        assert getfastq_search_term('SRR123456') == 'SRR123456'

    def test_with_additional_term(self):
        result = getfastq_search_term('SRR123456', 'Homo sapiens[Organism]')
        assert result == 'SRR123456 AND Homo sapiens[Organism]'

    def test_none_additional(self):
        result = getfastq_search_term('PRJNA1', None)
        assert result == 'PRJNA1'


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
        assert end == 5100

    def test_total_exceeds_max_offset_too_large(self):
        sra_stat = {'total_spot': 6000, 'num_read_per_sra': 5000}
        start, end = get_range(sra_stat, offset=2000, total_sra_bp=2000, max_bp=1000)
        # total_spot > num_read_per_sra but total_spot <= num_read_per_sra + offset
        assert start == 1000  # total_spot - num_read_per_sra
        assert end == 6000

    def test_total_spot_less_than_num_reads(self):
        sra_stat = {'total_spot': 3000, 'num_read_per_sra': 5000}
        start, end = get_range(sra_stat, offset=0, total_sra_bp=2000, max_bp=1000)
        assert start == 1
        assert end == 3000


# ---------------------------------------------------------------------------
# get_layout (wiki: auto-detects paired/single from metadata)
# ---------------------------------------------------------------------------

class TestGetLayout:
    def test_auto_prefers_paired(self):
        """Wiki: auto layout prefers paired when multiple layouts exist."""
        class Args:
            layout = 'auto'
        data = {'lib_layout': ['paired', 'single', 'paired']}
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


# ---------------------------------------------------------------------------
# check_metadata_validity (issues #96, #110: empty total_bases/total_spots)
# ---------------------------------------------------------------------------

class TestCheckMetadataValidity:
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


# ---------------------------------------------------------------------------
# rename_fastq (renames fastq files by extension)
# ---------------------------------------------------------------------------

class TestRenameFastq:
    def test_rename_single(self, tmp_path):
        """Renames single-end fastq file."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        rename_fastq(sra_stat, str(tmp_path), '.fastq.gz', '.amalgkit.fastq.gz')
        assert os.path.exists(str(tmp_path / 'SRR001.amalgkit.fastq.gz'))
        assert not os.path.exists(str(tmp_path / 'SRR001.fastq.gz'))

    def test_rename_paired(self, tmp_path):
        """Renames paired-end fastq files."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        (tmp_path / 'SRR001_1.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.fastq.gz').write_text('data')
        rename_fastq(sra_stat, str(tmp_path), '.fastq.gz', '.amalgkit.fastq.gz')
        assert os.path.exists(str(tmp_path / 'SRR001_1.amalgkit.fastq.gz'))
        assert os.path.exists(str(tmp_path / 'SRR001_2.amalgkit.fastq.gz'))


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
        assert is_getfastq_output_present(sra_stat)

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
        assert m.df.loc[0, 'bp_until_target_size'] == 1000000
        assert m.df.loc[0, 'bp_dumped'] == 0
        assert m.df.loc[0, 'layout_amalgkit'] == ''
