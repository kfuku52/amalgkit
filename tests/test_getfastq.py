import pytest
import pandas
import numpy
import time
import gzip
import subprocess
import xml.etree.ElementTree as ET

import os

from amalgkit.getfastq import (
    getfastq_search_term,
    getfastq_getxml,
    getfastq_metadata,
    get_range,
    get_layout,
    remove_experiment_without_run,
    check_metadata_validity,
    initialize_global_params,
    rename_fastq,
    rename_reads,
    remove_old_intermediate_files,
    remove_intermediate_files,
    is_getfastq_output_present,
    initialize_columns,
    calc_2nd_ranges,
    is_2nd_round_needed,
    get_identical_paired_ratio,
    maybe_treat_paired_as_single,
    parse_fastp_metrics,
    parse_fastp_summary_counts,
    update_fastp_metrics,
    write_fastp_stats,
    run_fastp,
    compress_fasterq_output_files,
    run_fasterq_dump,
    check_getfastq_dependency,
    remove_sra_path,
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


class TestGetfastqXmlRetrieval:
    class _DummyTree:
        def __init__(self, root):
            self._root = root

        def getroot(self):
            return self._root

    def test_returns_empty_root_when_no_records(self, monkeypatch):
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': []})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: (_ for _ in ()).throw(AssertionError('efetch should not be called')))

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE_SET'

    def test_batches_without_extra_request_on_exact_multiple(self, monkeypatch):
        id_list = ['ID{}'.format(i) for i in range(2000)]
        efetch_calls = []

        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': id_list})

        def fake_efetch(**kwargs):
            efetch_calls.append(list(kwargs['id']))
            return object()

        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', fake_efetch)
        monkeypatch.setattr(
            'amalgkit.getfastq.ET.parse',
            lambda handle: self._DummyTree(ET.Element('EXPERIMENT_PACKAGE'))
        )

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root is not None
        assert len(efetch_calls) == 2
        assert [len(c) for c in efetch_calls] == [1000, 1000]


class TestGetfastqMetadataIdListParsing:
    def test_ignores_blank_and_comment_lines(self, tmp_path, monkeypatch):
        id_list_path = tmp_path / 'ids.txt'
        id_list_path.write_text('\n'.join([
            '# comment',
            'SRR100',
            '',
            '   # indented comment',
            'SRR200  ',
            '   ',
        ]))

        called_terms = []

        def fake_getxml(search_term, retmax=1000):
            called_terms.append(search_term)
            root = ET.Element('ROOT')
            root.set('sid', search_term.split(' AND ')[0])
            return root

        def fake_from_xml(xml_root):
            sid = xml_root.get('sid')
            return Metadata.from_DataFrame(pandas.DataFrame({
                'run': [sid],
                'lib_layout': ['single'],
                'scientific_name': ['Testus species'],
                'total_bases': ['100'],
                'spot_length': ['50'],
                'exclusion': ['no'],
            }))

        class Args:
            id = None
            id_list = str(id_list_path)
            entrez_email = 'test@example.org'
            entrez_additional_search_term = None
            layout = 'single'
            sci_name = None
            metadata = 'unused'
            out_dir = str(tmp_path)

        monkeypatch.setattr('amalgkit.getfastq.getfastq_getxml', fake_getxml)
        monkeypatch.setattr('amalgkit.getfastq.Metadata.from_xml', fake_from_xml)

        metadata = getfastq_metadata(Args())

        assert called_terms == ['SRR100', 'SRR200']
        assert set(metadata.df['run'].tolist()) == {'SRR100', 'SRR200'}


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
        assert 'fastp_duplication_rate' in m.df.columns
        assert 'fastp_insert_size_peak' in m.df.columns
        assert m.df.loc[0, 'bp_until_target_size'] == 1000000
        assert m.df.loc[0, 'bp_dumped'] == 0
        assert m.df.loc[0, 'layout_amalgkit'] == ''
        assert numpy.isnan(m.df.loc[0, 'fastp_duplication_rate'])
        assert numpy.isnan(m.df.loc[0, 'fastp_insert_size_peak'])


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
        assert out.df.loc[10, 'spot_end_2nd'] == 108
        assert out.df.loc[20, 'spot_end_2nd'] == 205
        assert out.df.loc[30, 'spot_end_2nd'] == 305


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

    def test_write_fastp_stats_writes_tsv(self, tmp_path):
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
        out_path = tmp_path / 'fastp_stats.tsv'
        assert out_path.exists()
        out = pandas.read_csv(out_path, sep='\t')
        assert out.loc[0, 'run'] == 'SRR001'
        assert out.loc[0, 'fastp_duplication_rate'] == pytest.approx(22.9565)
        assert out.loc[0, 'fastp_insert_size_peak'] == pytest.approx(138.0)
        assert out.loc[0, 'num_fastp_in'] == 1000
        assert out.loc[0, 'num_fastp_out'] == 900
        assert out.loc[0, 'bp_fastp_in'] == 100000
        assert out.loc[0, 'bp_fastp_out'] == 90000


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

    def test_uses_configured_fastp_exe_and_shlex_parsing(self, tmp_path, monkeypatch):
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
            lambda sra_stat, work_dir: '.fastq.gz'
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
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=stderr_txt)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

        assert captured['cmd'][0] == 'fastp-custom'
        assert '--adapter_sequence' in captured['cmd']
        assert 'A B' in captured['cmd']
        assert metadata.df.loc[0, 'num_fastp_in'] == 10
        assert metadata.df.loc[0, 'num_fastp_out'] == 8
        assert metadata.df.loc[0, 'bp_fastp_in'] == 100
        assert metadata.df.loc[0, 'bp_fastp_out'] == 80

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
            lambda sra_stat, work_dir: '.fastq.gz'
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
            lambda sra_stat, work_dir: '.fastq.gz'
        )

        def fake_run(cmd, stdout=None, stderr=None):
            bad_stderr = '\n'.join([
                ' before filtering:',
                'total reads: 10',
            ]).encode('utf8')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=bad_stderr)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(RuntimeError, match='Unexpected fastp stderr format'):
            run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)


class TestIdenticalPairedReads:
    @staticmethod
    def _write_fastq_gz(path, seqs):
        with gzip.open(path, 'wt') as out:
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


class TestSraRecovery:
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
        return Args()

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

    def test_run_fasterq_dump_retries_once_after_redownload(self, tmp_path, monkeypatch):
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
                return subprocess.CompletedProcess(cmd, 1, stdout=b'', stderr=b'fasterq-dump failed')
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
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda _s: 4)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda x: x)
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda x: None)

        metadata, sra_stat_out = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 2
        assert redownload_calls == [True]
        assert sra_path.is_file()
        assert metadata.df.loc[0, 'num_written'] == 4
        assert metadata.df.loc[0, 'num_dumped'] == 4
        assert metadata.df.loc[0, 'num_rejected'] == 0
        assert metadata.df.loc[0, 'bp_written'] == 400
        assert metadata.df.loc[0, 'bp_dumped'] == 400
        assert metadata.df.loc[0, 'bp_rejected'] == 0
        assert sra_stat_out['layout'] == 'paired'

    def test_run_fasterq_dump_exits_when_retry_fails(self, tmp_path, monkeypatch):
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

        with pytest.raises(SystemExit):
            run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 2
        assert redownload_calls == [True]

    def test_run_fasterq_dump_no_redownload_when_first_attempt_succeeds(self, tmp_path, monkeypatch):
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
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda _s: 5)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda x: x)
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda x: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 1
        assert metadata.df.loc[0, 'num_written'] == 5
        assert metadata.df.loc[0, 'num_dumped'] == 5
        assert metadata.df.loc[0, 'num_rejected'] == 0
        assert metadata.df.loc[0, 'bp_written'] == 500
        assert metadata.df.loc[0, 'bp_dumped'] == 500

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
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda _s: 6667)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda x: x)
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda x: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert metadata.df.loc[0, 'num_written'] == 6667
        assert metadata.df.loc[0, 'num_dumped'] == 6667
        assert metadata.df.loc[0, 'num_rejected'] == 0


class TestGetfastqDependencyChecks:
    def test_uses_fasterq_dump_dependency(self, monkeypatch):
        class Args:
            fasterq_dump_exe = 'fasterq-dump'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'default'

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'][0][0] == 'fasterq-dump'

    def test_obsolete_flags_are_ignored_and_still_uses_fasterq_dump(self, monkeypatch):
        class Args:
            obsolete_pfd = True
            obsolete_pfd_exe = '/tmp/legacy_pfd_exe'
            obsolete_fastq_dump_exe = '/tmp/legacy_fastq_dump_exe'
            obsolete_prefetch_exe = '/tmp/legacy_prefetch_exe'
            fasterq_dump_exe = 'fasterq-dump'
            fastp = False
            fastp_exe = 'fastp'
            read_name = 'default'

        called = {'cmds': []}

        def fake_run(cmd, stdout=None, stderr=None):
            called['cmds'].append(cmd)
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
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        check_getfastq_dependency(Args())
        assert called['cmds'][0][0] == 'fasterq-dump'


class TestFasterqCompressionFallback:
    def _write_fastq(self, path, reads):
        with open(path, 'wt') as out:
            for i, seq in enumerate(reads):
                out.write('@r{}\n'.format(i))
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_uses_python_gzip_when_pigz_not_available(self, tmp_path, monkeypatch):
        sra_id = 'SRR999'
        fastq_path = tmp_path / '{}.fastq'.format(sra_id)
        self._write_fastq(str(fastq_path), ['ACGT', 'TGCA'])

        class Args:
            threads = 2
            dump_print = False

        sra_stat = {
            'sra_id': sra_id,
            'getfastq_sra_dir': str(tmp_path),
        }

        def fail_run(*_args, **_kwargs):
            raise AssertionError('subprocess.run should not be used when pigz is unavailable.')

        monkeypatch.setattr('amalgkit.getfastq.shutil.which', lambda _x: None)
        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fail_run)
        compress_fasterq_output_files(sra_stat=sra_stat, args=Args())

        gz_path = tmp_path / '{}.fastq.gz'.format(sra_id)
        assert not fastq_path.exists()
        assert gz_path.exists()
        with gzip.open(str(gz_path), 'rt') as f:
            content = f.read()
        assert '@r0' in content
        assert '@r1' in content


class TestRenameReadsPythonFallback:
    def _write_fastq_gz(self, path, headers_and_seqs):
        with gzip.open(path, 'wt') as out:
            for header, seq in headers_and_seqs:
                out.write(header + '\n')
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_single_end_trinity_rename(self, tmp_path):
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

        rename_reads(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path))
        out_path = tmp_path / '{}.rename.fastq.gz'.format(sra_id)
        assert out_path.exists()
        with gzip.open(str(out_path), 'rt') as f:
            lines = [f.readline().rstrip('\n') for _ in range(4)]
        assert lines[0] == '@r0/1'
