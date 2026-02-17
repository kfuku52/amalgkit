import pytest
import pandas
import numpy
import time
import gzip
import subprocess
import xml.etree.ElementTree as ET

import os
import urllib.error
from types import SimpleNamespace

from amalgkit.getfastq import (
    getfastq_search_term,
    getfastq_getxml,
    getfastq_metadata,
    get_range,
    get_layout,
    concat_fastq,
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
    sequence_extraction,
    estimate_num_written_spots_from_fastq,
    compress_fasterq_output_files,
    normalize_fasterq_size_check,
    download_sra,
    run_fasterq_dump,
    getfastq_main,
    check_getfastq_dependency,
    count_fastq_records,
    remove_sra_files,
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

class TestNormalizeFasterqSizeCheck:
    @pytest.mark.parametrize('raw, expected', [
        ('on', 'on'),
        ('off', 'off'),
        ('only', 'only'),
        (' yes ', 'on'),
        ('NO', 'off'),
        ('unexpected', 'on'),
        (True, 'on'),
        (False, 'off'),
        (1, 'on'),
        (0, 'off'),
    ])
    def test_normalize(self, raw, expected):
        assert normalize_fasterq_size_check(raw) == expected


class TestCountFastqRecords:
    @staticmethod
    def _write_fastq(path, reads):
        with open(path, 'wt') as out:
            for i, seq in enumerate(reads):
                out.write('@r{}\n'.format(i))
                out.write(seq + '\n')
                out.write('+\n')
                out.write('I' * len(seq) + '\n')

    def test_counts_records_plain_and_gz(self, tmp_path):
        plain = tmp_path / 'x.fastq'
        self._write_fastq(str(plain), ['AAAA', 'CCCC', 'GGGG'])
        gz = tmp_path / 'x.fastq.gz'
        with open(str(plain), 'rb') as fin, gzip.open(str(gz), 'wb') as fout:
            fout.write(fin.read())

        assert count_fastq_records(str(plain)) == 3
        assert count_fastq_records(str(gz)) == 3

    def test_warns_on_truncated_fastq(self, tmp_path, capsys):
        path = tmp_path / 'bad.fastq'
        with open(path, 'wt') as out:
            out.write('@r0\nAAAA\n+\nIIII\n')
            out.write('@r1\nCCCC\n+\n')  # truncated

        assert count_fastq_records(str(path)) == 1
        assert 'not divisible by 4' in capsys.readouterr().err


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

    def test_merges_package_set_chunks_without_nested_container(self, monkeypatch):
        id_list = ['ID{}'.format(i) for i in range(2000)]
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': id_list})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: object())

        def fake_parse(_handle):
            root = ET.Element('EXPERIMENT_PACKAGE_SET')
            ET.SubElement(root, 'EXPERIMENT_PACKAGE')
            return self._DummyTree(root)

        monkeypatch.setattr('amalgkit.getfastq.ET.parse', fake_parse)

        root = getfastq_getxml(search_term='SRR000000', retmax=1000)

        assert root.tag == 'EXPERIMENT_PACKAGE_SET'
        assert len(root.findall('./EXPERIMENT_PACKAGE')) == 2
        assert len(root.findall('./EXPERIMENT_PACKAGE_SET')) == 0

    def test_raises_when_error_tag_present(self, monkeypatch):
        monkeypatch.setattr('amalgkit.getfastq.Entrez.esearch', lambda **kwargs: object())
        monkeypatch.setattr('amalgkit.getfastq.Entrez.read', lambda handle: {'IdList': ['ID1']})
        monkeypatch.setattr('amalgkit.getfastq.Entrez.efetch', lambda **kwargs: object())
        err_root = ET.Element('EXPERIMENT_PACKAGE')
        err = ET.SubElement(err_root, 'Error')
        err.text = 'SRA error'
        monkeypatch.setattr(
            'amalgkit.getfastq.ET.parse',
            lambda handle: self._DummyTree(err_root)
        )

        with pytest.raises(Exception, match='<Error> found in the xml'):
            getfastq_getxml(search_term='SRR000000', retmax=1000)


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
            bad_stderr = '\n'.join([
                ' before filtering:',
                'total reads: 10',
            ]).encode('utf8')
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=bad_stderr)

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)

        with pytest.raises(RuntimeError, match='Unexpected fastp stderr format'):
            run_fastp(sra_stat=sra_stat, args=Args(), output_dir=str(tmp_path), metadata=metadata)

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
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 4)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

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
        monkeypatch.setattr('amalgkit.getfastq.estimate_num_written_spots_from_fastq', lambda *args, **kwargs: 5)
        monkeypatch.setattr('amalgkit.getfastq.detect_layout_from_file', lambda *args, **kwargs: args[0])
        monkeypatch.setattr('amalgkit.getfastq.remove_unpaired_files', lambda *args, **kwargs: None)

        metadata, _ = run_fasterq_dump(sra_stat, args, metadata, start=1, end=10)
        assert run_calls['count'] == 1
        assert metadata.df.loc[0, 'num_written'] == 5
        assert metadata.df.loc[0, 'num_dumped'] == 5
        assert metadata.df.loc[0, 'num_rejected'] == 0
        assert metadata.df.loc[0, 'bp_written'] == 500
        assert metadata.df.loc[0, 'bp_dumped'] == 500

    def test_run_fasterq_dump_reuses_trimmed_counts_without_recount(self, tmp_path, monkeypatch):
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
            return subprocess.CompletedProcess(cmd, 0, stdout=b'', stderr=b'')

        def fail_estimate(*_args, **_kwargs):
            raise AssertionError('estimate_num_written_spots_from_fastq should not be called after trim counts are available.')

        monkeypatch.setattr('amalgkit.getfastq.subprocess.run', fake_run)
        monkeypatch.setattr('amalgkit.getfastq.trim_fasterq_output_files', lambda *args, **kwargs: {'': 7, '_1': 0, '_2': 0})
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
        assert metadata.df.loc[0, 'num_dumped'] == 7
        assert metadata.df.loc[0, 'num_rejected'] == 0

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


class TestDownloadSraUrlSchemes:
    @staticmethod
    def _make_metadata(sra_id, aws_link, gcp_link, ncbi_link):
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': [sra_id],
            'AWS_Link': [aws_link],
            'GCP_Link': [gcp_link],
            'NCBI_Link': [ncbi_link],
        }))

    @staticmethod
    def _make_args(gcp_project=''):
        args = type('Args', (), {})()
        args.aws = True
        args.gcp = True
        args.ncbi = True
        args.gcp_project = gcp_project
        return args

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

        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlretrieve', fake_urlretrieve)

        with pytest.raises(AssertionError):
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

        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlretrieve', fake_urlretrieve)
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

        monkeypatch.setattr('amalgkit.getfastq.urllib.request.urlretrieve', fake_urlretrieve)

        with pytest.raises(AssertionError):
            download_sra(metadata=metadata, sra_stat=sra_stat, args=args, work_dir=str(tmp_path), overwrite=False)

        assert called_urls == ['https://storage.googleapis.com/bucket/path/to.sra?userProject=test-project']


class TestGetfastqMainJobs:
    def test_rejects_nonpositive_jobs(self):
        args = SimpleNamespace(jobs=0)
        with pytest.raises(ValueError, match='--jobs must be > 0'):
            getfastq_main(args)

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
            jobs=2,
            redo=False,
            out_dir=str(tmp_path),
            remove_sra=False,
            tol=1,
            fastp=False,
            read_name='default',
        )
        processed = []

        def fake_process_getfastq_run(args, row_index, sra_id, run_row_df, g):
            processed.append(sra_id)
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

        getfastq_main(args)

        assert set(processed) == {'SRR001', 'SRR002'}

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
            jobs=4,
            threads=4,
            cpu_budget=4,
            redo=False,
            out_dir=str(tmp_path),
            remove_sra=False,
            tol=1,
            fastp=False,
            read_name='default',
        )
        processed = []

        def fake_process_getfastq_run(args, row_index, sra_id, run_row_df, g):
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
            raise AssertionError('run_tasks_with_optional_threads should not be used when --cpu_budget caps jobs to 1.')

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
        monkeypatch.setattr('amalgkit.getfastq.run_tasks_with_optional_threads', fail_if_called)

        getfastq_main(args)

        assert set(processed) == {'SRR001', 'SRR002'}
        assert args.threads == 4

    def test_private_file_in_any_run_skips_second_round(self, tmp_path, monkeypatch):
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
            jobs=2,
            redo=False,
            out_dir=str(tmp_path),
            remove_sra=False,
            tol=1,
            fastp=False,
            read_name='default',
        )

        def fake_process_getfastq_run(args, row_index, sra_id, run_row_df, g):
            return {
                'row_index': row_index,
                'sra_id': sra_id,
                'row': run_row_df.iloc[0].copy(),
                'flag_any_output_file_present': False,
                'flag_private_file': (sra_id == 'SRR001'),
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
        monkeypatch.setattr('amalgkit.getfastq.print_read_stats', lambda *_args, **_kwargs: None)

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('Second-round extraction should be skipped when any run is private.')

        monkeypatch.setattr('amalgkit.getfastq.is_2nd_round_needed', fail_if_called)

        getfastq_main(args)


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
