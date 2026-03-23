import json
import os
import pytest
import pandas
import numpy as np
from types import SimpleNamespace

from amalgkit.exceptions import AmalgkitExit
from amalgkit.sanity import (
    list_duplicates,
    parse_metadata,
    check_quant_output,
    check_quant_index,
    check_getfastq_outputs,
    run_sanity_check_quant,
    sanity_main,
)
from amalgkit.util import Metadata


def make_sanity_args(tmp_path, **overrides):
    defaults = {
        'out_dir': str(tmp_path / 'out'),
        'metadata': 'metadata.tsv',
        'all': False,
        'getfastq': False,
        'index': False,
        'quant': False,
        'merge': False,
        'busco': False,
        'finalize': False,
        'check': None,
        'run': None,
        'species': None,
        'getfastq_dir': None,
        'index_dir': None,
        'quant_dir': None,
        'merge_dir': None,
        'busco_dir': None,
        'finalize_dir': None,
        'quiet': True,
        'verbose_runs': 0,
        'threads': 1,
        'strict': False,
        'strict_level': 'none',
    }
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


class TestListDuplicates:
    def test_finds_duplicates(self):
        assert sorted(list_duplicates([1, 2, 3, 2, 4, 3])) == [2, 3]

    def test_no_duplicates(self):
        assert list_duplicates([1, 2, 3, 4]) == []

    def test_empty_list(self):
        assert list_duplicates([]) == []

    def test_all_duplicates(self):
        assert sorted(list_duplicates(['a', 'a', 'b', 'b'])) == ['a', 'b']

    def test_single_element(self):
        assert list_duplicates([42]) == []


# ---------------------------------------------------------------------------
# parse_metadata
# ---------------------------------------------------------------------------

class TestParseMetadata:
    def test_valid_metadata(self, sample_metadata):
        class Args:
            metadata = 'metadata.tsv'
        uni_species, sra_ids = parse_metadata(Args(), sample_metadata)
        assert len(uni_species) == 2  # Homo sapiens, Mus musculus
        assert len(sra_ids) == 5

    def test_duplicate_sra_ids_raises(self):
        """Issue: Duplicate SRA IDs should raise ValueError."""
        data = {
            'scientific_name': ['Sp1', 'Sp1'],
            'run': ['SRR001', 'SRR001'],  # duplicate
            'exclusion': ['no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        class Args:
            metadata = 'metadata.tsv'
        with pytest.raises(ValueError):
            parse_metadata(Args(), m)

    def test_raises_on_missing_run_ids(self):
        data = {
            'scientific_name': ['Sp1', 'Sp1', 'Sp1'],
            'run': ['SRR001', np.nan, '  '],
            'exclusion': ['no', 'no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))

        class Args:
            metadata = 'metadata.tsv'

        with pytest.raises(ValueError, match='without Run ID'):
            parse_metadata(Args(), m)

    def test_ignores_blank_scientific_names(self):
        data = {
            'scientific_name': ['Sp1', '', np.nan],
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'exclusion': ['no', 'no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))

        class Args:
            metadata = 'metadata.tsv'

        uni_species, sra_ids = parse_metadata(Args(), m)
        assert uni_species.tolist() == ['Sp1']
        assert sra_ids.tolist() == ['SRR001', 'SRR002', 'SRR003']

    def test_raises_when_scientific_name_column_missing(self):
        data = {
            'run': ['SRR001'],
            'exclusion': ['no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.df = m.df.drop(columns=['scientific_name'])

        class Args:
            metadata = 'metadata.tsv'

        with pytest.raises(ValueError, match='Missing required metadata column\\(s\\): scientific_name'):
            parse_metadata(Args(), m)

    def test_raises_when_run_column_missing(self):
        data = {
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.df = m.df.drop(columns=['run'])

        class Args:
            metadata = 'metadata.tsv'

        with pytest.raises(ValueError, match='Missing required metadata column\\(s\\): run'):
            parse_metadata(Args(), m)


# ---------------------------------------------------------------------------
# check_quant_output (filesystem checks)
# ---------------------------------------------------------------------------

class TestCheckQuantOutput:
    def test_all_outputs_present(self, tmp_path, sample_metadata):
        """When all quant outputs exist, all SRA IDs should be available."""
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'
        quant_dir = tmp_path / 'quant'
        # Create output files for first two SRA IDs
        for sra_id in ['SRR001', 'SRR002']:
            sra_dir = quant_dir / sra_id
            sra_dir.mkdir(parents=True)
            (sra_dir / f'{sra_id}_abundance.tsv').write_text('data')
            (sra_dir / f'{sra_id}_run_info.json').write_text('{}')

        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001', 'SRR002'])
        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))
        assert avail == ['SRR001', 'SRR002']
        assert unavail == []

    def test_missing_outputs(self, tmp_path):
        """Missing quant output should be reported as unavailable."""
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'
        quant_dir = tmp_path / 'quant'
        quant_dir.mkdir()
        # SRR001 has no output directory
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))
        assert avail == []
        assert unavail == ['SRR001']

    def test_partial_outputs(self, tmp_path):
        """Missing abundance.tsv should mark SRA as unavailable."""
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        # Only run_info.json present, missing abundance.tsv
        (sra_dir / 'SRR001_run_info.json').write_text('{}')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))
        assert 'SRR001' in unavail

    def test_missing_run_info_marks_unavailable(self, tmp_path):
        """Missing run_info.json should mark SRA as unavailable."""
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_abundance.tsv').write_text('data')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))
        assert avail == []
        assert unavail == ['SRR001']

    def test_directory_named_quant_output_file_is_ignored(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_abundance.tsv').mkdir()
        (sra_dir / 'SRR001_run_info.json').mkdir()
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))
        assert avail == []
        assert unavail == ['SRR001']

    def test_sra_entry_that_is_not_directory_is_unavailable(self, tmp_path):
        """If quant/SRRxxx is a file, treat output as unavailable without crashing."""
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'
        quant_dir = tmp_path / 'quant'
        quant_dir.mkdir()
        (quant_dir / 'SRR001').write_text('not a directory')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))
        assert avail == []
        assert unavail == ['SRR001']

    def test_no_quant_directory(self, tmp_path):
        """No quant directory at all."""
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))
        assert avail == []
        assert unavail == ['SRR001']

    def test_respects_quant_dir_argument(self, tmp_path):
        """Custom --quant_dir should be used instead of out_dir/quant."""
        class Args:
            out_dir = str(tmp_path / 'out')
            quant_dir = str(tmp_path / 'custom_quant')
            metadata = 'metadata.tsv'

        custom_quant_dir = tmp_path / 'custom_quant'
        sra_dir = custom_quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_abundance.tsv').write_text('data')
        (sra_dir / 'SRR001_run_info.json').write_text('{}')

        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))
        assert avail == ['SRR001']
        assert unavail == []

    def test_quant_dir_file_path_marks_runs_unavailable(self, tmp_path):
        """If --quant_dir points to a file, mark all runs unavailable without crashing."""
        class Args:
            out_dir = str(tmp_path / 'out')
            quant_dir = str(tmp_path / 'quant_path')
            metadata = 'metadata.tsv'
        (tmp_path / 'quant_path').write_text('not a directory')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))
        assert avail == []
        assert unavail == ['SRR001']

    def test_suppresses_per_run_logs_when_verbose_runs_exceeded(self, tmp_path, capsys):
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'
            quiet = False
            verbose_runs = 1
        quant_dir = tmp_path / 'quant'
        quant_dir.mkdir()
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001', 'SRR002'])

        check_quant_output(Args(), sra_ids, str(output_dir))

        stdout = capsys.readouterr().out
        assert 'Per-run logs suppressed for 2 runs' in stdout
        assert 'Looking for SRR001' not in stdout

    def test_ignores_missing_run_ids(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'

        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001', np.nan, '  '])

        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))

        assert avail == []
        assert unavail == ['SRR001']

    def test_creates_missing_output_dir_for_unavailable_report(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'

        sra_ids = pandas.Series(['SRR001'])
        output_dir = tmp_path / 'sanity'

        avail, unavail = check_quant_output(Args(), sra_ids, str(output_dir))

        assert avail == []
        assert unavail == ['SRR001']
        assert output_dir.exists()
        assert (output_dir / 'SRA_IDs_without_quant.txt').exists()

    def test_uses_parallel_workers_when_logs_suppressed(self, tmp_path, monkeypatch):
        class Args:
            out_dir = str(tmp_path)
            metadata = 'metadata.tsv'
            quiet = True
            verbose_runs = 0
            threads = 2

        quant_dir = tmp_path / 'quant'
        for sra_id in ['SRR001', 'SRR002']:
            sra_dir = quant_dir / sra_id
            sra_dir.mkdir(parents=True)
            (sra_dir / f'{sra_id}_abundance.tsv').write_text('data')
            (sra_dir / f'{sra_id}_run_info.json').write_text('{}')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        observed = {'max_workers': None}

        def fake_run_tasks(task_items, task_fn, max_workers=1):
            observed['max_workers'] = max_workers
            results = {}
            failures = []
            for task_item in task_items:
                try:
                    results[task_item] = task_fn(task_item)
                except Exception as exc:
                    failures.append((task_item, exc))
            return results, failures

        monkeypatch.setattr('amalgkit.sanity.run_tasks_with_optional_threads', fake_run_tasks)
        avail, unavail = check_quant_output(Args(), pandas.Series(['SRR001', 'SRR002']), str(output_dir))

        assert observed['max_workers'] == 2
        assert avail == ['SRR001', 'SRR002']
        assert unavail == []


# ---------------------------------------------------------------------------
# check_getfastq_outputs (filesystem checks)
# ---------------------------------------------------------------------------

class TestCheckGetfastqOutputs:
    def test_outputs_present(self, tmp_path, sample_metadata):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = None
            metadata = 'metadata.tsv'

        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_1.fastq.gz').write_text('data')
        (sra_dir / 'SRR001_2.fastq.gz').write_text('data')

        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_getfastq_outputs(Args(), sra_ids, sample_metadata, str(output_dir))
        assert avail == ['SRR001']
        assert unavail == []

    def test_missing_fastq_files(self, tmp_path, sample_metadata):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = None
            metadata = 'metadata.tsv'

        (tmp_path / 'getfastq' / 'SRR001').mkdir(parents=True)
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_getfastq_outputs(Args(), sra_ids, sample_metadata, str(output_dir))
        assert avail == []
        assert unavail == ['SRR001']

    def test_directory_named_fastq_file_is_ignored(self, tmp_path, sample_metadata):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = None
            metadata = 'metadata.tsv'

        sra_dir = tmp_path / 'getfastq' / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_1.fastq.gz').mkdir()
        (sra_dir / 'SRR001_2.fastq.gz').mkdir()
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_getfastq_outputs(Args(), sra_ids, sample_metadata, str(output_dir))
        assert avail == []
        assert unavail == ['SRR001']

    def test_sra_entry_that_is_not_directory_is_unavailable(self, tmp_path, sample_metadata):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = None
            metadata = 'metadata.tsv'

        getfastq_root = tmp_path / 'getfastq'
        getfastq_root.mkdir(parents=True)
        (getfastq_root / 'SRR001').write_text('not a directory')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_getfastq_outputs(Args(), sra_ids, sample_metadata, str(output_dir))
        assert avail == []
        assert unavail == ['SRR001']

    def test_getfastq_dir_file_path_marks_runs_unavailable(self, tmp_path, sample_metadata):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = str(tmp_path / 'getfastq_path')
            metadata = 'metadata.tsv'

        (tmp_path / 'getfastq_path').write_text('not a directory')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001'])
        avail, unavail = check_getfastq_outputs(Args(), sra_ids, sample_metadata, str(output_dir))
        assert avail == []
        assert unavail == ['SRR001']

    def test_suppresses_per_run_logs_when_quiet(self, tmp_path, sample_metadata, capsys):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = None
            metadata = 'metadata.tsv'
            quiet = True
            verbose_runs = 99

        getfastq_dir = tmp_path / 'getfastq'
        getfastq_dir.mkdir()
        (getfastq_dir / 'SRR001').mkdir()
        (getfastq_dir / 'SRR002').mkdir()
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001', 'SRR002'])

        check_getfastq_outputs(Args(), sra_ids, sample_metadata, str(output_dir))

        stdout = capsys.readouterr().out
        assert 'Per-run logs suppressed for 2 runs' in stdout
        assert 'Looking for SRR001' not in stdout

    def test_ignores_missing_run_ids(self, tmp_path, sample_metadata):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = None
            metadata = 'metadata.tsv'

        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        sra_ids = pandas.Series(['SRR001', np.nan, ''])

        avail, unavail = check_getfastq_outputs(Args(), sra_ids, sample_metadata, str(output_dir))

        assert avail == []
        assert unavail == ['SRR001']

    def test_handles_metadata_run_whitespace_without_crashing(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = None
            metadata = 'metadata.tsv'
            quiet = True
            verbose_runs = 0

        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['sp1'],
            'run': [' SRR001 '],
            'lib_layout': ['single'],
            'total_spots': [1],
            'spot_length': [100],
            'total_bases': [100],
            'exclusion': ['no'],
        }))
        getfastq_dir = tmp_path / 'getfastq' / 'SRR001'
        getfastq_dir.mkdir(parents=True)
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()

        avail, unavail = check_getfastq_outputs(Args(), metadata.df['run'], metadata, str(output_dir))

        assert avail == []
        assert unavail == ['SRR001']

    def test_metadata_run_whitespace_detects_existing_output(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = None
            metadata = 'metadata.tsv'
            quiet = True
            verbose_runs = 0

        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['sp1'],
            'run': [' SRR001 '],
            'lib_layout': ['single'],
            'total_spots': [1],
            'spot_length': [100],
            'total_bases': [100],
            'exclusion': ['no'],
        }))
        getfastq_dir = tmp_path / 'getfastq' / 'SRR001'
        getfastq_dir.mkdir(parents=True)
        (getfastq_dir / 'SRR001.fastq.gz').write_text('data')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()

        avail, unavail = check_getfastq_outputs(Args(), metadata.df['run'], metadata, str(output_dir))

        assert avail == ['SRR001']
        assert unavail == []

    def test_creates_missing_output_dir_for_getfastq_unavailable_report(self, tmp_path, sample_metadata):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = None
            metadata = 'metadata.tsv'

        sra_ids = pandas.Series(['SRR001'])
        output_dir = tmp_path / 'sanity'

        avail, unavail = check_getfastq_outputs(Args(), sra_ids, sample_metadata, str(output_dir))

        assert avail == []
        assert unavail == ['SRR001']
        assert output_dir.exists()
        assert (output_dir / 'SRA_IDs_without_fastq.txt').exists()

    def test_uses_parallel_workers_when_logs_suppressed(self, tmp_path, sample_metadata, monkeypatch):
        class Args:
            out_dir = str(tmp_path)
            getfastq_dir = None
            metadata = 'metadata.tsv'
            quiet = True
            verbose_runs = 0
            threads = 2

        for sra_id in ['SRR001', 'SRR002']:
            sra_dir = tmp_path / 'getfastq' / sra_id
            sra_dir.mkdir(parents=True)
            (sra_dir / f'{sra_id}_1.fastq.gz').write_text('data')
            (sra_dir / f'{sra_id}_2.fastq.gz').write_text('data')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        observed = {'max_workers': None}

        def fake_run_tasks(task_items, task_fn, max_workers=1):
            observed['max_workers'] = max_workers
            results = {}
            failures = []
            for task_item in task_items:
                try:
                    results[task_item] = task_fn(task_item)
                except Exception as exc:
                    failures.append((task_item, exc))
            return results, failures

        monkeypatch.setattr('amalgkit.sanity.run_tasks_with_optional_threads', fake_run_tasks)
        avail, unavail = check_getfastq_outputs(
            Args(),
            pandas.Series(['SRR001', 'SRR002']),
            sample_metadata,
            str(output_dir),
        )

        assert observed['max_workers'] == 2
        assert avail == ['SRR001', 'SRR002']
        assert unavail == []


# ---------------------------------------------------------------------------
# check_quant_index (issue #72: subspecies fallback)
# ---------------------------------------------------------------------------

class TestCheckQuantIndex:
    def test_index_found(self, tmp_path):
        """Index file is found for a species."""
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
        index_dir = tmp_path / 'index'
        index_dir.mkdir()
        (index_dir / 'Homo_sapiens.idx').write_text('')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        avail, unavail = check_quant_index(
            Args(), np.array(['Homo sapiens']), str(output_dir))
        assert 'Homo sapiens' in avail
        assert unavail == []

    def test_subspecies_fallback(self, tmp_path):
        """Issue #72: If Gorilla_gorilla_gorilla.idx not found, try Gorilla_gorilla.idx."""
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
        index_dir = tmp_path / 'index'
        index_dir.mkdir()
        (index_dir / 'Gorilla_gorilla.idx').write_text('')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        avail, unavail = check_quant_index(
            Args(), np.array(['Gorilla gorilla gorilla']), str(output_dir))
        assert 'Gorilla gorilla gorilla' in avail

    def test_index_not_found(self, tmp_path):
        """Missing index for a species."""
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
        index_dir = tmp_path / 'index'
        index_dir.mkdir()
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        avail, unavail = check_quant_index(
            Args(), np.array(['Danio rerio']), str(output_dir))
        assert 'Danio rerio' in unavail

    def test_index_directory_named_idx_is_ignored(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
        index_dir = tmp_path / 'index'
        index_dir.mkdir()
        (index_dir / 'Homo_sapiens.idx').mkdir()
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        avail, unavail = check_quant_index(
            Args(), np.array(['Homo sapiens']), str(output_dir))
        assert avail == []
        assert unavail == ['Homo sapiens']

    def test_marks_species_unavailable_when_multiple_index_files_match_prefix(self, tmp_path):
        """Ambiguous prefix matches should be treated as unavailable."""
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
        index_dir = tmp_path / 'index'
        index_dir.mkdir()
        (index_dir / 'Homo_sapiens_k31.idx').write_text('')
        (index_dir / 'Homo_sapiens_k63.idx').write_text('')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        avail, unavail = check_quant_index(
            Args(), np.array(['Homo sapiens']), str(output_dir))
        assert avail == []
        assert unavail == ['Homo sapiens']

    def test_does_not_match_similar_species_prefix(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
        index_dir = tmp_path / 'index'
        index_dir.mkdir()
        (index_dir / 'Homo_sapiens2.idx').write_text('')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        avail, unavail = check_quant_index(
            Args(), np.array(['Homo sapiens']), str(output_dir))
        assert avail == []
        assert unavail == ['Homo sapiens']

    def test_index_dir_file_path_does_not_crash(self, tmp_path):
        """If --index_dir points to a file, check should not crash."""
        class Args:
            out_dir = str(tmp_path)
            index_dir = str(tmp_path / 'index_path')
            metadata = 'metadata.tsv'
        (tmp_path / 'index_path').write_text('not a directory')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        avail, unavail = check_quant_index(
            Args(), np.array(['Homo sapiens']), str(output_dir))
        assert avail == []
        assert unavail == ['Homo sapiens']

    def test_missing_index_dir_marks_all_species_unavailable(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        avail, unavail = check_quant_index(
            Args(), np.array(['Homo sapiens', 'Mus musculus']), str(output_dir))
        assert avail == []
        assert unavail == ['Homo sapiens', 'Mus musculus']
        report = output_dir / 'species_without_index.txt'
        assert report.exists()

    def test_species_with_dot_matches_quant_style_index_name(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
        index_dir = tmp_path / 'index'
        index_dir.mkdir()
        (index_dir / 'Candidatus_sp._X.idx').write_text('')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()

        avail, unavail = check_quant_index(
            Args(), np.array(['Candidatus sp. X']), str(output_dir))

        assert avail == ['Candidatus sp. X']
        assert unavail == []

    def test_species_with_redundant_spaces_uses_normalized_fallback_prefix(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
        index_dir = tmp_path / 'index'
        index_dir.mkdir()
        (index_dir / 'Canis_lupus.idx').write_text('')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()

        avail, unavail = check_quant_index(
            Args(), np.array(['Canis   lupus familiaris']), str(output_dir))

        assert avail == ['Canis   lupus familiaris']
        assert unavail == []

    def test_suppresses_per_species_logs_when_verbose_runs_exceeded(self, tmp_path, capsys):
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
            quiet = False
            verbose_runs = 1
        index_dir = tmp_path / 'index'
        index_dir.mkdir()
        (index_dir / 'Species_alpha.idx').write_text('')
        (index_dir / 'Species_beta.idx').write_text('')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()

        check_quant_index(Args(), np.array(['Species alpha', 'Species beta']), str(output_dir))

        stdout = capsys.readouterr().out
        assert 'Per-species logs suppressed for 2 species' in stdout
        assert 'Looking for index file' not in stdout

    def test_uses_parallel_workers_when_logs_suppressed(self, tmp_path, monkeypatch):
        class Args:
            out_dir = str(tmp_path)
            index_dir = None
            metadata = 'metadata.tsv'
            quiet = True
            verbose_runs = 0
            threads = 2

        index_dir = tmp_path / 'index'
        index_dir.mkdir()
        (index_dir / 'Species_alpha.idx').write_text('')
        (index_dir / 'Species_beta.idx').write_text('')
        output_dir = tmp_path / 'sanity'
        output_dir.mkdir()
        observed = {'max_workers': None}

        def fake_run_tasks(task_items, task_fn, max_workers=1):
            observed['max_workers'] = max_workers
            results = {}
            failures = []
            for task_item in task_items:
                try:
                    results[task_item] = task_fn(task_item)
                except Exception as exc:
                    failures.append((task_item, exc))
            return results, failures

        monkeypatch.setattr('amalgkit.sanity.run_tasks_with_optional_threads', fake_run_tasks)

        avail, unavail = check_quant_index(
            Args(),
            np.array(['Species alpha', 'Species beta']),
            str(output_dir),
        )

        assert observed['max_workers'] == 2
        assert avail == ['Species alpha', 'Species beta']
        assert unavail == []


class TestSanityMain:
    def test_rejects_out_dir_file_path(self, tmp_path):
        out_path = tmp_path / 'out_path'
        out_path.write_text('not a directory')

        class Args:
            out_dir = str(out_path)
            metadata = 'metadata.tsv'
            all = False
            getfastq = False
            index = False
            quant = False

        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            sanity_main(Args())

    def test_rejects_sanity_output_file_path(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        (out_dir / 'sanity').write_text('not a directory')
        args = make_sanity_args(tmp_path, out_dir=str(out_dir))

        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Sp1'],
            'run': ['SRR001'],
            'exclusion': ['no'],
        }))
        monkeypatch.setattr('amalgkit.sanity.load_metadata', lambda _args: metadata)

        with pytest.raises(NotADirectoryError, match='Sanity output path exists but is not a directory'):
            sanity_main(args)

    def test_defaults_to_all_checks_and_writes_summary(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        args = make_sanity_args(tmp_path, out_dir=str(out_dir))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Sp1'],
            'run': ['SRR001'],
            'exclusion': ['no'],
        }))
        call_order = []

        monkeypatch.setattr('amalgkit.sanity.load_metadata', lambda _args: metadata)

        def fake_runner(check_name, target_type):
            def _runner(args, metadata, uni_species, sra_ids, output_dir, metadata_path):
                _ = (args, metadata, uni_species, sra_ids, metadata_path)
                call_order.append(check_name)
                return {
                    'check': check_name,
                    'target_type': target_type,
                    'checked_count': 1,
                    'available_count': 1,
                    'unavailable_count': 0,
                    'warning_count': 0,
                    'issue_count': 0,
                    'status': 'ok',
                    'scanned_path': str(tmp_path / check_name),
                    'missing_report_path': '',
                    'rerun_manifest_path': '',
                    'example_targets': '',
                }, []
            return _runner

        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_getfastq', fake_runner('getfastq', 'runs'))
        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_index', fake_runner('index', 'species'))
        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_quant', fake_runner('quant', 'runs'))
        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_merge', fake_runner('merge', 'species'))
        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_busco', fake_runner('busco', 'species'))
        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_finalize', fake_runner('finalize', 'species'))

        sanity_main(args)

        assert call_order == ['getfastq', 'index', 'quant', 'merge', 'busco', 'finalize']
        summary_df = pandas.read_csv(out_dir / 'sanity' / 'sanity_summary.tsv', sep='\t')
        assert summary_df['check'].tolist() == ['getfastq', 'index', 'quant', 'merge', 'busco', 'finalize']
        assert summary_df['status'].tolist() == ['ok', 'ok', 'ok', 'ok', 'ok', 'ok']
        assert (out_dir / 'sanity' / 'sanity_issues.tsv').exists()
        with open(out_dir / 'sanity' / 'sanity_report.json', 'r', encoding='utf-8') as handle:
            report_payload = json.load(handle)
        assert report_payload['requested_checks'] == ['getfastq', 'index', 'quant', 'merge', 'busco', 'finalize']

    def test_strict_mode_raises_after_writing_summary(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        args = make_sanity_args(tmp_path, out_dir=str(out_dir), quant=True, strict=True)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Sp1'],
            'run': ['SRR001'],
            'exclusion': ['no'],
        }))
        monkeypatch.setattr('amalgkit.sanity.load_metadata', lambda _args: metadata)

        with pytest.raises(
            AmalgkitExit,
            match=r'Sanity checks reported 1 error\(s\) and 0 warning\(s\) at strict_level=error\.',
        ):
            sanity_main(args)

        assert (out_dir / 'sanity' / 'SRA_IDs_without_quant.txt').exists()
        assert (out_dir / 'sanity' / 'rerun_quant_run_ids.txt').exists()
        summary_df = pandas.read_csv(out_dir / 'sanity' / 'sanity_summary.tsv', sep='\t')
        assert summary_df.loc[0, 'check'] == 'quant'
        assert int(summary_df.loc[0, 'unavailable_count']) == 1
        assert summary_df.loc[0, 'status'] == 'error'

    def test_check_and_filter_options_limit_selected_runners_and_write_comparison(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        output_dir = out_dir / 'sanity'
        output_dir.mkdir(parents=True)
        previous_payload = {
            'summary': [
                {
                    'check': 'quant',
                    'status': 'error',
                    'unavailable_count': 1,
                    'warning_count': 0,
                },
            ],
        }
        (output_dir / 'sanity_report.json').write_text(json.dumps(previous_payload), encoding='utf-8')
        args = make_sanity_args(
            tmp_path,
            out_dir=str(out_dir),
            check='quant,merge',
            run='SRR001',
            species='Sp1',
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Sp1', 'Sp2'],
            'run': ['SRR001', 'SRR002'],
            'exclusion': ['no', 'no'],
        }))
        observed = {}

        monkeypatch.setattr('amalgkit.sanity.load_metadata', lambda _args: metadata)

        def fake_quant(args, metadata, uni_species, sra_ids, output_dir, metadata_path):
            observed['quant'] = {
                'species': list(uni_species),
                'runs': list(sra_ids),
                'metadata_runs': metadata.df['run'].tolist(),
                'metadata_species': metadata.df['scientific_name'].tolist(),
            }
            return {
                'check': 'quant',
                'target_type': 'runs',
                'checked_count': 1,
                'available_count': 1,
                'unavailable_count': 0,
                'warning_count': 0,
                'issue_count': 0,
                'status': 'ok',
                'scanned_path': str(tmp_path / 'quant'),
                'missing_report_path': '',
                'rerun_manifest_path': '',
                'example_targets': '',
            }, []

        def fake_merge(args, metadata, uni_species, sra_ids, output_dir, metadata_path):
            observed['merge'] = {
                'species': list(uni_species),
                'runs': list(sra_ids),
            }
            return {
                'check': 'merge',
                'target_type': 'species',
                'checked_count': 1,
                'available_count': 1,
                'unavailable_count': 0,
                'warning_count': 0,
                'issue_count': 0,
                'status': 'ok',
                'scanned_path': str(tmp_path / 'merge'),
                'missing_report_path': '',
                'rerun_manifest_path': '',
                'example_targets': '',
            }, []

        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_quant', fake_quant)
        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_merge', fake_merge)
        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_getfastq', lambda *args, **kwargs: pytest.fail('getfastq should not run'))
        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_index', lambda *args, **kwargs: pytest.fail('index should not run'))
        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_busco', lambda *args, **kwargs: pytest.fail('busco should not run'))
        monkeypatch.setattr('amalgkit.sanity.run_sanity_check_finalize', lambda *args, **kwargs: pytest.fail('finalize should not run'))

        sanity_main(args)

        assert observed['quant']['species'] == ['Sp1']
        assert observed['quant']['runs'] == ['SRR001']
        assert observed['quant']['metadata_runs'] == ['SRR001']
        assert observed['quant']['metadata_species'] == ['Sp1']
        assert observed['merge']['species'] == ['Sp1']
        comparison_df = pandas.read_csv(output_dir / 'sanity_comparison.tsv', sep='\t')
        assert comparison_df['check'].tolist() == ['quant', 'merge']
        assert comparison_df.loc[0, 'delta_unavailable_count'] == -1
        with open(output_dir / 'sanity_report.json', 'r', encoding='utf-8') as handle:
            report_payload = json.load(handle)
        assert report_payload['run_filters'] == ['SRR001']
        assert report_payload['species_filters'] == ['Sp1']
        assert report_payload['requested_checks'] == ['quant', 'merge']

    def test_quant_runner_reports_invalid_content_orphans_and_rerun_manifest(self, tmp_path):
        out_dir = tmp_path / 'out'
        quant_dir = out_dir / 'quant'
        output_dir = out_dir / 'sanity'
        quant_dir.mkdir(parents=True)
        output_dir.mkdir(parents=True)
        metadata_path = tmp_path / 'metadata.tsv'
        metadata_path.write_text(
            'scientific_name\trun\texclusion\n'
            'Sp1\tSRR001\tno\n',
            encoding='utf-8',
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Sp1'],
            'run': ['SRR001'],
            'exclusion': ['no'],
        }))
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir()
        (sra_dir / 'SRR001_abundance.tsv').write_text('bad_column\nx\n', encoding='utf-8')
        (sra_dir / 'SRR001_run_info.json').write_text('{}', encoding='utf-8')
        orphan_dir = quant_dir / 'SRR999'
        orphan_dir.mkdir()
        (orphan_dir / 'SRR999_abundance.tsv').write_text('target_id\tlength\teff_length\test_counts\ttpm\nx\t1\t1\t1\t1\n', encoding='utf-8')
        (orphan_dir / 'SRR999_run_info.json').write_text('{"p_pseudoaligned": 10}', encoding='utf-8')
        args = make_sanity_args(tmp_path, out_dir=str(out_dir), quant=True)

        row, issues = run_sanity_check_quant(
            args=args,
            metadata=metadata,
            uni_species=np.array(['Sp1']),
            sra_ids=pandas.Series(['SRR001']),
            output_dir=str(output_dir),
            metadata_path=str(metadata_path),
        )

        assert row['status'] == 'error'
        assert row['unavailable_count'] == 1
        assert row['warning_count'] == 1
        assert os.path.exists(row['rerun_manifest_path'])
        with open(row['rerun_manifest_path'], 'r', encoding='utf-8') as handle:
            assert handle.read().strip() == 'SRR001'
        issue_types = {(issue['issue_type'], issue['target_id']) for issue in issues}
        assert ('invalid_content', 'SRR001') in issue_types
        assert ('orphan_output', 'SRR999') in issue_types
