import os
import pytest
import pandas
import numpy as np

from amalgkit.sanity import list_duplicates, parse_metadata, check_quant_output, check_quant_index, check_getfastq_outputs
from amalgkit.util import Metadata


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
        # unavail is empty because the code only populates it when quant_path exists
        assert unavail == []

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

    def test_prefix_match_with_suffix_filename(self, tmp_path):
        """Prefix search should detect index files with extra suffixes."""
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
        assert 'Homo sapiens' in avail
        assert unavail == []
