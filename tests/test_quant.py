import os
import pytest

from amalgkit.quant import quant_output_exists, check_layout_mismatch


# ---------------------------------------------------------------------------
# quant_output_exists (checks for abundance.tsv in output directory)
# ---------------------------------------------------------------------------

class TestQuantOutputExists:
    def test_output_exists(self, tmp_path):
        """Returns True when abundance.tsv exists."""
        (tmp_path / 'SRR001_abundance.tsv').write_text('target_id\tdata\n')
        assert quant_output_exists('SRR001', str(tmp_path)) is True

    def test_output_missing(self, tmp_path):
        """Returns False when abundance.tsv is missing."""
        assert quant_output_exists('SRR001', str(tmp_path)) is False


# ---------------------------------------------------------------------------
# check_layout_mismatch (issue #80: corrects layout when files disagree)
# ---------------------------------------------------------------------------

class TestCheckLayoutMismatch:
    def test_paired_metadata_single_file(self, tmp_path):
        """Paired layout corrected to single when only one fastq found."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('data')
        result = check_layout_mismatch(sra_stat, str(tmp_path))
        assert result['layout'] == 'single'

    def test_paired_metadata_paired_files(self, tmp_path):
        """Paired layout unchanged when two fastq files found."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.amalgkit.fastq.gz').write_text('data')
        result = check_layout_mismatch(sra_stat, str(tmp_path))
        assert result['layout'] == 'paired'

    def test_single_metadata_no_change(self, tmp_path):
        """Single layout is not affected by mismatch check."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        result = check_layout_mismatch(sra_stat, str(tmp_path))
        assert result['layout'] == 'single'
