import os
import subprocess
import pytest
import pandas
from types import SimpleNamespace

from amalgkit.quant import quant_output_exists, check_layout_mismatch, get_index, call_kallisto
from amalgkit.util import Metadata


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


class TestQuantEdgeCases:
    def test_call_kallisto_uses_per_sra_nominal_length(self, tmp_path, monkeypatch):
        args = SimpleNamespace(threads=1)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'nominal_length': [150, 350],
            'exclusion': ['no', 'no'],
        }))
        sra_stat = {'sra_id': 'SRR002', 'layout': 'single'}
        captured = {}

        def fake_run(cmd, stdout, stderr):
            captured['cmd'] = cmd
            return SimpleNamespace(returncode=0, stdout=b'', stderr=b'')

        monkeypatch.setattr(subprocess, 'run', fake_run)
        call_kallisto(
            args=args,
            in_files=['dummy.fastq.gz'],
            metadata=metadata,
            sra_stat=sra_stat,
            output_dir=str(tmp_path),
            index='dummy.idx',
        )
        idx = captured['cmd'].index('-l')
        assert captured['cmd'][idx + 1] == '350'

    def test_get_index_raises_clear_error_for_multiple_fasta_candidates(self, tmp_path):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        (fasta_dir / 'Homo_sapiens_a.fa').write_text('>a\nAAAA\n')
        (fasta_dir / 'Homo_sapiens_b.fasta').write_text('>b\nCCCC\n')

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
        )
        with pytest.raises(ValueError, match='Found multiple reference fasta files'):
            get_index(args, 'Homo_sapiens')

    def test_call_kallisto_rejects_unsupported_layout(self, tmp_path):
        args = SimpleNamespace(threads=1)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'nominal_length': [200],
            'exclusion': ['no'],
        }))
        sra_stat = {'sra_id': 'SRR001', 'layout': 'unknown'}
        with pytest.raises(ValueError, match='Unsupported library layout'):
            call_kallisto(
                args=args,
                in_files=['dummy.fastq.gz'],
                metadata=metadata,
                sra_stat=sra_stat,
                output_dir=str(tmp_path),
                index='dummy.idx',
            )
