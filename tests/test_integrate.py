import gzip
from types import SimpleNamespace
import pytest

from amalgkit.integrate import count_fastq_lines, write_fastq_head, get_fastq_stats


def _write_fastq(path, reads):
    with open(path, 'wt') as fh:
        for i, seq in enumerate(reads):
            fh.write('@r{}\n'.format(i))
            fh.write(seq + '\n')
            fh.write('+\n')
            fh.write('I' * len(seq) + '\n')


class TestIntegrateFastqHelpers:
    def test_count_fastq_lines_plain_and_gz(self, tmp_path):
        plain = tmp_path / 'x.fastq'
        _write_fastq(str(plain), ['AAAA', 'CCCC', 'GGGG'])
        gz = tmp_path / 'x.fastq.gz'
        with open(str(plain), 'rb') as fin, gzip.open(str(gz), 'wb') as fout:
            fout.write(fin.read())

        assert count_fastq_lines(str(plain)) == 12
        assert count_fastq_lines(str(gz)) == 12

    def test_write_fastq_head_from_gz(self, tmp_path):
        plain = tmp_path / 'y.fastq'
        _write_fastq(str(plain), ['AAAA', 'CCCC', 'GGGG'])
        gz = tmp_path / 'y.fastq.gz'
        with open(str(plain), 'rb') as fin, gzip.open(str(gz), 'wb') as fout:
            fout.write(fin.read())

        out_head = tmp_path / 'head.fastq'
        write_fastq_head(path_fastq=str(gz), path_out=str(out_head), max_lines=5)
        with open(str(out_head), 'rt') as f:
            lines = f.readlines()
        assert len(lines) == 5
        assert lines[0].startswith('@r0')


class TestIntegrateGetFastqStats:
    @staticmethod
    def _write_one_read_fastq(path):
        with open(path, 'wt') as fh:
            fh.write('@r0\n')
            fh.write('AAAA\n')
            fh.write('+\n')
            fh.write('IIII\n')

    def test_keeps_dotted_sample_ids_distinct(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        self._write_one_read_fastq(str(fastq_dir / 'sample.v1_1.fastq'))
        self._write_one_read_fastq(str(fastq_dir / 'sample.v2_1.fastq'))
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_dir))

        df = get_fastq_stats(args)

        assert set(df['run'].tolist()) == {'sample.v1', 'sample.v2'}
        assert set(df['lib_layout'].tolist()) == {'single'}
        assert all(df['read2_path'] == 'unavailable')

    def test_paired_mates_are_assigned_to_read1_read2(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        # Intentionally create _2 first to ensure assignment is by suffix, not directory listing order.
        self._write_one_read_fastq(str(fastq_dir / 'alpha_2.fastq'))
        self._write_one_read_fastq(str(fastq_dir / 'alpha_1.fastq'))
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_dir))

        df = get_fastq_stats(args)

        assert df.shape[0] == 1
        assert df.loc[0, 'run'] == 'alpha'
        assert df.loc[0, 'lib_layout'] == 'paired'
        assert df.loc[0, 'read1_path'].endswith('alpha_1.fastq')
        assert df.loc[0, 'read2_path'].endswith('alpha_2.fastq')

    def test_raises_on_paired_read_count_mismatch(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        # read1: 1 read
        self._write_one_read_fastq(str(fastq_dir / 'beta_1.fastq'))
        # read2: 2 reads
        with open(str(fastq_dir / 'beta_2.fastq'), 'wt') as fh:
            fh.write('@r0\nAAAA\n+\nIIII\n')
            fh.write('@r1\nCCCC\n+\nIIII\n')
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_dir))

        with pytest.raises(ValueError, match='Mismatched paired-end read counts'):
            get_fastq_stats(args)

    def test_raises_on_malformed_gz_fastq_line_count_in_quick_mode(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        malformed = fastq_dir / 'gamma_1.fastq.gz'
        with gzip.open(str(malformed), 'wt') as fh:
            fh.write('@r0\nAAAA\n+\nIIII\n')
            fh.write('orphan_line\n')  # total lines = 5, not divisible by 4
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=False, out_dir=str(out_dir))

        with pytest.raises(ValueError, match='Malformed FASTQ'):
            get_fastq_stats(args)

    def test_parallel_scan_with_threads(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        self._write_one_read_fastq(str(fastq_dir / 'runA.fastq'))
        self._write_one_read_fastq(str(fastq_dir / 'runB.fastq'))
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_dir), threads=2)

        df = get_fastq_stats(args)

        assert set(df['run'].tolist()) == {'runA', 'runB'}
        assert set(df['lib_layout'].tolist()) == {'single'}
