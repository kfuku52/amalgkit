import gzip
from types import SimpleNamespace
import pytest
import pandas
import numpy

from amalgkit.integrate import (
    count_fastq_lines,
    write_fastq_head,
    get_fastq_stats,
    integrate_main,
    get_gzip_isize,
    estimate_total_spots_from_gzip_sample,
)
from amalgkit.util import Metadata


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

    def test_get_gzip_isize_raises_clear_error_for_too_small_file(self, tmp_path):
        bad = tmp_path / 'broken.fastq.gz'
        bad.write_bytes(b'\x1f\x8b')
        with pytest.raises(ValueError, match='too small for ISIZE footer'):
            get_gzip_isize(str(bad))

    def test_estimate_total_spots_from_gzip_sample_rejects_large_gzip_for_quick_mode(self, tmp_path, monkeypatch):
        gz_path = tmp_path / 'x.fastq.gz'
        gz_path.write_bytes(b'\x1f\x8b\x08\x00')
        monkeypatch.setattr('amalgkit.integrate.os.path.getsize', lambda _p: (1 << 32))
        monkeypatch.setattr('amalgkit.integrate.get_gzip_isize', lambda _p: 1024)

        with pytest.raises(ValueError, match='ISIZE footer may overflow'):
            estimate_total_spots_from_gzip_sample(
                path_fastq=str(gz_path),
                sampled_reads=100,
                sampled_record_chars=1000,
            )


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
        expected_size = (fastq_dir / 'alpha_1.fastq').stat().st_size + (fastq_dir / 'alpha_2.fastq').stat().st_size
        assert int(df.loc[0, 'size']) == int(expected_size)

    def test_paired_total_bases_uses_both_read_lengths(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        with open(str(fastq_dir / 'lenmix_1.fastq'), 'wt') as fh:
            fh.write('@r0\nAAAA\n+\nIIII\n')
        with open(str(fastq_dir / 'lenmix_2.fastq'), 'wt') as fh:
            fh.write('@r0\nAA\n+\nII\n')
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_dir))

        df = get_fastq_stats(args)

        assert df.shape[0] == 1
        assert int(df.loc[0, 'total_spots']) == 1
        assert int(df.loc[0, 'spot_length']) == 4
        assert int(df.loc[0, 'total_bases']) == 6

    def test_raises_when_only_read2_file_is_present(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        self._write_one_read_fastq(str(fastq_dir / 'alpha_2.fastq'))
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_dir))

        with pytest.raises(ValueError, match='Only one paired-end mate file'):
            get_fastq_stats(args)

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

    def test_quick_mode_does_not_use_full_line_counter(self, tmp_path, monkeypatch):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()

        plain = fastq_dir / 'delta_1.fastq'
        _write_fastq(str(plain), ['AAAA', 'CCCC', 'GGGG', 'TTTT'])
        gz = fastq_dir / 'delta_1.fastq.gz'
        with open(str(plain), 'rb') as fin, gzip.open(str(gz), 'wb') as fout:
            fout.write(fin.read())
        plain.unlink()

        monkeypatch.setattr(
            'amalgkit.integrate.estimate_total_spots_from_line_count',
            lambda _path: (_ for _ in ()).throw(AssertionError('full line counter should not be used in quick mode')),
        )
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=False, out_dir=str(out_dir))

        df = get_fastq_stats(args)

        assert df.shape[0] == 1
        assert df.loc[0, 'run'] == 'delta'
        assert int(df.loc[0, 'total_spots']) > 0

    def test_quick_mode_counts_reads_exactly_when_gzip_is_fully_sampled(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()

        part1 = fastq_dir / 'epsilon.part1.fastq.gz'
        part2 = fastq_dir / 'epsilon.part2.fastq.gz'
        with gzip.open(str(part1), 'wt') as fh:
            fh.write('@r0\nAAAA\n+\nIIII\n')
            fh.write('@r1\nCCCC\n+\nIIII\n')
            fh.write('@r2\nGGGG\n+\nIIII\n')
        with gzip.open(str(part2), 'wt') as fh:
            fh.write('@r3\nTTTT\n+\nIIII\n')
            fh.write('@r4\nNNNN\n+\nIIII\n')

        merged = fastq_dir / 'epsilon.fastq.gz'
        with open(str(merged), 'wb') as fout:
            fout.write(part1.read_bytes())
            fout.write(part2.read_bytes())
        part1.unlink()
        part2.unlink()

        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=False, out_dir=str(out_dir))
        df = get_fastq_stats(args)

        assert df.shape[0] == 1
        assert df.loc[0, 'run'] == 'epsilon'
        assert int(df.loc[0, 'total_spots']) == 5

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

    def test_ignores_fastq_files_without_run_id(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        self._write_one_read_fastq(str(fastq_dir / 'alpha.fastq'))
        self._write_one_read_fastq(str(fastq_dir / '.fastq'))
        self._write_one_read_fastq(str(fastq_dir / '_1.fastq'))
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_dir))

        df = get_fastq_stats(args)

        assert df['run'].tolist() == ['alpha']

    def test_detects_uppercase_fastq_extension(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        self._write_one_read_fastq(str(fastq_dir / 'alpha.FASTQ'))
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_dir))

        df = get_fastq_stats(args)

        assert df['run'].tolist() == ['alpha']

    def test_raises_clear_error_when_fastq_dir_is_missing(self, tmp_path):
        args = SimpleNamespace(fastq_dir=None, accurate_size=True, out_dir=str(tmp_path))
        with pytest.raises(ValueError, match='--fastq_dir is required'):
            get_fastq_stats(args)

    def test_raises_clear_error_when_fastq_dir_is_file(self, tmp_path):
        fastq_file = tmp_path / 'fq_path'
        fastq_file.write_text('not a directory')
        args = SimpleNamespace(fastq_dir=str(fastq_file), accurate_size=True, out_dir=str(tmp_path))
        with pytest.raises(NotADirectoryError, match='not a directory'):
            get_fastq_stats(args)

    def test_raises_when_out_dir_is_file(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_path = tmp_path / 'out_path'
        fastq_dir.mkdir()
        out_path.write_text('not a directory')
        self._write_one_read_fastq(str(fastq_dir / 'alpha.fastq'))
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_path))

        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            get_fastq_stats(args)

    def test_raises_when_metadata_output_path_is_file(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        self._write_one_read_fastq(str(fastq_dir / 'alpha.fastq'))
        (out_dir / 'metadata').write_text('not a directory')
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_dir))

        with pytest.raises(NotADirectoryError, match='Metadata output path exists but is not a directory'):
            get_fastq_stats(args)

    def test_raises_when_single_fastq_has_no_reads_in_accurate_mode(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        (fastq_dir / 'empty.fastq').write_text('')
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=True, out_dir=str(out_dir))

        with pytest.raises(ValueError, match='No reads detected in FASTQ file'):
            get_fastq_stats(args)

    def test_raises_when_single_fastq_has_no_reads_in_quick_mode(self, tmp_path):
        fastq_dir = tmp_path / 'fq'
        out_dir = tmp_path / 'out'
        fastq_dir.mkdir()
        out_dir.mkdir()
        with gzip.open(str(fastq_dir / 'empty.fastq.gz'), 'wt') as fh:
            fh.write('')
        args = SimpleNamespace(fastq_dir=str(fastq_dir), accurate_size=False, out_dir=str(out_dir))

        with pytest.raises(ValueError, match='No reads detected in FASTQ file'):
            get_fastq_stats(args)


class TestIntegrateMain:
    def test_marks_data_available_with_whitespace_run_ids(self, tmp_path, monkeypatch):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        (metadata_dir / 'metadata.tsv').write_text('dummy\n')
        args = SimpleNamespace(
            metadata='inferred',
            out_dir=str(tmp_path),
            fastq_dir=str(tmp_path / 'fq'),
            accurate_size=True,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [' SRR001 '],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))

        monkeypatch.setattr('amalgkit.integrate.load_metadata', lambda _args: metadata)

        def fake_check_getfastq_outputs(_args, sra_ids, _metadata, _output_dir):
            assert list(sra_ids) == ['SRR001']
            return ['SRR001'], []

        monkeypatch.setattr('amalgkit.integrate.check_getfastq_outputs', fake_check_getfastq_outputs)
        monkeypatch.setattr(
            'amalgkit.integrate.get_fastq_stats',
            lambda _args: pandas.DataFrame(columns=['run', 'scientific_name', 'exclusion']),
        )

        integrate_main(args)

        out_path = metadata_dir / 'metadata_updated_for_private_fastq.tsv'
        out_df = pandas.read_csv(str(out_path), sep='\t')
        assert out_df.loc[0, 'run'] == 'SRR001'
        assert out_df.loc[0, 'data_available'] == 'yes'

    def test_rejects_missing_run_ids_in_existing_metadata(self, tmp_path, monkeypatch):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        (metadata_dir / 'metadata.tsv').write_text('dummy\n')
        args = SimpleNamespace(
            metadata='inferred',
            out_dir=str(tmp_path),
            fastq_dir=str(tmp_path / 'fq'),
            accurate_size=True,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [numpy.nan],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))

        monkeypatch.setattr('amalgkit.integrate.load_metadata', lambda _args: metadata)
        monkeypatch.setattr(
            'amalgkit.integrate.check_getfastq_outputs',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('check_getfastq_outputs should not be called')),
        )

        with pytest.raises(ValueError, match='Missing Run ID\\(s\\) were detected'):
            integrate_main(args)

    def test_rejects_existing_metadata_without_run_column(self, tmp_path, monkeypatch):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        (metadata_dir / 'metadata.tsv').write_text('dummy\n')
        args = SimpleNamespace(
            metadata='inferred',
            out_dir=str(tmp_path),
            fastq_dir=str(tmp_path / 'fq'),
            accurate_size=True,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))
        metadata.df = metadata.df.drop(columns=['run'])

        monkeypatch.setattr('amalgkit.integrate.load_metadata', lambda _args: metadata)
        monkeypatch.setattr(
            'amalgkit.integrate.check_getfastq_outputs',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('check_getfastq_outputs should not be called')),
        )

        with pytest.raises(ValueError, match='Column \"run\" is required in metadata'):
            integrate_main(args)

    def test_rejects_duplicate_run_ids_in_existing_metadata(self, tmp_path, monkeypatch):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        (metadata_dir / 'metadata.tsv').write_text('dummy\n')
        args = SimpleNamespace(
            metadata='inferred',
            out_dir=str(tmp_path),
            fastq_dir=str(tmp_path / 'fq'),
            accurate_size=True,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', ' SRR001 '],
            'scientific_name': ['Sp1', 'Sp1'],
            'exclusion': ['no', 'no'],
        }))

        monkeypatch.setattr('amalgkit.integrate.load_metadata', lambda _args: metadata)
        monkeypatch.setattr(
            'amalgkit.integrate.check_getfastq_outputs',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('check_getfastq_outputs should not be called')),
        )

        with pytest.raises(ValueError, match='Duplicate Run ID\\(s\\) were detected in metadata: SRR001'):
            integrate_main(args)

    def test_rejects_duplicate_run_ids_after_merging_existing_and_private_metadata(self, tmp_path, monkeypatch):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        (metadata_dir / 'metadata.tsv').write_text('dummy\n')
        args = SimpleNamespace(
            metadata='inferred',
            out_dir=str(tmp_path),
            fastq_dir=str(tmp_path / 'fq'),
            accurate_size=True,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Sp1'],
            'exclusion': ['no'],
        }))

        monkeypatch.setattr('amalgkit.integrate.load_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.integrate.check_getfastq_outputs', lambda *_args, **_kwargs: (['SRR001'], []))
        monkeypatch.setattr(
            'amalgkit.integrate.get_fastq_stats',
            lambda _args: pandas.DataFrame({
                'run': [' SRR001 '],
                'scientific_name': ['Sp1'],
                'exclusion': ['no'],
            }),
        )

        with pytest.raises(ValueError, match='Duplicate Run ID\\(s\\) were detected after merging existing metadata and private fastq info: SRR001'):
            integrate_main(args)

    def test_rejects_metadata_path_that_is_directory(self, tmp_path):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        (metadata_dir / 'metadata.tsv').mkdir()
        args = SimpleNamespace(
            metadata='inferred',
            out_dir=str(tmp_path),
            fastq_dir=str(tmp_path / 'fq'),
            accurate_size=True,
        )

        with pytest.raises(IsADirectoryError, match='Metadata path exists but is not a file'):
            integrate_main(args)

    def test_rejects_explicit_missing_metadata_path(self, tmp_path):
        args = SimpleNamespace(
            metadata=str(tmp_path / 'missing.tsv'),
            out_dir=str(tmp_path),
            fastq_dir=str(tmp_path / 'fq'),
            accurate_size=True,
        )

        with pytest.raises(FileNotFoundError, match='Metadata file not found'):
            integrate_main(args)
