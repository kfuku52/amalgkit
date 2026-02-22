import os
import subprocess
import pytest
import pandas
from types import SimpleNamespace

from amalgkit.quant import (
    quant_output_exists,
    purge_existing_quant_outputs,
    check_layout_mismatch,
    resolve_input_fastq_files,
    find_species_index_files,
    find_species_fasta_files,
    check_kallisto_dependency,
    get_index,
    call_kallisto,
    run_quant,
    pre_resolve_species_indices,
    prefetch_getfastq_run_files,
    build_quant_tasks,
    quant_main,
)
from amalgkit.util import Metadata


# ---------------------------------------------------------------------------
# quant_output_exists (checks for abundance.tsv in output directory)
# ---------------------------------------------------------------------------

class TestQuantOutputExists:
    def test_output_exists(self, tmp_path):
        """Returns True when both abundance.tsv and run_info.json exist."""
        (tmp_path / 'SRR001_abundance.tsv').write_text('target_id\tdata\n')
        (tmp_path / 'SRR001_run_info.json').write_text('{}')
        assert quant_output_exists('SRR001', str(tmp_path)) is True

    def test_output_missing(self, tmp_path):
        """Returns False when required output files are missing."""
        assert quant_output_exists('SRR001', str(tmp_path)) is False

    def test_output_missing_run_info(self, tmp_path):
        (tmp_path / 'SRR001_abundance.tsv').write_text('target_id\tdata\n')
        assert quant_output_exists('SRR001', str(tmp_path)) is False

    def test_ignores_directory_named_as_abundance_file(self, tmp_path):
        (tmp_path / 'SRR001_abundance.tsv').mkdir()
        assert quant_output_exists('SRR001', str(tmp_path)) is False


class TestPurgeExistingQuantOutputs:
    def test_removes_known_stale_outputs(self, tmp_path):
        for name in [
            'SRR001_abundance.tsv',
            'SRR001_run_info.json',
            'SRR001_abundance.h5',
            'abundance.tsv',
            'run_info.json',
            'abundance.h5',
        ]:
            (tmp_path / name).write_text('x')

        purge_existing_quant_outputs('SRR001', str(tmp_path))

        for name in [
            'SRR001_abundance.tsv',
            'SRR001_run_info.json',
            'SRR001_abundance.h5',
            'abundance.tsv',
            'run_info.json',
            'abundance.h5',
        ]:
            assert not (tmp_path / name).exists()

    def test_raises_when_stale_output_is_directory(self, tmp_path):
        (tmp_path / 'SRR001_abundance.tsv').mkdir()
        with pytest.raises(IsADirectoryError, match='not a file'):
            purge_existing_quant_outputs('SRR001', str(tmp_path))


class TestCheckKallistoDependency:
    def test_raises_clear_error_when_kallisto_missing(self, monkeypatch):
        def fake_run(_cmd, stdout=None, stderr=None):
            raise FileNotFoundError('kallisto')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        with pytest.raises(FileNotFoundError, match='kallisto executable not found'):
            check_kallisto_dependency()

    def test_raises_when_kallisto_probe_returns_nonzero(self, monkeypatch):
        def fake_run(cmd, stdout=None, stderr=None):
            return subprocess.CompletedProcess(cmd, 127, stdout=b'', stderr=b'')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        with pytest.raises(RuntimeError, match='kallisto dependency probe failed'):
            check_kallisto_dependency()


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

    def test_uses_prefetched_files_without_rescan(self, tmp_path, monkeypatch):
        """When files are given, directory listing helper should not be called."""
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        files = {'SRR001.fastq.gz'}

        def fail_if_called(_output_dir):
            raise AssertionError('list_getfastq_run_files should not be called when files are provided.')

        monkeypatch.setattr('amalgkit.quant.list_getfastq_run_files', fail_if_called)
        result = check_layout_mismatch(sra_stat, str(tmp_path), files=files)
        assert result['layout'] == 'single'

    def test_ignores_similar_prefix_run_ids_when_counting_fastq_files(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        files = {'SRR0010.fastq.gz'}
        result = check_layout_mismatch(sra_stat, str(tmp_path), files=files)
        assert result['layout'] == 'paired'


class TestResolveInputFastqFiles:
    def test_paired_returns_deterministic_order(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'paired'}
        files = {'SRR001_2.fastq.gz', 'SRR001_1.fastq.gz'}
        result = resolve_input_fastq_files(
            sra_stat=sra_stat,
            output_dir_getfastq=str(tmp_path),
            ext='.fastq.gz',
            files=files,
        )
        assert result == [
            str(tmp_path / 'SRR001_1.fastq.gz'),
            str(tmp_path / 'SRR001_2.fastq.gz'),
        ]

    def test_single_returns_single_file(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        files = {'SRR001.fastq.gz'}
        result = resolve_input_fastq_files(
            sra_stat=sra_stat,
            output_dir_getfastq=str(tmp_path),
            ext='.fastq.gz',
            files=files,
        )
        assert result == [str(tmp_path / 'SRR001.fastq.gz')]

    def test_fallback_ignores_similar_prefix_run_ids(self, tmp_path):
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        files = {'SRR0010.fastq.gz'}
        result = resolve_input_fastq_files(
            sra_stat=sra_stat,
            output_dir_getfastq=str(tmp_path),
            ext='.fastq.gz',
            files=files,
        )
        assert result == []


class TestIndexDiscoveryHelpers:
    def test_find_species_index_files_filters_prefix(self, tmp_path):
        (tmp_path / 'Homo_sapiens.idx').write_text('idx')
        (tmp_path / 'Homo_sapiens.v2.idx').write_text('idx')
        (tmp_path / 'Mus_musculus.idx').write_text('idx')
        result = find_species_index_files(str(tmp_path), 'Homo_sapiens')
        assert result == [
            str(tmp_path / 'Homo_sapiens.idx'),
            str(tmp_path / 'Homo_sapiens.v2.idx'),
        ]

    def test_find_species_fasta_files_filters_suffix(self, tmp_path):
        (tmp_path / 'Homo_sapiens.fa').write_text('>a\nAAAA\n')
        (tmp_path / 'Homo_sapiens.fasta.gz').write_text('gz')
        (tmp_path / 'Homo_sapiens.txt').write_text('no')
        (tmp_path / 'Mus_musculus.fa').write_text('>a\nAAAA\n')
        result = find_species_fasta_files(str(tmp_path), 'Homo_sapiens')
        assert result == [
            str(tmp_path / 'Homo_sapiens.fa'),
            str(tmp_path / 'Homo_sapiens.fasta.gz'),
        ]

    def test_find_species_fasta_files_accepts_uppercase_suffix(self, tmp_path):
        (tmp_path / 'Homo_sapiens.FA').write_text('>a\nAAAA\n')
        (tmp_path / 'Homo_sapiens.FASTA.GZ').write_text('gz')
        result = find_species_fasta_files(str(tmp_path), 'Homo_sapiens')
        assert result == [
            str(tmp_path / 'Homo_sapiens.FA'),
            str(tmp_path / 'Homo_sapiens.FASTA.GZ'),
        ]

    def test_find_species_index_files_rejects_similar_species_prefix(self, tmp_path):
        (tmp_path / 'Homo_sapiens.idx').write_text('idx')
        (tmp_path / 'Homo_sapiens2.idx').write_text('idx')
        result = find_species_index_files(str(tmp_path), 'Homo_sapiens')
        assert result == [str(tmp_path / 'Homo_sapiens.idx')]

    def test_find_species_index_files_ignores_directories(self, tmp_path):
        (tmp_path / 'Homo_sapiens.idx').mkdir()
        (tmp_path / 'Homo_sapiens.v2.idx').write_text('idx')
        result = find_species_index_files(str(tmp_path), 'Homo_sapiens')
        assert result == [str(tmp_path / 'Homo_sapiens.v2.idx')]

    def test_find_species_fasta_files_raises_for_file_path(self, tmp_path):
        fasta_path = tmp_path / 'fasta_path'
        fasta_path.write_text('not a directory')
        with pytest.raises(NotADirectoryError, match='Fasta path exists but is not a directory'):
            find_species_fasta_files(str(fasta_path), 'Homo_sapiens')

    def test_find_species_fasta_files_ignores_directories(self, tmp_path):
        (tmp_path / 'Homo_sapiens.fa').mkdir()
        (tmp_path / 'Homo_sapiens.fasta').write_text('>a\nAAAA\n')
        result = find_species_fasta_files(str(tmp_path), 'Homo_sapiens')
        assert result == [str(tmp_path / 'Homo_sapiens.fasta')]


class TestGetfastqPrefetch:
    def test_prefetch_getfastq_run_files_scans_only_targets(self, tmp_path):
        out_dir = tmp_path / 'out'
        getfastq_root = out_dir / 'getfastq'
        (getfastq_root / 'SRR001').mkdir(parents=True)
        (getfastq_root / 'SRR002').mkdir(parents=True)
        (getfastq_root / 'SRR003').mkdir(parents=True)
        (getfastq_root / 'SRR001' / 'SRR001.fastq.gz').write_text('x')
        (getfastq_root / 'SRR002' / 'SRR002.fastq.gz').write_text('y')
        (getfastq_root / 'SRR003' / 'SRR003.fastq.gz').write_text('z')
        tasks = [('SRR001', 'Species A'), ('SRR002', 'Species B')]
        args = SimpleNamespace(out_dir=str(out_dir))

        out = prefetch_getfastq_run_files(args, tasks)

        assert set(out.keys()) == {'SRR001', 'SRR002'}
        assert out['SRR001'] == {'SRR001.fastq.gz'}
        assert out['SRR002'] == {'SRR002.fastq.gz'}

    def test_prefetch_getfastq_run_files_avoids_root_directory_scan(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        getfastq_root = out_dir / 'getfastq'
        (getfastq_root / 'SRR001').mkdir(parents=True)
        (getfastq_root / 'SRR001' / 'SRR001.fastq.gz').write_text('x')
        args = SimpleNamespace(out_dir=str(out_dir))
        tasks = [('SRR001', 'Species A')]

        real_scandir = os.scandir
        root_realpath = os.path.realpath(str(getfastq_root))

        def fail_on_root_scan(path):
            if os.path.realpath(path) == root_realpath:
                raise AssertionError('Root getfastq directory should not be scanned.')
            return real_scandir(path)

        monkeypatch.setattr('amalgkit.quant.os.scandir', fail_on_root_scan)

        out = prefetch_getfastq_run_files(args, tasks)

        assert out == {'SRR001': {'SRR001.fastq.gz'}}

    def test_prefetch_getfastq_run_files_ignores_non_file_entries(self, tmp_path):
        out_dir = tmp_path / 'out'
        getfastq_root = out_dir / 'getfastq'
        run_dir = getfastq_root / 'SRR001'
        run_dir.mkdir(parents=True)
        (run_dir / 'SRR001.fastq.gz').write_text('x')
        (run_dir / 'tmp_dir.fastq.gz').mkdir()
        args = SimpleNamespace(out_dir=str(out_dir))
        tasks = [('SRR001', 'Species A')]

        out = prefetch_getfastq_run_files(args, tasks)

        assert out == {'SRR001': {'SRR001.fastq.gz'}}

    def test_run_quant_uses_prefetched_getfastq_files(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
            _prefetched_getfastq_run_files={'SRR001': {'SRR001.fastq.gz'}},
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'spot_length': [100],
            'total_bases': [1000],
            'nominal_length': [200],
            'exclusion': ['no'],
        }))
        called = {'kallisto': 0}

        def fail_if_listdir_used(_path):
            raise AssertionError('list_getfastq_run_files should not be called when prefetched files are available.')

        def fake_call_kallisto(_args, in_files, _metadata, _sra_stat, _output_dir, _index):
            called['kallisto'] += 1
            assert len(in_files) == 1
            assert in_files[0].endswith('SRR001.fastq.gz')
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr('amalgkit.quant.list_getfastq_run_files', fail_if_listdir_used)
        monkeypatch.setattr('amalgkit.quant.call_kallisto', fake_call_kallisto)

        run_quant(args, metadata, 'SRR001', 'dummy.idx')

        assert called['kallisto'] == 1

    def test_run_quant_redo_purges_stale_outputs_before_call(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        quant_run_dir = out_dir / 'quant' / 'SRR001'
        quant_run_dir.mkdir(parents=True)
        (quant_run_dir / 'SRR001_abundance.tsv').write_text('stale')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=True,
            clean_fastq=False,
            threads=1,
            _prefetched_getfastq_run_files={'SRR001': {'SRR001.fastq.gz'}},
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'spot_length': [100],
            'total_bases': [1000],
            'nominal_length': [200],
            'exclusion': ['no'],
        }))
        called = {'kallisto': 0}

        def fake_call_kallisto(_args, in_files, _metadata, _sra_stat, _output_dir, _index):
            called['kallisto'] += 1
            assert not (quant_run_dir / 'SRR001_abundance.tsv').exists()
            assert len(in_files) == 1
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr('amalgkit.quant.call_kallisto', fake_call_kallisto)

        run_quant(args, metadata, 'SRR001', 'dummy.idx')
        assert called['kallisto'] == 1

    def test_run_quant_raises_when_fastq_is_missing(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
            _prefetched_getfastq_run_files={'SRR001': set()},
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'spot_length': [100],
            'total_bases': [1000],
            'nominal_length': [200],
            'exclusion': ['no'],
        }))

        monkeypatch.setattr(
            'amalgkit.quant.get_newest_intermediate_file_extension',
            lambda *_args, **_kwargs: '.fastq.gz',
        )
        monkeypatch.setattr('amalgkit.quant.list_getfastq_run_files', lambda _path: set())
        monkeypatch.setattr(
            'amalgkit.quant.call_kallisto',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('call_kallisto should not be called')),
        )

        with pytest.raises(FileNotFoundError, match=r'SRR001: Fastq file not found\. Check .*getfastq.*SRR001'):
            run_quant(args, metadata, 'SRR001', 'dummy.idx')

    def test_run_quant_rejects_getfastq_run_path_that_is_file(self, tmp_path):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        getfastq_root = out_dir / 'getfastq'
        getfastq_root.mkdir()
        (getfastq_root / 'SRR001').write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'spot_length': [100],
            'total_bases': [1000],
            'nominal_length': [200],
            'exclusion': ['no'],
        }))

        with pytest.raises(NotADirectoryError, match='getfastq run path exists but is not a directory'):
            run_quant(args, metadata, 'SRR001', 'dummy.idx')

    def test_run_quant_rejects_quant_run_output_path_that_is_file(self, tmp_path):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        quant_root = out_dir / 'quant'
        quant_root.mkdir()
        (quant_root / 'SRR001').write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'spot_length': [100],
            'total_bases': [1000],
            'nominal_length': [200],
            'exclusion': ['no'],
        }))

        with pytest.raises(NotADirectoryError, match='Quant run output path exists but is not a directory'):
            run_quant(args, metadata, 'SRR001', 'dummy.idx')


class TestQuantEdgeCases:
    def test_build_quant_tasks_rejects_missing_required_columns(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'exclusion': ['no'],
        }))
        metadata.df = metadata.df.drop(columns=['scientific_name'])
        with pytest.raises(ValueError, match='Missing required metadata column\\(s\\) for quant: scientific_name'):
            build_quant_tasks(metadata)

    def test_build_quant_tasks_rejects_missing_scientific_name(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': [float('nan')],
            'exclusion': ['no'],
            'nominal_length': [200],
        }))
        with pytest.raises(ValueError, match='Missing scientific_name in metadata for run\\(s\\): SRR001'):
            build_quant_tasks(metadata)

    def test_build_quant_tasks_rejects_missing_run_id(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [float('nan')],
            'scientific_name': ['Species A'],
            'exclusion': ['no'],
            'nominal_length': [200],
        }))
        with pytest.raises(ValueError, match='Missing run ID in metadata'):
            build_quant_tasks(metadata)

    def test_build_quant_tasks_strips_run_and_species_whitespace(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [' SRR001 '],
            'scientific_name': [' Species A '],
            'exclusion': ['no'],
            'nominal_length': [200],
        }))
        tasks = build_quant_tasks(metadata)
        assert tasks == [('SRR001', 'Species A')]
        assert metadata.df.loc[0, 'run'] == 'SRR001'
        assert metadata.df.loc[0, 'scientific_name'] == 'Species A'

    def test_build_quant_tasks_rejects_duplicate_run_ids(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', ' SRR001 '],
            'scientific_name': ['Species A', 'Species A'],
            'exclusion': ['no', 'no'],
            'nominal_length': [200, 200],
        }))
        with pytest.raises(ValueError, match='Duplicate run ID in metadata for run\\(s\\): SRR001'):
            build_quant_tasks(metadata)

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
            (tmp_path / 'run_info.json').write_text('{}')
            (tmp_path / 'abundance.tsv').write_text('target_id\tlength\teff_length\test_counts\ttpm\n')
            (tmp_path / 'abundance.h5').write_bytes(b'h5')
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

    def test_call_kallisto_tolerates_non_utf8_process_output(self, tmp_path, monkeypatch):
        args = SimpleNamespace(threads=1)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'nominal_length': [250],
            'exclusion': ['no'],
        }))
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        def fake_run(cmd, stdout, stderr):
            (tmp_path / 'run_info.json').write_text('{}')
            (tmp_path / 'abundance.tsv').write_text('target_id\tlength\teff_length\test_counts\ttpm\n')
            return SimpleNamespace(returncode=0, stdout=b'\xff\xfe', stderr=b'\xff')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        observed = call_kallisto(
            args=args,
            in_files=['dummy.fastq.gz'],
            metadata=metadata,
            sra_stat=sra_stat,
            output_dir=str(tmp_path),
            index='dummy.idx',
        )

        assert observed.returncode == 0

    def test_call_kallisto_raises_when_required_output_file_missing(self, tmp_path, monkeypatch):
        args = SimpleNamespace(threads=1)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'nominal_length': [250],
            'exclusion': ['no'],
        }))
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}

        def fake_run(cmd, stdout, stderr):
            (tmp_path / 'run_info.json').write_text('{}')
            return SimpleNamespace(returncode=0, stdout=b'', stderr=b'')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        with pytest.raises(FileNotFoundError, match='kallisto output file was not generated'):
            call_kallisto(
                args=args,
                in_files=['dummy.fastq.gz'],
                metadata=metadata,
                sra_stat=sra_stat,
                output_dir=str(tmp_path),
                index='dummy.idx',
            )

    def test_call_kallisto_raises_on_nonzero_return_and_skips_rename(self, tmp_path, monkeypatch):
        args = SimpleNamespace(threads=1)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'nominal_length': [250],
            'exclusion': ['no'],
        }))
        sra_stat = {'sra_id': 'SRR001', 'layout': 'single'}
        renamed = {'called': False}

        def fake_run(cmd, stdout, stderr):
            return SimpleNamespace(returncode=1, stdout=b'', stderr=b'failure')

        def fake_rename(*_args, **_kwargs):
            renamed['called'] = True

        monkeypatch.setattr(subprocess, 'run', fake_run)
        monkeypatch.setattr('amalgkit.quant.rename_kallisto_outputs', fake_rename)

        with pytest.raises(RuntimeError, match='kallisto quant failed with exit code 1'):
            call_kallisto(
                args=args,
                in_files=['dummy.fastq.gz'],
                metadata=metadata,
                sra_stat=sra_stat,
                output_dir=str(tmp_path),
                index='dummy.idx',
            )
        assert renamed['called'] is False

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

    def test_get_index_waits_for_existing_builder_lock(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        index_dir.mkdir()
        lock_path = index_dir / '.Homo_sapiens.idx.lock'
        lock_path.write_text('12345\n')

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir='inferred',
        )

        index_path = index_dir / 'Homo_sapiens.idx'

        def fake_sleep(_seconds):
            index_path.write_text('index_ready')
            lock_path.unlink()

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('subprocess.run should not be called when waiting for another builder.')

        monkeypatch.setattr('amalgkit.quant.time.sleep', fake_sleep)
        monkeypatch.setattr(subprocess, 'run', fail_if_called)
        observed = get_index(args, 'Homo_sapiens')
        assert observed == str(index_path)

    def test_get_index_rejects_lock_path_that_is_directory(self, tmp_path):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        (fasta_dir / 'Homo_sapiens.fa').write_text('>a\nAAAA\n')
        lock_path = index_dir / '.Homo_sapiens.idx.lock'
        lock_path.mkdir()

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
            index_lock_poll=1,
            index_lock_timeout=1,
        )

        with pytest.raises(IsADirectoryError, match='Index lock path exists but is not a file'):
            get_index(args, 'Homo_sapiens')

    def test_get_index_rejects_lock_path_that_is_symlink(self, tmp_path):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        (fasta_dir / 'Homo_sapiens.fa').write_text('>a\nAAAA\n')
        lock_path = index_dir / '.Homo_sapiens.idx.lock'
        os.symlink(index_dir / 'dangling.lock.target', lock_path)

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
            index_lock_poll=1,
            index_lock_timeout=1,
        )

        with pytest.raises(IsADirectoryError, match='Index lock path exists but is not a file'):
            get_index(args, 'Homo_sapiens')

    def test_get_index_releases_lock_after_build(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        (fasta_dir / 'Homo_sapiens.fa').write_text('>a\nAAAA\n')
        index_path = index_dir / 'Homo_sapiens.idx'
        lock_path = index_dir / '.Homo_sapiens.idx.lock'

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
        )

        def fake_run(cmd, stdout, stderr):
            assert cmd[:3] == ['kallisto', 'index', '-i']
            index_path.write_text('built')
            return SimpleNamespace(returncode=0, stdout=b'', stderr=b'')

        monkeypatch.setattr(subprocess, 'run', fake_run)
        observed = get_index(args, 'Homo_sapiens')
        assert observed == str(index_path)
        assert not lock_path.exists()

    def test_get_index_uses_prefetched_fasta_entries(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        index_path = index_dir / 'Homo_sapiens.idx'
        prefetched_entries = {'Homo_sapiens.fa'}
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
            _prefetched_fasta_entries=prefetched_entries,
            _prefetched_fasta_dir=os.path.realpath(str(fasta_dir)),
        )
        seen = {'entries': None}

        def fake_find_species_fasta_files(path_fasta_dir, sci_name, entries=None):
            seen['entries'] = entries
            return [os.path.join(path_fasta_dir, sci_name + '.fa')]

        def fake_run(cmd, stdout, stderr):
            index_path.write_text('built')
            return SimpleNamespace(returncode=0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.quant.find_species_fasta_files', fake_find_species_fasta_files)
        monkeypatch.setattr(subprocess, 'run', fake_run)

        observed = get_index(args, 'Homo_sapiens')

        assert observed == str(index_path)
        assert seen['entries'] is prefetched_entries

    def test_get_index_uses_prefetched_index_entries_first(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        index_dir.mkdir()
        index_path = index_dir / 'Homo_sapiens.idx'
        index_path.write_text('ready')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=False,
            fasta_dir='inferred',
            _prefetched_index_entries={'Homo_sapiens.idx'},
            _prefetched_index_dir=os.path.realpath(str(index_dir)),
        )

        def fail_if_called(_path_dir):
            raise AssertionError('list_dir_entries should not be used when prefetched index entries are available.')

        monkeypatch.setattr('amalgkit.quant.list_dir_entries', fail_if_called)

        observed = get_index(args, 'Homo_sapiens')

        assert observed == str(index_path)

    def test_get_index_uses_subspecies_fallback_prefix(self, tmp_path):
        out_dir = tmp_path / 'out'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        index_dir.mkdir()
        index_path = index_dir / 'Gorilla_gorilla.idx'
        index_path.write_text('ready')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=False,
            fasta_dir='inferred',
        )

        observed = get_index(args, 'Gorilla_gorilla_gorilla')

        assert observed == str(index_path)

    def test_get_index_normalizes_redundant_underscores_for_subspecies_fallback(self, tmp_path):
        out_dir = tmp_path / 'out'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        index_dir.mkdir()
        index_path = index_dir / 'Canis_lupus.idx'
        index_path.write_text('ready')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=False,
            fasta_dir='inferred',
        )

        observed = get_index(args, 'Canis__lupus_familiaris')

        assert observed == str(index_path)

    def test_get_index_rejects_nonpositive_lock_poll(self, tmp_path):
        out_dir = tmp_path / 'out'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        index_dir.mkdir()
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=False,
            fasta_dir='inferred',
            index_lock_poll=0,
        )
        with pytest.raises(ValueError, match='--index_lock_poll must be > 0'):
            get_index(args, 'Homo_sapiens')

    def test_get_index_rejects_nonpositive_lock_timeout(self, tmp_path):
        out_dir = tmp_path / 'out'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        index_dir.mkdir()
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=False,
            fasta_dir='inferred',
            index_lock_timeout=0,
        )
        with pytest.raises(ValueError, match='--index_lock_timeout must be > 0'):
            get_index(args, 'Homo_sapiens')

    def test_get_index_rejects_index_path_that_is_file(self, tmp_path):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        index_file = tmp_path / 'index'
        index_file.write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_file),
            build_index=False,
            fasta_dir='inferred',
        )
        with pytest.raises(NotADirectoryError, match='not a directory'):
            get_index(args, 'Homo_sapiens')

    def test_get_index_ignores_directory_named_idx_file(self, tmp_path):
        out_dir = tmp_path / 'out'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        index_dir.mkdir()
        (index_dir / 'Homo_sapiens.idx').mkdir()
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=False,
            fasta_dir='inferred',
            index_lock_poll=1,
            index_lock_timeout=10,
        )
        with pytest.raises(FileNotFoundError, match='Could not find index file'):
            get_index(args, 'Homo_sapiens')

    def test_get_index_rejects_directory_when_build_reports_success(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        (fasta_dir / 'Homo_sapiens.fa').write_text('>a\nAAAA\n')
        (index_dir / 'Homo_sapiens.idx').mkdir()

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
            index_lock_poll=1,
            index_lock_timeout=10,
        )

        monkeypatch.setattr(
            subprocess,
            'run',
            lambda cmd, stdout, stderr: SimpleNamespace(returncode=0, stdout=b'', stderr=b''),
        )

        with pytest.raises(RuntimeError, match='Index file was not generated'):
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

    def test_quant_main_rejects_nonpositive_jobs(self):
        args = SimpleNamespace(internal_jobs=0)
        with pytest.raises(ValueError, match='--internal_jobs must be > 0'):
            quant_main(args)

    def test_quant_main_rejects_out_dir_file_path(self, tmp_path, monkeypatch):
        out_path = tmp_path / 'out_path'
        out_path.write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_path),
            internal_jobs=1,
            threads=1,
            redo=False,
            metadata='inferred',
            index_dir=None,
            build_index=False,
            fasta_dir='inferred',
            clean_fastq=False,
            index_lock_poll=1,
            index_lock_timeout=10,
        )
        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda: None)
        monkeypatch.setattr(
            'amalgkit.quant.load_metadata',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('load_metadata should not be called')),
        )

        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            quant_main(args)

    def test_quant_main_rejects_quant_path_file(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        (out_dir / 'quant').write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            internal_jobs=1,
            threads=1,
            redo=False,
            metadata='inferred',
            index_dir=None,
            build_index=False,
            fasta_dir='inferred',
            clean_fastq=False,
            index_lock_poll=1,
            index_lock_timeout=10,
        )
        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda: None)
        monkeypatch.setattr(
            'amalgkit.quant.load_metadata',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('load_metadata should not be called')),
        )

        with pytest.raises(NotADirectoryError, match='Quant path exists but is not a directory'):
            quant_main(args)

    def test_quant_main_parallel_dispatch(self, tmp_path, monkeypatch):
        args = SimpleNamespace(
            out_dir=str(tmp_path),
            internal_jobs=2,
            threads=2,
            redo=False,
            metadata='inferred',
            index_dir=None,
            build_index=False,
            fasta_dir='inferred',
            clean_fastq=False,
            index_lock_poll=1,
            index_lock_timeout=10,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Species A', 'Species B'],
            'exclusion': ['no', 'no'],
            'nominal_length': [200, 200],
        }))
        dispatched = []

        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda: None)
        monkeypatch.setattr('amalgkit.quant.load_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.quant.pre_resolve_species_indices', lambda _args, _tasks: {})
        monkeypatch.setattr(
            'amalgkit.quant.run_quant_for_sra',
            lambda _args, _metadata, sra_id, sci_name: dispatched.append((sra_id, sci_name)),
        )

        quant_main(args)

        assert set(dispatched) == {('SRR001', 'Species A'), ('SRR002', 'Species B')}

    def test_quant_main_parallel_worker_system_exit_is_reported_as_runtime_error(self, tmp_path, monkeypatch):
        args = SimpleNamespace(
            out_dir=str(tmp_path),
            internal_jobs=2,
            threads=2,
            redo=False,
            metadata='inferred',
            index_dir=None,
            build_index=False,
            fasta_dir='inferred',
            clean_fastq=False,
            index_lock_poll=1,
            index_lock_timeout=10,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Species A', 'Species B'],
            'exclusion': ['no', 'no'],
            'nominal_length': [200, 200],
        }))

        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda: None)
        monkeypatch.setattr('amalgkit.quant.load_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.quant.pre_resolve_species_indices', lambda _args, _tasks: {})

        def fake_run_quant_for_sra(_args, _metadata, sra_id, _sci_name):
            if sra_id == 'SRR002':
                raise SystemExit(1)
            return None

        monkeypatch.setattr('amalgkit.quant.run_quant_for_sra', fake_run_quant_for_sra)

        with pytest.raises(RuntimeError, match='quant failed for 1/2 SRA runs'):
            quant_main(args)

    def test_quant_main_cpu_budget_caps_jobs_to_serial(self, tmp_path, monkeypatch):
        args = SimpleNamespace(
            out_dir=str(tmp_path),
            internal_jobs=4,
            threads=4,
            internal_cpu_budget=1,
            redo=False,
            metadata='inferred',
            index_dir=None,
            build_index=False,
            fasta_dir='inferred',
            clean_fastq=False,
            index_lock_poll=1,
            index_lock_timeout=10,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Species A', 'Species B'],
            'exclusion': ['no', 'no'],
            'nominal_length': [200, 200],
        }))
        dispatched = []

        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda: None)
        monkeypatch.setattr('amalgkit.quant.load_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.quant.pre_resolve_species_indices', lambda _args, _tasks: {})
        monkeypatch.setattr(
            'amalgkit.quant.run_quant_for_sra',
            lambda _args, _metadata, sra_id, sci_name: dispatched.append((sra_id, sci_name)),
        )

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('run_tasks_with_optional_threads should not be used when --internal_cpu_budget caps internal_jobs to 1.')

        monkeypatch.setattr('amalgkit.quant.run_tasks_with_optional_threads', fail_if_called)

        quant_main(args)

        assert set(dispatched) == {('SRR001', 'Species A'), ('SRR002', 'Species B')}
        assert args.threads == 1

    def test_quant_main_pre_resolves_index_once_per_species(self, tmp_path, monkeypatch):
        args = SimpleNamespace(
            out_dir=str(tmp_path),
            internal_jobs=2,
            threads=1,
            redo=False,
            metadata='inferred',
            index_dir=None,
            build_index=False,
            fasta_dir='inferred',
            clean_fastq=False,
            index_lock_poll=1,
            index_lock_timeout=10,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'scientific_name': ['Species A', 'Species A', 'Species B'],
            'exclusion': ['no', 'no', 'no'],
            'nominal_length': [200, 200, 200],
        }))
        dispatched = []
        resolved_species = []

        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda: None)
        monkeypatch.setattr('amalgkit.quant.load_metadata', lambda _args: metadata)
        monkeypatch.setattr(
            'amalgkit.quant.get_index',
            lambda _args, sci_name: resolved_species.append(sci_name) or '{}.idx'.format(sci_name),
        )
        monkeypatch.setattr(
            'amalgkit.quant.run_quant_for_sra',
            lambda _args, _metadata, sra_id, sci_name: dispatched.append((sra_id, sci_name)),
        )

        quant_main(args)

        assert resolved_species.count('Species_A') == 1
        assert resolved_species.count('Species_B') == 1
        assert set(dispatched) == {('SRR001', 'Species A'), ('SRR002', 'Species A'), ('SRR003', 'Species B')}

    def test_pre_resolve_species_indices_prefetches_fasta_entries(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        fasta_dir = out_dir / 'fasta'
        index_dir = out_dir / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir(parents=True)
        index_dir.mkdir(parents=True)
        (fasta_dir / 'Species_A.fa').write_text('>a\nAAAA\n')
        (fasta_dir / 'Species_B.fa').write_text('>b\nCCCC\n')
        (index_dir / 'Species_A.idx').write_text('idx')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            build_index=True,
            fasta_dir='inferred',
            index_dir=None,
        )
        tasks = [('SRR001', 'Species A'), ('SRR002', 'Species B')]
        resolved = []

        monkeypatch.setattr('amalgkit.quant.get_index', lambda _args, sci_name: resolved.append(sci_name) or sci_name + '.idx')
        out = pre_resolve_species_indices(args, tasks)

        assert out == {'Species_A': 'Species_A.idx', 'Species_B': 'Species_B.idx'}
        assert resolved == ['Species_A', 'Species_B']
        assert getattr(args, '_prefetched_fasta_entries') == {'Species_A.fa', 'Species_B.fa'}
        assert getattr(args, '_prefetched_fasta_dir') == os.path.realpath(str(fasta_dir))
        assert getattr(args, '_prefetched_index_entries') == {'Species_A.idx'}
        assert getattr(args, '_prefetched_index_dir') == os.path.realpath(str(index_dir))

    def test_pre_resolve_species_indices_normalizes_redundant_species_whitespace(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        index_dir = out_dir / 'index'
        out_dir.mkdir()
        index_dir.mkdir(parents=True)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            build_index=False,
            fasta_dir='inferred',
            index_dir=None,
            internal_jobs=1,
        )
        tasks = [('SRR001', 'Species  A'), ('SRR002', 'Species A')]
        resolved = []

        monkeypatch.setattr('amalgkit.quant.get_index', lambda _args, sci_name: resolved.append(sci_name) or sci_name + '.idx')

        out = pre_resolve_species_indices(args, tasks)

        assert resolved == ['Species_A']
        assert out == {'Species_A': 'Species_A.idx'}

    def test_pre_resolve_species_indices_uses_parallel_jobs(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        index_dir = out_dir / 'index'
        out_dir.mkdir()
        index_dir.mkdir(parents=True)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            build_index=False,
            fasta_dir='inferred',
            index_dir=None,
            internal_jobs=4,
        )
        tasks = [('SRR001', 'Species A'), ('SRR002', 'Species B')]
        seen = {'max_workers': None}

        def fake_run_tasks(task_items, task_fn, max_workers):
            items = list(task_items)
            seen['max_workers'] = max_workers
            return {item: task_fn(item) for item in items}, []

        monkeypatch.setattr('amalgkit.quant.run_tasks_with_optional_threads', fake_run_tasks)
        monkeypatch.setattr('amalgkit.quant.get_index', lambda _args, sci_name: sci_name + '.idx')

        out = pre_resolve_species_indices(args, tasks)

        assert seen['max_workers'] == 2
        assert out == {'Species_A': 'Species_A.idx', 'Species_B': 'Species_B.idx'}

    def test_pre_resolve_species_indices_rejects_index_path_that_is_file(self, tmp_path):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        index_file = tmp_path / 'index_path'
        index_file.write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            build_index=False,
            fasta_dir='inferred',
            index_dir=str(index_file),
        )
        tasks = [('SRR001', 'Species A')]

        with pytest.raises(NotADirectoryError, match='not a directory'):
            pre_resolve_species_indices(args, tasks)

    def test_pre_resolve_species_indices_rejects_fasta_path_that_is_file(self, tmp_path):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        fasta_file = tmp_path / 'fasta_path'
        fasta_file.write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            build_index=True,
            fasta_dir=str(fasta_file),
            index_dir=None,
        )
        tasks = [('SRR001', 'Species A')]

        with pytest.raises(NotADirectoryError, match='Fasta path exists but is not a directory'):
            pre_resolve_species_indices(args, tasks)
