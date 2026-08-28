import os
import json
import pathlib

import pytest
import numpy
import pandas
from types import SimpleNamespace

from amalgkit.command_context import QuantRuntimeContext
from amalgkit.quant import (
    _metadata_with_quant_input_sra_stats,
    run_quant,
    run_quant_for_sra,
    prefetch_getfastq_run_files,
)
from amalgkit.util import Metadata


def write_valid_quant_outputs(output_dir, sra_id='SRR001', target_id='tx1'):
    os.makedirs(output_dir, exist_ok=True)
    pandas.DataFrame([{
        'target_id': target_id,
        'length': 100,
        'eff_length': 90,
        'est_counts': 5,
        'tpm': 1.0,
    }]).to_csv(os.path.join(output_dir, sra_id + '_abundance.tsv'), sep='\t', index=False)
    with open(os.path.join(output_dir, sra_id + '_run_info.json'), 'w', encoding='utf-8') as handle:
        json.dump({'p_pseudoaligned': 50.0}, handle)


def write_safely_removed_fastq_marker(out_dir, sra_id='SRR001'):
    run_dir = pathlib.Path(out_dir) / 'getfastq' / sra_id
    run_dir.mkdir(parents=True, exist_ok=True)
    marker = run_dir / (sra_id + '.fastq.gz.safely_removed')
    marker.write_text(
        'This fastq file was safely removed after `amalgkit quant`.',
        encoding='utf-8',
    )
    return marker


def make_single_run_quant_metadata():
    return Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['SRR001'],
        'scientific_name': ['Species A'],
        'lib_layout': ['single'],
        'total_spots': [10],
        'spot_length': [100],
        'total_bases': [1000],
        'nominal_length': [200],
        'exclusion': ['no'],
    }))


# ---------------------------------------------------------------------------
# quant_output_exists (checks for abundance.tsv in output directory)
# ---------------------------------------------------------------------------


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

    def test_prefetch_getfastq_run_files_accepts_symlinked_input_root(self, tmp_path):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        shared_getfastq = tmp_path / 'shared_getfastq'
        run_dir = shared_getfastq / 'SRR001'
        run_dir.mkdir(parents=True)
        (run_dir / 'SRR001.fastq.gz').write_text('reads', encoding='utf-8')
        (out_dir / 'getfastq').symlink_to(shared_getfastq, target_is_directory=True)
        args = SimpleNamespace(out_dir=str(out_dir))

        observed = prefetch_getfastq_run_files(args, [('SRR001', 'Species A')])

        assert observed == {'SRR001': {'SRR001.fastq.gz'}}

    def test_prefetch_getfastq_run_files_rejects_symlinked_run_directory(self, tmp_path):
        out_dir = tmp_path / 'out'
        getfastq_root = out_dir / 'getfastq'
        getfastq_root.mkdir(parents=True)
        outside = tmp_path / 'outside'
        outside.mkdir()
        (getfastq_root / 'SRR001').symlink_to(outside, target_is_directory=True)
        args = SimpleNamespace(out_dir=str(out_dir))

        with pytest.raises(ValueError, match='symbolic-link run ID'):
            prefetch_getfastq_run_files(args, [('SRR001', 'Species A')])

    def test_run_quant_uses_prefetched_getfastq_files(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(run_files_by_run={'SRR001': {'SRR001.fastq.gz'}})
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
            write_valid_quant_outputs(_output_dir)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr('amalgkit.quant.list_getfastq_run_files', fail_if_listdir_used)
        monkeypatch.setattr('amalgkit.quant.call_kallisto', fake_call_kallisto)

        run_quant(args, metadata, 'SRR001', 'dummy.idx', runtime_context=runtime_context)

        assert called['kallisto'] == 1

    def test_run_quant_accepts_symlinked_read_only_getfastq_root(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        shared_getfastq = tmp_path / 'shared_getfastq'
        run_dir = shared_getfastq / 'SRR001'
        run_dir.mkdir(parents=True)
        (run_dir / 'SRR001.fastq.gz').write_text('reads', encoding='utf-8')
        (out_dir / 'getfastq').symlink_to(shared_getfastq, target_is_directory=True)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        metadata = make_single_run_quant_metadata()
        observed = {}

        def fake_call_kallisto(_args, in_files, _metadata, _sra_stat, output_dir, _index):
            observed['in_files'] = in_files
            write_valid_quant_outputs(output_dir)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr('amalgkit.quant.call_kallisto', fake_call_kallisto)

        run_quant(args, metadata, 'SRR001', 'dummy.idx', backend='kallisto')

        assert observed['in_files'] == [str(run_dir / 'SRR001.fastq.gz')]

    def test_run_quant_uses_getfastq_stats_without_mutating_source_metadata(self, tmp_path, monkeypatch, capsys):
        out_dir = tmp_path / 'out'
        run_dir = out_dir / 'getfastq' / 'SRR001'
        run_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'num_written': 12,
            'bp_fastp_in': 2400,
        }]).to_csv(run_dir / 'getfastq_stats.tsv', sep='\t', index=False)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(run_files_by_run={'SRR001': {'SRR001_1.fastq.gz', 'SRR001_2.fastq.gz'}})
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'platform': ['ILLUMINA'],
            'lib_layout': ['paired'],
            'total_spots': [numpy.nan],
            'spot_length': [numpy.nan],
            'total_bases': [numpy.nan],
            'nominal_length': [100],
            'exclusion': ['no'],
        }))
        observed = {}

        def fake_call_kallisto(_args, in_files, _metadata, sra_stat, _output_dir, _index):
            observed['in_files'] = in_files
            observed['sra_stat'] = sra_stat
            write_valid_quant_outputs(_output_dir)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr('amalgkit.quant.call_kallisto', fake_call_kallisto)

        run_quant(args, metadata, 'SRR001', 'dummy.idx', runtime_context=runtime_context)

        assert observed['in_files'] == [
            str(run_dir / 'SRR001_1.fastq.gz'),
            str(run_dir / 'SRR001_2.fastq.gz'),
        ]
        assert observed['sra_stat']['total_spot'] == 12
        assert observed['sra_stat']['spot_length'] == 200
        assert pandas.isna(metadata.df.loc[0, 'total_spots'])
        assert pandas.isna(metadata.df.loc[0, 'total_bases'])
        assert pandas.isna(metadata.df.loc[0, 'spot_length'])
        stderr = capsys.readouterr().err
        assert 'Using extracted-read statistics as a quant-only fallback' in stderr
        assert 'The source metadata is unchanged.' in stderr

    def test_run_quant_auto_infers_backend_from_getfastq_stats(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        run_dir = out_dir / 'getfastq' / 'SRR001'
        run_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'num_written': 10,
            'bp_fastp_in': 1000,
        }]).to_csv(run_dir / 'getfastq_stats.tsv', sep='\t', index=False)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
            quant_backend='auto',
            oarfish_seq_tech='auto',
        )
        runtime_context = QuantRuntimeContext(
            run_files_by_run={'SRR001': {'SRR001.fastq.gz'}},
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'lib_layout': ['single'],
            'total_spots': [numpy.nan],
            'spot_length': [numpy.nan],
            'total_bases': [numpy.nan],
            'nominal_length': [200],
            'exclusion': ['no'],
        }))
        observed = {}

        def fake_call_kallisto(_args, _in_files, quant_metadata, sra_stat, output_dir, _index):
            observed['spot_length'] = quant_metadata.df.loc[0, 'spot_length']
            observed['total_spot'] = sra_stat['total_spot']
            write_valid_quant_outputs(output_dir)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr('amalgkit.quant.call_kallisto', fake_call_kallisto)

        run_quant(args, metadata, 'SRR001', 'dummy.idx', runtime_context=runtime_context)

        assert observed == {'spot_length': 100.0, 'total_spot': 10}
        assert pandas.isna(metadata.df.loc[0, 'total_spots'])
        assert pandas.isna(metadata.df.loc[0, 'total_bases'])
        assert pandas.isna(metadata.df.loc[0, 'spot_length'])

    def test_run_quant_for_sra_auto_recovers_stats_before_backend_inference(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        run_dir = out_dir / 'getfastq' / 'SRR001'
        run_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'num_written': 10,
            'bp_fastp_in': 1000,
        }]).to_csv(run_dir / 'getfastq_stats.tsv', sep='\t', index=False)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            quant_backend='auto',
            oarfish_seq_tech='auto',
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'lib_layout': ['single'],
            'total_spots': [numpy.nan],
            'spot_length': [numpy.nan],
            'total_bases': [numpy.nan],
            'exclusion': ['no'],
        }))
        observed = {}

        monkeypatch.setattr(
            'amalgkit.quant.get_index',
            lambda *_args, **_kwargs: 'dummy.idx',
        )

        def fake_run_quant(_args, quant_metadata, _sra_id, _index, **kwargs):
            observed['backend'] = kwargs['backend']
            observed['spot_length'] = quant_metadata.df.loc[0, 'spot_length']

        monkeypatch.setattr('amalgkit.quant.run_quant', fake_run_quant)

        run_quant_for_sra(args, metadata, 'SRR001', 'Species A')

        assert observed == {'backend': 'kallisto', 'spot_length': 100.0}
        assert pandas.isna(metadata.df.loc[0, 'spot_length'])

    def test_quant_stats_fallback_supports_single_end_runs(self, tmp_path):
        pandas.DataFrame([{
            'run': 'SRR001',
            'num_written': 10,
            'bp_fastp_in': 1000,
        }]).to_csv(tmp_path / 'getfastq_stats.tsv', sep='\t', index=False)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'lib_layout': ['single'],
            'total_spots': [numpy.nan],
            'spot_length': [numpy.nan],
            'total_bases': [numpy.nan],
            'exclusion': ['no'],
        }))

        quant_metadata = _metadata_with_quant_input_sra_stats(metadata, 'SRR001', str(tmp_path))

        assert int(quant_metadata.df.loc[0, 'total_spots']) == 10
        assert int(quant_metadata.df.loc[0, 'total_bases']) == 1000
        assert int(quant_metadata.df.loc[0, 'spot_length']) == 100
        assert pandas.isna(metadata.df.loc[0, 'total_spots'])
        assert pandas.isna(metadata.df.loc[0, 'total_bases'])
        assert pandas.isna(metadata.df.loc[0, 'spot_length'])

    def test_quant_stats_fallback_normalizes_nullable_string_numeric_columns(self, tmp_path):
        pandas.DataFrame([{
            'run': 'SRR001',
            'num_written': 10,
            'bp_fastp_in': 1497,
        }]).to_csv(tmp_path / 'getfastq_stats.tsv', sep='\t', index=False)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'lib_layout': ['single'],
            'total_spots': pandas.Series([pandas.NA], dtype='string'),
            'spot_length': pandas.Series([pandas.NA], dtype='string'),
            'total_bases': pandas.Series([pandas.NA], dtype='string'),
            'exclusion': ['no'],
        }))

        quant_metadata = _metadata_with_quant_input_sra_stats(
            metadata,
            'SRR001',
            str(tmp_path),
        )

        assert quant_metadata.df.loc[0, 'total_spots'] == 10
        assert quant_metadata.df.loc[0, 'total_bases'] == 1497
        assert quant_metadata.df.loc[0, 'spot_length'] == 149.7
        assert str(quant_metadata.df['spot_length'].dtype) == 'float64'
        assert pandas.isna(metadata.df.loc[0, 'spot_length'])
        assert str(metadata.df['spot_length'].dtype) == 'string'

    def test_run_quant_preserves_valid_sra_stats_when_getfastq_stats_differ(self, tmp_path, monkeypatch, capsys):
        out_dir = tmp_path / 'out'
        run_dir = out_dir / 'getfastq' / 'SRR001'
        run_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'num_written': 12,
            'bp_fastp_in': 2400,
        }]).to_csv(run_dir / 'getfastq_stats.tsv', sep='\t', index=False)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(run_files_by_run={'SRR001': {'SRR001_1.fastq.gz', 'SRR001_2.fastq.gz'}})
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'platform': ['ILLUMINA'],
            'lib_layout': ['paired'],
            'total_spots': [100],
            'spot_length': [150],
            'total_bases': [15000],
            'nominal_length': [200],
            'exclusion': ['no'],
        }))
        observed = {}

        def fake_call_kallisto(_args, _in_files, _metadata, sra_stat, _output_dir, _index):
            observed['sra_stat'] = sra_stat
            write_valid_quant_outputs(_output_dir)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr('amalgkit.quant.call_kallisto', fake_call_kallisto)

        run_quant(args, metadata, 'SRR001', 'dummy.idx', runtime_context=runtime_context)

        assert observed['sra_stat']['total_spot'] == 100
        assert observed['sra_stat']['spot_length'] == 150
        assert 'quant-only fallback' not in capsys.readouterr().err

    def test_run_quant_combines_raw_total_spots_with_extracted_spot_length(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        run_dir = out_dir / 'getfastq' / 'SRR001'
        run_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'num_written': 12,
            'bp_fastp_in': 1800,
            'bp_rrna_in': 2400,
        }]).to_csv(run_dir / 'getfastq_stats.tsv', sep='\t', index=False)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(run_files_by_run={'SRR001': {'SRR001_1.fastq.gz', 'SRR001_2.fastq.gz'}})
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'platform': ['ILLUMINA'],
            'lib_layout': ['paired'],
            'total_spots': [100],
            'spot_length': [numpy.nan],
            'total_bases': [numpy.nan],
            'nominal_length': [200],
            'exclusion': ['no'],
        }))
        observed = {}

        def fake_call_kallisto(_args, _in_files, _metadata, sra_stat, _output_dir, _index):
            observed['sra_stat'] = sra_stat
            write_valid_quant_outputs(_output_dir)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr('amalgkit.quant.call_kallisto', fake_call_kallisto)

        run_quant(args, metadata, 'SRR001', 'dummy.idx', runtime_context=runtime_context)

        assert observed['sra_stat']['total_spot'] == 100
        assert observed['sra_stat']['spot_length'] == 200
        assert pandas.isna(metadata.df.loc[0, 'spot_length'])
        assert pandas.isna(metadata.df.loc[0, 'total_bases'])

    def test_run_quant_rejects_mismatched_getfastq_stats_run(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        run_dir = out_dir / 'getfastq' / 'SRR001'
        run_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR999',
            'num_written': 12,
            'bp_fastp_in': 2400,
        }]).to_csv(run_dir / 'getfastq_stats.tsv', sep='\t', index=False)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(run_files_by_run={'SRR001': {'SRR001_1.fastq.gz', 'SRR001_2.fastq.gz'}})
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'platform': ['ILLUMINA'],
            'lib_layout': ['paired'],
            'total_spots': [numpy.nan],
            'spot_length': [numpy.nan],
            'total_bases': [numpy.nan],
            'nominal_length': [200],
            'exclusion': ['no'],
        }))
        monkeypatch.setattr(
            'amalgkit.quant.call_kallisto',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('call_kallisto should not be called')),
        )

        with pytest.raises(ValueError, match='total_spots must be > 0'):
            run_quant(args, metadata, 'SRR001', 'dummy.idx', runtime_context=runtime_context)

    def test_run_quant_redo_failure_preserves_previous_valid_outputs(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        quant_run_dir = out_dir / 'quant' / 'SRR001'
        quant_run_dir.mkdir(parents=True)
        write_valid_quant_outputs(quant_run_dir, target_id='old')
        old_abundance = (quant_run_dir / 'SRR001_abundance.tsv').read_text()
        old_run_info = (quant_run_dir / 'SRR001_run_info.json').read_text()
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=True,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(run_files_by_run={'SRR001': {'SRR001.fastq.gz'}})
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
        def fake_call_kallisto(_args, in_files, _metadata, _sra_stat, _output_dir, _index):
            assert (quant_run_dir / 'SRR001_abundance.tsv').exists()
            assert len(in_files) == 1
            raise RuntimeError('simulated backend failure')

        monkeypatch.setattr('amalgkit.quant.call_kallisto', fake_call_kallisto)

        with pytest.raises(RuntimeError, match='simulated backend failure'):
            run_quant(args, metadata, 'SRR001', 'dummy.idx', runtime_context=runtime_context)
        assert (quant_run_dir / 'SRR001_abundance.tsv').read_text() == old_abundance
        assert (quant_run_dir / 'SRR001_run_info.json').read_text() == old_run_info

    def test_run_quant_rejects_safely_removed_fastq_when_outputs_are_missing(self, tmp_path):
        out_dir = tmp_path / 'out'
        marker = write_safely_removed_fastq_marker(out_dir)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(
            run_files_by_run={'SRR001': {marker.name}},
        )

        with pytest.raises(
            FileNotFoundError,
            match='valid quant outputs are missing or invalid.*amalgkit getfastq',
        ):
            run_quant(
                args,
                make_single_run_quant_metadata(),
                'SRR001',
                'dummy.idx',
                runtime_context=runtime_context,
            )

    def test_run_quant_rejects_safely_removed_fastq_when_outputs_are_invalid(self, tmp_path):
        out_dir = tmp_path / 'out'
        marker = write_safely_removed_fastq_marker(out_dir)
        quant_run_dir = out_dir / 'quant' / 'SRR001'
        quant_run_dir.mkdir(parents=True)
        (quant_run_dir / 'SRR001_abundance.tsv').write_text('invalid\n', encoding='utf-8')
        (quant_run_dir / 'SRR001_run_info.json').write_text('invalid\n', encoding='utf-8')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(
            run_files_by_run={'SRR001': {marker.name}},
        )

        with pytest.raises(
            FileNotFoundError,
            match='valid quant outputs are missing or invalid',
        ):
            run_quant(
                args,
                make_single_run_quant_metadata(),
                'SRR001',
                'dummy.idx',
                runtime_context=runtime_context,
            )

    def test_run_quant_rejects_redo_with_safely_removed_fastq_and_preserves_outputs(self, tmp_path):
        out_dir = tmp_path / 'out'
        marker = write_safely_removed_fastq_marker(out_dir)
        quant_run_dir = out_dir / 'quant' / 'SRR001'
        write_valid_quant_outputs(quant_run_dir, target_id='old')
        old_abundance = (quant_run_dir / 'SRR001_abundance.tsv').read_text(encoding='utf-8')
        old_run_info = (quant_run_dir / 'SRR001_run_info.json').read_text(encoding='utf-8')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=True,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(
            run_files_by_run={'SRR001': {marker.name}},
        )

        with pytest.raises(FileNotFoundError, match='--redo requires re-quantification'):
            run_quant(
                args,
                make_single_run_quant_metadata(),
                'SRR001',
                'dummy.idx',
                runtime_context=runtime_context,
            )

        assert (quant_run_dir / 'SRR001_abundance.tsv').read_text(encoding='utf-8') == old_abundance
        assert (quant_run_dir / 'SRR001_run_info.json').read_text(encoding='utf-8') == old_run_info

    def test_run_quant_reuses_valid_outputs_with_safely_removed_fastq(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        marker = write_safely_removed_fastq_marker(out_dir)
        write_valid_quant_outputs(out_dir / 'quant' / 'SRR001')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(
            run_files_by_run={'SRR001': {marker.name}},
        )
        monkeypatch.setattr(
            'amalgkit.quant.get_newest_intermediate_file_extension',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError('FASTQ marker should not be inspected when valid output is reused')
            ),
        )

        run_quant(
            args,
            make_single_run_quant_metadata(),
            'SRR001',
            'dummy.idx',
            runtime_context=runtime_context,
        )

    def test_run_quant_redo_preserves_only_unknown_regular_sidecars(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        quant_run_dir = out_dir / 'quant' / 'SRR001'
        quant_run_dir.mkdir(parents=True)
        write_valid_quant_outputs(quant_run_dir, target_id='old')
        (quant_run_dir / 'diagnostic.log').write_text('keep\n', encoding='utf-8')
        (quant_run_dir / 'SRR001.prob').write_text('stale-managed\n', encoding='utf-8')
        (quant_run_dir / 'debug-directory').mkdir()
        outside_log = tmp_path / 'outside.log'
        outside_log.write_text('outside\n', encoding='utf-8')
        (quant_run_dir / 'linked.log').symlink_to(outside_log)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=True,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(
            run_files_by_run={'SRR001': {'SRR001.fastq.gz'}}
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

        def fake_call_kallisto(_args, _in_files, _metadata, _sra_stat, stage_dir, _index):
            assert not (pathlib.Path(stage_dir) / 'diagnostic.log').exists()
            assert not (pathlib.Path(stage_dir) / 'SRR001.prob').exists()
            assert not (pathlib.Path(stage_dir) / 'debug-directory').exists()
            assert not (pathlib.Path(stage_dir) / 'linked.log').exists()
            write_valid_quant_outputs(stage_dir, target_id='new')
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr('amalgkit.quant.call_kallisto', fake_call_kallisto)

        run_quant(
            args,
            metadata,
            'SRR001',
            'dummy.idx',
            runtime_context=runtime_context,
        )

        abundance = pandas.read_csv(
            quant_run_dir / 'SRR001_abundance.tsv',
            sep='\t',
        )
        assert abundance['target_id'].tolist() == ['new']
        assert (quant_run_dir / 'diagnostic.log').read_text(encoding='utf-8') == 'keep\n'
        assert not (quant_run_dir / 'SRR001.prob').exists()
        assert not (quant_run_dir / 'debug-directory').exists()
        assert not (quant_run_dir / 'linked.log').exists()
        assert outside_log.read_text(encoding='utf-8') == 'outside\n'

    def test_run_quant_raises_when_fastq_is_missing(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
        )
        runtime_context = QuantRuntimeContext(run_files_by_run={'SRR001': set()})
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
            run_quant(args, metadata, 'SRR001', 'dummy.idx', runtime_context=runtime_context)

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
