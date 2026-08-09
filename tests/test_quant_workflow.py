import os
import json
import pathlib
import subprocess
from contextlib import contextmanager

import pytest
import numpy
import pandas
from types import SimpleNamespace

from amalgkit.command_context import PrefetchedDirEntries, QuantRuntimeContext
from amalgkit.quant import (
    _build_kallisto_index,
    _write_index_ready_marker,
    get_index,
    call_kallisto,
    call_oarfish,
    run_quant,
    pre_resolve_species_indices,
    build_quant_tasks,
    quant_main,
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


class TestQuantEdgeCases:
    def test_build_quant_tasks_uses_only_selected_nonexcluded_rows(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'scientific_name': ['Species A', 'Species A', 'Species A'],
            'exclusion': ['no', 'low_nspots', 'no'],
            'is_sampled': ['yes', 'yes', 'no'],
        }))

        assert build_quant_tasks(metadata) == [('SRR001', 'Species A')]
        assert metadata.df['run'].tolist() == ['SRR001']

    def test_build_quant_tasks_rejects_invalid_sampling_flag(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'exclusion': ['no'],
            'is_sampled': ['maybe'],
        }))

        with pytest.raises(ValueError, match='invalid flag'):
            build_quant_tasks(metadata)

    def test_build_quant_tasks_rejects_unsafe_run_identifier(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['../escape'],
            'scientific_name': ['Species A'],
            'exclusion': ['no'],
        }))

        with pytest.raises(ValueError, match='run ID'):
            build_quant_tasks(metadata)

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

    def test_build_quant_tasks_rejects_placeholder_scientific_name(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Please add in format: Genus species'],
            'exclusion': ['no'],
            'nominal_length': [200],
        }))
        with pytest.raises(
            ValueError,
            match='Placeholder scientific_name from amalgkit integrate was found for run\\(s\\): SRR001',
        ):
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

    def test_run_quant_ignores_scientific_name_original_for_index_lookup(self, tmp_path):
        out_dir = tmp_path / 'out'
        index_dir = out_dir / 'index'
        out_dir.mkdir()
        index_dir.mkdir(parents=True)
        (index_dir / 'Pygsuia_biforma.idx').write_text('idx')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
            build_index=False,
            fasta_dir='inferred',
            index_dir=None,
            internal_jobs=1,
            quant_backend='kallisto',
            index_lock_poll=1,
            index_lock_timeout=1,
        )
        runtime_context = QuantRuntimeContext(
            run_files_by_run={'SRR001': {'SRR001.fastq.gz'}},
            quant_backend_by_run={'SRR001': 'kallisto'},
            species_identifier_values_by_run={
                'SRR001': ['unclassified eukaryotes', 'Pygsuia biforma'],
            },
        )

        with pytest.raises(FileNotFoundError, match='Could not find index file'):
            pre_resolve_species_indices(
                args,
                [('SRR001', 'unclassified eukaryotes')],
                runtime_context=runtime_context,
            )

    def test_get_index_raises_clear_error_for_multiple_fasta_candidates(self, tmp_path):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        (fasta_dir / 'Homo_sapiens.fa').write_text('>a\nAAAA\n')
        (fasta_dir / 'Homo_sapiens_for_kallisto_index.fasta').write_text('>b\nCCCC\n')

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
        )
        with pytest.raises(ValueError, match='Found multiple reference fasta files'):
            get_index(args, 'Homo_sapiens')

    def test_get_index_prefers_exact_fasta_stem_before_broad_prefix_matches(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        exact_fasta = fasta_dir / 'Mayorella_for_kallisto_index.fasta'
        broad_fasta = fasta_dir / 'Mayorella_sp_for_kallisto_index.fasta'
        exact_fasta.write_text('>a\nAAAA\n')
        broad_fasta.write_text('>b\nCCCC\n')
        observed = {'fasta_file': None}

        def fake_build(index_path, fasta_file, sci_name):
            observed['fasta_file'] = fasta_file
            open(index_path, 'w').write('idx')

        monkeypatch.setattr('amalgkit.quant._build_kallisto_index', fake_build)

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
        )

        observed_index = get_index(args, 'Mayorella')

        assert observed_index == str(index_dir / 'Mayorella.idx')
        assert observed['fasta_file'] == str(exact_fasta)

    def test_get_index_ignores_alias_exact_match(self, tmp_path):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        alias_exact_fasta = fasta_dir / 'Mayorella_sp_for_kallisto_index.fasta'
        ambiguous_broad_fasta = fasta_dir / 'Mayorella_marianaensis_for_kallisto_index.fasta'
        alias_exact_fasta.write_text('>a\nAAAA\n')
        ambiguous_broad_fasta.write_text('>b\nCCCC\n')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
        )

        with pytest.raises(FileNotFoundError, match='Could not find reference fasta file'):
            get_index(args, 'Mayorella', alias_names=['Mayorella sp'])

    def test_get_index_does_not_use_variety_fasta_as_base_species_fallback(self, tmp_path):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        (fasta_dir / 'Asimitellaria_furusei_var._furusei_for_kallisto_index.fasta').write_text('>a\nAAAA\n')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
        )

        with pytest.raises(FileNotFoundError, match='Could not find reference fasta file'):
            get_index(args, 'Asimitellaria_furusei')

    def test_get_index_prefers_uncompressed_fasta_over_matching_gz_copy(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        fasta_plain = fasta_dir / 'Homo_sapiens.fa'
        fasta_gz = fasta_dir / 'Homo_sapiens.fa.gz'
        fasta_plain.write_text('>a\nAAAA\n')
        fasta_gz.write_text('gz')
        observed = {'fasta_file': None}

        def fake_build(index_path, fasta_file, sci_name):
            observed['fasta_file'] = fasta_file
            assert sci_name == 'Homo_sapiens'
            open(index_path, 'w').write('idx')

        monkeypatch.setattr('amalgkit.quant._build_kallisto_index', fake_build)

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
        )

        observed_index = get_index(args, 'Homo_sapiens')

        assert observed_index == str(index_dir / 'Homo_sapiens.idx')
        assert observed['fasta_file'] == str(fasta_plain)

    def test_get_index_requires_exact_fasta_stem_when_building_index(self, tmp_path):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        fasta_path = fasta_dir / 'Gorilla_gorilla.fa'
        fasta_path.write_text('>a\nAAAA\n')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
        )

        with pytest.raises(FileNotFoundError, match='Could not find reference fasta file'):
            get_index(args, 'Gorilla_gorilla_gorilla')

    def test_get_index_uses_shared_lock_and_rechecks_after_acquire(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        index_dir.mkdir()
        lock_path = index_dir / '.Homo_sapiens.idx.lock'

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir='inferred',
            index_lock_poll=7,
            index_lock_timeout=11,
        )

        index_path = index_dir / 'Homo_sapiens.idx'
        captured = {}

        @contextmanager
        def fake_acquire_exclusive_lock(lock_path, lock_label='Lock', poll_seconds=5, timeout_seconds=3600):
            captured['lock_path'] = lock_path
            captured['lock_label'] = lock_label
            captured['poll_seconds'] = poll_seconds
            captured['timeout_seconds'] = timeout_seconds
            index_path.write_text('index_ready')
            _write_index_ready_marker(str(index_path), backend='kallisto')
            yield

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('subprocess.run should not be called when waiting for another builder.')

        monkeypatch.setattr('amalgkit.quant.acquire_exclusive_lock', fake_acquire_exclusive_lock)
        monkeypatch.setattr(subprocess, 'run', fail_if_called)
        observed = get_index(args, 'Homo_sapiens')

        assert observed == str(index_path)
        assert captured == {
            'lock_path': str(lock_path),
            'lock_label': 'Index lock',
            'poll_seconds': 7,
            'timeout_seconds': 11,
        }

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
            with open(cmd[3], 'w', encoding='utf-8') as handle:
                handle.write('built')
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
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
        )
        runtime_context = QuantRuntimeContext(
            prefetched_fasta=PrefetchedDirEntries.from_entries(
                entries=['Homo_sapiens.fa'],
                path_dir=str(fasta_dir),
            ),
        )
        (fasta_dir / 'Homo_sapiens.fa').write_text('>a\nAAAA\n', encoding='utf-8')
        seen = {'entries': None}

        def fake_find_species_fasta_files(path_fasta_dir, sci_name, entries=None):
            seen['entries'] = entries
            return sci_name, [os.path.join(path_fasta_dir, sci_name + '.fa')]

        def fake_run(cmd, stdout, stderr):
            with open(cmd[cmd.index('-i') + 1], 'w', encoding='utf-8') as handle:
                handle.write('built')
            return SimpleNamespace(returncode=0, stdout=b'', stderr=b'')

        monkeypatch.setattr('amalgkit.quant._find_species_fasta_files', fake_find_species_fasta_files)
        monkeypatch.setattr(subprocess, 'run', fake_run)

        observed = get_index(args, 'Homo_sapiens', runtime_context=runtime_context)

        assert observed == str(index_path)
        assert seen['entries'] == ['Homo_sapiens.fa']

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
        )
        runtime_context = QuantRuntimeContext(
            prefetched_index=PrefetchedDirEntries.from_entries(
                entries=['Homo_sapiens.idx'],
                path_dir=str(index_dir),
            ),
        )

        def fail_if_called(_path_dir):
            raise AssertionError('list_dir_entries should not be used when prefetched index entries are available.')

        monkeypatch.setattr('amalgkit.quant.list_dir_entries', fail_if_called)

        observed = get_index(args, 'Homo_sapiens', runtime_context=runtime_context)

        assert observed == str(index_path)

    def test_get_index_requires_exact_index_stem(self, tmp_path):
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

        with pytest.raises(FileNotFoundError, match='Could not find index file'):
            get_index(args, 'Gorilla_gorilla_gorilla')

    def test_get_index_normalizes_redundant_underscores_without_shortening(self, tmp_path):
        out_dir = tmp_path / 'out'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        index_dir.mkdir()
        index_path = index_dir / 'Canis_lupus_familiaris.idx'
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

    def test_get_index_rebuilds_empty_index_file(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        (fasta_dir / 'Homo_sapiens.fa').write_text('>a\nAAAA\n')
        stale_index = index_dir / 'Homo_sapiens.idx'
        stale_index.touch()
        built = []

        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
            index_lock_poll=1,
            index_lock_timeout=10,
        )

        def fake_build(index_path, fasta_file, sci_name):
            built.append((index_path, fasta_file, sci_name))
            with open(index_path, 'w') as handle:
                handle.write('idx')

        monkeypatch.setattr('amalgkit.quant._build_kallisto_index', fake_build)

        observed = get_index(args, 'Homo_sapiens')

        assert observed == str(stale_index)
        assert built == [(str(stale_index), str(fasta_dir / 'Homo_sapiens.fa'), 'Homo_sapiens')]
        assert stale_index.read_text() == 'idx'

    def test_get_index_rebuilds_nonempty_index_without_ready_marker(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        fasta_dir = tmp_path / 'fasta'
        index_dir = tmp_path / 'index'
        out_dir.mkdir()
        fasta_dir.mkdir()
        index_dir.mkdir()
        (fasta_dir / 'Homo_sapiens.fa').write_text('>a\nAAAA\n', encoding='utf-8')
        partial_index = index_dir / 'Homo_sapiens.idx'
        partial_index.write_text('partial-but-nonempty', encoding='utf-8')
        built = []
        args = SimpleNamespace(
            out_dir=str(out_dir),
            index_dir=str(index_dir),
            build_index=True,
            fasta_dir=str(fasta_dir),
            index_lock_poll=1,
            index_lock_timeout=10,
        )

        def fake_build(index_path, fasta_file, sci_name):
            built.append((index_path, fasta_file, sci_name))
            with open(index_path, 'w', encoding='utf-8') as handle:
                handle.write('complete-index')
            _write_index_ready_marker(index_path, backend='kallisto')

        monkeypatch.setattr('amalgkit.quant._build_kallisto_index', fake_build)

        observed = get_index(args, 'Homo_sapiens')

        assert observed == str(partial_index)
        assert len(built) == 1
        assert partial_index.read_text(encoding='utf-8') == 'complete-index'

    def test_build_kallisto_index_uses_isolated_workdir(self, tmp_path, monkeypatch):
        index_dir = tmp_path / 'index'
        fasta_dir = tmp_path / 'fasta'
        index_dir.mkdir()
        fasta_dir.mkdir()
        index_path = index_dir / 'Homo_sapiens.idx'
        fasta_file = fasta_dir / 'Homo_sapiens.fa'
        fasta_file.write_text('>a\nAAAA\n')
        observed_cwds = []

        def fake_run(cmd, stdout, stderr, cwd=None):
            observed_cwds.append(cwd)
            assert cwd is not None
            assert os.path.isdir(cwd)
            assert os.path.dirname(cwd) == str(index_dir)
            assert cmd[:3] == ['kallisto', 'index', '-i']
            with open(cmd[3], 'w') as handle:
                handle.write('idx')
            return SimpleNamespace(returncode=0, stdout=b'', stderr=b'')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        _build_kallisto_index(str(index_path), str(fasta_file), 'Homo_sapiens')

        assert index_path.read_text() == 'idx'
        assert len(observed_cwds) == 1
        assert observed_cwds[0] != os.getcwd()
        assert not os.path.exists(observed_cwds[0])

    def test_build_kallisto_index_does_not_retry_internal_type_error(self, tmp_path, monkeypatch):
        index_path = tmp_path / 'Homo_sapiens.idx'
        fasta_file = tmp_path / 'Homo_sapiens.fa'
        fasta_file.write_text('>a\nAAAA\n', encoding='utf-8')
        calls = []

        def fake_run(cmd, stdout, stderr, cwd=None, timeout=None):
            calls.append((cwd, timeout))
            raise TypeError('timeout calculation failed inside runner')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        with pytest.raises(TypeError, match='timeout calculation failed'):
            _build_kallisto_index(
                str(index_path),
                str(fasta_file),
                'Homo_sapiens',
                args=SimpleNamespace(tool_timeout_seconds=30),
            )

        assert len(calls) == 1

    def test_failed_kallisto_index_build_preserves_previous_index(self, tmp_path, monkeypatch):
        index_path = tmp_path / 'Homo_sapiens.idx'
        fasta_file = tmp_path / 'Homo_sapiens.fa'
        index_path.write_text('previous-good-index', encoding='utf-8')
        fasta_file.write_text('>a\nAAAA\n', encoding='utf-8')

        def fake_run(cmd, stdout, stderr, cwd=None):
            _ = (stdout, stderr, cwd)
            with open(cmd[3], 'w', encoding='utf-8') as handle:
                handle.write('partial-index')
            return SimpleNamespace(returncode=1, stdout=b'', stderr=b'failed')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        with pytest.raises(RuntimeError, match='kallisto index failed'):
            _build_kallisto_index(str(index_path), str(fasta_file), 'Homo_sapiens')

        assert index_path.read_text(encoding='utf-8') == 'previous-good-index'

    def test_build_kallisto_index_rejects_symbolic_link_output(self, tmp_path):
        outside = tmp_path / 'outside.idx'
        outside.write_text('outside-index', encoding='utf-8')
        index_path = tmp_path / 'Homo_sapiens.idx'
        index_path.symlink_to(outside)
        fasta_file = tmp_path / 'Homo_sapiens.fa'
        fasta_file.write_text('>a\nAAAA\n', encoding='utf-8')

        with pytest.raises(ValueError, match='symbolic-link index output'):
            _build_kallisto_index(str(index_path), str(fasta_file), 'Homo_sapiens')

        assert outside.read_text(encoding='utf-8') == 'outside-index'

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

    def test_call_oarfish_writes_compatibility_outputs(self, tmp_path, monkeypatch):
        args = SimpleNamespace(threads=1)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'platform': ['OXFORD_NANOPORE'],
            'instrument': ['PromethION'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'spot_length': [1000],
            'total_bases': [10000],
            'exclusion': ['no'],
        }))
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'total_spot': 10,
        }

        def fake_run(cmd, stdout, stderr):
            prefix = tmp_path / 'SRR001'
            (prefix.with_suffix('.quant')).write_text(
                'tname\tlen\tnum_reads\n'
                'tx1\t1000\t4\n'
                'tx2\t500\t1\n'
            )
            (prefix.with_suffix('.meta_info.json')).write_text(json.dumps({'alignment_source': 'from_raw_reads'}))
            return SimpleNamespace(returncode=0, stdout=b'', stderr=b'')

        monkeypatch.setattr(subprocess, 'run', fake_run)

        observed = call_oarfish(
            args=args,
            in_files=['reads.fastq.gz'],
            metadata=metadata,
            sra_stat=sra_stat,
            output_dir=str(tmp_path),
            index='dummy.mmi',
            seq_tech='ont-cdna',
        )

        assert observed.returncode == 0
        abundance_df = pandas.read_csv(tmp_path / 'SRR001_abundance.tsv', sep='\t')
        assert list(abundance_df.columns) == ['target_id', 'length', 'eff_length', 'est_counts', 'tpm']
        assert abundance_df['target_id'].tolist() == ['tx1', 'tx2']
        with open(tmp_path / 'SRR001_run_info.json') as handle:
            run_info = json.load(handle)
        assert run_info['quant_backend'] == 'oarfish'
        assert run_info['oarfish_seq_tech'] == 'ont-cdna'
        assert run_info['p_pseudoaligned'] == 50.0

    def test_run_quant_auto_uses_oarfish_for_long_reads(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        args = SimpleNamespace(
            out_dir=str(out_dir),
            redo=False,
            clean_fastq=False,
            threads=1,
            quant_backend='auto',
            oarfish_seq_tech='auto',
        )
        runtime_context = QuantRuntimeContext(run_files_by_run={'SRR001': {'SRR001.fastq.gz'}})
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'platform': ['PACBIO_SMRT'],
            'instrument': ['Sequel IIe'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'spot_length': [2000],
            'total_bases': [20000],
            'nominal_length': [numpy.nan],
            'exclusion': ['no'],
        }))
        observed = {'kallisto': 0, 'oarfish': 0}

        monkeypatch.setattr(
            'amalgkit.quant.call_kallisto',
            lambda *_args, **_kwargs: observed.__setitem__('kallisto', observed['kallisto'] + 1),
        )

        def fake_call_oarfish(_args, in_files, _metadata, _sra_stat, _output_dir, _index, seq_tech):
            observed['oarfish'] += 1
            assert in_files == [str(out_dir / 'getfastq' / 'SRR001' / 'SRR001.fastq.gz')]
            assert seq_tech == 'pac-bio'
            write_valid_quant_outputs(_output_dir)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr('amalgkit.quant.call_oarfish', fake_call_oarfish)

        run_quant(args, metadata, 'SRR001', 'dummy.mmi', runtime_context=runtime_context)

        assert observed == {'kallisto': 0, 'oarfish': 1}

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
        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda args=None: None)
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
        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda args=None: None)
        monkeypatch.setattr(
            'amalgkit.quant.load_metadata',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('load_metadata should not be called')),
        )

        with pytest.raises(NotADirectoryError, match='Quant path exists but is not a directory'):
            quant_main(args)

    def test_quant_main_recovers_stats_before_auto_backend_inference(self, tmp_path, monkeypatch):
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
            internal_jobs=1,
            threads=1,
            redo=False,
            metadata='inferred',
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

        monkeypatch.setattr('amalgkit.quant.load_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.quant.check_quant_dependencies', lambda _backends, _args=None: None)
        monkeypatch.setattr(
            'amalgkit.quant.pre_resolve_species_indices',
            lambda _args, _tasks, runtime_context=None: {},
        )

        def fake_run_quant_for_sra(_args, quant_metadata, sra_id, _sci_name, runtime_context=None):
            observed['run'] = sra_id
            observed['backend'] = runtime_context.quant_backend_by_run[sra_id]
            observed['spot_length'] = quant_metadata.df.loc[0, 'spot_length']

        monkeypatch.setattr('amalgkit.quant.run_quant_for_sra', fake_run_quant_for_sra)

        quant_main(args)

        assert observed == {
            'run': 'SRR001',
            'backend': 'kallisto',
            'spot_length': 100.0,
        }
        assert pandas.isna(metadata.df.loc[0, 'total_spots'])
        assert pandas.isna(metadata.df.loc[0, 'total_bases'])
        assert pandas.isna(metadata.df.loc[0, 'spot_length'])

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
        observed = {}

        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda args=None: None)
        def fake_load_metadata(runtime_args):
            observed['load_threads'] = runtime_args.threads
            return metadata

        monkeypatch.setattr('amalgkit.quant.load_metadata', fake_load_metadata)
        monkeypatch.setattr('amalgkit.quant.pre_resolve_species_indices', lambda _args, _tasks, runtime_context=None: {})
        monkeypatch.setattr(
            'amalgkit.quant.run_quant_for_sra',
            lambda runtime_args, _metadata, sra_id, sci_name, runtime_context=None: (
                observed.setdefault('run_threads', runtime_args.threads),
                observed.setdefault('run_jobs', runtime_args.internal_jobs),
                dispatched.append((sra_id, sci_name)),
            )[-1],
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

        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda args=None: None)
        monkeypatch.setattr('amalgkit.quant.load_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.quant.pre_resolve_species_indices', lambda _args, _tasks, runtime_context=None: {})

        def fake_run_quant_for_sra(_args, _metadata, sra_id, _sci_name, runtime_context=None):
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
        observed = {}

        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda args=None: None)
        def fake_load_metadata(runtime_args):
            observed['load_threads'] = runtime_args.threads
            return metadata

        monkeypatch.setattr('amalgkit.quant.load_metadata', fake_load_metadata)
        monkeypatch.setattr('amalgkit.quant.pre_resolve_species_indices', lambda _args, _tasks, runtime_context=None: {})
        monkeypatch.setattr(
            'amalgkit.quant.run_quant_for_sra',
            lambda runtime_args, _metadata, sra_id, sci_name, runtime_context=None: (
                observed.setdefault('run_threads', runtime_args.threads),
                observed.setdefault('run_jobs', runtime_args.internal_jobs),
                dispatched.append((sra_id, sci_name)),
            )[-1],
        )

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('run_tasks_with_optional_threads should not be used when --internal_cpu_budget caps internal_jobs to 1.')

        monkeypatch.setattr('amalgkit.quant.run_tasks_with_optional_threads', fail_if_called)

        quant_main(args)

        assert set(dispatched) == {('SRR001', 'Species A'), ('SRR002', 'Species B')}
        assert observed['load_threads'] == 1
        assert observed['run_threads'] == 1
        assert observed['run_jobs'] == 1
        assert args.threads == 4
        assert args.internal_jobs == 4

    def test_quant_main_uses_runtime_copy_without_mutating_caller_namespace(self, tmp_path, monkeypatch):
        raw_out_dir = str(tmp_path / 'nested' / '..' / 'out')
        args = SimpleNamespace(
            out_dir=raw_out_dir,
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
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'exclusion': ['no'],
            'nominal_length': [200],
        }))
        observed = {}

        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda args=None: None)

        def fake_load_metadata(runtime_args):
            observed['load_out_dir'] = runtime_args.out_dir
            return metadata

        def fake_prefetch_getfastq_run_files(runtime_args, _tasks):
            observed['prefetch_out_dir'] = runtime_args.out_dir
            return {'SRR001': {'SRR001_1.fastq.gz'}}

        def fake_pre_resolve_species_indices(runtime_args, _tasks, runtime_context=None):
            observed['resolve_out_dir'] = runtime_args.out_dir
            observed['prefetched_context_run_files'] = runtime_context.run_files_by_run
            return {'Species_A': 'Species_A.idx'}

        def fake_run_quant_for_sra(runtime_args, _metadata, sra_id, sci_name, runtime_context=None):
            observed['runtime_has_run_files'] = runtime_context.run_files_by_run
            observed['runtime_has_index_cache'] = runtime_context.resolved_index_cache
            observed['runtime_is_caller'] = (runtime_args is args)
            observed['dispatched'] = (sra_id, sci_name)

        monkeypatch.setattr('amalgkit.quant.load_metadata', fake_load_metadata)
        monkeypatch.setattr('amalgkit.quant.prefetch_getfastq_run_files', fake_prefetch_getfastq_run_files)
        monkeypatch.setattr('amalgkit.quant.pre_resolve_species_indices', fake_pre_resolve_species_indices)
        monkeypatch.setattr('amalgkit.quant.run_quant_for_sra', fake_run_quant_for_sra)

        quant_main(args)

        normalized_out_dir = os.path.realpath(raw_out_dir)
        assert observed['load_out_dir'] == normalized_out_dir
        assert observed['prefetch_out_dir'] == normalized_out_dir
        assert observed['resolve_out_dir'] == normalized_out_dir
        assert observed['prefetched_context_run_files'] == {'SRR001': {'SRR001_1.fastq.gz'}}
        assert observed['runtime_has_run_files'] == {'SRR001': {'SRR001_1.fastq.gz'}}
        assert observed['runtime_has_index_cache'] == {'Species_A': 'Species_A.idx'}
        assert observed['runtime_is_caller'] is False
        assert observed['dispatched'] == ('SRR001', 'Species A')
        assert args.out_dir == raw_out_dir

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

        monkeypatch.setattr('amalgkit.quant.check_kallisto_dependency', lambda args=None: None)
        monkeypatch.setattr('amalgkit.quant.load_metadata', lambda _args: metadata)
        monkeypatch.setattr(
            'amalgkit.quant.get_index',
            lambda _args, sci_name, runtime_context=None: resolved_species.append(sci_name) or '{}.idx'.format(sci_name),
        )
        monkeypatch.setattr(
            'amalgkit.quant.run_quant_for_sra',
            lambda _args, _metadata, sra_id, sci_name, runtime_context=None: dispatched.append((sra_id, sci_name)),
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
        runtime_context = QuantRuntimeContext()

        monkeypatch.setattr(
            'amalgkit.quant.get_index',
            lambda _args, sci_name, runtime_context=None: resolved.append(sci_name) or sci_name + '.idx',
        )
        out = pre_resolve_species_indices(args, tasks, runtime_context=runtime_context)

        assert out == {'Species_A': 'Species_A.idx', 'Species_B': 'Species_B.idx'}
        assert resolved == ['Species_A', 'Species_B']
        assert runtime_context.prefetched_fasta.resolve_entries(str(fasta_dir)) == ['Species_A.fa', 'Species_B.fa']
        assert runtime_context.prefetched_index.resolve_entries(str(index_dir)) == ['Species_A.idx']

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

        monkeypatch.setattr(
            'amalgkit.quant.get_index',
            lambda _args, sci_name, runtime_context=None: resolved.append(sci_name) or sci_name + '.idx',
        )

        out = pre_resolve_species_indices(args, tasks)

        assert resolved == ['Species_A']
        assert out == {'Species_A': 'Species_A.idx'}

    def test_pre_resolve_species_indices_ignores_runtime_alias_values(self, tmp_path, monkeypatch):
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
        tasks = [
            ('SRR001', 'unclassified eukaryotes'),
            ('SRR002', 'unclassified eukaryotes'),
        ]
        runtime_context = QuantRuntimeContext(
            species_identifier_values_by_run={
                'SRR001': ['unclassified eukaryotes', 'Pygsuia biforma'],
                'SRR002': ['unclassified eukaryotes', 'Thecamonas trahens'],
            }
        )
        resolved = []

        def fake_get_index(_args, sci_name, runtime_context=None):
            _ = runtime_context
            resolved.append(sci_name)
            return '{}.idx'.format(sci_name)

        monkeypatch.setattr('amalgkit.quant.get_index', fake_get_index)

        out = pre_resolve_species_indices(args, tasks, runtime_context=runtime_context)

        assert resolved == ['unclassified_eukaryotes']
        assert out == {'unclassified_eukaryotes': 'unclassified_eukaryotes.idx'}

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
        monkeypatch.setattr(
            'amalgkit.quant.get_index',
            lambda _args, sci_name, runtime_context=None: sci_name + '.idx',
        )

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
