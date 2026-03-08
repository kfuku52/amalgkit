import os
import sys
import pytest
import pandas
from types import SimpleNamespace

from amalgkit.command_context import PerSpeciesTableContext
from amalgkit.per_species_tables import (
    get_sample_group,
    resolve_per_species_input,
    collect_pending_species_for_tables,
    write_completion_flag,
    generate_per_species_tables,
    list_selected_species,
)
from amalgkit.util import Metadata


# ---------------------------------------------------------------------------
# get_sample_group (wiki: curate extracts sample groups from args or metadata)
# ---------------------------------------------------------------------------

class TestGetSampleGroup:
    def test_from_args(self):
        """When --sample_group is specified, parse it."""
        class Args:
            sample_group = 'brain,liver,heart'
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'], 'exclusion': ['no'],
            'sample_group': ['brain'],
        }))
        result = get_sample_group(Args(), m)
        assert result == 'brain|liver|heart'

    def test_from_args_keeps_hyphen_and_trims(self):
        class Args:
            sample_group = 'non-treated, treated '
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'], 'exclusion': ['no'],
            'sample_group': ['non-treated'],
        }))
        result = get_sample_group(Args(), m)
        assert result == 'non-treated|treated'

    def test_from_metadata(self):
        """When --sample_group is None, extract from metadata."""
        class Args:
            sample_group = None
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2', 'R3'],
            'exclusion': ['no', 'no', 'no'],
            'sample_group': ['brain', 'liver', 'brain'],
        }))
        result = get_sample_group(Args(), m)
        # Should contain brain and liver separated by pipe
        groups = result.split('|')
        assert 'brain' in groups
        assert 'liver' in groups

    def test_from_metadata_trims_and_deduplicates(self):
        class Args:
            sample_group = None
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2', 'R3', 'R4'],
            'exclusion': ['no', 'no', 'no', 'no'],
            'sample_group': [' brain ', 'brain', '', float('nan')],
        }))
        result = get_sample_group(Args(), m)
        assert result == 'brain'

    def test_empty_sample_group_exits(self):
        """Wiki: curate exits with error if sample_group is empty."""
        class Args:
            sample_group = None
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['no'],
            'sample_group': [float('nan')],
        }))
        with pytest.raises(ValueError, match='sample_group'):
            get_sample_group(Args(), m)

    def test_missing_sample_group_column_exits_with_clear_error(self, capsys):
        class Args:
            sample_group = None
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['no'],
        }))
        with pytest.raises(ValueError, match='sample_group'):
            get_sample_group(Args(), m)


class TestResolveCurateInput:
    @staticmethod
    def _args(out_dir, input_dir):
        class Args:
            pass
        args = Args()
        args.out_dir = str(out_dir)
        args.input_dir = str(input_dir)
        args.metadata = 'inferred'
        return args

    def test_explicit_input_dir_with_inferred_metadata_reads_input_metadata(self, tmp_path):
        out_dir = tmp_path / 'out'
        input_dir = tmp_path / 'custom_input'
        input_dir.mkdir(parents=True)
        pandas.DataFrame({
            'run': ['R1'],
            'scientific_name': ['Species A'],
            'sample_group': ['g1'],
            'exclusion': ['no'],
        }).to_csv(input_dir / 'metadata.tsv', sep='\t', index=False)

        args = self._args(out_dir=out_dir, input_dir=input_dir)
        metadata, resolved_input_dir = resolve_per_species_input(args)

        assert resolved_input_dir == str(input_dir)
        assert metadata.df['run'].tolist() == ['R1']

    def test_explicit_input_dir_with_missing_inferred_metadata_raises(self, tmp_path):
        out_dir = tmp_path / 'out'
        input_dir = tmp_path / 'custom_input'
        input_dir.mkdir(parents=True)
        args = self._args(out_dir=out_dir, input_dir=input_dir)

        with pytest.raises(FileNotFoundError, match='metadata.tsv not found in --input_dir'):
            resolve_per_species_input(args)

    def test_explicit_input_dir_file_path_raises(self, tmp_path):
        out_dir = tmp_path / 'out'
        input_path = tmp_path / 'custom_input'
        input_path.write_text('not a directory')
        args = self._args(out_dir=out_dir, input_dir=input_path)

        with pytest.raises(NotADirectoryError, match='Input path exists but is not a directory'):
            resolve_per_species_input(args)

    def test_inferred_cstmm_path_file_raises(self, tmp_path):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        (out_dir / 'cstmm').write_text('not a directory')
        args = self._args(out_dir=out_dir, input_dir='inferred')

        with pytest.raises(NotADirectoryError, match='cstmm input path exists but is not a directory'):
            resolve_per_species_input(args)


class TestCollectPendingSpeciesForCurate:
    def test_redo_replaces_species_symlink(self, tmp_path):
        curate_dir = tmp_path / 'per_species'
        curate_dir.mkdir()
        species = 'Species_A'
        real_species_dir = tmp_path / 'real_species_dir'
        real_species_dir.mkdir()
        (real_species_dir / 'keep.txt').write_text('keep')
        os.symlink(real_species_dir, curate_dir / species)
        write_completion_flag(str(curate_dir), species)

        args = SimpleNamespace(redo=True)
        pending = collect_pending_species_for_tables(args, str(curate_dir), [species])

        assert pending == [species]
        assert not os.path.lexists(str(curate_dir / species))
        assert os.path.exists(str(real_species_dir / 'keep.txt'))


class TestCurateMain:
    @staticmethod
    def _args(out_dir):
        class Args:
            pass
        args = Args()
        args.input_dir = 'inferred'
        args.out_dir = str(out_dir)
        args.norm = 'log2p1-fpkm'
        args.redo = False
        args.sample_group = None
        args.metadata = 'metadata.tsv'
        args.dist_method = 'pearson'
        args.mapping_rate = 0.2
        args.correlation_threshold = 0.3
        args.plot_intermediate = False
        args.sample_group_color = 'DEFAULT'
        args.one_outlier_per_iter = False
        args.batch_effect_alg = 'no'
        args.clip_negative = True
        args.maintain_zero = True
        args.skip_curation = False
        args.internal_jobs = 1
        return args

    def test_creates_nested_out_dir_and_completion_flags(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A', 'Species B'],
            'run': ['R1', 'R2'],
            'sample_group': ['g1', 'g2'],
            'exclusion': ['no', 'no'],
        }))

        monkeypatch.setattr('amalgkit.per_species_tables.load_metadata', lambda *_args, **_kwargs: metadata)

        def fake_run_curate(args, metadata, sp, input_dir):
            os.makedirs(os.path.join(args.out_dir, 'per_species', sp), exist_ok=True)
            return 0

        monkeypatch.setattr('amalgkit.per_species_tables.run_per_species_job', fake_run_curate)
        generate_per_species_tables(args)
        assert (out_dir / 'per_species' / 'Species_A' / 'per_species_completion_flag.txt').exists()
        assert (out_dir / 'per_species' / 'Species_B' / 'per_species_completion_flag.txt').exists()

    def test_uses_explicit_context_without_resolving_input(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        input_dir = tmp_path / 'custom_input'
        args = self._args(out_dir)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A'],
            'run': ['R1'],
            'sample_group': ['g1'],
            'exclusion': ['no'],
        }))

        monkeypatch.setattr(
            'amalgkit.per_species_tables.resolve_per_species_input',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('resolve_per_species_input should not be called')),
        )

        def fake_run_curate(args, metadata, sp, input_dir):
            os.makedirs(os.path.join(args.out_dir, 'per_species', sp), exist_ok=True)
            return 0

        monkeypatch.setattr('amalgkit.per_species_tables.run_per_species_job', fake_run_curate)

        generate_per_species_tables(
            args,
            context=PerSpeciesTableContext(metadata=metadata, input_dir=str(input_dir)),
        )

        assert (out_dir / 'per_species' / 'Species_A' / 'per_species_completion_flag.txt').exists()

    def test_exits_nonzero_when_any_species_fails(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A', 'Species B'],
            'run': ['R1', 'R2'],
            'sample_group': ['g1', 'g2'],
            'exclusion': ['no', 'no'],
        }))

        monkeypatch.setattr('amalgkit.per_species_tables.load_metadata', lambda *_args, **_kwargs: metadata)

        def fake_run_curate(args, metadata, sp, input_dir):
            os.makedirs(os.path.join(args.out_dir, 'per_species', sp), exist_ok=True)
            return 1 if sp == 'Species_A' else 0

        monkeypatch.setattr('amalgkit.per_species_tables.run_per_species_job', fake_run_curate)
        with pytest.raises(RuntimeError, match='Per-species table generation failed for 1/2 species'):
            generate_per_species_tables(args)
        assert not (out_dir / 'per_species' / 'Species_A' / 'per_species_completion_flag.txt').exists()
        assert (out_dir / 'per_species' / 'Species_B' / 'per_species_completion_flag.txt').exists()

    def test_raises_for_tpm_with_cstmm_counts_input(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        args.norm = 'tpm'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A'],
            'run': ['R1'],
            'sample_group': ['g1'],
            'exclusion': ['no'],
        }))
        input_dir = tmp_path / 'input_custom'
        (input_dir / 'Species_A').mkdir(parents=True)
        (input_dir / 'Species_A' / 'Species_A_cstmm_counts.tsv').write_text('target_id\tR1\nG1\t1\n')

        monkeypatch.setattr('amalgkit.per_species_tables.resolve_per_species_input', lambda _args: (metadata, str(input_dir)))
        monkeypatch.setattr(
            'amalgkit.per_species_tables.run_per_species_job',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('run_per_species_job should not be called')),
        )

        with pytest.raises(ValueError, match='TPM and TMM are incompatible'):
            generate_per_species_tables(args)

    def test_tpm_does_not_false_positive_on_input_dir_name(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        args.norm = 'tpm'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A'],
            'run': ['R1'],
            'sample_group': ['g1'],
            'exclusion': ['no'],
        }))
        input_dir = tmp_path / 'custom_path_with_cstmm_in_name'
        (input_dir / 'Species_A').mkdir(parents=True)
        (input_dir / 'Species_A' / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')

        monkeypatch.setattr('amalgkit.per_species_tables.resolve_per_species_input', lambda _args: (metadata, str(input_dir)))

        def fake_run_curate(args, metadata, sp, input_dir):
            os.makedirs(os.path.join(args.out_dir, 'per_species', sp), exist_ok=True)
            return 0

        monkeypatch.setattr('amalgkit.per_species_tables.run_per_species_job', fake_run_curate)

        generate_per_species_tables(args)
        assert (out_dir / 'per_species' / 'Species_A' / 'per_species_completion_flag.txt').exists()

    def test_tpm_does_not_false_positive_on_directory_named_cstmm_counts_file(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        args.norm = 'tpm'
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A'],
            'run': ['R1'],
            'sample_group': ['g1'],
            'exclusion': ['no'],
        }))
        input_dir = tmp_path / 'input_custom'
        (input_dir / 'Species_A').mkdir(parents=True)
        (input_dir / 'Species_A' / 'Species_A_cstmm_counts.tsv').mkdir()
        (input_dir / 'Species_A' / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')

        monkeypatch.setattr('amalgkit.per_species_tables.resolve_per_species_input', lambda _args: (metadata, str(input_dir)))

        def fake_run_curate(args, metadata, sp, input_dir):
            os.makedirs(os.path.join(args.out_dir, 'per_species', sp), exist_ok=True)
            return 0

        monkeypatch.setattr('amalgkit.per_species_tables.run_per_species_job', fake_run_curate)

        generate_per_species_tables(args)
        assert (out_dir / 'per_species' / 'Species_A' / 'per_species_completion_flag.txt').exists()

    def test_parallel_species_jobs_creates_completion_flags(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        args.internal_jobs = 2
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A', 'Species B'],
            'run': ['R1', 'R2'],
            'sample_group': ['g1', 'g2'],
            'exclusion': ['no', 'no'],
        }))

        monkeypatch.setattr('amalgkit.per_species_tables.load_metadata', lambda *_args, **_kwargs: metadata)

        def fake_run_curate(args, metadata, sp, input_dir):
            os.makedirs(os.path.join(args.out_dir, 'per_species', sp), exist_ok=True)
            return 0

        monkeypatch.setattr('amalgkit.per_species_tables.run_per_species_job', fake_run_curate)
        generate_per_species_tables(args)

        assert (out_dir / 'per_species' / 'Species_A' / 'per_species_completion_flag.txt').exists()
        assert (out_dir / 'per_species' / 'Species_B' / 'per_species_completion_flag.txt').exists()

    def test_cpu_budget_caps_species_jobs_to_serial(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        args.internal_jobs = 4
        args.internal_cpu_budget = 1
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A', 'Species B'],
            'run': ['R1', 'R2'],
            'sample_group': ['g1', 'g2'],
            'exclusion': ['no', 'no'],
        }))
        processed = []

        monkeypatch.setattr('amalgkit.per_species_tables.load_metadata', lambda *_args, **_kwargs: metadata)

        def fake_run_curate(args, metadata, sp, input_dir):
            processed.append(sp)
            os.makedirs(os.path.join(args.out_dir, 'per_species', sp), exist_ok=True)
            return 0

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError('run_tasks_with_optional_threads should not be used when --internal_cpu_budget caps internal_jobs to 1.')

        monkeypatch.setattr('amalgkit.per_species_tables.run_per_species_job', fake_run_curate)
        monkeypatch.setattr('amalgkit.per_species_tables.run_tasks_with_optional_threads', fail_if_called)

        generate_per_species_tables(args)

        assert set(processed) == {'Species_A', 'Species_B'}

    def test_rejects_nonpositive_species_jobs(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'nested' / 'output'
        args = self._args(out_dir)
        args.internal_jobs = 0
        with pytest.raises(ValueError, match='--internal_jobs must be > 0'):
            generate_per_species_tables(args)

    def test_rejects_out_dir_file_path(self, tmp_path, monkeypatch):
        out_path = tmp_path / 'out_path'
        out_path.write_text('not a directory')
        args = self._args(out_path)
        monkeypatch.setattr(
            'amalgkit.per_species_tables.resolve_per_species_input',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('resolve_per_species_input should not be called')),
        )

        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            generate_per_species_tables(args)

    def test_rejects_curate_path_file(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        (out_dir / 'per_species').write_text('not a directory')
        args = self._args(out_dir)
        monkeypatch.setattr(
            'amalgkit.per_species_tables.resolve_per_species_input',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('resolve_per_species_input should not be called')),
        )

        with pytest.raises(NotADirectoryError, match='Per-species path exists but is not a directory'):
            generate_per_species_tables(args)

    def test_list_selected_species_ignores_missing_scientific_name(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A', float('nan')],
            'run': ['R1', 'R2'],
            'sample_group': ['g1', 'g2'],
            'exclusion': ['no', 'no'],
        }))
        assert list_selected_species(metadata) == ['Species_A']

    def test_list_selected_species_ignores_blank_scientific_name(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A', '', '   '],
            'run': ['R1', 'R2', 'R3'],
            'sample_group': ['g1', 'g2', 'g3'],
            'exclusion': ['no', 'no', 'no'],
        }))
        assert list_selected_species(metadata) == ['Species_A']

    def test_list_selected_species_normalizes_exclusion_flags(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A', 'Species B', 'Species C'],
            'run': ['R1', 'R2', 'R3'],
            'sample_group': ['g1', 'g2', 'g3'],
            'exclusion': [' NO ', 'No', 'yes'],
        }))
        assert list_selected_species(metadata) == ['Species_A', 'Species_B']

    def test_list_selected_species_rejects_missing_exclusion_column(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Species A'],
            'run': ['R1'],
            'sample_group': ['g1'],
        }))
        metadata.df = metadata.df.drop(columns=['exclusion'])

        with pytest.raises(ValueError, match='Missing required metadata column\\(s\\) for per-species table generation: exclusion'):
            list_selected_species(metadata)

    def test_list_selected_species_rejects_missing_scientific_name_column(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'sample_group': ['g1'],
            'exclusion': ['no'],
        }))
        metadata.df = metadata.df.drop(columns=['scientific_name'])

        with pytest.raises(ValueError, match='Missing required metadata column\\(s\\) for per-species table generation: scientific_name'):
            list_selected_species(metadata)
