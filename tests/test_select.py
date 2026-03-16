import os
import subprocess
import sys
import pandas
import pytest

from types import SimpleNamespace
from pathlib import Path

from amalgkit.select import (
    apply_select_filters,
    write_select_outputs,
    resolve_select_config_dir,
    filter_metadata_by_sample_group,
    select_main,
)
from amalgkit.util import Metadata


# ---------------------------------------------------------------------------
# write_select_outputs (writes metadata and pivot tables)
# ---------------------------------------------------------------------------

class TestWriteSelectOutputs:
    def test_writes_metadata_and_pivots(self, tmp_path, sample_metadata):
        """Writes metadata.tsv, pivot_qualified.tsv, and pivot_selected.tsv."""
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = str(tmp_path / 'metadata' / 'metadata_original.tsv')
        path_table = str(tmp_path / 'metadata' / 'metadata.tsv')
        # Write a source metadata table to copy
        sample_metadata.df.to_csv(path_table, sep='\t', index=False)
        sample_metadata.label_sampled_data(max_sample=10)
        write_select_outputs(path_original, path_table, str(metadata_dir), sample_metadata)
        assert os.path.exists(path_original)
        assert os.path.exists(path_table)
        assert os.path.exists(str(metadata_dir / 'pivot_qualified.tsv'))
        assert os.path.exists(str(metadata_dir / 'pivot_selected.tsv'))

    def test_refreshes_original_from_current_input(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = str(tmp_path / 'metadata' / 'metadata_original.tsv')
        path_table = str(tmp_path / 'metadata' / 'metadata.tsv')
        original_df = sample_metadata.df.copy(deep=True)
        with open(path_original, 'w') as f:
            f.write('marker_content\n')
        sample_metadata.df.to_csv(path_table, sep='\t', index=False)
        sample_metadata.label_sampled_data(max_sample=10)
        write_select_outputs(
            path_original,
            path_table,
            str(metadata_dir),
            sample_metadata,
            metadata_original_df=original_df,
        )
        loaded = pandas.read_csv(path_original, sep='\t')
        assert loaded.shape[0] == original_df.shape[0]
        assert set(loaded['run'].tolist()) == set(original_df['run'].tolist())

    def test_writes_original_without_existing_metadata_table(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        path_original = str(tmp_path / 'metadata' / 'metadata_original.tsv')
        path_table = str(tmp_path / 'metadata' / 'metadata.tsv')
        sample_metadata.label_sampled_data(max_sample=10)
        original_df = sample_metadata.df.copy(deep=True)

        write_select_outputs(
            path_metadata_original=path_original,
            path_metadata_table=path_table,
            metadata_dir=str(metadata_dir),
            metadata=sample_metadata,
            metadata_original_df=original_df,
        )

        assert os.path.exists(path_original)
        assert os.path.exists(path_table)
        original_loaded = pandas.read_csv(path_original, sep='\t')
        assert original_loaded.shape[0] == original_df.shape[0]
        assert set(original_loaded['run'].tolist()) == set(original_df['run'].tolist())

    def test_uses_current_input_for_original_even_when_metadata_table_exists(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = str(metadata_dir / 'metadata_original.tsv')
        path_table = str(metadata_dir / 'metadata.tsv')
        existing_filtered = sample_metadata.df.iloc[[0], :].copy()
        existing_filtered.to_csv(path_table, sep='\t', index=False)
        original_df = sample_metadata.df.copy(deep=True)
        sample_metadata.label_sampled_data(max_sample=10)

        write_select_outputs(
            path_metadata_original=path_original,
            path_metadata_table=path_table,
            metadata_dir=str(metadata_dir),
            metadata=sample_metadata,
            metadata_original_df=original_df,
        )

        original_loaded = pandas.read_csv(path_original, sep='\t')
        assert original_loaded.shape[0] == original_df.shape[0]
        assert set(original_loaded['run'].tolist()) == set(original_df['run'].tolist())
        assert original_loaded.shape[0] != existing_filtered.shape[0]

    def test_rejects_metadata_original_directory_path(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = metadata_dir / 'metadata_original.tsv'
        path_original.mkdir()
        path_table = str(metadata_dir / 'metadata.tsv')
        sample_metadata.label_sampled_data(max_sample=10)

        with pytest.raises(NotADirectoryError, match='not a file'):
            write_select_outputs(str(path_original), path_table, str(metadata_dir), sample_metadata)

    def test_rejects_metadata_table_directory_path(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = str(metadata_dir / 'metadata_original.tsv')
        path_table = metadata_dir / 'metadata.tsv'
        path_table.mkdir()
        sample_metadata.label_sampled_data(max_sample=10)

        with pytest.raises(NotADirectoryError, match='not a file'):
            write_select_outputs(path_original, str(path_table), str(metadata_dir), sample_metadata)

    def test_rejects_metadata_dir_file_path(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.write_text('not a directory')
        path_original = str(tmp_path / 'metadata_original.tsv')
        path_table = str(tmp_path / 'metadata.tsv')
        sample_metadata.label_sampled_data(max_sample=10)

        with pytest.raises(NotADirectoryError, match='not a directory'):
            write_select_outputs(path_original, path_table, str(metadata_dir), sample_metadata)


class TestSelectHelpers:
    def test_resolve_select_config_dir_inferred(self):
        args = SimpleNamespace(out_dir='/tmp/out', config_dir='inferred')
        assert resolve_select_config_dir(args) == '/tmp/out/config'

    def test_resolve_select_config_dir_explicit(self):
        args = SimpleNamespace(out_dir='/tmp/out', config_dir='/tmp/custom')
        assert resolve_select_config_dir(args) == '/tmp/custom'

    def test_filter_metadata_by_sample_group(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'sample_group': ['brain', 'liver', 'brain'],
            'scientific_name': ['sp1', 'sp1', 'sp2'],
            'exclusion': ['no', 'no', 'no'],
        }))
        out = filter_metadata_by_sample_group(metadata, 'brain')
        assert set(out.df['run']) == {'SRR001', 'SRR003'}

    def test_filter_metadata_by_sample_group_strips_tokens(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'sample_group': ['brain', 'liver', 'heart'],
            'scientific_name': ['sp1', 'sp1', 'sp2'],
            'exclusion': ['no', 'no', 'no'],
        }))
        out = filter_metadata_by_sample_group(metadata, 'brain, liver')
        assert set(out.df['run']) == {'SRR001', 'SRR002'}

    def test_filter_metadata_by_sample_group_supports_pipe_separator(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'sample_group': ['brain', 'liver', 'heart'],
            'scientific_name': ['sp1', 'sp1', 'sp2'],
            'exclusion': ['no', 'no', 'no'],
        }))
        out = filter_metadata_by_sample_group(metadata, 'brain|heart')
        assert set(out.df['run']) == {'SRR001', 'SRR003'}

    def test_filter_metadata_by_sample_group_strips_metadata_values(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'sample_group': [' brain ', 'liver'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        out = filter_metadata_by_sample_group(metadata, 'brain')
        assert set(out.df['run']) == {'SRR001'}

    def test_filter_metadata_by_sample_group_none_keeps_all(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'sample_group': ['brain', 'liver'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        out = filter_metadata_by_sample_group(metadata, None)
        assert out.df.shape[0] == 2

    def test_filter_metadata_by_sample_group_empty_argument_raises(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'sample_group': ['brain', 'liver'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        with pytest.raises(ValueError, match='No sample_group was selected'):
            filter_metadata_by_sample_group(metadata, '  ')

    def test_filter_metadata_by_sample_group_requires_column(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        metadata.df = metadata.df.drop(columns=['sample_group'])
        with pytest.raises(ValueError, match='Column \"sample_group\" is required'):
            filter_metadata_by_sample_group(metadata, 'brain')


class TestApplySelectFilters:
    def test_preserves_grouped_source_columns_until_keyword_filters_finish(self, tmp_path):
        (tmp_path / 'group_attribute.config').write_text('sample_group\tsample_attribute\n')
        (tmp_path / 'exclude_keyword.config').write_text('sample_attribute\tbad_tissue\tcancer\n')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'sample_group': [''],
            'sample_attribute': ['cancer'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
        }))
        args = SimpleNamespace(
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=10,
        )

        out = apply_select_filters(metadata, args, str(tmp_path))

        assert out.df.loc[0, 'exclusion'] == 'bad_tissue'
        assert out.df.loc[0, 'sample_group'] == 'cancer[sample_attribute]'
        assert 'sample_attribute' not in out.df.columns

    def test_preserves_taxid_and_extra_columns_after_filtering(self, tmp_path):
        (tmp_path / 'group_attribute.config').write_text('')
        (tmp_path / 'exclude_keyword.config').write_text('')
        (tmp_path / 'control_term.config').write_text('')
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'sample_group': ['leaf'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
            'taxid_species': [12345],
            'sample_attribute_tissue': ['leaf'],
            'sample_group_normalization_source': ['sample_attribute_tissue'],
        }))
        args = SimpleNamespace(
            min_nspots=0,
            mark_missing_rank='species',
            mark_redundant_biosamples=False,
            max_sample=10,
        )

        out = apply_select_filters(metadata, args, str(tmp_path))

        assert 'taxid_species' in out.df.columns
        assert out.df.loc[0, 'taxid_species'] == 12345
        assert 'sample_attribute_tissue' in out.df.columns
        assert out.df.loc[0, 'sample_attribute_tissue'] == 'leaf'
        assert 'sample_group_normalization_source' in out.df.columns
        assert out.df.loc[0, 'sample_group_normalization_source'] == 'sample_attribute_tissue'


class TestSelectMain:
    def test_rejects_out_dir_file_path(self, tmp_path, monkeypatch):
        out_path = tmp_path / 'out_path'
        out_path.write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_path),
            config_dir='inferred',
            metadata='inferred',
            sample_group=None,
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=1,
        )
        monkeypatch.setattr(
            'amalgkit.select.check_config_dir',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('check_config_dir should not be called')),
        )

        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            select_main(args)

    def test_uses_runtime_copy_without_mutating_caller_args(self, tmp_path, monkeypatch):
        raw_out_dir = str(tmp_path / 'nested' / '..' / 'out')
        args = SimpleNamespace(
            out_dir=raw_out_dir,
            config_dir='inferred',
            metadata='inferred',
            sample_group=None,
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=1,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'sample_group': ['brain'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
        }))
        observed = {}

        def fake_check_config_dir(dir_path, mode):
            observed['config_dir'] = dir_path
            observed['config_mode'] = mode

        def fake_load_metadata(runtime_args):
            observed['load_out_dir'] = runtime_args.out_dir
            return metadata

        def fake_apply_select_filters(current_metadata, runtime_args, dir_config):
            observed['filter_out_dir'] = runtime_args.out_dir
            observed['filter_config_dir'] = dir_config
            return current_metadata

        def fake_write_select_outputs(**kwargs):
            observed['write_table'] = kwargs['path_metadata_table']

        monkeypatch.setattr('amalgkit.select.check_config_dir', fake_check_config_dir)
        monkeypatch.setattr('amalgkit.select.load_metadata', fake_load_metadata)
        monkeypatch.setattr('amalgkit.select.apply_select_filters', fake_apply_select_filters)
        monkeypatch.setattr('amalgkit.select.write_select_outputs', fake_write_select_outputs)

        select_main(args)

        normalized_out_dir = os.path.realpath(raw_out_dir)
        assert observed['config_dir'] == os.path.join(normalized_out_dir, 'config')
        assert observed['config_mode'] == 'select'
        assert observed['load_out_dir'] == normalized_out_dir
        assert observed['filter_out_dir'] == normalized_out_dir
        assert observed['filter_config_dir'] == os.path.join(normalized_out_dir, 'config')
        assert observed['write_table'] == os.path.join(normalized_out_dir, 'metadata', 'metadata.tsv')
        assert args.out_dir == raw_out_dir


class TestSelectBatchMain:
    def test_batch_mode_writes_summary_queue_and_manifest(self, tmp_path):
        metadata_specieswise_dir = tmp_path / 'metadata_specieswise'
        metadata_specieswise_dir.mkdir()
        species_tsv = tmp_path / 'species.tsv'
        config_dir = tmp_path / 'config'
        config_dir.mkdir()
        for filename in ['group_attribute.config', 'exclude_keyword.config', 'control_term.config']:
            (config_dir / filename).write_text('')

        species_rows = pandas.DataFrame({
            'scientific_name': ['Species alpha', 'Species beta'],
            'species_token': ['Species_alpha', 'Species_beta'],
        })
        species_rows.to_csv(species_tsv, sep='\t', index=False)

        def write_merged_metadata(species_token, scientific_name, organ_counts):
            species_dir = metadata_specieswise_dir / species_token
            species_dir.mkdir()
            rows = []
            counter = 1
            for organ, count in organ_counts.items():
                for _ in range(count):
                    rows.append({
                        'run': 'SRR{:04d}'.format(counter),
                        'scientific_name': scientific_name,
                        'sample_group': '',
                        'sample_attribute_tissue': organ,
                        'bioproject': 'PRJ{}'.format(counter),
                        'biosample': 'SAM{}'.format(counter),
                        'total_spots': 100,
                        'exclusion': 'no',
                    })
                    counter += 1
            pandas.DataFrame(rows).to_csv(
                species_dir / '{}.metadata.tsv'.format(species_token),
                sep='\t',
                index=False,
            )

        write_merged_metadata(
            species_token='Species_alpha',
            scientific_name='Species alpha',
            organ_counts={'flower': 2, 'leaf': 2, 'root': 2},
        )
        write_merged_metadata(
            species_token='Species_beta',
            scientific_name='Species beta',
            organ_counts={'flower': 1, 'leaf': 2, 'root': 2},
        )

        out_dir = tmp_path / 'select_batch'
        args = SimpleNamespace(
            out_dir=str(out_dir),
            config_dir=str(config_dir),
            metadata='inferred',
            sample_group='flower,leaf,root',
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=2,
            species_tsv=str(species_tsv),
            metadata_specieswise_dir=str(metadata_specieswise_dir),
            summary_tsv='inferred',
            queue_tsv='inferred',
            manifest_tsv='inferred',
            batch_label='inferred',
        )

        select_main(args)

        normalization_summary = tmp_path / 'select_batch' / 'normalization_summary.tsv'
        select_summary = tmp_path / 'select_batch' / 'select_summary.tsv'
        select_queue = tmp_path / 'select_batch' / 'select_queue.tsv'
        manifest = tmp_path / 'select_batch' / 'external_manifest.tsv'
        manifest_strict = tmp_path / 'select_batch' / 'external_manifest_strict.tsv'
        manifest_relaxed = tmp_path / 'select_batch' / 'external_manifest_relaxed.tsv'

        assert normalization_summary.exists()
        assert select_summary.exists()
        assert select_queue.exists()
        assert manifest.exists()
        assert manifest_strict.exists()
        assert manifest_relaxed.exists()

        summary_df = pandas.read_csv(select_summary, sep='\t')
        assert set(summary_df['species_token'].tolist()) == {'Species_alpha', 'Species_beta'}
        queue_by_species = summary_df.set_index('species_token')['queue_tier'].to_dict()
        assert queue_by_species['Species_alpha'] == 'strict'
        assert queue_by_species['Species_beta'] == 'defer'

        manifest_df = pandas.read_csv(manifest, sep='\t')
        assert set(manifest_df['queue_tier'].tolist()) == {'strict', 'defer'}
        assert 'selected_metadata_path' in manifest_df.columns


class TestCliEntry:
    def test_python_module_entrypoint_supports_help(self):
        repo_root = Path(__file__).resolve().parents[1]
        completed = subprocess.run(
            [sys.executable, '-m', 'amalgkit', 'select', '--help'],
            cwd=str(repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert completed.returncode == 0
        assert 'usage: amalgkit select' in completed.stdout
