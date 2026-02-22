import os
import pandas
import pytest

from types import SimpleNamespace

from amalgkit.select import (
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

    def test_does_not_overwrite_original(self, tmp_path, sample_metadata):
        """If metadata_original.tsv already exists, it is not overwritten."""
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = str(tmp_path / 'metadata' / 'metadata_original.tsv')
        path_table = str(tmp_path / 'metadata' / 'metadata.tsv')
        # Pre-create original with marker content
        with open(path_original, 'w') as f:
            f.write('marker_content\n')
        sample_metadata.df.to_csv(path_table, sep='\t', index=False)
        sample_metadata.label_sampled_data(max_sample=10)
        write_select_outputs(path_original, path_table, str(metadata_dir), sample_metadata)
        with open(path_original) as f:
            content = f.read()
        assert 'marker_content' in content

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
