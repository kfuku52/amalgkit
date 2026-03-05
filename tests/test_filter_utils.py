import os

import pandas
import pytest

from amalgkit.filter_utils import (
    build_species_prefixed_filename,
    infer_latest_filter_metadata,
    merge_metadata_by_run,
)


def test_merge_metadata_by_run_updates_existing_rows_and_adds_new_columns():
    source = pandas.DataFrame({
        'run': ['R1', 'R2', 'R3'],
        'exclusion': ['no', 'no', 'no'],
        'sample_group': ['a', 'a', 'b'],
    })
    update = pandas.DataFrame({
        'run': ['R2', 'R3'],
        'exclusion': ['low_within_sample_group_correlation', 'no'],
        'ws_margin': [-0.12, 0.05],
    })
    merged = merge_metadata_by_run(source_df=source, update_df=update)
    assert list(merged['run']) == ['R1', 'R2', 'R3']
    assert merged.loc[merged['run'] == 'R1', 'exclusion'].iloc[0] == 'no'
    assert merged.loc[merged['run'] == 'R2', 'exclusion'].iloc[0] == 'low_within_sample_group_correlation'
    assert 'ws_margin' in merged.columns
    assert pandas.isna(merged.loc[merged['run'] == 'R1', 'ws_margin'].iloc[0])


def test_merge_metadata_by_run_requires_run_column():
    source = pandas.DataFrame({'run': ['R1']})
    update = pandas.DataFrame({'x': [1]})
    with pytest.raises(ValueError, match='run'):
        merge_metadata_by_run(source_df=source, update_df=update)


def test_merge_metadata_by_run_handles_dtype_change_numeric_to_string():
    source = pandas.DataFrame({
        'run': ['R1', 'R2'],
        'pca_axis': [0.1, 0.2],
    })
    update = pandas.DataFrame({
        'run': ['R2'],
        'pca_axis': ['not_provided'],
    })
    merged = merge_metadata_by_run(source_df=source, update_df=update)
    value = merged.loc[merged['run'] == 'R2', 'pca_axis'].iloc[0]
    assert str(value) == 'not_provided'


def test_infer_latest_filter_metadata_returns_none_when_no_filter_outputs(tmp_path):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    assert infer_latest_filter_metadata(str(out_dir)) is None


def test_infer_latest_filter_metadata_prefers_latest_timestamp(tmp_path):
    out_dir = tmp_path / 'out'
    ws_meta = out_dir / 'wsfilter' / 'metadata.tsv'
    cs_meta = out_dir / 'csfilter' / 'metadata.tsv'
    ws_meta.parent.mkdir(parents=True)
    cs_meta.parent.mkdir(parents=True)
    ws_meta.write_text('run\texclusion\nR1\tno\n')
    cs_meta.write_text('run\texclusion\nR1\tno\n')
    os.utime(ws_meta, (1_700_000_000, 1_700_000_000))
    os.utime(cs_meta, (1_700_000_100, 1_700_000_100))
    latest = infer_latest_filter_metadata(str(out_dir))
    assert latest == os.path.realpath(str(cs_meta))


def test_build_species_prefixed_filename_rewrites_legacy_plot_name():
    out = build_species_prefixed_filename('Species_A', 'Species_A.2.correlation_cutoff.no.pdf')
    assert out == 'Species_A_2_correlation_cutoff_no.pdf'


def test_build_species_prefixed_filename_adds_species_prefix_when_missing():
    out = build_species_prefixed_filename('Species_A', 'csca_overview.pdf')
    assert out == 'Species_A_csca_overview.pdf'
