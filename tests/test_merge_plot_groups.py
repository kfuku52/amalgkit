import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes

from amalgkit import merge_plots


@pytest.mark.parametrize('kind', ['boxplot', 'histogram'])
@pytest.mark.parametrize('name_dtype', ['object', 'string', 'category'])
def test_species_plots_keep_normalized_sorted_groups_and_row_order(tmp_path, monkeypatch, kind, name_dtype):
    metadata = pd.DataFrame(
        {
            'scientific_name': [' Species B ', 'Species A', 'Species B', None, '', 'NA', '0001', 1, '1', 'Species A', 'Species A', 'Species A', 'Species B'],
            'metric': ['8', '3', '2', '9', '10', '6', '5', '7', '4', 'bad', np.inf, np.nan, '100'],
            'exclusion': [' NO ', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'no', 'yes'],
        },
        index=[3, 0, 3, 1, 1, 0, 4, 4, 4, 7, 7, 8, 9],
    )
    metadata['scientific_name'] = metadata['scientific_name'].astype(name_dtype)
    original = metadata.copy(deep=True)
    captured = {}
    method_name = 'boxplot' if kind == 'boxplot' else 'hist'
    original_plot = getattr(Axes, method_name)

    def capture_plot(axis, values, *args, **kwargs):
        captured['values'] = [value.copy() for value in values]
        captured['kwargs'] = kwargs
        return original_plot(axis, values, *args, **kwargs)

    def capture_figure(figure, _path):
        captured['axis'] = figure.axes[0]

    monkeypatch.setattr(Axes, method_name, capture_plot)
    monkeypatch.setattr(merge_plots, '_save_figure', capture_figure)
    output = str(tmp_path / 'plot.pdf')
    if kind == 'boxplot':
        merge_plots._save_species_boxplot_pdf(metadata, 'metric', output, 'Metric', y_limits=(0, 12))
    else:
        merge_plots._save_species_histogram_pdf(
            metadata, 'metric', output, 'Metric', bin_breaks=[0, 4, 8, 12], x_breaks=[0, 6, 12], x_limits=(0, 12),
        )
    assert len(captured['values']) == 5
    for actual, expected in zip(captured['values'], [[5], [7, 4], [6], [3], [8, 2]], strict=True):
        np.testing.assert_array_equal(actual, expected)
    axis = captured['axis']
    if kind == 'boxplot':
        assert [label.get_text() for label in axis.get_xticklabels()] == ['0001', '1', 'NA', 'Species\nA', 'Species\nB']
        assert axis.get_ylabel() == 'Metric'
        assert axis.get_ylim() == (0, 12)
    else:
        assert captured['kwargs']['label'] == ['0001', '1', 'NA', 'Species A', 'Species B']
        np.testing.assert_array_equal(captured['kwargs']['bins'], [0, 4, 8, 12])
        np.testing.assert_array_equal(axis.get_xticks(), [0, 6, 12])
        assert axis.get_xlabel() == 'Metric'
        assert axis.get_xlim() == (0, 12)
        assert [label.get_text() for label in axis.get_legend().get_texts()] == ['0001', '1', 'NA', 'Species A', 'Species B']
    pd.testing.assert_frame_equal(metadata, original)


@pytest.mark.parametrize('kind', ['boxplot', 'histogram'])
@pytest.mark.parametrize('case', ['missing_metric', 'missing_exclusion', 'excluded', 'empty_species', 'nonfinite'])
def test_species_plots_skip_when_no_valid_values(tmp_path, monkeypatch, kind, case, capsys):
    metadata = pd.DataFrame({'scientific_name': ['Species A'], 'metric': [1.0], 'exclusion': ['no']})
    if case == 'missing_metric':
        metadata = metadata.drop(columns='metric')
    elif case == 'missing_exclusion':
        metadata = metadata.drop(columns='exclusion')
    elif case == 'excluded':
        metadata['exclusion'] = pd.NA
    elif case == 'empty_species':
        metadata['scientific_name'] = '  '
    else:
        metadata['metric'] = np.inf

    def unexpected_pyplot():
        raise AssertionError('Skipped plots must not create figures')

    monkeypatch.setattr(merge_plots, '_load_pyplot', unexpected_pyplot)
    helper = getattr(merge_plots, f'_save_species_{kind}_pdf')
    helper(metadata, 'metric', str(tmp_path / 'plot.pdf'), 'Metric')
    assert 'Skipping plot generation.' in capsys.readouterr().out
    assert not (tmp_path / 'plot.pdf').exists()
