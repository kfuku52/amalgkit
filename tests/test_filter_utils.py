import pandas
import pytest

from amalgkit.filter_utils import save_exclusion_plot_pdf


def test_save_exclusion_plot_pdf_writes_pdf(tmp_path):
    metadata_df = pandas.DataFrame(
        {
            'scientific_name': [
                'Homo_sapiens',
                'Homo_sapiens',
                'Mus_musculus',
                'Mus_musculus',
                'Danio_rerio',
            ],
            'exclusion': ['no', 'yes', 'no', 'no', 'yes'],
        }
    )
    out_pdf_path = tmp_path / 'exclusion.pdf'

    save_exclusion_plot_pdf(
        df_metadata=metadata_df,
        out_pdf_path=str(out_pdf_path),
        y_label='Sample count',
        font_size=8,
    )

    assert out_pdf_path.exists()
    assert out_pdf_path.stat().st_size > 0


def test_save_exclusion_plot_pdf_warns_and_skips_for_missing_columns(tmp_path):
    metadata_df = pandas.DataFrame({'scientific_name': ['Homo_sapiens']})
    out_pdf_path = tmp_path / 'missing.pdf'

    with pytest.warns(UserWarning, match='Missing scientific_name/exclusion columns'):
        save_exclusion_plot_pdf(
            df_metadata=metadata_df,
            out_pdf_path=str(out_pdf_path),
            y_label='Sample count',
            font_size=8,
        )

    assert not out_pdf_path.exists()
