import os

import pandas
import pytest

from amalgkit.filter_utils import save_exclusion_plot_pdf, staged_output_dir
from amalgkit.output_utils import atomic_output_path


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


def test_atomic_outputs_use_shareable_default_permissions(tmp_path):
    previous_umask = os.umask(0o022)
    try:
        output_path = tmp_path / 'result.tsv'
        with atomic_output_path(str(output_path)) as temporary_path:
            with open(temporary_path, 'w', encoding='utf-8') as handle:
                handle.write('ready\n')
    finally:
        os.umask(previous_umask)

    assert os.stat(output_path).st_mode & 0o777 == 0o644


def test_staged_output_directory_uses_shareable_default_permissions(tmp_path):
    previous_umask = os.umask(0o022)
    try:
        output_dir = tmp_path / 'result'
        with staged_output_dir(str(output_dir), redo=True) as stage_dir:
            with open(os.path.join(stage_dir, 'ready.txt'), 'w', encoding='utf-8') as handle:
                handle.write('ready\n')
    finally:
        os.umask(previous_umask)

    assert os.stat(output_dir).st_mode & 0o777 == 0o755


def test_atomic_and_staged_outputs_respect_restrictive_umask(tmp_path):
    previous_umask = os.umask(0o077)
    try:
        output_path = tmp_path / 'private.tsv'
        with atomic_output_path(str(output_path)) as temporary_path:
            with open(temporary_path, 'w', encoding='utf-8') as handle:
                handle.write('private\n')
        output_dir = tmp_path / 'private-dir'
        with staged_output_dir(str(output_dir), redo=True):
            pass
    finally:
        os.umask(previous_umask)

    assert os.stat(output_path).st_mode & 0o777 == 0o600
    assert os.stat(output_dir).st_mode & 0o777 == 0o700


def test_atomic_output_can_replace_read_only_file_and_preserves_mode(tmp_path):
    output_path = tmp_path / 'read-only.tsv'
    output_path.write_text('old\n', encoding='utf-8')
    os.chmod(output_path, 0o444)

    with atomic_output_path(str(output_path)) as temporary_path:
        with open(temporary_path, 'w', encoding='utf-8') as handle:
            handle.write('new\n')

    assert output_path.read_text(encoding='utf-8') == 'new\n'
    assert os.stat(output_path).st_mode & 0o777 == 0o444
