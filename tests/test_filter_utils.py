import json
import os

import pandas
import pytest

from amalgkit.filter_utils import (
    FILTER_METADATA_STATE_FILENAME,
    get_filter_metadata_state_path,
    infer_latest_filter_metadata,
    merge_metadata_by_run,
    save_exclusion_plot_pdf,
    staged_output_dir,
    write_filter_metadata_state,
)
from amalgkit.output_utils import atomic_output_path


def _write_filter_metadata(out_dir, command):
    metadata_path = out_dir / command / 'metadata.tsv'
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text('run\texclusion\nR1\tno\n', encoding='utf-8')
    return metadata_path


def test_filter_metadata_state_tracks_last_completed_filter_independent_of_mtime(tmp_path):
    ws_metadata = _write_filter_metadata(tmp_path, 'wsfilter')
    cs_metadata = _write_filter_metadata(tmp_path, 'csfilter')
    write_filter_metadata_state(tmp_path, 'csfilter')
    os.utime(ws_metadata, (2_000_000_000, 2_000_000_000))
    os.utime(cs_metadata, (1_000_000_000, 1_000_000_000))

    assert infer_latest_filter_metadata(tmp_path) == str(cs_metadata)

    write_filter_metadata_state(tmp_path, 'wsfilter')
    os.utime(cs_metadata, (2_000_000_000, 2_000_000_000))
    os.utime(ws_metadata, (1_000_000_000, 1_000_000_000))

    assert infer_latest_filter_metadata(tmp_path) == str(ws_metadata)
    state = json.loads((tmp_path / FILTER_METADATA_STATE_FILENAME).read_text(encoding='utf-8'))
    assert state == {
        'command': 'wsfilter',
        'metadata_path': 'wsfilter/metadata.tsv',
        'schema_version': 1,
    }


def test_filter_metadata_legacy_fallback_is_deterministic_and_warns(tmp_path):
    ws_metadata = _write_filter_metadata(tmp_path, 'wsfilter')
    cs_metadata = _write_filter_metadata(tmp_path, 'csfilter')
    os.utime(ws_metadata, (2_000_000_000, 2_000_000_000))
    os.utime(cs_metadata, (1_000_000_000, 1_000_000_000))

    with pytest.warns(UserWarning, match='legacy deterministic fallback.*csfilter/metadata.tsv'):
        selected = infer_latest_filter_metadata(tmp_path)

    assert selected == str(cs_metadata)


def test_filter_metadata_invalid_state_cannot_escape_expected_stage_path(tmp_path):
    ws_metadata = _write_filter_metadata(tmp_path, 'wsfilter')
    state_path = get_filter_metadata_state_path(tmp_path)
    with open(state_path, 'w', encoding='utf-8') as handle:
        json.dump(
            {
                'schema_version': 1,
                'command': 'wsfilter',
                'metadata_path': '../outside.tsv',
            },
            handle,
        )

    with pytest.warns(UserWarning, match='metadata_path does not match.*legacy deterministic fallback'):
        selected = infer_latest_filter_metadata(tmp_path)

    assert selected == str(ws_metadata)


def test_write_filter_metadata_state_rejects_unknown_command_and_missing_output(tmp_path):
    with pytest.raises(ValueError, match='Unknown filter command'):
        write_filter_metadata_state(tmp_path, 'unknown')
    with pytest.raises(FileNotFoundError, match='Filter metadata was not generated'):
        write_filter_metadata_state(tmp_path, 'wsfilter')


def test_merge_metadata_by_run_rejects_update_runs_absent_from_source():
    source_df = pandas.DataFrame(
        {
            'run': ['R1', 'R2'],
            'mapping_rate': [0.1, 0.2],
        }
    )
    update_df = pandas.DataFrame(
        {
            'run': ['R1', 'R3'],
            'mapping_rate': [0.5, 0.7],
        }
    )

    with pytest.raises(ValueError, match='absent from source metadata: R3'):
        merge_metadata_by_run(source_df, update_df)


def test_merge_metadata_by_run_preserves_numeric_dtype_for_numeric_strings():
    source_df = pandas.DataFrame(
        {
            'run': ['R1', 'R2'],
            'mapping_rate': [0.1, 0.2],
        }
    )
    update_df = pandas.DataFrame(
        {
            'run': ['R1'],
            'mapping_rate': pandas.Series(['0.5'], dtype='object'),
        }
    )

    observed = merge_metadata_by_run(source_df, update_df)

    assert pandas.api.types.is_float_dtype(observed['mapping_rate'].dtype)
    assert observed['mapping_rate'].tolist() == [0.5, 0.2]


def test_merge_metadata_by_run_warns_before_promoting_incompatible_numeric_column():
    source_df = pandas.DataFrame(
        {
            'run': ['R1', 'R2'],
            'mapping_rate': [0.1, 0.2],
        }
    )
    update_df = pandas.DataFrame(
        {
            'run': ['R1'],
            'mapping_rate': ['not_available'],
        }
    )

    with pytest.warns(UserWarning, match='mapping_rate.*float64.*object.*R1'):
        observed = merge_metadata_by_run(source_df, update_df)

    assert pandas.api.types.is_object_dtype(observed['mapping_rate'].dtype)
    assert observed['mapping_rate'].tolist() == ['not_available', 0.2]


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
