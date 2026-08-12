import json
import os
import re
import shutil
import tempfile
import warnings
from contextlib import contextmanager

import pandas

from amalgkit.download_utils import acquire_exclusive_lock
from amalgkit.output_utils import atomic_output_path, get_default_creation_mode
from amalgkit.runtime_utils import format_species_label as _format_genus_species_label


FILTER_METADATA_STATE_FILENAME = 'filter_metadata_state.json'
FILTER_METADATA_STATE_SCHEMA_VERSION = 1
FILTER_COMMANDS = ('wsfilter', 'csfilter')


def _require_unique_nonempty_run_ids(run_series, context):
    normalized = run_series.fillna('').astype(str).str.strip()
    if bool(normalized.eq('').any()):
        raise ValueError(
            '{} metadata contains empty run IDs; every row must have a run identifier.'.format(
                context
            )
        )
    duplicated = set(normalized.loc[normalized.duplicated()].tolist())
    if duplicated:
        raise ValueError(
            '{} metadata contains duplicate run IDs: {}.'.format(
                context,
                ', '.join(sorted(duplicated)[:5]),
            )
        )


def merge_metadata_by_run(source_df, update_df):
    if 'run' not in source_df.columns:
        raise ValueError('Source metadata is missing required "run" column.')
    if 'run' not in update_df.columns:
        raise ValueError('Updated metadata is missing required "run" column.')
    source = source_df.copy()
    update = update_df.copy()
    source['run'] = source['run'].fillna('').astype(str).str.strip()
    update['run'] = update['run'].fillna('').astype(str).str.strip()
    _require_unique_nonempty_run_ids(source['run'], 'Source')
    update = update.loc[update['run'] != '', :]
    _require_unique_nonempty_run_ids(update['run'], 'Updated')
    source = source.set_index('run', drop=False)
    update = update.set_index('run', drop=False)
    source_run_ids = set(source.index.tolist())
    missing_run_ids = [run_id for run_id in update.index.tolist() if run_id not in source_run_ids]
    if missing_run_ids:
        raise ValueError(
            'Updated metadata contains run IDs absent from source metadata: {}.'.format(
                ', '.join(missing_run_ids)
            )
        )
    for col in update.columns:
        if col not in source.columns:
            source[col] = pandas.NA
        update_values = update[col].dropna()
        if update_values.shape[0] == 0:
            continue

        source_dtype = source[col].dtype
        if pandas.api.types.is_numeric_dtype(source_dtype):
            numeric_values = pandas.to_numeric(update_values, errors='coerce')
            incompatible = numeric_values.isna() & update_values.notna()
            if not bool(incompatible.any()):
                combined_values = pandas.concat(
                    [
                        source[col].reset_index(drop=True),
                        numeric_values.reset_index(drop=True),
                    ],
                    ignore_index=True,
                )
                target_dtype = combined_values.dtype
                source[col] = source[col].astype(target_dtype)
                for run_id, value in numeric_values.items():
                    source.loc[run_id, col] = value
                continue

            warnings.warn(
                'Metadata column "{}" contains values incompatible with dtype {}; '
                'promoting the column to object while merging run(s): {}.'.format(
                    col,
                    source_dtype,
                    ', '.join(update_values.index[incompatible].astype(str).tolist()),
                ),
                UserWarning,
                stacklevel=2,
            )
            source[col] = source[col].astype('object')
            for run_id, value in update_values.items():
                source.loc[run_id, col] = value
            continue

        if pandas.api.types.is_string_dtype(source_dtype) and not pandas.api.types.is_object_dtype(source_dtype):
            for run_id, value in update_values.astype(str).items():
                source.loc[run_id, col] = value
            continue

        candidate = source[col].copy()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('error', FutureWarning)
                for run_id, value in update_values.items():
                    candidate.loc[run_id] = value
        except (FutureWarning, TypeError, ValueError):
            warnings.warn(
                'Metadata column "{}" contains values incompatible with dtype {}; '
                'promoting the column to object while merging run(s): {}.'.format(
                    col,
                    source_dtype,
                    ', '.join(update_values.index.astype(str).tolist()),
                ),
                UserWarning,
                stacklevel=2,
            )
            candidate = source[col].astype('object')
            for run_id, value in update_values.items():
                candidate.loc[run_id] = value
        source[col] = candidate
    return source.reset_index(drop=True)


def collect_per_species_metadata_tables(per_species_dir):
    tables = []
    if not os.path.isdir(per_species_dir):
        return tables
    for species in sorted(os.listdir(per_species_dir)):
        species_dir = os.path.join(per_species_dir, species)
        tables_dir = os.path.join(species_dir, 'tables')
        if not os.path.isdir(tables_dir):
            continue
        for name in sorted(os.listdir(tables_dir)):
            if not name.endswith('.metadata.tsv'):
                continue
            path = os.path.join(tables_dir, name)
            if os.path.isfile(path):
                tables.append(path)
    return tables


def load_merged_per_species_metadata(per_species_dir):
    metadata_tables = collect_per_species_metadata_tables(per_species_dir=per_species_dir)
    if len(metadata_tables) == 0:
        raise FileNotFoundError('No per-species metadata table was found under: {}'.format(per_species_dir))
    frames = [
        pandas.read_csv(path, sep='\t', low_memory=False)
        for path in metadata_tables
    ]
    return pandas.concat(frames, axis=0, ignore_index=True, sort=False)


def prepare_output_dir(path_dir, redo=False):
    if os.path.lexists(path_dir):
        if not redo:
            raise FileExistsError('Output already exists. Use --redo yes to overwrite: {}'.format(path_dir))
        if os.path.islink(path_dir) or os.path.isfile(path_dir):
            os.remove(path_dir)
        else:
            shutil.rmtree(path_dir)
    os.makedirs(path_dir, exist_ok=True)


@contextmanager
def staged_output_dir(target_dir, redo=False, prefix='amalgkit_stage_'):
    absolute_target_dir = os.path.abspath(target_dir)
    if os.path.lexists(absolute_target_dir) and os.path.islink(absolute_target_dir):
        raise ValueError('Refusing to replace symbolic-link output directory: {}'.format(absolute_target_dir))
    parent_dir = os.path.realpath(os.path.dirname(absolute_target_dir))
    target_dir = os.path.join(parent_dir, os.path.basename(absolute_target_dir))
    if parent_dir != '':
        if os.path.exists(parent_dir) and (not os.path.isdir(parent_dir)):
            raise NotADirectoryError('Output parent path exists but is not a directory: {}'.format(parent_dir))
        os.makedirs(parent_dir, exist_ok=True)
    if os.path.lexists(target_dir) and (not redo):
        raise FileExistsError('Output already exists. Use --redo yes to overwrite: {}'.format(target_dir))
    stage_dir = tempfile.mkdtemp(prefix=prefix, dir=parent_dir if parent_dir != '' else None)
    lock_path = os.path.join(
        parent_dir if parent_dir != '' else '.',
        '.{}.amalgkit-output.lock'.format(os.path.basename(target_dir)),
    )
    backup_path = None
    committed = False
    try:
        yield stage_dir
        with acquire_exclusive_lock(
            lock_path=lock_path,
            lock_label='output directory commit',
            poll_seconds=1,
        ):
            # The target may have appeared while work was being written to the
            # stage directory. Re-check under the commit lock so redo=False is
            # a durable no-overwrite guarantee rather than a start-time hint.
            if os.path.lexists(target_dir) and os.path.islink(target_dir):
                raise ValueError('Refusing to replace symbolic-link output directory: {}'.format(target_dir))
            if os.path.lexists(target_dir) and (not redo):
                raise FileExistsError('Output already exists. Use --redo yes to overwrite: {}'.format(target_dir))
            target_mode = get_default_creation_mode(parent_dir, is_directory=True)
            if os.path.isdir(target_dir):
                target_mode = os.stat(target_dir, follow_symlinks=False).st_mode & 0o777
            os.chmod(stage_dir, target_mode)
            if os.path.lexists(target_dir):
                backup_path = tempfile.mkdtemp(prefix=prefix + 'backup_', dir=parent_dir if parent_dir != '' else None)
                os.rmdir(backup_path)
                os.rename(target_dir, backup_path)
            try:
                os.rename(stage_dir, target_dir)
            except BaseException:
                if (backup_path is not None) and (not os.path.lexists(target_dir)) and os.path.lexists(backup_path):
                    os.rename(backup_path, target_dir)
                    backup_path = None
                raise
            committed = True
            if backup_path is not None:
                if os.path.islink(backup_path) or os.path.isfile(backup_path):
                    os.remove(backup_path)
                elif os.path.isdir(backup_path):
                    shutil.rmtree(backup_path)
                backup_path = None
    except Exception:
        raise
    finally:
        if (not committed) and os.path.isdir(stage_dir):
            shutil.rmtree(stage_dir, ignore_errors=True)
        if (not committed) and (backup_path is not None) and os.path.lexists(backup_path):
            if os.path.islink(backup_path) or os.path.isfile(backup_path):
                os.remove(backup_path)
            elif os.path.isdir(backup_path):
                shutil.rmtree(backup_path, ignore_errors=True)


def copy_per_species_plots(per_species_dir, dst_plot_dir):
    os.makedirs(dst_plot_dir, exist_ok=True)
    if not os.path.isdir(per_species_dir):
        return
    for species in sorted(os.listdir(per_species_dir)):
        src_plot_dir = os.path.join(per_species_dir, species, 'plots')
        if not os.path.isdir(src_plot_dir):
            continue
        dst_species_plot_dir = os.path.join(dst_plot_dir, species)
        if os.path.lexists(dst_species_plot_dir):
            if os.path.islink(dst_species_plot_dir) or os.path.isfile(dst_species_plot_dir):
                os.remove(dst_species_plot_dir)
            else:
                shutil.rmtree(dst_species_plot_dir)
        shutil.copytree(src_plot_dir, dst_species_plot_dir)


def copy_root_pdf_plots(src_dir, dst_plot_dir):
    os.makedirs(dst_plot_dir, exist_ok=True)
    if not os.path.isdir(src_dir):
        return
    for name in sorted(os.listdir(src_dir)):
        if not name.lower().endswith('.pdf'):
            continue
        src_path = os.path.join(src_dir, name)
        if not os.path.isfile(src_path):
            continue
        shutil.copy2(src_path, os.path.join(dst_plot_dir, name))


def build_species_prefixed_filename(species, filename):
    stem, ext = os.path.splitext(filename)
    species = str(species)
    stem = str(stem)
    if stem.startswith(species + '.'):
        suffix = stem[len(species) + 1:]
    elif stem.startswith(species + '_'):
        suffix = stem[len(species) + 1:]
    elif stem == species:
        suffix = ''
    else:
        suffix = stem
    suffix = re.sub(r'[^A-Za-z0-9]+', '_', suffix).strip('_')
    if suffix == '':
        return species + ext
    return '{}_{}{}'.format(species, suffix, ext)


def copy_per_species_pdfs(per_species_dir, dst_dir, species_subset=None):
    os.makedirs(dst_dir, exist_ok=True)
    wanted = None
    if species_subset is not None:
        wanted = set(species_subset)
    if not os.path.isdir(per_species_dir):
        return
    for species in sorted(os.listdir(per_species_dir)):
        if (wanted is not None) and (species not in wanted):
            continue
        src_plot_dir = os.path.join(per_species_dir, species, 'plots')
        if not os.path.isdir(src_plot_dir):
            continue
        dst_species_dir = os.path.join(dst_dir, species)
        os.makedirs(dst_species_dir, exist_ok=True)
        for name in sorted(os.listdir(src_plot_dir)):
            if not name.lower().endswith('.pdf'):
                continue
            src_path = os.path.join(src_plot_dir, name)
            if not os.path.isfile(src_path):
                continue
            dst_name = build_species_prefixed_filename(species=species, filename=name)
            shutil.copy2(src_path, os.path.join(dst_species_dir, dst_name))


def copy_root_pdfs_to_species_dirs(src_dir, dst_dir, species_list):
    if not os.path.isdir(src_dir):
        return
    root_pdfs = []
    for name in sorted(os.listdir(src_dir)):
        if not name.lower().endswith('.pdf'):
            continue
        src_path = os.path.join(src_dir, name)
        if not os.path.isfile(src_path):
            continue
        root_pdfs.append((name, src_path))
    if len(root_pdfs) == 0:
        return
    for species in species_list:
        dst_species_dir = os.path.join(dst_dir, species)
        os.makedirs(dst_species_dir, exist_ok=True)
        for name, src_path in root_pdfs:
            dst_name = build_species_prefixed_filename(species=species, filename=name)
            shutil.copy2(src_path, os.path.join(dst_species_dir, dst_name))


def save_exclusion_plot_pdf(df_metadata, out_pdf_path, y_label='Sample count', font_size=8):
    os.makedirs(os.path.dirname(os.path.realpath(out_pdf_path)), exist_ok=True)
    if ('scientific_name' not in df_metadata.columns) or ('exclusion' not in df_metadata.columns):
        warnings.warn('Missing scientific_name/exclusion columns. Skipping exclusion plot: {}'.format(out_pdf_path))
        return
    df_plot = df_metadata.loc[:, ['scientific_name', 'exclusion']].copy()
    df_plot['scientific_name'] = df_plot['scientific_name'].fillna('').astype(str).str.strip()
    df_plot['exclusion'] = df_plot['exclusion'].fillna('').astype(str).str.strip()
    df_plot = df_plot.loc[
        (df_plot['scientific_name'] != '') &
        (df_plot['exclusion'] != ''),
        :
    ]
    if df_plot.empty:
        warnings.warn('No valid scientific_name/exclusion rows. Skipping exclusion plot: {}'.format(out_pdf_path))
        return

    try:
        import matplotlib
        matplotlib.use('Agg', force=True)
        from matplotlib import pyplot
    except ImportError as exc:
        warnings.warn('matplotlib is required to generate exclusion plot {}: {}'.format(out_pdf_path, exc))
        return

    species_order = df_plot.loc[:, 'scientific_name'].drop_duplicates().tolist()
    exclusion_order = df_plot.loc[:, 'exclusion'].drop_duplicates().tolist()
    summary = (
        df_plot.groupby(['scientific_name', 'exclusion'], sort=False)
        .size()
        .rename('count')
        .reset_index()
    )
    pivot = (
        summary.pivot(index='scientific_name', columns='exclusion', values='count')
        .reindex(index=species_order, columns=exclusion_order)
        .fillna(0.0)
    )

    plot_width = max(3.6, 0.11 * len(species_order))
    figure, axis = pyplot.subplots(figsize=(plot_width, 3.6))
    x_positions = list(range(len(species_order)))
    bottoms = [0.0] * len(species_order)
    cmap = pyplot.get_cmap('tab20')
    for idx, exclusion_value in enumerate(exclusion_order):
        heights = pivot.loc[:, exclusion_value].astype(float).tolist()
        axis.bar(
            x_positions,
            heights,
            bottom=bottoms,
            label=str(exclusion_value),
            color=cmap(idx % max(1, cmap.N)),
            width=0.8,
        )
        bottoms = [bottom + height for bottom, height in zip(bottoms, heights)]

    axis.set_xlim(-0.5, max(len(species_order) - 0.5, 0.5))
    axis.set_ylabel(str(y_label), fontsize=float(font_size))
    axis.set_xlabel('')
    axis.set_xticks(x_positions)
    axis.set_xticklabels(
        [_format_genus_species_label(value) for value in species_order],
        rotation=90,
        ha='center',
        fontsize=float(font_size),
    )
    axis.tick_params(axis='y', labelsize=float(font_size))
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.set_facecolor('white')
    axis.grid(axis='y', color='#d0d0d0', linewidth=0.6)
    if exclusion_order:
        axis.legend(
            title='',
            loc='upper center',
            bbox_to_anchor=(0.5, -0.22),
            ncol=max(1, min(4, len(exclusion_order))),
            frameon=False,
            fontsize=float(font_size),
        )
    figure.tight_layout()
    figure.savefig(os.path.realpath(out_pdf_path), format='pdf', bbox_inches='tight')
    pyplot.close(figure)


def get_filter_metadata_state_path(out_dir):
    return os.path.join(os.path.realpath(out_dir), FILTER_METADATA_STATE_FILENAME)


def _get_filter_metadata_paths(out_dir):
    out_root = os.path.realpath(out_dir)
    return {
        command: os.path.join(out_root, command, 'metadata.tsv')
        for command in FILTER_COMMANDS
    }


def write_filter_metadata_state(out_dir, command):
    command = str(command).strip()
    if command not in FILTER_COMMANDS:
        raise ValueError('Unknown filter command for metadata state: {}'.format(command))
    metadata_paths = _get_filter_metadata_paths(out_dir)
    metadata_path = metadata_paths[command]
    if not os.path.isfile(metadata_path):
        raise FileNotFoundError('Filter metadata was not generated: {}'.format(metadata_path))
    state_path = get_filter_metadata_state_path(out_dir)
    payload = {
        'schema_version': FILTER_METADATA_STATE_SCHEMA_VERSION,
        'command': command,
        'metadata_path': '{}/metadata.tsv'.format(command),
    }
    with atomic_output_path(state_path, suffix='.json') as temporary_path:
        with open(temporary_path, 'w', encoding='utf-8') as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write('\n')
    return state_path


def _read_filter_metadata_state(out_dir):
    state_path = get_filter_metadata_state_path(out_dir)
    if not os.path.lexists(state_path):
        return None, 'state file is missing'
    if os.path.islink(state_path):
        return None, 'state file is a symbolic link'
    if not os.path.isfile(state_path):
        return None, 'state path is not a file'
    try:
        with open(state_path, 'r', encoding='utf-8') as handle:
            payload = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return None, 'state file could not be read: {}'.format(exc)
    if not isinstance(payload, dict):
        return None, 'state payload is not an object'
    if payload.get('schema_version') != FILTER_METADATA_STATE_SCHEMA_VERSION:
        return None, 'unsupported state schema_version: {}'.format(payload.get('schema_version'))
    command = str(payload.get('command', '')).strip()
    if command not in FILTER_COMMANDS:
        return None, 'unknown filter command: {}'.format(command)
    expected_relative_path = '{}/metadata.tsv'.format(command)
    recorded_relative_path = str(payload.get('metadata_path', '')).strip()
    if os.path.normpath(recorded_relative_path) != os.path.normpath(expected_relative_path):
        return None, 'metadata_path does not match command {}'.format(command)
    metadata_path = _get_filter_metadata_paths(out_dir)[command]
    if not os.path.isfile(metadata_path):
        return None, 'recorded metadata file is missing: {}'.format(metadata_path)
    return metadata_path, None


def infer_latest_filter_metadata(out_dir):
    state_metadata_path, state_error = _read_filter_metadata_state(out_dir)
    if state_metadata_path is not None:
        return state_metadata_path

    metadata_paths = _get_filter_metadata_paths(out_dir)
    existing = [
        metadata_paths[command]
        for command in ('csfilter', 'wsfilter')
        if os.path.isfile(metadata_paths[command])
    ]
    if len(existing) == 0:
        if state_error != 'state file is missing':
            warnings.warn(
                'Ignoring invalid filter metadata state ({}); no completed filter metadata was found.'.format(
                    state_error
                ),
                UserWarning,
                stacklevel=2,
            )
        return None
    selected = existing[0]
    warnings.warn(
        'Filter metadata state is unavailable ({}). Using legacy deterministic fallback: {}. '
        'Run wsfilter or csfilter to refresh {}.'.format(
            state_error,
            selected,
            FILTER_METADATA_STATE_FILENAME,
        ),
        UserWarning,
        stacklevel=2,
    )
    return selected
