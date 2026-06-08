import os
import re
import shutil
import tempfile
import warnings
from contextlib import contextmanager

import pandas


def merge_metadata_by_run(source_df, update_df):
    if 'run' not in source_df.columns:
        raise ValueError('Source metadata is missing required "run" column.')
    if 'run' not in update_df.columns:
        raise ValueError('Updated metadata is missing required "run" column.')
    source = source_df.copy()
    update = update_df.copy()
    source['run'] = source['run'].fillna('').astype(str).str.strip()
    update['run'] = update['run'].fillna('').astype(str).str.strip()
    update = update.loc[update['run'] != '', :].drop_duplicates(subset=['run'], keep='last')
    source = source.set_index('run', drop=False)
    update = update.set_index('run', drop=False)
    for col in update.columns:
        if col not in source.columns:
            source[col] = pandas.NA
        update_col = update[col]
        if update_col.isna().all():
            continue
        for run_id, value in update_col.dropna().items():
            if run_id not in source.index:
                continue
            try:
                source.loc[run_id, col] = value
            except (TypeError, ValueError):
                # Some metadata columns can change type across steps (e.g., numeric -> string labels).
                source[col] = source[col].astype('object')
                source.loc[run_id, col] = value
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
    target_dir = os.path.realpath(target_dir)
    parent_dir = os.path.dirname(target_dir)
    if parent_dir != '':
        if os.path.exists(parent_dir) and (not os.path.isdir(parent_dir)):
            raise NotADirectoryError('Output parent path exists but is not a directory: {}'.format(parent_dir))
        os.makedirs(parent_dir, exist_ok=True)
    if os.path.lexists(target_dir) and (not redo):
        raise FileExistsError('Output already exists. Use --redo yes to overwrite: {}'.format(target_dir))
    stage_dir = tempfile.mkdtemp(prefix=prefix, dir=parent_dir if parent_dir != '' else None)
    backup_path = None
    committed = False
    try:
        yield stage_dir
        if os.path.lexists(target_dir):
            backup_path = tempfile.mkdtemp(prefix=prefix + 'backup_', dir=parent_dir if parent_dir != '' else None)
            os.rmdir(backup_path)
            os.rename(target_dir, backup_path)
        os.rename(stage_dir, target_dir)
        committed = True
        if backup_path is not None:
            if os.path.islink(backup_path) or os.path.isfile(backup_path):
                os.remove(backup_path)
            elif os.path.isdir(backup_path):
                shutil.rmtree(backup_path)
    except Exception:
        if (backup_path is not None) and (not os.path.lexists(target_dir)) and os.path.lexists(backup_path):
            os.rename(backup_path, target_dir)
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


def _format_genus_species_label(value):
    tokens = str(value).replace('_', ' ').split()
    if len(tokens) >= 2:
        return '{}\n{}'.format(tokens[0], ' '.join(tokens[1:]))
    if len(tokens) == 1:
        return tokens[0]
    return ''


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


def infer_latest_filter_metadata(out_dir):
    out_root = os.path.realpath(out_dir)
    candidates = [
        os.path.join(out_root, 'wsfilter', 'metadata.tsv'),
        os.path.join(out_root, 'csfilter', 'metadata.tsv'),
    ]
    existing = [path for path in candidates if os.path.isfile(path)]
    if len(existing) == 0:
        return None
    return max(existing, key=lambda path: os.path.getmtime(path))
