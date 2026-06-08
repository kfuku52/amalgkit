import math
import os
import warnings

import numpy
import pandas

from amalgkit.filter_utils import save_exclusion_plot_pdf


def _load_pyplot():
    import matplotlib

    matplotlib.use('Agg', force=True)
    from matplotlib import pyplot

    return pyplot


def _normalize_exclusion_values(exclusion_values):
    normalized = pandas.Series(exclusion_values).astype('string')
    normalized = normalized.fillna(pandas.NA).str.strip().str.lower()
    return normalized


def _is_non_excluded(df):
    if 'exclusion' not in df.columns:
        return pandas.Series([False] * df.shape[0], index=df.index, dtype=bool)
    normalized = _normalize_exclusion_values(df.loc[:, 'exclusion'])
    return normalized.eq('no').fillna(False)


def _format_species_label(value):
    tokens = str(value).replace('_', ' ').split()
    if len(tokens) >= 2:
        return '{}\n{}'.format(tokens[0], ' '.join(tokens[1:]))
    if len(tokens) == 1:
        return tokens[0]
    return ''


def _ordered_unique(values):
    seen = set()
    ordered = []
    for value in values:
        if pandas.isna(value):
            continue
        text = str(value).strip()
        if text == '':
            continue
        if text in seen:
            continue
        seen.add(text)
        ordered.append(text)
    return ordered


def _resolve_plot_width(num_species):
    return max(3.6, 0.11 * max(1, int(num_species)))


def _save_figure(figure, out_path):
    figure.tight_layout()
    figure.savefig(os.path.realpath(out_path), format='pdf', bbox_inches='tight')


def _normalize_hist_breaks(values, bin_breaks=None, bins=40):
    if bin_breaks is not None:
        breaks = numpy.asarray(bin_breaks, dtype=float)
        breaks = numpy.unique(breaks[numpy.isfinite(breaks)])
        if breaks.size >= 2:
            return breaks
    values = numpy.asarray(values, dtype=float)
    values = values[numpy.isfinite(values)]
    if values.size == 0:
        return None
    value_min = float(values.min())
    value_max = float(values.max())
    if value_min == value_max:
        pad = max(1.0, abs(value_min) * 0.1)
        value_min -= pad
        value_max += pad
    try:
        breaks = numpy.histogram_bin_edges(values, bins=max(3, int(bins)))
    except Exception:
        breaks = numpy.linspace(value_min, value_max, num=max(4, int(bins) + 1))
    breaks = numpy.unique(breaks[numpy.isfinite(breaks)])
    if breaks.size < 2:
        step = max(1.0, math.ceil(abs(value_min) * 0.05))
        center = float((value_min + value_max) / 2.0)
        breaks = numpy.arange(center - 3.0 * step, center + 3.0 * step + step, step)
    return breaks


def _insert_axis_breaks(values):
    values = numpy.asarray(values, dtype=float)
    values = values[numpy.isfinite(values)]
    if values.size == 0:
        return None
    value_min = math.floor(float(values.min()))
    value_max = math.ceil(float(values.max()))
    value_range = value_max - value_min
    if value_range <= 0:
        step = 10 if value_max <= 250 else 25
    elif value_range <= 200:
        step = 25
    elif value_range <= 500:
        step = 50
    elif value_range <= 1000:
        step = 100
    else:
        step = 200
    min_tick = math.floor(value_min / step) * step
    max_tick = math.ceil(value_max / step) * step
    if min_tick == max_tick:
        min_tick -= 2 * step
        max_tick += 2 * step
    ticks = numpy.arange(min_tick, max_tick + step, step, dtype=float)
    if ticks.size < 3:
        ticks = numpy.arange(min_tick - step, max_tick + 2 * step, step, dtype=float)
    return ticks


def _save_species_boxplot_pdf(df, value_col, out_path, y_label, font_size=8, y_limits=None):
    if value_col not in df.columns:
        print('{} column not found. Skipping plot generation.'.format(value_col), flush=True)
        return
    df_plot = df.loc[_is_non_excluded(df), ['scientific_name', value_col]].copy()
    df_plot.loc[:, 'scientific_name'] = df_plot.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
    df_plot.loc[:, value_col] = pandas.to_numeric(df_plot.loc[:, value_col], errors='coerce')
    df_plot = df_plot.loc[
        (df_plot.loc[:, 'scientific_name'] != '') &
        numpy.isfinite(df_plot.loc[:, value_col].to_numpy(dtype=float)),
        :
    ]
    print(
        'Number of SRA samples for {} plotting: {:,}'.format(value_col, int(df_plot.shape[0])),
        flush=True,
    )
    if df_plot.empty:
        print('No data available for {}. Skipping plot generation.'.format(value_col), flush=True)
        return
    species_order = sorted(_ordered_unique(df_plot.loc[:, 'scientific_name']))
    grouped = [
        df_plot.loc[df_plot.loc[:, 'scientific_name'] == species, value_col].to_numpy(dtype=float)
        for species in species_order
    ]
    grouped = [values for values in grouped if values.size > 0]
    species_order = [species for species in species_order if (df_plot.loc[df_plot.loc[:, 'scientific_name'] == species, :].shape[0] > 0)]
    if len(grouped) == 0:
        return
    pyplot = _load_pyplot()
    figure, axis = pyplot.subplots(figsize=(_resolve_plot_width(len(species_order)), 3.6))
    axis.boxplot(
        grouped,
        tick_labels=[_format_species_label(species) for species in species_order],
        patch_artist=False,
    )
    axis.set_xlabel('')
    axis.set_ylabel(str(y_label), fontsize=float(font_size))
    axis.tick_params(axis='x', labelrotation=90, labelsize=float(font_size))
    axis.tick_params(axis='y', labelsize=float(font_size))
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.grid(axis='y', color='#d0d0d0', linewidth=0.6)
    if y_limits is not None:
        axis.set_ylim(y_limits)
    else:
        ymin, ymax = axis.get_ylim()
        axis.set_ylim(bottom=0.0, top=ymax)
    _save_figure(figure, out_path)
    pyplot.close(figure)


def _save_stacked_bar_pdf(summary_df, x_col, fill_col, y_col, out_path, y_label, legend_label, font_size=8):
    if summary_df.empty:
        return
    x_order = _ordered_unique(summary_df.loc[:, x_col])
    fill_order = _ordered_unique(summary_df.loc[:, fill_col])
    pivot = (
        summary_df.pivot(index=x_col, columns=fill_col, values=y_col)
        .reindex(index=x_order, columns=fill_order)
        .fillna(0.0)
    )
    pyplot = _load_pyplot()
    figure, axis = pyplot.subplots(figsize=(_resolve_plot_width(len(x_order)), 3.6))
    positions = numpy.arange(len(x_order), dtype=float)
    bottoms = numpy.zeros((len(x_order),), dtype=float)
    cmap = pyplot.get_cmap('tab20')
    for idx, category in enumerate(fill_order):
        heights = pivot.loc[:, category].to_numpy(dtype=float)
        axis.bar(
            positions,
            heights,
            bottom=bottoms,
            label=str(category),
            color=cmap(idx % max(1, cmap.N)),
            width=0.8,
        )
        bottoms = bottoms + heights
    axis.set_xticks(positions)
    axis.set_xticklabels([_format_species_label(value) for value in x_order], rotation=90, ha='center', fontsize=float(font_size))
    axis.set_xlabel('')
    axis.set_ylabel(str(y_label), fontsize=float(font_size))
    axis.tick_params(axis='y', labelsize=float(font_size))
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.grid(axis='y', color='#d0d0d0', linewidth=0.6)
    axis.legend(
        title='',
        loc='upper center',
        bbox_to_anchor=(0.5, -0.22),
        ncol=max(1, min(4, len(fill_order))),
        frameon=False,
        fontsize=float(font_size),
    )
    _ = legend_label
    _save_figure(figure, out_path)
    pyplot.close(figure)


def _save_species_histogram_pdf(
    df,
    value_col,
    out_path,
    x_label,
    font_size=8,
    bins=40,
    bin_breaks=None,
    x_breaks=None,
    x_limits=None,
):
    if value_col not in df.columns:
        print('{} column not found. Skipping plot generation.'.format(value_col), flush=True)
        return
    df_plot = df.loc[_is_non_excluded(df), ['scientific_name', value_col]].copy()
    df_plot.loc[:, 'scientific_name'] = df_plot.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
    df_plot.loc[:, value_col] = pandas.to_numeric(df_plot.loc[:, value_col], errors='coerce')
    df_plot = df_plot.loc[
        (df_plot.loc[:, 'scientific_name'] != '') &
        numpy.isfinite(df_plot.loc[:, value_col].to_numpy(dtype=float)),
        :
    ]
    print(
        'Number of SRA samples for {} plotting: {:,}'.format(value_col, int(df_plot.shape[0])),
        flush=True,
    )
    if df_plot.empty:
        print('No data available for {}. Skipping plot generation.'.format(value_col), flush=True)
        return
    species_order = sorted(_ordered_unique(df_plot.loc[:, 'scientific_name']))
    hist_values = [
        df_plot.loc[df_plot.loc[:, 'scientific_name'] == species, value_col].to_numpy(dtype=float)
        for species in species_order
    ]
    hist_values = [values for values in hist_values if values.size > 0]
    species_order = [species for species in species_order if (df_plot.loc[df_plot.loc[:, 'scientific_name'] == species, :].shape[0] > 0)]
    if len(hist_values) == 0:
        return
    breaks = _normalize_hist_breaks(df_plot.loc[:, value_col].to_numpy(dtype=float), bin_breaks=bin_breaks, bins=bins)
    if breaks is None:
        return
    pyplot = _load_pyplot()
    figure, axis = pyplot.subplots(figsize=(3.6, 3.6))
    cmap = pyplot.get_cmap('tab20')
    colors = [cmap(idx % max(1, cmap.N)) for idx in range(len(hist_values))]
    axis.hist(
        hist_values,
        bins=breaks,
        stacked=True,
        label=species_order,
        color=colors,
        edgecolor='black',
        linewidth=0.2,
    )
    axis.set_xlabel(str(x_label), fontsize=float(font_size))
    axis.set_ylabel('Sample count', fontsize=float(font_size))
    axis.tick_params(axis='x', labelsize=float(font_size))
    axis.tick_params(axis='y', labelsize=float(font_size))
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.grid(axis='y', color='#d0d0d0', linewidth=0.6)
    if x_breaks is not None:
        axis.set_xticks(list(x_breaks))
    if x_limits is not None:
        axis.set_xlim(x_limits)
    axis.legend(
        title='',
        loc='upper center',
        bbox_to_anchor=(0.5, -0.22),
        ncol=max(1, min(3, len(species_order))),
        frameon=False,
        fontsize=float(font_size),
    )
    _save_figure(figure, out_path)
    pyplot.close(figure)


def _read_est_count_summary(file_path):
    try:
        dat = pandas.read_csv(file_path, sep='\t', low_memory=False)
    except Exception as exc:
        warnings.warn('Failed to read est_counts file. Skipping {}: {}'.format(file_path, exc))
        return None
    sample_cols = [column for column in dat.columns if column != 'target_id']
    if len(sample_cols) == 0:
        warnings.warn('No sample columns found in est_counts file. Skipping {}'.format(file_path))
        return None
    dat_numeric = dat.loc[:, sample_cols].apply(pandas.to_numeric, errors='coerce')
    mean_before = dat_numeric.mean(axis=0, skipna=True)
    lib_size = dat_numeric.sum(axis=0, skipna=True)
    species_tag = os.path.basename(os.path.dirname(file_path))
    species_name = species_tag.replace('_', ' ')
    return pandas.DataFrame(
        {
            'scientific_name': [species_name] * len(sample_cols),
            'run': sample_cols,
            'mean_before': mean_before.to_numpy(dtype=float),
            'lib_size': lib_size.to_numpy(dtype=float),
        }
    )


def _collect_est_count_summaries(merge_dir):
    est_count_files = []
    for root, _dirs, files in os.walk(merge_dir):
        for name in files:
            if name.endswith('_est_counts.tsv'):
                est_count_files.append(os.path.join(root, name))
    est_count_files.sort()
    summaries = []
    for file_path in est_count_files:
        summary_df = _read_est_count_summary(file_path)
        if (summary_df is None) or summary_df.empty:
            continue
        summaries.append(summary_df)
    if len(summaries) == 0:
        return pandas.DataFrame(columns=['scientific_name', 'run', 'mean_before', 'lib_size'])
    out = pandas.concat(summaries, axis=0, ignore_index=True)
    row_key = out.loc[:, 'scientific_name'].astype(str) + '|||' + out.loc[:, 'run'].astype(str)
    duplicated = row_key.duplicated(keep='first')
    if duplicated.any():
        duplicated_keys = sorted(row_key.loc[duplicated].unique().tolist())
        warnings.warn(
            'Detected duplicated species/run rows across est_counts tables; keeping first entries: {}'.format(
                ', '.join(duplicated_keys)
            )
        )
        out = out.loc[~duplicated, :].reset_index(drop=True)
    return out


def _save_mean_expression_boxplot(metadata_df, merge_dir, out_path, font_size=8):
    summary_df = _collect_est_count_summaries(merge_dir=merge_dir)
    if summary_df.empty:
        print(
            'No est_counts tables were found for mean expression plotting. Skipping {}'.format(
                os.path.basename(out_path)
            ),
            flush=True,
        )
        return
    if {'scientific_name', 'run'}.issubset(metadata_df.columns):
        metadata_cols = [column for column in ['scientific_name', 'run', 'exclusion'] if column in metadata_df.columns]
        metadata_subset = metadata_df.loc[:, metadata_cols].drop_duplicates()
        summary_df = summary_df.merge(
            metadata_subset,
            on=['scientific_name', 'run'],
            how='left',
            sort=False,
        )
    elif 'run' in metadata_df.columns:
        metadata_cols = [column for column in ['run', 'exclusion'] if column in metadata_df.columns]
        metadata_subset = metadata_df.loc[:, metadata_cols].drop_duplicates()
        summary_df = summary_df.merge(metadata_subset, on='run', how='left', sort=False)
    else:
        summary_df.loc[:, 'exclusion'] = pandas.NA
    if 'exclusion' in summary_df.columns:
        has_metadata_match = summary_df.loc[:, 'exclusion'].notna()
        if int(has_metadata_match.sum()) == 0:
            print(
                'No est_counts runs matched metadata rows. Skipping {}'.format(os.path.basename(out_path)),
                flush=True,
            )
            return
        exclusion_norm = _normalize_exclusion_values(summary_df.loc[:, 'exclusion'])
        summary_df = summary_df.loc[has_metadata_match & exclusion_norm.eq('no').fillna(False), :].copy()
    print(
        'Number of SRA samples for mean_expression plotting: {:,}'.format(int(summary_df.shape[0])),
        flush=True,
    )
    if summary_df.empty:
        print(
            'No non-excluded samples available for mean expression plotting. Skipping {}'.format(
                os.path.basename(out_path)
            ),
            flush=True,
        )
        return
    valid_libsize = pandas.to_numeric(summary_df.loc[:, 'lib_size'], errors='coerce')
    valid_libsize = valid_libsize[numpy.isfinite(valid_libsize.to_numpy(dtype=float)) & (valid_libsize.to_numpy(dtype=float) > 0)]
    ref_libsize = float(valid_libsize.median()) if valid_libsize.shape[0] > 0 else 1.0
    if (not numpy.isfinite(ref_libsize)) or (ref_libsize <= 0):
        ref_libsize = 1.0
    scale_factor = pandas.to_numeric(summary_df.loc[:, 'lib_size'], errors='coerce') / ref_libsize
    scale_factor = scale_factor.where(numpy.isfinite(scale_factor.to_numpy(dtype=float)) & (scale_factor.to_numpy(dtype=float) > 0), 1.0)
    mean_before = pandas.to_numeric(summary_df.loc[:, 'mean_before'], errors='coerce').to_numpy(dtype=float)
    mean_after = mean_before / scale_factor.to_numpy(dtype=float)
    pyplot = _load_pyplot()
    figure, axis = pyplot.subplots(figsize=(2.6, 3.6))
    axis.boxplot(
        [mean_before[numpy.isfinite(mean_before)], mean_after[numpy.isfinite(mean_after)]],
        tick_labels=['Raw\ncounts', 'Library-size\ncorrected\ncounts'],
    )
    finite_values = numpy.concatenate([
        mean_before[numpy.isfinite(mean_before)],
        mean_after[numpy.isfinite(mean_after)],
    ])
    ymax = 1.0 if finite_values.size == 0 else float(finite_values.max())
    if (not numpy.isfinite(ymax)) or (ymax <= 0):
        ymax = 1.0
    else:
        ymax = ymax * 1.05
    axis.set_ylim(0.0, ymax)
    axis.set_xlabel('')
    axis.set_ylabel('Mean count of all genes', fontsize=float(font_size))
    axis.tick_params(axis='x', labelsize=float(font_size))
    axis.tick_params(axis='y', labelsize=float(font_size))
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.grid(axis='y', color='#d0d0d0', linewidth=0.6)
    _save_figure(figure, out_path)
    pyplot.close(figure)


def generate_merge_plots(merge_dir, metadata_path, font_size=8):
    metadata_df = pandas.read_csv(metadata_path, sep='\t', low_memory=False)
    if 'scientific_name' in metadata_df.columns:
        metadata_df.loc[:, 'scientific_name'] = metadata_df.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
    _save_species_boxplot_pdf(
        df=metadata_df,
        value_col='mapping_rate',
        out_path=os.path.join(merge_dir, 'merge_mapping_rate.pdf'),
        y_label='Mapping rate',
        font_size=font_size,
        y_limits=(0.0, 100.0),
    )
    _save_species_boxplot_pdf(
        df=metadata_df,
        value_col='total_spots',
        out_path=os.path.join(merge_dir, 'merge_total_spots.pdf'),
        y_label='Total spots',
        font_size=font_size,
    )
    _save_species_boxplot_pdf(
        df=metadata_df,
        value_col='total_bases',
        out_path=os.path.join(merge_dir, 'merge_total_bases.pdf'),
        y_label='Total bases',
        font_size=font_size,
    )

    layout_df = metadata_df.loc[_is_non_excluded(metadata_df), :].copy()
    if {'scientific_name', 'lib_layout'}.issubset(layout_df.columns):
        layout_df.loc[:, 'scientific_name'] = layout_df.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
        layout_df.loc[:, 'lib_layout'] = layout_df.loc[:, 'lib_layout'].fillna('').astype(str).str.strip()
        layout_df = layout_df.loc[(layout_df.loc[:, 'scientific_name'] != '') & (layout_df.loc[:, 'lib_layout'] != ''), :]
        if not layout_df.empty:
            summary_df = (
                layout_df.groupby(['scientific_name', 'lib_layout'], sort=False)
                .size()
                .rename('count')
                .reset_index()
            )
            _save_stacked_bar_pdf(
                summary_df=summary_df,
                x_col='scientific_name',
                fill_col='lib_layout',
                y_col='count',
                out_path=os.path.join(merge_dir, 'merge_library_layout.pdf'),
                y_label='Count',
                legend_label='Library layout',
                font_size=font_size,
            )
        else:
            print('No non-excluded samples with scientific_name/lib_layout values. Skipping merge_library_layout.pdf', flush=True)

    _save_mean_expression_boxplot(
        metadata_df=metadata_df,
        merge_dir=merge_dir,
        out_path=os.path.join(merge_dir, 'merge_mean_expression_boxplot.pdf'),
        font_size=font_size,
    )
    _save_species_histogram_pdf(
        df=metadata_df,
        value_col='fastp_duplication_rate',
        out_path=os.path.join(merge_dir, 'merge_fastp_duplication_rate_histogram.pdf'),
        x_label='Fastp duplication rate (%)',
        font_size=font_size,
        bin_breaks=numpy.arange(0, 101, 10, dtype=float),
        x_breaks=numpy.arange(0, 101, 10, dtype=float),
        x_limits=(0.0, 100.0),
    )
    if 'fastp_insert_size_peak' in metadata_df.columns:
        insert_df = metadata_df.loc[
            _is_non_excluded(metadata_df) & pandas.to_numeric(metadata_df.loc[:, 'fastp_insert_size_peak'], errors='coerce').notna(),
            :,
        ]
        insert_breaks = _insert_axis_breaks(pandas.to_numeric(insert_df.loc[:, 'fastp_insert_size_peak'], errors='coerce').to_numpy(dtype=float))
        insert_limits = None if insert_breaks is None else (float(insert_breaks.min()), float(insert_breaks.max()))
    else:
        insert_breaks = None
        insert_limits = None
    _save_species_histogram_pdf(
        df=metadata_df,
        value_col='fastp_insert_size_peak',
        out_path=os.path.join(merge_dir, 'merge_fastp_insert_size_peak_histogram.pdf'),
        x_label='Fastp insert size peak',
        font_size=font_size,
        bin_breaks=insert_breaks,
        x_breaks=insert_breaks,
        x_limits=insert_limits,
    )
    print('Number of SRA samples for exclusion plotting: {:,}'.format(int(metadata_df.shape[0])), flush=True)
    save_exclusion_plot_pdf(
        df_metadata=metadata_df,
        out_pdf_path=os.path.join(merge_dir, 'merge_exclusion.pdf'),
        y_label='Sample count',
        font_size=font_size,
    )


__all__ = [
    'generate_merge_plots',
]
