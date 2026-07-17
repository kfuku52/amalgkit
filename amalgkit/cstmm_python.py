import os
import re
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import pandas

from amalgkit.normalization_tmm import run_tmm_rounds_for_cstmm


def _normalize_species_prefix(scientific_name):
    prefix = str(scientific_name).strip()
    prefix = re.sub(r'\s+', '_', prefix)
    prefix = re.sub(r'_+', '_', prefix)
    return prefix


def _list_matching_files(path_dir, pattern):
    regex = re.compile(pattern)
    if not os.path.isdir(path_dir):
        return []
    items = []
    for name in sorted(os.listdir(path_dir)):
        path = os.path.join(path_dir, name)
        if os.path.isfile(path) and regex.match(name):
            items.append(name)
    return items


def _read_est_counts(dir_count, species_name):
    species_dir = os.path.join(dir_count, species_name)
    infile = _list_matching_files(species_dir, r'.*est_counts\.tsv$')
    if len(infile) > 1:
        raise ValueError('Multiple *est_counts.tsv files found: {}'.format(species_name))
    if len(infile) == 0:
        raise FileNotFoundError('No *est_counts.tsv files found: {}'.format(species_name))
    infile_path = os.path.join(species_dir, infile[0])
    dat = pandas.read_csv(infile_path, sep='\t', index_col=0)
    if 'length' in dat.columns:
        dat = dat.drop(columns=['length'])
    dat.columns = ['{}_{}'.format(species_name, col) for col in dat.columns]
    return dat


def _get_spp_filled(dir_count, df_gc=None):
    species_names = []
    for species_name in sorted(os.listdir(dir_count)):
        species_dir = os.path.join(dir_count, species_name)
        if (not os.path.isdir(species_dir)) or species_name.startswith('.') or species_name.startswith('tmp.'):
            continue
        count_files = _list_matching_files(species_dir, r'.*est_counts\.tsv$')
        if len(count_files) == 1:
            species_names.append(species_name)
    if df_gc is not None:
        species_names = [species_name for species_name in species_names if species_name in df_gc.columns]
    if len(species_names) == 0:
        raise FileNotFoundError('No species with valid *est_counts.tsv files were found under: {}'.format(dir_count))
    return species_names


def _read_genecount_table(file_genecount):
    df_gc = pandas.read_csv(file_genecount, sep='\t')
    if 'orthogroup_id' not in df_gc.columns:
        raise ValueError('Column "orthogroup_id" is required in genecount table: {}'.format(file_genecount))
    df_gc = df_gc.set_index('orthogroup_id')
    return df_gc


def _read_orthogroup_table(file_orthogroup_table):
    return pandas.read_csv(file_orthogroup_table, sep='\t', index_col=0)


def _load_uncorrected(dir_count, species_names):
    return {
        species_name: _read_est_counts(dir_count=dir_count, species_name=species_name)
        for species_name in species_names
    }


def _copy_eff_length_file(dir_count, dir_cstmm, species_name):
    species_dir = os.path.join(dir_count, species_name)
    eff_length_files = _list_matching_files(species_dir, r'.*eff_length\.tsv$')
    if len(eff_length_files) != 1:
        return
    src = os.path.join(species_dir, eff_length_files[0])
    dst = os.path.join(dir_cstmm, species_name, eff_length_files[0])
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)


def _fill_missing_with_row_mean(df):
    if not df.isna().any().any():
        return df
    row_means = df.mean(axis=1, skipna=True).fillna(0.0)
    return df.T.fillna(row_means).T


def _get_df_nonzero(df_counts):
    is_zero_col = (df_counts.sum(axis=0, skipna=True) == 0)
    df_nonzero = df_counts.loc[:, ~is_zero_col].copy()
    return _fill_missing_with_row_mean(df_nonzero)


def _get_singlecopy_bool_index(df_gc, spp_filled, percent_singlecopy_threshold=50.0):
    num_species = float(len(spp_filled))
    df_species = df_gc.loc[:, spp_filled]
    num_singlecopy_species = df_species.eq(1).sum(axis=1)
    percent_singlecopy_species = (num_singlecopy_species / num_species) * 100.0
    return percent_singlecopy_species.ge(float(percent_singlecopy_threshold))


def _get_df_exp_single_copy_ortholog(file_genecount, file_orthogroup_table, dir_count, uncorrected_by_species):
    df_gc = _read_genecount_table(file_genecount=file_genecount)
    df_og = _read_orthogroup_table(file_orthogroup_table=file_orthogroup_table)
    spp_filled = _get_spp_filled(dir_count=dir_count, df_gc=df_gc)
    is_singlecopy = _get_singlecopy_bool_index(df_gc=df_gc, spp_filled=spp_filled)
    df_singleog = df_og.loc[is_singlecopy, spp_filled].copy()
    df_sog = df_singleog.copy()
    for species_name in spp_filled:
        if species_name not in uncorrected_by_species:
            continue
        df_sog = df_sog.merge(
            uncorrected_by_species[species_name],
            left_on=species_name,
            right_index=True,
            how='left',
            sort=False,
        )
    if len(spp_filled) > 0:
        df_sog = df_sog.iloc[:, len(spp_filled):].copy()
    df_sog.index = df_singleog.index
    df_sog = df_sog.apply(pandas.to_numeric, errors='coerce')
    return df_sog


def _get_library_sizes(df_nonzero, uncorrected_by_species):
    library_sizes = pandas.Series(index=df_nonzero.columns, dtype=float)
    for species_name, counts_df in uncorrected_by_species.items():
        del species_name
        for sample_id in counts_df.columns:
            if sample_id in library_sizes.index:
                library_sizes.loc[sample_id] = float(counts_df.loc[:, sample_id].sum(skipna=True))
    return library_sizes


def _build_norm_factor_frame(metadata_df, roundtrip):
    df_nf = pandas.DataFrame(
        {
            'sample_id': list(roundtrip.round2_factors.index),
            'tmm_library_size': roundtrip.library_sizes.reindex(roundtrip.round2_factors.index).to_numpy(dtype=float),
            'tmm_normalization_factor': roundtrip.round2_factors.to_numpy(dtype=float),
        }
    )
    return df_nf


def append_tmm_stats_to_metadata_python(metadata_df, roundtrip):
    required_cols = {'scientific_name', 'run', 'exclusion'}
    missing_cols = sorted(required_cols.difference(metadata_df.columns))
    if missing_cols:
        raise ValueError('Required metadata columns are missing for cstmm: {}'.format(', '.join(missing_cols)))
    df_metadata = metadata_df.copy()
    df_metadata['scientific_name'] = df_metadata['scientific_name'].astype(str).str.strip()
    df_metadata['run'] = df_metadata['run'].astype(str).str.strip()
    df_metadata['sample_id'] = [
        '{}_{}'.format(_normalize_species_prefix(scientific_name), run_id)
        for scientific_name, run_id in zip(df_metadata['scientific_name'], df_metadata['run'])
    ]
    df_nf = _build_norm_factor_frame(metadata_df=df_metadata, roundtrip=roundtrip)
    df_nf_keys = set(df_nf['sample_id'].astype(str))
    out_cols = list(df_metadata.columns) + ['tmm_library_size', 'tmm_normalization_factor']
    merged = df_metadata.merge(df_nf, on='sample_id', how='left', sort=False)
    exclusion_norm = merged['exclusion'].astype(str).str.strip().str.lower()
    is_retained = exclusion_norm.eq('no')
    if 'mapping_rate' in merged.columns:
        mapping_rate = pandas.to_numeric(merged['mapping_rate'], errors='coerce').fillna(-999)
        merged.loc[is_retained & mapping_rate.eq(0), 'exclusion'] = 'no_mapping'
    exclusion_norm = merged['exclusion'].astype(str).str.strip().str.lower()
    is_retained = exclusion_norm.eq('no')
    merged.loc[is_retained & ~merged['sample_id'].astype(str).isin(df_nf_keys), 'exclusion'] = 'no_cstmm_output'
    exclusion_norm = merged['exclusion'].astype(str).str.strip().str.lower()
    is_retained = exclusion_norm.eq('no')
    merged.loc[is_retained & merged['tmm_normalization_factor'].isna(), 'exclusion'] = 'cstmm_failed'
    merged = merged[out_cols].copy()
    merged = merged.drop(columns=['sample_id'])
    return merged


def save_corrected_output_files_python(uncorrected_by_species, roundtrip, dir_count, dir_cstmm):
    corrected = {}
    factors = roundtrip.round2_factors
    for species_name, dat in uncorrected_by_species.items():
        corrected_dat = dat.astype(float).copy()
        for sample_id, factor in factors.items():
            if sample_id in corrected_dat.columns:
                corrected_dat.loc[:, sample_id] = corrected_dat.loc[:, sample_id].to_numpy(dtype=float) / float(factor)
        dat_out = corrected_dat.copy()
        dat_out.insert(0, 'target_id', dat_out.index)
        dat_out.columns = [
            col if col == 'target_id' else re.sub(r'^{}_'.format(re.escape(species_name)), '', str(col))
            for col in dat_out.columns
        ]
        species_dir = os.path.join(dir_cstmm, species_name)
        os.makedirs(species_dir, exist_ok=True)
        _copy_eff_length_file(dir_count=dir_count, dir_cstmm=dir_cstmm, species_name=species_name)
        out_path = os.path.join(species_dir, '{}_cstmm_counts.tsv'.format(species_name))
        dat_out.to_csv(out_path, sep='\t', index=False)
        corrected[species_name] = dat_out.set_index('target_id')
    return corrected


def _category_colors(values):
    unique_values = [str(v) for v in pandas.Series(values).fillna('').astype(str).tolist()]
    unique_values = [v for idx, v in enumerate(unique_values) if v != '' and v not in unique_values[:idx]]
    cmap = plt.get_cmap('tab20')
    return {value: cmap(idx % 20) for idx, value in enumerate(unique_values)}


def save_norm_factor_histograms(df_metadata, dir_cstmm):
    tmp = df_metadata.loc[df_metadata['tmm_normalization_factor'].notna(), :].copy()
    if tmp.empty:
        return
    values = numpy.log2(tmp['tmm_normalization_factor'].to_numpy(dtype=float))
    values = values[numpy.isfinite(values)]
    if values.size == 0:
        return
    bins = numpy.histogram_bin_edges(values, bins=40)
    for fill_by in ['scientific_name', 'sample_group']:
        if fill_by not in tmp.columns:
            continue
        fig, ax = plt.subplots(figsize=(4.8, 2.4))
        colors = _category_colors(tmp[fill_by])
        bottom = numpy.zeros((len(bins) - 1,), dtype=float)
        for value, color in colors.items():
            group_values = numpy.log2(
                pandas.to_numeric(
                    tmp.loc[tmp[fill_by].astype(str) == value, 'tmm_normalization_factor'],
                    errors='coerce',
                ).dropna().to_numpy(dtype=float)
            )
            group_values = group_values[numpy.isfinite(group_values)]
            if group_values.size == 0:
                continue
            hist, edges = numpy.histogram(group_values, bins=bins)
            ax.bar(edges[:-1], hist, width=numpy.diff(edges), align='edge', bottom=bottom, label=value, color=color, alpha=0.75)
            bottom = bottom + hist
        ax.set_xlabel('log2(TMM normalization factor)')
        ax.set_ylabel('Count')
        ax.legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(dir_cstmm, 'cstmm_normalization_factor_histogram.{}.pdf'.format(fill_by)))
        plt.close(fig)


def save_norm_factor_scatter(df_metadata, dir_cstmm):
    tmp = df_metadata.loc[df_metadata['tmm_normalization_factor'].notna(), :].copy()
    if tmp.empty:
        return
    x = numpy.log10(pandas.to_numeric(tmp['tmm_library_size'], errors='coerce').to_numpy(dtype=float))
    y = numpy.log2(pandas.to_numeric(tmp['tmm_normalization_factor'], errors='coerce').to_numpy(dtype=float))
    finite = numpy.isfinite(x) & numpy.isfinite(y)
    if not finite.any():
        return
    tmp = tmp.loc[finite, :].copy()
    x = x[finite]
    y = y[finite]
    species_colors = _category_colors(tmp['scientific_name'])
    group_colors = _category_colors(tmp['sample_group']) if 'sample_group' in tmp.columns else {}
    fig, ax = plt.subplots(figsize=(4.8, 2.0))
    for idx, row in enumerate(tmp.itertuples(index=False)):
        face = species_colors.get(str(row.scientific_name), '#4c72b0')
        edge = group_colors.get(str(getattr(row, 'sample_group', '')), face)
        ax.scatter(x[idx], y[idx], s=36, facecolor=face, edgecolor=edge, linewidth=0.8, alpha=0.8)
    ax.set_xlabel('log10(Library size)')
    ax.set_ylabel('log2(TMM normalization factor)')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_cstmm, 'cstmm_normalization_factor_scatter.pdf'))
    plt.close(fig)


def save_mean_expression_boxplot_python(df_nonzero, roundtrip, uncorrected_by_species, corrected_by_species, dir_cstmm):
    df_nonzero_tmm = df_nonzero.astype(float).copy()
    for sample_id, factor in roundtrip.round2_factors.items():
        if sample_id in df_nonzero_tmm.columns:
            df_nonzero_tmm.loc[:, sample_id] = df_nonzero_tmm.loc[:, sample_id].to_numpy(dtype=float) / float(factor)
    mean_before = df_nonzero.mean(axis=0, skipna=True).to_numpy(dtype=float)
    mean_after = df_nonzero_tmm.mean(axis=0, skipna=True).to_numpy(dtype=float)
    mean_sra_before = []
    mean_sra_after = []
    for species_name in corrected_by_species:
        mean_sra_before.extend(uncorrected_by_species[species_name].mean(axis=0, skipna=True).tolist())
        mean_sra_after.extend(corrected_by_species[species_name].mean(axis=0, skipna=True).tolist())
    fig, axes = plt.subplots(1, 2, figsize=(3.6, 3.6))
    boxplot_labels = ['Raw\ncounts', 'TMM-\ncorrected\ncounts']
    axes[0].boxplot([mean_before, mean_after])
    axes[0].set_xticks([1, 2])
    axes[0].set_xticklabels(boxplot_labels)
    axes[0].set_ylabel('Mean count of single-copy genes')
    axes[1].boxplot([mean_sra_before, mean_sra_after])
    axes[1].set_xticks([1, 2])
    axes[1].set_xticklabels(boxplot_labels)
    axes[1].set_ylabel('Mean count of all genes')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_cstmm, 'cstmm_mean_expression_boxplot.pdf'))
    plt.close(fig)


def save_exclusion_plot_python(df_metadata, out_path, y_label='Sample count'):
    required_cols = {'scientific_name', 'exclusion'}
    if not required_cols.issubset(df_metadata.columns):
        return
    tmp = df_metadata.loc[:, ['scientific_name', 'exclusion']].copy()
    tmp['scientific_name'] = tmp['scientific_name'].astype(str).str.strip()
    tmp['exclusion'] = tmp['exclusion'].astype(str).str.strip()
    tmp = tmp.loc[(tmp['scientific_name'] != '') & (tmp['exclusion'] != ''), :]
    if tmp.empty:
        return
    summary = (
        tmp.groupby(['scientific_name', 'exclusion'], sort=False)
        .size()
        .rename('count')
        .reset_index()
    )
    pivot = summary.pivot(index='scientific_name', columns='exclusion', values='count').fillna(0)
    fig_width = max(3.6, 0.11 * pivot.shape[0])
    fig, ax = plt.subplots(figsize=(fig_width, 3.6))
    bottom = numpy.zeros((pivot.shape[0],), dtype=float)
    colors = _category_colors(pivot.columns)
    for exclusion_value in pivot.columns:
        heights = pivot.loc[:, exclusion_value].to_numpy(dtype=float)
        ax.bar(pivot.index.tolist(), heights, bottom=bottom, label=str(exclusion_value), color=colors[str(exclusion_value)])
        bottom = bottom + heights
    ax.set_ylabel(y_label)
    ax.set_xlabel('')
    ax.tick_params(axis='x', rotation=90)
    ax.legend(frameon=False, fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _run_cstmm_python(uncorrected_by_species, df_sog, dir_count, dir_cstmm):
    df_nonzero = _get_df_nonzero(df_sog)
    library_sizes = _get_library_sizes(df_nonzero=df_nonzero, uncorrected_by_species=uncorrected_by_species)
    roundtrip = run_tmm_rounds_for_cstmm(counts=df_nonzero, lib_size=library_sizes.reindex(df_nonzero.columns))
    metadata_path = os.path.join(dir_count, 'metadata.tsv')
    df_metadata = pandas.read_csv(metadata_path, sep='\t')
    df_metadata = append_tmm_stats_to_metadata_python(metadata_df=df_metadata, roundtrip=roundtrip)
    df_metadata = df_metadata.loc[:, ~pandas.Index(df_metadata.columns).astype(str).str.startswith('Unnamed')]
    os.makedirs(dir_cstmm, exist_ok=True)
    df_metadata.to_csv(os.path.join(dir_cstmm, 'metadata.tsv'), sep='\t', index=False)
    corrected = save_corrected_output_files_python(
        uncorrected_by_species=uncorrected_by_species,
        roundtrip=roundtrip,
        dir_count=dir_count,
        dir_cstmm=dir_cstmm,
    )
    save_norm_factor_histograms(df_metadata=df_metadata, dir_cstmm=dir_cstmm)
    save_norm_factor_scatter(df_metadata=df_metadata, dir_cstmm=dir_cstmm)
    save_mean_expression_boxplot_python(
        df_nonzero=df_nonzero,
        roundtrip=roundtrip,
        uncorrected_by_species=uncorrected_by_species,
        corrected_by_species=corrected,
        dir_cstmm=dir_cstmm,
    )
    save_exclusion_plot_python(df_metadata=df_metadata, out_path=os.path.join(dir_cstmm, 'cstmm_exclusion.pdf'))
    return roundtrip


def run_cstmm_python_single_species(dir_count, dir_cstmm, species_name):
    uncorrected = {species_name: _read_est_counts(dir_count=dir_count, species_name=species_name)}
    return _run_cstmm_python(
        uncorrected_by_species=uncorrected,
        df_sog=uncorrected[species_name],
        dir_count=dir_count,
        dir_cstmm=dir_cstmm,
    )


def run_cstmm_python_multi_species(dir_count, dir_cstmm, file_genecount, file_orthogroup_table):
    df_gc = _read_genecount_table(file_genecount=file_genecount)
    species_names = _get_spp_filled(dir_count=dir_count, df_gc=df_gc)
    uncorrected = _load_uncorrected(dir_count=dir_count, species_names=species_names)
    df_sog = _get_df_exp_single_copy_ortholog(
        file_genecount=file_genecount,
        file_orthogroup_table=file_orthogroup_table,
        dir_count=dir_count,
        uncorrected_by_species=uncorrected,
    )
    return _run_cstmm_python(
        uncorrected_by_species=uncorrected,
        df_sog=df_sog,
        dir_count=dir_count,
        dir_cstmm=dir_cstmm,
    )
