import os
import re
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import pandas
from amalgkit.table_io import read_annotation_tsv, read_identifier_tsv

from amalgkit.imputation import impute_expression
from amalgkit.normalization_tmm import run_tmm_rounds_for_cstmm
from amalgkit.orthology_utils import (
    DEFAULT_SINGLE_COPY_THRESHOLD,
    get_single_copy_orthogroup_mask,
    validate_single_copy_threshold,
)
from amalgkit.runtime_utils import build_species_token_map


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
    dat = read_identifier_tsv(infile_path, index_col=0)
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
    df_gc = read_identifier_tsv(file_genecount, identifier_columns=('orthogroup_id',))
    if 'orthogroup_id' not in df_gc.columns:
        raise ValueError('Column "orthogroup_id" is required in genecount table: {}'.format(file_genecount))
    df_gc = df_gc.set_index('orthogroup_id')
    return df_gc


def _read_orthogroup_table(file_orthogroup_table):
    frame = read_annotation_tsv(file_orthogroup_table)
    return frame.set_index(frame.columns[0])


def _load_uncorrected(dir_count, species_names):
    return {
        species_name: _read_est_counts(dir_count=dir_count, species_name=species_name)
        for species_name in species_names
    }


def _index_cstmm_metadata(metadata_df):
    required_cols = {'scientific_name', 'run', 'exclusion'}
    missing_cols = sorted(required_cols.difference(metadata_df.columns))
    if missing_cols:
        raise ValueError('Required metadata columns are missing for cstmm: {}'.format(', '.join(missing_cols)))
    metadata = metadata_df.copy()
    for column in ('scientific_name', 'run'):
        metadata[column] = metadata[column].fillna('').astype(str).str.strip()
        if metadata[column].eq('').any():
            raise ValueError('cstmm metadata contains an empty {}.'.format(column))
    token_by_species = build_species_token_map(
        scientific_names=metadata['scientific_name'].tolist(),
        explicit_tokens=(metadata['species_token'].fillna('').astype(str).str.strip().tolist()
                         if 'species_token' in metadata.columns else None),
        context='cstmm metadata',
    )
    metadata['sample_id'] = [
        '{}_{}'.format(token_by_species[species], run)
        for species, run in zip(metadata['scientific_name'], metadata['run'])
    ]
    duplicated = metadata['sample_id'].duplicated(keep=False)
    if duplicated.any():
        raise ValueError('Duplicate cstmm species/run identifiers: {}'.format(
            ', '.join(metadata.loc[duplicated, 'sample_id'].drop_duplicates())))
    return metadata


def _select_cstmm_inputs(uncorrected_by_species, metadata_path):
    metadata = _index_cstmm_metadata(read_annotation_tsv(metadata_path))
    count_ids = {column for counts in uncorrected_by_species.values() for column in counts.columns}
    unknown = sorted(count_ids.difference(metadata['sample_id']))
    if unknown:
        raise ValueError(
            'Count columns have no matching species/run in cstmm metadata: {}. '
            'Keep their metadata rows and set exclusion to exclude runs.'.format(', '.join(unknown)))
    eligible = metadata['exclusion'].fillna('').astype(str).str.strip().str.lower().eq('no')
    if 'is_sampled' in metadata.columns:
        sampled = metadata['is_sampled'].fillna('').astype(str).str.strip().str.lower()
        if sampled.ne('').any():
            invalid = sorted(set(sampled).difference({'', 'yes', 'no'}))
            if invalid:
                raise ValueError('Column "is_sampled" contains invalid flag(s): {}'.format(', '.join(invalid)))
            eligible &= sampled.eq('yes')
    selected_ids = set(metadata.loc[eligible, 'sample_id'])
    selected = {}
    for species, counts in uncorrected_by_species.items():
        columns = [column for column in counts.columns if column in selected_ids]
        if columns:
            selected[species] = counts if len(columns) == counts.shape[1] else counts.loc[:, columns].copy()
    if not selected:
        raise ValueError('No eligible cstmm count columns remained after exclusion/is_sampled filtering.')
    print('cstmm metadata: {} ({} eligible count columns)'.format(
        metadata_path, sum(counts.shape[1] for counts in selected.values())), flush=True)
    return selected, metadata.drop(columns='sample_id')


def _copy_eff_length_file(dir_count, dir_cstmm, species_name):
    species_dir = os.path.join(dir_count, species_name)
    eff_length_files = _list_matching_files(species_dir, r'.*eff_length\.tsv$')
    if len(eff_length_files) != 1:
        return
    src = os.path.join(species_dir, eff_length_files[0])
    dst = os.path.join(dir_cstmm, species_name, eff_length_files[0])
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)


def _copy_quant_model_file(dir_count, dir_cstmm, species_name):
    species_dir = os.path.join(dir_count, species_name)
    quant_model_files = _list_matching_files(species_dir, r'.*quant_model\.tsv$')
    if len(quant_model_files) == 0:
        return None
    if len(quant_model_files) > 1:
        raise ValueError('Multiple *quant_model.tsv files found: {}'.format(species_name))
    src = os.path.join(species_dir, quant_model_files[0])
    dst = os.path.join(dir_cstmm, species_name, '{}_quant_model.tsv'.format(species_name))
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)
    return dst


def _get_df_nonzero(df_counts):
    is_zero_col = (df_counts.sum(axis=0, skipna=True) == 0)
    df_nonzero = df_counts.loc[:, ~is_zero_col].copy()
    return impute_expression(
        df_nonzero,
        strategy='em_pca',
        minimum_imputed_value=0.0,
    )


def _get_singlecopy_bool_index(df_gc, spp_filled, percent_singlecopy_threshold=50.0):
    return get_single_copy_orthogroup_mask(
        df_genecount=df_gc,
        species=spp_filled,
        single_copy_threshold=percent_singlecopy_threshold,
    )


def _get_df_exp_single_copy_ortholog(
    file_genecount,
    file_orthogroup_table,
    dir_count,
    uncorrected_by_species,
    single_copy_threshold=DEFAULT_SINGLE_COPY_THRESHOLD,
):
    df_gc = _read_genecount_table(file_genecount=file_genecount)
    df_og = _read_orthogroup_table(file_orthogroup_table=file_orthogroup_table)
    spp_filled = [species for species in uncorrected_by_species if species in df_gc.columns]
    is_singlecopy = _get_singlecopy_bool_index(
        df_gc=df_gc,
        spp_filled=spp_filled,
        percent_singlecopy_threshold=single_copy_threshold,
    )
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
    sample_ids = list(roundtrip.round2_factors.index)
    library_sizes = roundtrip.library_sizes.reindex(sample_ids).to_numpy(dtype=float)
    normalization_factors = roundtrip.round2_factors.reindex(sample_ids).to_numpy(dtype=float)
    df_nf = pandas.DataFrame(
        {
            'sample_id': sample_ids,
            'tmm_library_size': library_sizes,
            'tmm_normalization_factor': normalization_factors,
            'tmm_effective_library_size': library_sizes * normalization_factors,
        }
    )
    return df_nf


def append_tmm_stats_to_metadata_python(metadata_df, roundtrip):
    stat_columns = ['tmm_library_size', 'tmm_normalization_factor', 'tmm_effective_library_size']
    df_metadata = _index_cstmm_metadata(metadata_df.drop(columns=stat_columns, errors='ignore'))
    df_nf = _build_norm_factor_frame(metadata_df=df_metadata, roundtrip=roundtrip)
    df_nf_keys = set(df_nf['sample_id'].astype(str))
    out_cols = list(df_metadata.columns) + stat_columns
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
        _copy_quant_model_file(dir_count=dir_count, dir_cstmm=dir_cstmm, species_name=species_name)
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


def _run_cstmm_python(
    uncorrected_by_species,
    df_sog,
    dir_count,
    dir_cstmm,
    metadata_df,
    single_copy_threshold=None,
):
    df_nonzero = _get_df_nonzero(df_sog)
    library_sizes = _get_library_sizes(df_nonzero=df_nonzero, uncorrected_by_species=uncorrected_by_species)
    roundtrip = run_tmm_rounds_for_cstmm(counts=df_nonzero, lib_size=library_sizes.reindex(df_nonzero.columns))
    df_metadata = append_tmm_stats_to_metadata_python(metadata_df=metadata_df, roundtrip=roundtrip)
    if single_copy_threshold is not None:
        df_metadata['single_copy_threshold'] = validate_single_copy_threshold(single_copy_threshold)
    else:
        df_metadata = df_metadata.drop(columns='single_copy_threshold', errors='ignore')
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


def run_cstmm_python_single_species(dir_count, dir_cstmm, species_name, metadata_path=None):
    uncorrected = {species_name: _read_est_counts(dir_count=dir_count, species_name=species_name)}
    uncorrected, metadata = _select_cstmm_inputs(
        uncorrected, metadata_path if metadata_path is not None else os.path.join(dir_count, 'metadata.tsv'))
    return _run_cstmm_python(
        uncorrected_by_species=uncorrected,
        df_sog=uncorrected[species_name],
        dir_count=dir_count,
        dir_cstmm=dir_cstmm,
        metadata_df=metadata,
    )


def run_cstmm_python_multi_species(
    dir_count,
    dir_cstmm,
    file_genecount,
    file_orthogroup_table,
    single_copy_threshold=DEFAULT_SINGLE_COPY_THRESHOLD,
    metadata_path=None,
):
    single_copy_threshold = validate_single_copy_threshold(single_copy_threshold)
    df_gc = _read_genecount_table(file_genecount=file_genecount)
    species_names = _get_spp_filled(dir_count=dir_count, df_gc=df_gc)
    uncorrected = _load_uncorrected(dir_count=dir_count, species_names=species_names)
    uncorrected, metadata = _select_cstmm_inputs(
        uncorrected, metadata_path if metadata_path is not None else os.path.join(dir_count, 'metadata.tsv'))
    if len(uncorrected) == 1:
        print('Only one species remains after metadata selection; applying standard TMM.', flush=True)
        df_sog = next(iter(uncorrected.values()))
        single_copy_threshold = None
    else:
        df_sog = _get_df_exp_single_copy_ortholog(
            file_genecount=file_genecount,
            file_orthogroup_table=file_orthogroup_table,
            dir_count=dir_count,
            uncorrected_by_species=uncorrected,
            single_copy_threshold=single_copy_threshold,
        )
    return _run_cstmm_python(
        uncorrected_by_species=uncorrected,
        df_sog=df_sog,
        dir_count=dir_count,
        dir_cstmm=dir_cstmm,
        metadata_df=metadata,
        single_copy_threshold=single_copy_threshold,
    )
