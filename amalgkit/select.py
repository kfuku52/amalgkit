import datetime
import os
import re
import shutil

import pandas

from amalgkit.arg_utils import clone_namespace
from amalgkit.filter_utils import staged_output_dir
from amalgkit.metadata_utils import Metadata, load_metadata
from amalgkit.output_utils import atomic_write_dataframe
from amalgkit.runtime_utils import check_config_dir

SELECT_PRIORITY_COLUMNS = [
    'sample_attribute_tissue',
    'tissue',
    'sample_group',
    'sample_title',
    'source_name',
    'sample_description',
    'exp_title',
    'design',
]

SELECT_IGNORE_VALUES = {
    '',
    'na',
    'n/a',
    'none',
    'not applicable',
    'not_applicable',
    'unknown',
}

SELECT_ORGAN_PATTERNS = {
    'flower': re.compile(
        r"\b(?:flower|flowers|floral|petal|petals|corolla|corollas|anther|anthers|ovary|ovaries|nectary|nectaries|"
        r"style|styles|stigma|stigmata|pistil|pistils|pistillate|stamen|stamens|catkin|catkins|"
        r"inflorescence|inflorescences|flower[\s_-]?bud|flower[\s_-]?buds)\b",
        re.IGNORECASE,
    ),
    'leaf': re.compile(
        r"\b(?:leaf|leaves|foliar|foliage|petiole|petioles|blade|lamina|trifoliate|leave)\b",
        re.IGNORECASE,
    ),
    'root': re.compile(
        r"\b(?:root|roots|radicle|radicles|root[\s_-]?tip|root[\s_-]?tips|primary root|root hair)\b",
        re.IGNORECASE,
    ),
}

SELECT_REVIEW_PATTERN = re.compile(
    r"\b(?:hairy root|hairy roots|hairy root culture|nodule|nodules|shoot|shoots|stem|stems|meristem|"
    r"whole plant|whole seedling|seedling|pseudostem|axillary bud|vegetative bud|hypocotyl|internode|"
    r"callus|embryo|embryos|cotyledon|cotyledons|pith|aerial part|aerial parts|young shoot|young shoots|"
    r"rootstock|inflorescence meristem|inflorescence bud)\b",
    re.IGNORECASE,
)

SELECT_NONTARGET_PATTERN = re.compile(
    r"\b(?:fruit|fruits|receptacle|receptacles|flesh|pulp|seed|seeds|testa|pod|pods|silique|ovule|ovules)\b",
    re.IGNORECASE,
)


def resolve_select_config_dir(args):
    if args.config_dir == 'inferred':
        return os.path.join(args.out_dir, 'config')
    return args.config_dir


def resolve_select_metadata_specieswise_dir(args, out_dir):
    if getattr(args, 'metadata_specieswise_dir', 'inferred') != 'inferred':
        return os.path.realpath(args.metadata_specieswise_dir)
    return os.path.join(os.path.dirname(out_dir), 'metadata_specieswise')


def resolve_select_batch_output_path(args, out_dir, arg_name, default_filename):
    value = getattr(args, arg_name, 'inferred')
    if value == 'inferred':
        return os.path.join(out_dir, default_filename)
    return os.path.realpath(value)


def resolve_select_batch_label(args, out_dir):
    value = getattr(args, 'batch_label', 'inferred')
    if value != 'inferred':
        return value
    return os.path.basename(out_dir.rstrip(os.sep)) or 'select_batch'


def filter_metadata_by_sample_group(metadata, sample_group_arg):
    if sample_group_arg is None:
        return metadata
    txt = '{}: Extracting pre-selected sample_group entries: {}'
    print(txt.format(datetime.datetime.now(), sample_group_arg), flush=True)
    if 'sample_group' not in metadata.df.columns:
        raise ValueError(
            'Column "sample_group" is required in metadata when --sample_group is specified.'
        )
    selected_sample_groups = [token.strip() for token in re.split(r'[,\|]+', str(sample_group_arg))]
    selected_sample_groups = [token for token in selected_sample_groups if token != '']
    if len(selected_sample_groups) == 0:
        raise ValueError('No sample_group was selected. Provide non-empty --sample_group value.')
    sample_groups = metadata.df['sample_group'].fillna('').astype(str).str.strip()
    metadata.df = metadata.df.loc[sample_groups.isin(selected_sample_groups), :].reset_index(drop=True)
    return metadata


def apply_select_filters(metadata, args, dir_config):
    metadata.nspot_cutoff(args.min_nspots)
    metadata.mark_missing_rank(args.mark_missing_rank)
    metadata.mark_redundant_biosample(args.mark_redundant_biosamples)
    metadata.remove_specialchars()
    dropped_group_columns = metadata.group_attributes(dir_config, drop_source_columns=False)
    metadata.mark_exclude_keywords(dir_config)
    metadata.mark_treatment_terms(dir_config)
    if len(dropped_group_columns) > 0:
        metadata.df = metadata.df.drop(columns=[col for col in dropped_group_columns if col in metadata.df.columns])
    metadata.label_sampled_data(args.max_sample)
    metadata.reorder(omit_misc=False)
    return metadata


def write_select_outputs(path_metadata_original, path_metadata_table, metadata_dir, metadata, metadata_original_df=None):
    if os.path.exists(metadata_dir) and (not os.path.isdir(metadata_dir)):
        raise NotADirectoryError('Output metadata path exists but is not a directory: {}'.format(metadata_dir))
    for output_path in [path_metadata_original, path_metadata_table]:
        if os.path.exists(output_path) and (not os.path.isfile(output_path)):
            raise NotADirectoryError('Output path exists but is not a file: {}'.format(output_path))
    if metadata_original_df is None:
        metadata_original_df = metadata.df
    if os.path.exists(path_metadata_original):
        print('Refreshing original metadata copy at: {}'.format(path_metadata_original), flush=True)
    else:
        print('Creating original metadata copy at: {}'.format(path_metadata_original), flush=True)
    print('Updating metadata table at: {}'.format(path_metadata_table), flush=True)
    with staged_output_dir(metadata_dir, redo=True, prefix='amalgkit_select_stage_') as stage_dir:
        if os.path.isdir(metadata_dir):
            shutil.copytree(metadata_dir, stage_dir, dirs_exist_ok=True)
        stage_original = os.path.join(stage_dir, os.path.basename(path_metadata_original))
        stage_table = os.path.join(stage_dir, os.path.basename(path_metadata_table))
        atomic_write_dataframe(metadata_original_df, stage_original, sep='\t', index=False)
        atomic_write_dataframe(metadata.df, stage_table, sep='\t', index=False)
        pivot_specs = [
            ('pivot_qualified.tsv', False),
            ('pivot_selected.tsv', True),
        ]
        for filename, sampled_only in pivot_specs:
            pivot_df = metadata.pivot(n_sp_cutoff=0, qualified_only=True, sampled_only=sampled_only)
            atomic_write_dataframe(pivot_df, os.path.join(stage_dir, filename), sep='\t')


def normalize_select_value(value):
    text = '' if pandas.isna(value) else str(value)
    text = text.replace('_', ' ')
    text = re.sub(r'\s+', ' ', text).strip()
    return text


def classify_select_text(text):
    normalized = normalize_select_value(text)
    if normalized.lower() in SELECT_IGNORE_VALUES:
        return {'status': 'empty', 'organ': '', 'text': normalized}
    organ_hits = [
        organ
        for organ, pattern in SELECT_ORGAN_PATTERNS.items()
        if pattern.search(normalized)
    ]
    if len(organ_hits) > 1:
        return {'status': 'mixed', 'organ': '', 'text': normalized}
    if SELECT_REVIEW_PATTERN.search(normalized):
        return {'status': 'review', 'organ': '', 'text': normalized}
    if len(organ_hits) == 1:
        return {'status': 'organ', 'organ': organ_hits[0], 'text': normalized}
    if SELECT_NONTARGET_PATTERN.search(normalized):
        return {'status': 'non_target', 'organ': '', 'text': normalized}
    return {'status': 'unknown', 'organ': '', 'text': normalized}


def normalize_select_row(row):
    original_sample_group = normalize_select_value(row.get('sample_group', ''))
    result = {
        'sample_group_original': original_sample_group,
        'sample_group_normalization_status': 'unchanged',
        'sample_group_normalization_source': '',
        'sample_group_normalization_text': '',
    }
    updated_sample_group = original_sample_group
    fallback_classification = None
    fallback_source = ''
    for column in SELECT_PRIORITY_COLUMNS:
        if column not in row.index:
            continue
        classification = classify_select_text(row[column])
        status = classification['status']
        if status == 'empty':
            continue
        if status == 'unknown':
            if fallback_classification is None:
                fallback_classification = classification
                fallback_source = column
            continue
        result['sample_group_normalization_status'] = status
        result['sample_group_normalization_source'] = column
        result['sample_group_normalization_text'] = classification['text']
        if status == 'organ':
            updated_sample_group = classification['organ']
        return updated_sample_group, result
    if fallback_classification is not None:
        result['sample_group_normalization_status'] = fallback_classification['status']
        result['sample_group_normalization_source'] = fallback_source
        result['sample_group_normalization_text'] = fallback_classification['text']
    return updated_sample_group, result


def normalize_select_metadata_frame(df):
    normalized_groups = []
    statuses = []
    sources = []
    source_texts = []
    original_groups = []
    for _, row in df.iterrows():
        updated_sample_group, result = normalize_select_row(row)
        normalized_groups.append(updated_sample_group)
        statuses.append(result['sample_group_normalization_status'])
        sources.append(result['sample_group_normalization_source'])
        source_texts.append(result['sample_group_normalization_text'])
        original_groups.append(result['sample_group_original'])
    out_df = df.copy(deep=True)
    out_df['sample_group_original'] = original_groups
    out_df['sample_group'] = normalized_groups
    out_df['sample_group_normalization_status'] = statuses
    out_df['sample_group_normalization_source'] = sources
    out_df['sample_group_normalization_text'] = source_texts
    return out_df


def summarize_select_normalization(species_name, species_token, df, workspace_dir, metadata_path):
    sample_group_series = df['sample_group'].fillna('').astype(str).str.strip()
    status_series = df['sample_group_normalization_status'].fillna('').astype(str).str.strip()
    return {
        'scientific_name': species_name,
        'species_token': species_token,
        'rows_total': int(df.shape[0]),
        'rows_flower': int((sample_group_series == 'flower').sum()),
        'rows_leaf': int((sample_group_series == 'leaf').sum()),
        'rows_root': int((sample_group_series == 'root').sum()),
        'rows_review': int((status_series == 'review').sum()),
        'rows_mixed': int((status_series == 'mixed').sum()),
        'rows_non_target': int((status_series == 'non_target').sum()),
        'rows_unknown': int((status_series == 'unknown').sum()),
        'rows_unchanged': int((status_series == 'unchanged').sum()),
        'workspace_dir': workspace_dir,
        'metadata_path': metadata_path,
    }


def summarize_selected_metadata(species_name, species_token, metadata, selected_metadata_path, max_sample):
    df = metadata.df
    exclusion_series = df['exclusion'].fillna('').astype(str).str.strip() if 'exclusion' in df.columns else pandas.Series('', index=df.index)
    sample_group_series = df['sample_group'].fillna('').astype(str).str.strip() if 'sample_group' in df.columns else pandas.Series('', index=df.index)
    qualified_series = df['is_qualified'].fillna('').astype(str).str.strip() if 'is_qualified' in df.columns else pandas.Series('', index=df.index)
    sampled_series = df['is_sampled'].fillna('').astype(str).str.strip() if 'is_sampled' in df.columns else pandas.Series('', index=df.index)
    row = {
        'species_token': species_token,
        'scientific_name': species_name,
        'rows_after_select': int(df.shape[0]),
        'excluded_non_control': int((exclusion_series == 'non_control').sum()),
        'excluded_low_nspots': int((exclusion_series == 'low_nspots').sum()),
        'excluded_misc': int(((exclusion_series != '') & (exclusion_series != 'no') & (exclusion_series != 'non_control') & (exclusion_series != 'low_nspots')).sum()),
    }
    for organ in ['flower', 'leaf', 'root']:
        organ_mask = (sample_group_series == organ)
        row['{}_rows'.format(organ)] = int(organ_mask.sum())
        row['{}_qualified_yes'.format(organ)] = int((organ_mask & (qualified_series == 'yes')).sum())
        row['{}_sampled_yes'.format(organ)] = int((organ_mask & (sampled_series == 'yes')).sum())
    row['qualified_yes_total'] = int((qualified_series == 'yes').sum())
    row['sampled_yes_total'] = int((sampled_series == 'yes').sum())
    strict_threshold, relaxed_threshold = resolve_select_queue_thresholds(max_sample)
    row['strict_ready'] = all(row['{}_sampled_yes'.format(organ)] >= strict_threshold for organ in ['flower', 'leaf', 'root'])
    row['relaxed_ready'] = all(row['{}_sampled_yes'.format(organ)] >= relaxed_threshold for organ in ['flower', 'leaf', 'root'])
    row['queue_tier'] = 'strict' if row['strict_ready'] else 'relaxed' if row['relaxed_ready'] else 'defer'
    row['selected_metadata_path'] = selected_metadata_path
    return row


def resolve_select_queue_thresholds(max_sample):
    max_sample_int = int(max_sample)
    return min(max_sample_int, 30), min(max_sample_int, 20)


def build_select_queue_df(summary_df):
    queue_df = summary_df.loc[:, [
        'species_token',
        'scientific_name',
        'rows_after_select',
        'flower_sampled_yes',
        'leaf_sampled_yes',
        'root_sampled_yes',
        'qualified_yes_total',
        'sampled_yes_total',
        'excluded_non_control',
        'excluded_low_nspots',
        'excluded_misc',
        'strict_ready',
        'relaxed_ready',
        'queue_tier',
        'selected_metadata_path',
    ]].copy(deep=True)
    queue_df['queue_sort'] = queue_df['queue_tier'].map({'strict': 0, 'relaxed': 1, 'defer': 2}).fillna(9)
    queue_df = queue_df.sort_values(
        ['queue_sort', 'sampled_yes_total', 'qualified_yes_total', 'rows_after_select'],
        ascending=[True, False, False, False],
    ).drop(columns=['queue_sort'])
    return queue_df.reset_index(drop=True)


def build_select_manifest_df(queue_df, batch_label):
    manifest_df = queue_df.copy(deep=True)
    manifest_df.insert(
        0,
        'recommended_action',
        manifest_df['queue_tier'].map({
            'strict': 'external_run_now',
            'relaxed': 'external_run_if_capacity',
            'defer': 'defer_rule_refinement',
        }),
    )
    manifest_df.insert(1, 'pilot_batch', batch_label)
    manifest_df.insert(0, 'queue_order', range(1, manifest_df.shape[0] + 1))
    return manifest_df


def manifest_sidecar_paths(manifest_tsv):
    base, ext = os.path.splitext(manifest_tsv)
    if ext == '':
        ext = '.tsv'
    return (
        base + '_strict' + ext,
        base + '_relaxed' + ext,
    )


def load_select_species_table(path_species_tsv):
    species_df = pandas.read_csv(path_species_tsv, sep='\t', dtype=str, keep_default_na=False).fillna('')
    if 'scientific_name' not in species_df.columns:
        raise ValueError('--species_tsv must contain a scientific_name column.')
    species_df['scientific_name'] = species_df['scientific_name'].astype(str).str.strip()
    species_df = species_df.loc[species_df['scientific_name'] != '', :].copy()
    if 'species_token' in species_df.columns:
        species_df['species_token'] = species_df['species_token'].astype(str).str.strip()
    else:
        species_df['species_token'] = ''
    inferred_tokens = species_df['scientific_name'].str.replace(r'\s+', '_', regex=True)
    species_df.loc[species_df['species_token'] == '', 'species_token'] = inferred_tokens.loc[species_df['species_token'] == '']
    if species_df['species_token'].duplicated().any():
        duplicated = sorted(species_df.loc[species_df['species_token'].duplicated(), 'species_token'].unique().tolist())
        raise ValueError('Duplicate species_token detected in --species_tsv: {}'.format(', '.join(duplicated)))
    return species_df.reset_index(drop=True)


def select_batch_main(args):
    out_dir = os.path.realpath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    dir_config = resolve_select_config_dir(args)
    check_config_dir(dir_path=dir_config, mode='select')
    metadata_specieswise_dir = resolve_select_metadata_specieswise_dir(args, out_dir)
    if not os.path.isdir(metadata_specieswise_dir):
        raise FileNotFoundError(
            'metadata_specieswise directory not found: {}'.format(metadata_specieswise_dir)
        )
    species_df = load_select_species_table(os.path.realpath(args.species_tsv))
    summary_tsv = resolve_select_batch_output_path(args, out_dir, 'summary_tsv', 'select_summary.tsv')
    queue_tsv = resolve_select_batch_output_path(args, out_dir, 'queue_tsv', 'select_queue.tsv')
    manifest_tsv = resolve_select_batch_output_path(args, out_dir, 'manifest_tsv', 'external_manifest.tsv')
    normalization_summary_tsv = os.path.join(out_dir, 'normalization_summary.tsv')
    batch_label = resolve_select_batch_label(args, out_dir)
    normalization_rows = []
    summary_rows = []
    print('select batch: processing {} species'.format(species_df.shape[0]), flush=True)
    for index, row in species_df.iterrows():
        species_name = row['scientific_name']
        species_token = row['species_token']
        print('[{}/{}] {}'.format(index + 1, species_df.shape[0], species_token), flush=True)
        merged_metadata_path = os.path.join(
            metadata_specieswise_dir,
            species_token,
            '{}.metadata.tsv'.format(species_token),
        )
        if not os.path.exists(merged_metadata_path):
            raise FileNotFoundError('Merged metadata not found: {}'.format(merged_metadata_path))
        merged_df = pandas.read_csv(
            merged_metadata_path,
            sep='\t',
            dtype=str,
            keep_default_na=False,
            low_memory=False,
        ).fillna('')
        normalized_df = normalize_select_metadata_frame(merged_df)
        species_out_dir = os.path.join(out_dir, species_token)
        metadata_dir = os.path.join(species_out_dir, 'metadata')
        normalized_metadata_path = os.path.join(metadata_dir, 'metadata.tsv')
        normalization_rows.append(
            summarize_select_normalization(
                species_name=species_name,
                species_token=species_token,
                df=normalized_df,
                workspace_dir=species_out_dir,
                metadata_path=normalized_metadata_path,
            )
        )
        runtime_args = clone_namespace(args, out_dir=species_out_dir)
        metadata = Metadata.from_DataFrame(normalized_df)
        metadata_original_df = metadata.df.copy(deep=True)
        metadata = filter_metadata_by_sample_group(metadata, runtime_args.sample_group)
        metadata = apply_select_filters(metadata, runtime_args, dir_config)
        path_metadata_original = os.path.join(metadata_dir, 'metadata_original.tsv')
        path_metadata_table = os.path.join(metadata_dir, 'metadata.tsv')
        write_select_outputs(
            path_metadata_original=path_metadata_original,
            path_metadata_table=path_metadata_table,
            metadata_dir=metadata_dir,
            metadata=metadata,
            metadata_original_df=metadata_original_df,
        )
        summary_rows.append(
            summarize_selected_metadata(
                species_name=species_name,
                species_token=species_token,
                metadata=metadata,
                selected_metadata_path=path_metadata_table,
                max_sample=runtime_args.max_sample,
            )
        )
    normalization_df = pandas.DataFrame(normalization_rows)
    summary_df = pandas.DataFrame(summary_rows)
    if not summary_df.empty:
        summary_df = summary_df.sort_values(
            ['sampled_yes_total', 'qualified_yes_total', 'rows_after_select'],
            ascending=[False, False, False],
        ).reset_index(drop=True)
    atomic_write_dataframe(normalization_df, normalization_summary_tsv, sep='\t', index=False)
    atomic_write_dataframe(
        summary_df.loc[:, [
            'species_token',
            'scientific_name',
            'rows_after_select',
            'excluded_non_control',
            'excluded_low_nspots',
            'excluded_misc',
            'flower_rows',
            'flower_qualified_yes',
            'flower_sampled_yes',
            'leaf_rows',
            'leaf_qualified_yes',
            'leaf_sampled_yes',
            'root_rows',
            'root_qualified_yes',
            'root_sampled_yes',
            'qualified_yes_total',
            'sampled_yes_total',
            'strict_ready',
            'relaxed_ready',
            'queue_tier',
        ]],
        summary_tsv,
        sep='\t',
        index=False,
    )
    queue_df = build_select_queue_df(summary_df)
    atomic_write_dataframe(queue_df, queue_tsv, sep='\t', index=False)
    manifest_df = build_select_manifest_df(queue_df, batch_label=batch_label)
    atomic_write_dataframe(manifest_df, manifest_tsv, sep='\t', index=False)
    strict_manifest_tsv, relaxed_manifest_tsv = manifest_sidecar_paths(manifest_tsv)
    atomic_write_dataframe(manifest_df.loc[manifest_df['queue_tier'] == 'strict', :], strict_manifest_tsv, sep='\t', index=False)
    atomic_write_dataframe(manifest_df.loc[manifest_df['queue_tier'] == 'relaxed', :], relaxed_manifest_tsv, sep='\t', index=False)
    print('select batch complete: species_processed={}'.format(species_df.shape[0]), flush=True)
    print('normalization_summary={}'.format(normalization_summary_tsv), flush=True)
    print('select_summary={}'.format(summary_tsv), flush=True)
    print('select_queue={}'.format(queue_tsv), flush=True)
    print('external_manifest={}'.format(manifest_tsv), flush=True)


def select_main(args):
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    runtime_args = clone_namespace(args, out_dir=out_dir)
    if getattr(runtime_args, 'species_tsv', None):
        select_batch_main(runtime_args)
        return
    metadata_dir = os.path.join(out_dir, 'metadata')
    dir_config = resolve_select_config_dir(runtime_args)
    check_config_dir(dir_path=dir_config, mode='select')
    path_metadata_table = os.path.join(metadata_dir, 'metadata.tsv')
    path_metadata_original = os.path.join(metadata_dir, 'metadata_original.tsv')

    metadata = load_metadata(runtime_args)
    metadata_original_df = metadata.df.copy(deep=True)
    metadata = filter_metadata_by_sample_group(metadata, runtime_args.sample_group)
    metadata = apply_select_filters(metadata, runtime_args, dir_config)
    write_select_outputs(
        path_metadata_original=path_metadata_original,
        path_metadata_table=path_metadata_table,
        metadata_dir=metadata_dir,
        metadata=metadata,
        metadata_original_df=metadata_original_df,
    )
