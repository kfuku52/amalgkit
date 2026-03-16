import datetime
import os
import re
import shutil

import pandas

from amalgkit.arg_utils import clone_namespace
from amalgkit.filter_utils import staged_output_dir
from amalgkit.metadata_utils import Metadata, load_metadata
from amalgkit.output_utils import atomic_write_dataframe

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

SELECT_RULES_FILENAME = 'select_rules.tsv'
SELECT_RULES_REQUIRED_COLUMNS = [
    'rule_id',
    'enabled',
    'stage',
    'priority',
    'columns',
    'pattern',
    'action',
    'target_column',
    'outcome',
    'scope_column',
    'scope_mode',
    'stop_on_match',
    'note',
]
SELECT_STAGE_ORDER = {
    'aggregate': 0,
    'normalize': 1,
    'exclude': 2,
    'control': 3,
}
SELECT_NORMALIZE_ORGANS = {'flower', 'leaf', 'root'}
SELECT_NORMALIZE_STATUSES = {'review', 'mixed', 'non_target'}
SELECT_DEFAULT_SCOPE_COLUMN = 'bioproject'
SELECT_DEFAULT_SCOPE_MODE = 'mark_other_rows_in_scope'
SELECT_DEFAULT_TARGET_COLUMN = 'exclusion'


def resolve_select_rules_tsv(args, out_dir=None):
    base_out_dir = os.path.realpath(args.out_dir if out_dir is None else out_dir)
    if getattr(args, 'select_rules_tsv', 'inferred') == 'inferred':
        return os.path.join(base_out_dir, SELECT_RULES_FILENAME)
    return os.path.realpath(args.select_rules_tsv)


def parse_select_rule_bool(value, field_name, rule_id):
    normalized = str(value).strip().lower()
    if normalized in {'yes', 'true', '1'}:
        return True
    if normalized in {'no', 'false', '0'}:
        return False
    raise ValueError(
        'Invalid {} value in select rule "{}": {} (expected yes/no).'.format(
            field_name,
            rule_id,
            value,
        )
    )


def parse_select_rule_priority(value, rule_id):
    normalized = str(value).strip()
    if normalized == '':
        return 0
    try:
        return int(normalized)
    except ValueError as exc:
        raise ValueError(
            'Invalid priority in select rule "{}": {}'.format(rule_id, value)
        ) from exc


def parse_select_rule_columns(value):
    return [token.strip() for token in str(value).split(',') if token.strip() != '']


def compile_select_rule_pattern(pattern, rule_id):
    try:
        return re.compile(str(pattern), re.IGNORECASE)
    except re.error as exc:
        raise ValueError(
            'Invalid regex pattern in select rule "{}": {}'.format(rule_id, pattern)
        ) from exc


def read_select_rules(select_rules_tsv):
    if not os.path.exists(select_rules_tsv):
        raise FileNotFoundError('select rules file not found: {}'.format(select_rules_tsv))
    if not os.path.isfile(select_rules_tsv):
        raise IsADirectoryError('select rules path exists but is not a file: {}'.format(select_rules_tsv))
    rules_df = pandas.read_csv(
        select_rules_tsv,
        sep='\t',
        dtype=str,
        keep_default_na=False,
        comment='#',
    ).fillna('')
    missing_columns = [col for col in SELECT_RULES_REQUIRED_COLUMNS if col not in rules_df.columns]
    if len(missing_columns) > 0:
        raise ValueError(
            'select_rules.tsv is missing required column(s): {}'.format(', '.join(missing_columns))
        )
    rules_df = rules_df.loc[:, SELECT_RULES_REQUIRED_COLUMNS].copy(deep=True)
    for col in SELECT_RULES_REQUIRED_COLUMNS:
        rules_df[col] = rules_df[col].astype(str).str.strip()
    if rules_df.shape[0] == 0:
        raise ValueError('select_rules.tsv does not contain any rules.')
    empty_rule_ids = rules_df['rule_id'] == ''
    if empty_rule_ids.any():
        empty_rows = (rules_df.index[empty_rule_ids] + 2).astype(str).tolist()
        raise ValueError(
            'select_rules.tsv contains empty rule_id values on row(s): {}'.format(', '.join(empty_rows))
        )
    if rules_df['rule_id'].duplicated().any():
        duplicated = sorted(rules_df.loc[rules_df['rule_id'].duplicated(), 'rule_id'].unique().tolist())
        raise ValueError('Duplicate rule_id detected in select_rules.tsv: {}'.format(', '.join(duplicated)))
    rules = []
    for _, row in rules_df.iterrows():
        rule_id = row['rule_id']
        enabled = parse_select_rule_bool(row['enabled'] if row['enabled'] != '' else 'yes', 'enabled', rule_id)
        if not enabled:
            continue
        stage = row['stage'].lower()
        if stage not in SELECT_STAGE_ORDER:
            raise ValueError(
                'Unsupported stage in select rule "{}": {}'.format(rule_id, row['stage'])
            )
        priority = parse_select_rule_priority(row['priority'], rule_id)
        columns = parse_select_rule_columns(row['columns'])
        action = row['action'].lower()
        target_column = row['target_column']
        outcome = row['outcome']
        scope_column = row['scope_column'] if row['scope_column'] != '' else SELECT_DEFAULT_SCOPE_COLUMN
        scope_mode = row['scope_mode'] if row['scope_mode'] != '' else SELECT_DEFAULT_SCOPE_MODE
        stop_on_match = parse_select_rule_bool(
            row['stop_on_match'] if row['stop_on_match'] != '' else 'yes',
            'stop_on_match',
            rule_id,
        )
        if stage == 'aggregate':
            if action != 'append':
                raise ValueError(
                    'Aggregate select rule "{}" must use action=append.'.format(rule_id)
                )
            if len(columns) == 0:
                raise ValueError(
                    'Aggregate select rule "{}" must define at least one source column.'.format(rule_id)
                )
            if target_column == '':
                raise ValueError(
                    'Aggregate select rule "{}" must define target_column.'.format(rule_id)
                )
            compiled_pattern = None
        else:
            if len(columns) == 0:
                raise ValueError(
                    'Select rule "{}" must define at least one column.'.format(rule_id)
                )
            if row['pattern'] == '':
                raise ValueError(
                    'Select rule "{}" must define pattern.'.format(rule_id)
                )
            compiled_pattern = compile_select_rule_pattern(row['pattern'], rule_id)
        if stage == 'normalize':
            if action != 'assign':
                raise ValueError(
                    'Normalize select rule "{}" must use action=assign.'.format(rule_id)
                )
            if outcome not in (SELECT_NORMALIZE_ORGANS | SELECT_NORMALIZE_STATUSES):
                raise ValueError(
                    'Normalize select rule "{}" has unsupported outcome: {}'.format(rule_id, outcome)
                )
        elif stage == 'exclude':
            if action != 'exclude':
                raise ValueError(
                    'Exclude select rule "{}" must use action=exclude.'.format(rule_id)
                )
            if outcome == '':
                raise ValueError(
                    'Exclude select rule "{}" must define outcome.'.format(rule_id)
                )
            if target_column == '':
                target_column = SELECT_DEFAULT_TARGET_COLUMN
        elif stage == 'control':
            if action != 'mark_non_control':
                raise ValueError(
                    'Control select rule "{}" must use action=mark_non_control.'.format(rule_id)
                )
            if target_column == '':
                target_column = SELECT_DEFAULT_TARGET_COLUMN
            if outcome == '':
                outcome = 'non_control'
            if scope_mode != SELECT_DEFAULT_SCOPE_MODE:
                raise ValueError(
                    'Control select rule "{}" has unsupported scope_mode: {}'.format(rule_id, scope_mode)
                )
        rules.append({
            'rule_id': rule_id,
            'stage': stage,
            'priority': priority,
            'columns': columns,
            'pattern': row['pattern'],
            'regex': compiled_pattern,
            'action': action,
            'target_column': target_column,
            'outcome': outcome,
            'scope_column': scope_column,
            'scope_mode': scope_mode,
            'stop_on_match': stop_on_match,
            'note': row['note'],
        })
    if len(rules) == 0:
        raise ValueError('select_rules.tsv does not contain any enabled rules.')
    return sorted(
        rules,
        key=lambda rule: (SELECT_STAGE_ORDER[rule['stage']], rule['priority'], rule['rule_id']),
    )


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


def ensure_select_target_column(df, target_column):
    if target_column not in df.columns:
        df[target_column] = 'no' if target_column == 'exclusion' else ''
    return df


def apply_select_aggregate_rules(df, select_rules):
    aggregate_rules = [rule for rule in select_rules if rule['stage'] == 'aggregate']
    out_df = df.copy(deep=True)
    for rule in aggregate_rules:
        target_column = rule['target_column']
        out_df = ensure_select_target_column(out_df, target_column)
        for source_column in rule['columns']:
            if source_column not in out_df.columns:
                continue
            source_series = out_df[source_column].fillna('').astype(str).str.strip()
            if not (source_series != '').any():
                continue
            target_series = out_df[target_column].fillna('').astype(str).str.strip()
            is_source_nonempty = source_series != ''
            is_target_empty = target_series == ''
            fill_mask = is_source_nonempty & is_target_empty
            append_mask = is_source_nonempty & (~is_target_empty)
            out_df.loc[fill_mask, target_column] = source_series.loc[fill_mask]
            out_df.loc[append_mask, target_column] = (
                target_series.loc[append_mask] + '; ' + source_series.loc[append_mask]
            )
    return out_df


def build_select_normalization_columns(normalize_rules):
    columns = []
    for rule in normalize_rules:
        for column in rule['columns']:
            if column not in columns:
                columns.append(column)
    if 'sample_group' not in columns:
        columns.insert(0, 'sample_group')
    return columns


def classify_select_text(text, normalize_rules):
    normalized = normalize_select_value(text)
    if normalized.lower() in SELECT_IGNORE_VALUES:
        return {'status': 'empty', 'organ': '', 'text': normalized, 'rule_id': ''}
    for rule in normalize_rules:
        if not rule['regex'].search(normalized):
            continue
        if rule['outcome'] in SELECT_NORMALIZE_ORGANS:
            return {
                'status': 'organ',
                'organ': rule['outcome'],
                'text': normalized,
                'rule_id': rule['rule_id'],
            }
        return {
            'status': rule['outcome'],
            'organ': '',
            'text': normalized,
            'rule_id': rule['rule_id'],
        }
    return {'status': 'unknown', 'organ': '', 'text': normalized, 'rule_id': ''}


def normalize_select_row(row, normalize_rules, normalization_columns):
    original_sample_group = normalize_select_value(row.get('sample_group', ''))
    result = {
        'sample_group_original': original_sample_group,
        'sample_group_normalization_status': 'unchanged',
        'sample_group_normalization_source': '',
        'sample_group_normalization_text': '',
        'sample_group_normalization_rule_id': '',
    }
    updated_sample_group = original_sample_group
    fallback_source = ''
    fallback_text = ''
    for rule in normalize_rules:
        for column in rule['columns']:
            if column not in row.index:
                continue
            classification = classify_select_text(row[column], [rule])
            status = classification['status']
            if status == 'empty':
                continue
            if (fallback_text == '') and (status == 'unknown'):
                fallback_source = column
                fallback_text = classification['text']
            if status == 'unknown':
                continue
            result['sample_group_normalization_status'] = status
            result['sample_group_normalization_source'] = column
            result['sample_group_normalization_text'] = classification['text']
            result['sample_group_normalization_rule_id'] = classification['rule_id']
            if status == 'organ':
                updated_sample_group = classification['organ']
            if rule['stop_on_match']:
                return updated_sample_group, result
    if fallback_text == '':
        for column in normalization_columns:
            if column not in row.index:
                continue
            fallback_value = normalize_select_value(row[column])
            if fallback_value.lower() in SELECT_IGNORE_VALUES:
                continue
            fallback_source = column
            fallback_text = fallback_value
            break
    if fallback_text != '':
        result['sample_group_normalization_status'] = 'unknown'
        result['sample_group_normalization_source'] = fallback_source
        result['sample_group_normalization_text'] = fallback_text
    return updated_sample_group, result


def normalize_select_metadata_frame(df, select_rules):
    normalize_rules = [rule for rule in select_rules if rule['stage'] == 'normalize']
    normalization_columns = build_select_normalization_columns(normalize_rules)
    out_df = df.copy(deep=True)
    if 'sample_group' not in out_df.columns:
        out_df['sample_group'] = ''
    normalized_groups = []
    statuses = []
    sources = []
    source_texts = []
    original_groups = []
    rule_ids = []
    for _, row in out_df.iterrows():
        updated_sample_group, result = normalize_select_row(
            row=row,
            normalize_rules=normalize_rules,
            normalization_columns=normalization_columns,
        )
        normalized_groups.append(updated_sample_group)
        statuses.append(result['sample_group_normalization_status'])
        sources.append(result['sample_group_normalization_source'])
        source_texts.append(result['sample_group_normalization_text'])
        original_groups.append(result['sample_group_original'])
        rule_ids.append(result['sample_group_normalization_rule_id'])
    out_df['sample_group_original'] = original_groups
    out_df['sample_group'] = normalized_groups
    out_df['sample_group_normalization_status'] = statuses
    out_df['sample_group_normalization_source'] = sources
    out_df['sample_group_normalization_text'] = source_texts
    out_df['sample_group_normalization_rule_id'] = rule_ids
    return out_df


def prepare_select_metadata(metadata, select_rules):
    metadata.remove_specialchars()
    metadata.df = apply_select_aggregate_rules(metadata.df, select_rules)
    metadata.df = normalize_select_metadata_frame(metadata.df, select_rules)
    return metadata


def build_select_rule_match_mask(df, rule):
    matched = pandas.Series(False, index=df.index)
    for column in rule['columns']:
        if column not in df.columns:
            continue
        text_col = df[column].fillna('').astype(str)
        matched = matched | text_col.str.contains(rule['pattern'], regex=True, case=False).fillna(False)
    return matched


def resolve_select_rule_write_mask(df, target_column, candidate_mask, stop_on_match):
    if not stop_on_match:
        return candidate_mask
    current_values = df[target_column].fillna('').astype(str).str.strip().str.lower()
    return candidate_mask & current_values.isin({'', 'no'})


def apply_select_exclude_rules(metadata, select_rules):
    exclude_rules = [rule for rule in select_rules if rule['stage'] == 'exclude']
    if len(exclude_rules) == 0:
        return metadata
    print('{}: Marking SRAs with select exclude rules'.format(datetime.datetime.now()), flush=True)
    metadata.df = ensure_select_target_column(metadata.df, SELECT_DEFAULT_TARGET_COLUMN)
    for rule in exclude_rules:
        target_column = rule['target_column'] if rule['target_column'] != '' else SELECT_DEFAULT_TARGET_COLUMN
        metadata.df = ensure_select_target_column(metadata.df, target_column)
        matched = build_select_rule_match_mask(metadata.df, rule)
        write_mask = resolve_select_rule_write_mask(metadata.df, target_column, matched, rule['stop_on_match'])
        metadata.df.loc[write_mask, target_column] = rule['outcome']
        print(
            '{}: Applying exclude rule "{}": matched {:,}, wrote {:,}'.format(
                datetime.datetime.now(),
                rule['rule_id'],
                int(matched.sum()),
                int(write_mask.sum()),
            ),
            flush=True,
        )
    return metadata


def apply_select_control_rules(metadata, select_rules):
    control_rules = [rule for rule in select_rules if rule['stage'] == 'control']
    if len(control_rules) == 0:
        return metadata
    print('{}: Marking SRAs with select control rules'.format(datetime.datetime.now()), flush=True)
    metadata.df = ensure_select_target_column(metadata.df, SELECT_DEFAULT_TARGET_COLUMN)
    for rule in control_rules:
        scope_column = rule['scope_column'] if rule['scope_column'] != '' else SELECT_DEFAULT_SCOPE_COLUMN
        if scope_column not in metadata.df.columns:
            print(
                '{}: Skipping control rule "{}" because scope column "{}" is missing'.format(
                    datetime.datetime.now(),
                    rule['rule_id'],
                    scope_column,
                ),
                flush=True,
            )
            continue
        target_column = rule['target_column'] if rule['target_column'] != '' else SELECT_DEFAULT_TARGET_COLUMN
        metadata.df = ensure_select_target_column(metadata.df, target_column)
        matched = build_select_rule_match_mask(metadata.df, rule)
        scope_values = metadata.df.loc[matched, scope_column].fillna('').astype(str).str.strip()
        num_control = 0
        num_treatment = 0
        for scope_value in scope_values.unique().tolist():
            if scope_value == '':
                continue
            is_scope = metadata.df.loc[:, scope_column].fillna('').astype(str).str.strip() == scope_value
            candidate_mask = is_scope & (~matched)
            write_mask = resolve_select_rule_write_mask(metadata.df, target_column, candidate_mask, rule['stop_on_match'])
            metadata.df.loc[write_mask, target_column] = rule['outcome']
            num_control += int((is_scope & matched).sum())
            num_treatment += int(write_mask.sum())
        print(
            '{}: Applying control rule "{}" within "{}": control {:,}, marked {:,}'.format(
                datetime.datetime.now(),
                rule['rule_id'],
                scope_column,
                num_control,
                num_treatment,
            ),
            flush=True,
        )
    return metadata


def apply_select_filters(metadata, args, select_rules):
    metadata.nspot_cutoff(args.min_nspots)
    metadata.mark_missing_rank(args.mark_missing_rank)
    metadata.mark_redundant_biosample(args.mark_redundant_biosamples)
    metadata = apply_select_exclude_rules(metadata, select_rules)
    metadata = apply_select_control_rules(metadata, select_rules)
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
    select_rules_tsv = resolve_select_rules_tsv(args, out_dir=out_dir)
    select_rules = read_select_rules(select_rules_tsv)
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
        metadata = Metadata.from_DataFrame(merged_df)
        metadata = prepare_select_metadata(metadata, select_rules)
        normalized_df = metadata.df.copy(deep=True)
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
        metadata_original_df = metadata.df.copy(deep=True)
        metadata = filter_metadata_by_sample_group(metadata, runtime_args.sample_group)
        metadata = apply_select_filters(metadata, runtime_args, select_rules)
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
    select_rules_tsv = resolve_select_rules_tsv(runtime_args, out_dir=out_dir)
    select_rules = read_select_rules(select_rules_tsv)
    path_metadata_table = os.path.join(metadata_dir, 'metadata.tsv')
    path_metadata_original = os.path.join(metadata_dir, 'metadata_original.tsv')

    metadata = load_metadata(runtime_args)
    metadata = prepare_select_metadata(metadata, select_rules)
    metadata_original_df = metadata.df.copy(deep=True)
    metadata = filter_metadata_by_sample_group(metadata, runtime_args.sample_group)
    metadata = apply_select_filters(metadata, runtime_args, select_rules)
    write_select_outputs(
        path_metadata_original=path_metadata_original,
        path_metadata_table=path_metadata_table,
        metadata_dir=metadata_dir,
        metadata=metadata,
        metadata_original_df=metadata_original_df,
    )
