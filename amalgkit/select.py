import datetime
import os
import re
import shlex
import shutil

import pandas

from amalgkit.arg_utils import clone_namespace
from amalgkit.filter_utils import staged_output_dir
from amalgkit.metadata_utils import Metadata, load_metadata, SELECT_SAMPLING_STRATEGIES
from amalgkit.output_utils import atomic_write_dataframe
from amalgkit.parallel_utils import resolve_worker_allocation, run_tasks_with_optional_threads

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
    'parameter_name',
    'parameter_value',
    'note',
]
SELECT_STAGE_ORDER = {
    'aggregate': 0,
    'validate': 1,
    'normalize': 2,
    'exclude': 3,
    'control': 4,
    'filter': 5,
    'dedup': 6,
}
SELECT_NORMALIZE_ORGANS = {'flower', 'leaf', 'root'}
SELECT_NORMALIZE_STATUSES = {'review', 'mixed', 'non_target'}
SELECT_VALIDATION_STRONG_SPLIT_PATTERN = re.compile(r'\s*(?:/|&|\b(?:and|or|plus)\b)\s*', re.IGNORECASE)
SELECT_VALIDATION_LIST_SPLIT_PATTERN = re.compile(r'\s*[;,]\s*')
SELECT_CROSS_FIELD_TISSUE_CONFLICT_RULE_ID = 'normalize_review_cross_field_tissue_conflict'
SELECT_CROSS_FIELD_STRONG_TISSUE_COLUMNS = [
    'sample_attribute_tissue',
    'tissue',
    'organism_part',
    'organ',
    'plant_structure',
    'plant_anatomical_entity',
    'tissue_type',
    'tissue_types',
    'tissues',
    'tissue_source',
]
SELECT_CROSS_FIELD_ORGAN_PATTERNS = {
    'flower': re.compile(
        r'\b(?:flower|flowers|inflorescence|inflorescences|catkin|catkins|'
        r'panicle|panicles|spikelet|spikelets|floret|florets)\b',
        re.IGNORECASE,
    ),
    'leaf': re.compile(
        r'\b(?:flag leaf|flag leaves|leaf|leaves|foliar|foliage|trifoliate|'
        r'lamina|laminae|leaflet|leaflets)\b',
        re.IGNORECASE,
    ),
    'root': re.compile(
        r'\b(?:root|roots|taproot|taproots|primary root|primary roots|'
        r'lateral root|lateral roots|radicle|radicles|root[\s_-]?tip|'
        r'root[\s_-]?tips|root hair|root hairs)\b',
        re.IGNORECASE,
    ),
}
SELECT_DEFAULT_SCOPE_COLUMN = 'bioproject'
SELECT_DEFAULT_SCOPE_MODE = 'mark_other_rows_in_scope'
SELECT_DEFAULT_TARGET_COLUMN = 'exclusion'
SELECT_PARAMETER_DEFINITIONS = {
    'min_nspots': {
        'kind': 'int',
        'minimum': 0,
        'required': True,
    },
    'max_sample': {
        'kind': 'int',
        'minimum': 1,
        'required': True,
    },
    'mark_missing_rank': {
        'kind': 'choice',
        'choices': ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'domain', 'none'],
        'required': True,
    },
    'mark_redundant_biosamples': {
        'kind': 'bool',
        'required': True,
    },
    'sample_group': {
        'kind': 'string',
        'required': False,
    },
    'sampling_strategy': {
        'kind': 'choice',
        'choices': list(SELECT_SAMPLING_STRATEGIES),
        'required': False,
    },
}


def resolve_select_rules_tsv(args, out_dir=None):
    base_out_dir = os.path.realpath(args.out_dir if out_dir is None else out_dir)
    if getattr(args, 'select_rules_tsv', 'inferred') == 'inferred':
        return os.path.join(base_out_dir, SELECT_RULES_FILENAME)
    return os.path.realpath(args.select_rules_tsv)


def build_missing_select_rules_message(select_rules_tsv):
    out_dir = os.path.dirname(os.path.realpath(select_rules_tsv))
    quoted_out_dir = shlex.quote(out_dir)
    return '\n'.join([
        'select rules file not found: {}'.format(select_rules_tsv),
        'Create it with one of:',
        '  amalgkit dataset --name init --out_dir {}'.format(quoted_out_dir),
        '  amalgkit dataset --rule_set base --out_dir {}'.format(quoted_out_dir),
    ])


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


def build_select_scope_key_series(df, scope_columns):
    normalized_columns = []
    for column in scope_columns:
        if column not in df.columns:
            raise KeyError(column)
        normalized_columns.append(df[column].fillna('').astype(str).str.strip())
    return pandas.Series(list(zip(*normalized_columns)), index=df.index)


def is_valid_select_scope_key(scope_key):
    return all(str(value).strip() != '' for value in scope_key)


def compile_select_rule_pattern(pattern, rule_id):
    try:
        return re.compile(str(pattern), re.IGNORECASE)
    except re.error as exc:
        raise ValueError(
            'Invalid regex pattern in select rule "{}": {}'.format(rule_id, pattern)
        ) from exc


def parse_select_parameter_value(parameter_name, parameter_value, rule_id):
    definition = SELECT_PARAMETER_DEFINITIONS[parameter_name]
    raw_value = str(parameter_value).strip()
    if raw_value == '':
        raise ValueError(
            'Parameter select rule "{}" must define parameter_value.'.format(rule_id)
        )
    if definition['kind'] == 'int':
        try:
            value = int(raw_value)
        except ValueError as exc:
            raise ValueError(
                'Invalid integer parameter_value in select rule "{}": {}'.format(rule_id, parameter_value)
            ) from exc
        minimum = definition.get('minimum')
        if (minimum is not None) and (value < minimum):
            raise ValueError(
                'Parameter select rule "{}" must be >= {}: {}'.format(rule_id, minimum, parameter_value)
            )
        return value
    if definition['kind'] == 'choice':
        if raw_value not in definition['choices']:
            raise ValueError(
                'Invalid parameter_value in select rule "{}": {} (expected one of {}).'.format(
                    rule_id,
                    parameter_value,
                    ', '.join(definition['choices']),
                )
            )
        return raw_value
    if definition['kind'] == 'bool':
        return parse_select_rule_bool(raw_value, 'parameter_value', rule_id)
    if definition['kind'] == 'string':
        return raw_value
    raise ValueError(
        'Unsupported parameter kind for "{}": {}'.format(parameter_name, definition['kind'])
    )


def resolve_select_runtime_parameter(args, parameter_name, rule_id):
    if parameter_name == '':
        raise ValueError(
            'Select rule "{}" must define parameter_name.'.format(rule_id)
        )
    if parameter_name not in SELECT_PARAMETER_DEFINITIONS:
        raise ValueError(
            'Unsupported parameter_name in select rule "{}": {}'.format(rule_id, parameter_name)
        )
    if not hasattr(args, parameter_name):
        raise ValueError(
            'Runtime args do not define parameter "{}" required by select rule "{}".'.format(
                parameter_name,
                rule_id,
            )
        )
    return getattr(args, parameter_name)


def load_select_rules_table(select_rules_tsv):
    if not os.path.exists(select_rules_tsv):
        raise FileNotFoundError(build_missing_select_rules_message(select_rules_tsv))
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
    return rules_df


def read_select_config(select_rules_tsv):
    rules_df = load_select_rules_table(select_rules_tsv)
    rules = []
    parameters = {}
    for _, row in rules_df.iterrows():
        rule_id = row['rule_id']
        enabled = parse_select_rule_bool(row['enabled'] if row['enabled'] != '' else 'yes', 'enabled', rule_id)
        if not enabled:
            continue
        stage = row['stage'].lower()
        if stage == 'parameter':
            parameter_name = row['parameter_name']
            if parameter_name == '':
                raise ValueError(
                    'Parameter select rule "{}" must define parameter_name.'.format(rule_id)
                )
            if parameter_name not in SELECT_PARAMETER_DEFINITIONS:
                raise ValueError(
                    'Unsupported parameter_name in select rule "{}": {}'.format(rule_id, parameter_name)
                )
            if parameter_name in parameters:
                raise ValueError(
                    'Duplicate parameter_name detected in select_rules.tsv: {}'.format(parameter_name)
                )
            parameters[parameter_name] = parse_select_parameter_value(
                parameter_name,
                row['parameter_value'],
                rule_id,
            )
            continue
        if stage not in SELECT_STAGE_ORDER:
            raise ValueError(
                'Unsupported stage in select rule "{}": {}'.format(rule_id, row['stage'])
            )
        priority = parse_select_rule_priority(row['priority'], rule_id)
        columns = parse_select_rule_columns(row['columns'])
        action = row['action'].lower()
        target_column = row['target_column']
        outcome = row['outcome']
        parameter_name = row['parameter_name']
        scope_column = row['scope_column'] if row['scope_column'] != '' else SELECT_DEFAULT_SCOPE_COLUMN
        scope_columns = parse_select_rule_columns(scope_column)
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
        elif stage in {'normalize', 'exclude', 'control', 'validate'}:
            if len(columns) == 0:
                if stage != 'validate':
                    raise ValueError(
                        'Select rule "{}" must define at least one column.'.format(rule_id)
                    )
            if row['pattern'] == '':
                raise ValueError(
                    'Select rule "{}" must define pattern.'.format(rule_id)
                )
            compiled_pattern = compile_select_rule_pattern(row['pattern'], rule_id)
        else:
            compiled_pattern = None
        if stage == 'normalize':
            if action not in {'assign', 'assign_safe'}:
                raise ValueError(
                    'Normalize select rule "{}" must use action=assign or action=assign_safe.'.format(rule_id)
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
            if len(scope_columns) == 0:
                raise ValueError(
                    'Control select rule "{}" must define at least one scope column.'.format(rule_id)
                )
            if target_column == '':
                target_column = SELECT_DEFAULT_TARGET_COLUMN
            if outcome == '':
                outcome = 'non_control'
            if scope_mode != SELECT_DEFAULT_SCOPE_MODE:
                raise ValueError(
                    'Control select rule "{}" has unsupported scope_mode: {}'.format(rule_id, scope_mode)
                )
        elif stage == 'filter':
            if action == 'exclude_if_lt_parameter':
                if len(columns) != 1:
                    raise ValueError(
                        'Filter select rule "{}" must define exactly one column for action=exclude_if_lt_parameter.'.format(rule_id)
                    )
                if target_column == '':
                    target_column = SELECT_DEFAULT_TARGET_COLUMN
                if outcome == '':
                    raise ValueError(
                        'Filter select rule "{}" must define outcome.'.format(rule_id)
                    )
                if parameter_name == '':
                    raise ValueError(
                        'Filter select rule "{}" must define parameter_name.'.format(rule_id)
                    )
            elif action == 'exclude_if_empty':
                if len(columns) == 0:
                    raise ValueError(
                        'Filter select rule "{}" must define at least one column for action=exclude_if_empty.'.format(rule_id)
                    )
                if target_column == '':
                    target_column = SELECT_DEFAULT_TARGET_COLUMN
                if outcome == '':
                    raise ValueError(
                        'Filter select rule "{}" must define outcome.'.format(rule_id)
                    )
            elif action == 'exclude_if_missing_selected_rank':
                if len(columns) != 1:
                    raise ValueError(
                        'Filter select rule "{}" must define exactly one column template for action=exclude_if_missing_selected_rank.'.format(rule_id)
                    )
                if target_column == '':
                    target_column = SELECT_DEFAULT_TARGET_COLUMN
                if outcome == '':
                    raise ValueError(
                        'Filter select rule "{}" must define outcome.'.format(rule_id)
                    )
                if parameter_name == '':
                    raise ValueError(
                        'Filter select rule "{}" must define parameter_name.'.format(rule_id)
                    )
            else:
                raise ValueError(
                    'Unsupported filter action in select rule "{}": {}'.format(rule_id, action)
                )
        elif stage == 'dedup':
            if action != 'exclude_redundant_by_best_numeric':
                raise ValueError(
                    'Dedup select rule "{}" must use action=exclude_redundant_by_best_numeric.'.format(rule_id)
                )
            if len(columns) < 2:
                raise ValueError(
                    'Dedup select rule "{}" must define at least two key columns.'.format(rule_id)
                )
            if target_column == '':
                raise ValueError(
                    'Dedup select rule "{}" must define target_column.'.format(rule_id)
                )
            if outcome == '':
                raise ValueError(
                    'Dedup select rule "{}" must define outcome.'.format(rule_id)
                )
        elif stage == 'validate':
            if action == 'hint_organ':
                if outcome not in SELECT_NORMALIZE_ORGANS:
                    raise ValueError(
                        'Validate select rule "{}" has unsupported organ outcome: {}'.format(rule_id, outcome)
                    )
            elif action == 'hint_review':
                if outcome == '':
                    outcome = 'review'
                if outcome != 'review':
                    raise ValueError(
                        'Validate select rule "{}" with action=hint_review must use outcome=review.'.format(rule_id)
                    )
            elif action == 'ignore_segment':
                if outcome == '':
                    outcome = 'ignore'
                if outcome != 'ignore':
                    raise ValueError(
                        'Validate select rule "{}" with action=ignore_segment must use outcome=ignore.'.format(rule_id)
                    )
            else:
                raise ValueError(
                    'Unsupported validate action in select rule "{}": {}'.format(rule_id, action)
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
            'parameter_name': parameter_name,
            'scope_column': scope_column,
            'scope_columns': scope_columns,
            'scope_mode': scope_mode,
            'stop_on_match': stop_on_match,
            'skip_validation': (stage == 'normalize') and (action == 'assign_safe'),
            'note': row['note'],
        })
    if len(rules) == 0:
        raise ValueError('select_rules.tsv does not contain any enabled rules.')
    return {
        'rules': sorted(
            rules,
            key=lambda rule: (SELECT_STAGE_ORDER[rule['stage']], rule['priority'], rule['rule_id']),
        ),
        'parameters': parameters,
    }


def read_select_rules(select_rules_tsv):
    return read_select_config(select_rules_tsv)['rules']


def read_select_parameters(select_rules_tsv):
    return read_select_config(select_rules_tsv)['parameters']


def apply_select_config_parameters(runtime_args, select_parameters):
    for parameter_name, definition in SELECT_PARAMETER_DEFINITIONS.items():
        current_value = getattr(runtime_args, parameter_name, None)
        if current_value is not None:
            continue
        if parameter_name not in select_parameters:
            if definition.get('required', False):
                raise ValueError(
                    'Select parameter "{}" must be provided via select_rules.tsv or runtime args.'.format(
                        parameter_name
                    )
                )
            continue
        setattr(runtime_args, parameter_name, select_parameters[parameter_name])
    return runtime_args


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


def classify_select_text_raw(text, normalize_rules):
    normalized = normalize_select_value(text)
    if normalized.lower() in SELECT_IGNORE_VALUES:
        return {'status': 'empty', 'organ': '', 'text': normalized, 'rule_id': '', 'skip_validation': False}
    for rule in normalize_rules:
        if not rule['regex'].search(normalized):
            continue
        if rule['outcome'] in SELECT_NORMALIZE_ORGANS:
            return {
                'status': 'organ',
                'organ': rule['outcome'],
                'text': normalized,
                'rule_id': rule['rule_id'],
                'skip_validation': rule.get('skip_validation', False),
            }
        return {
            'status': rule['outcome'],
            'organ': '',
            'text': normalized,
            'rule_id': rule['rule_id'],
            'skip_validation': False,
        }
    return {'status': 'unknown', 'organ': '', 'text': normalized, 'rule_id': '', 'skip_validation': False}


def split_select_validation_segments(normalized, split_pattern):
    return [
        segment.strip(" \t\r\n;:,.()[]{}")
        for segment in split_pattern.split(normalized)
        if segment.strip(" \t\r\n;:,.()[]{}") != ''
    ]


def classify_select_validation_segment(segment, normalize_rules, validate_rules=None):
    classification = classify_select_text_raw(segment, normalize_rules)
    if classification['status'] != 'unknown':
        return classification
    normalized = normalize_select_value(segment)
    for rule in (validate_rules or []):
        if not rule['regex'].search(normalized):
            continue
        if rule['action'] == 'hint_organ':
            return {
                'status': 'organ',
                'organ': rule['outcome'],
                'text': normalized,
                'rule_id': rule['rule_id'],
            }
        if rule['action'] == 'hint_review':
            return {
                'status': 'review',
                'organ': '',
                'text': normalized,
                'rule_id': rule['rule_id'],
            }
        if rule['action'] == 'ignore_segment':
            return {
                'status': 'ignore',
                'organ': '',
                'text': normalized,
                'rule_id': rule['rule_id'],
            }
    return classification


def summarize_select_segment_classifications(segments, organ, normalize_rules, validate_rules, unknown_policy):
    seen_target = False
    seen_other_target = False
    seen_review_like = False
    for segment in segments:
        classification = classify_select_validation_segment(segment, normalize_rules, validate_rules=validate_rules)
        status = classification['status']
        if status == 'empty':
            continue
        if status == 'ignore':
            continue
        if status == 'unknown':
            if unknown_policy == 'review':
                return {
                    'status': 'review',
                    'organ': '',
                    'text': '',
                    'rule_id': 'validate_review_unknown_segment',
                }
            continue
        if status == 'organ':
            if classification['organ'] == organ:
                seen_target = True
            else:
                seen_other_target = True
            continue
        if status == 'mixed':
            seen_other_target = True
            continue
        seen_review_like = True
    if seen_other_target:
        return {
            'status': 'mixed',
            'organ': '',
            'text': '',
            'rule_id': 'validate_mixed_segment',
        }
    if seen_review_like and seen_target:
        return {
            'status': 'review',
            'organ': '',
            'text': '',
            'rule_id': 'validate_review_segment',
        }
    return None


def validate_select_organ_text(text, organ, normalize_rules, validate_rules=None):
    normalized = normalize_select_value(text)
    if normalized.lower() in SELECT_IGNORE_VALUES:
        return None
    if SELECT_VALIDATION_STRONG_SPLIT_PATTERN.search(normalized):
        segments = split_select_validation_segments(normalized, SELECT_VALIDATION_STRONG_SPLIT_PATTERN)
        if len(segments) > 1:
            summary = summarize_select_segment_classifications(
                segments=segments,
                organ=organ,
                normalize_rules=normalize_rules,
                validate_rules=validate_rules,
                unknown_policy='review',
            )
            if summary is not None:
                summary['text'] = normalized
                return summary
    if SELECT_VALIDATION_LIST_SPLIT_PATTERN.search(normalized):
        segments = split_select_validation_segments(normalized, SELECT_VALIDATION_LIST_SPLIT_PATTERN)
        if len(segments) > 1:
            summary = summarize_select_segment_classifications(
                segments=segments,
                organ=organ,
                normalize_rules=normalize_rules,
                validate_rules=validate_rules,
                unknown_policy='ignore',
            )
            if summary is not None:
                summary['text'] = normalized
                return summary
    return None


def classify_select_text(text, normalize_rules, validate_rules=None):
    classification = classify_select_text_raw(text, normalize_rules)
    if classification['status'] != 'organ':
        return classification
    if classification.get('skip_validation', False):
        return classification
    validation = validate_select_organ_text(
        text=classification['text'],
        organ=classification['organ'],
        normalize_rules=normalize_rules,
        validate_rules=validate_rules,
    )
    if validation is not None:
        return validation
    return classification


def normalize_select_row(row, normalize_rules, normalization_columns, validate_rules=None):
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
            classification = classify_select_text_raw(row[column], [rule])
            status = classification['status']
            if status == 'empty':
                continue
            if (fallback_text == '') and (status == 'unknown'):
                fallback_source = column
                fallback_text = classification['text']
            if status == 'unknown':
                continue
            if status == 'organ':
                if not rule.get('skip_validation', False):
                    validation = validate_select_organ_text(
                        text=classification['text'],
                        organ=classification['organ'],
                        normalize_rules=normalize_rules,
                        validate_rules=validate_rules,
                    )
                    if validation is not None:
                        classification = validation
                        status = classification['status']
            result['sample_group_normalization_status'] = status
            result['sample_group_normalization_source'] = column
            result['sample_group_normalization_text'] = classification['text']
            result['sample_group_normalization_rule_id'] = classification['rule_id']
            if status == 'organ':
                updated_sample_group = classification['organ']
            else:
                updated_sample_group = status
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
    validate_rules = [rule for rule in select_rules if rule['stage'] == 'validate']
    normalize_rules = [rule for rule in select_rules if rule['stage'] == 'normalize']
    normalization_columns = build_select_normalization_columns(normalize_rules)
    out_df = df.copy(deep=True)
    if 'sample_group' not in out_df.columns:
        out_df['sample_group'] = ''
    normalized_columns = {}
    ignored_columns = {}
    for column in normalization_columns:
        if column not in out_df.columns:
            continue
        normalized_series = normalize_select_series(out_df[column])
        normalized_columns[column] = normalized_series
        ignored_columns[column] = normalized_series.str.lower().isin(SELECT_IGNORE_VALUES)
    original_groups = normalized_columns.get(
        'sample_group',
        pandas.Series('', index=out_df.index, dtype=str),
    ).copy()
    updated_groups = original_groups.copy()
    statuses = pandas.Series('unchanged', index=out_df.index, dtype=str)
    sources = pandas.Series('', index=out_df.index, dtype=str)
    source_texts = pandas.Series('', index=out_df.index, dtype=str)
    rule_ids = pandas.Series('', index=out_df.index, dtype=str)
    fallback_sources = pandas.Series('', index=out_df.index, dtype=str)
    fallback_texts = pandas.Series('', index=out_df.index, dtype=str)
    resolved = pandas.Series(False, index=out_df.index)
    for rule in normalize_rules:
        compiled_pattern = rule['regex']
        for column in rule['columns']:
            if column not in normalized_columns:
                continue
            text_series = normalized_columns[column]
            ignored_mask = ignored_columns[column]
            nonempty_mask = ~ignored_mask
            fallback_mask = (fallback_texts == '') & nonempty_mask
            if fallback_mask.any():
                fallback_sources.loc[fallback_mask] = column
                fallback_texts.loc[fallback_mask] = text_series.loc[fallback_mask]
            matched_mask = text_series.str.contains(compiled_pattern, regex=True).fillna(False)
            candidate_mask = nonempty_mask & matched_mask
            if rule['stop_on_match']:
                candidate_mask = candidate_mask & (~resolved)
            if not candidate_mask.any():
                continue
            if rule['outcome'] in SELECT_NORMALIZE_ORGANS:
                if not rule.get('skip_validation', False):
                    validation_series = text_series.loc[candidate_mask].apply(
                        lambda text: validate_select_organ_text(
                            text=text,
                            organ=rule['outcome'],
                            normalize_rules=normalize_rules,
                            validate_rules=validate_rules,
                        )
                    )
                    validation_mask = validation_series.notna()
                    if validation_mask.any():
                        validation_records = pandas.DataFrame(
                            list(validation_series.loc[validation_mask]),
                            index=validation_series.loc[validation_mask].index,
                        )
                        if not validation_records.empty:
                            mixed_mask = validation_records['status'] == 'mixed'
                            review_mask = validation_records['status'] == 'review'
                            if mixed_mask.any():
                                mixed_index = validation_records.index[mixed_mask]
                                updated_groups.loc[mixed_index] = 'mixed'
                                statuses.loc[mixed_index] = 'mixed'
                                sources.loc[mixed_index] = column
                                source_texts.loc[mixed_index] = validation_records.loc[mixed_index, 'text']
                                rule_ids.loc[mixed_index] = validation_records.loc[mixed_index, 'rule_id']
                            if review_mask.any():
                                review_index = validation_records.index[review_mask]
                                updated_groups.loc[review_index] = 'review'
                                statuses.loc[review_index] = 'review'
                                sources.loc[review_index] = column
                                source_texts.loc[review_index] = validation_records.loc[review_index, 'text']
                                rule_ids.loc[review_index] = validation_records.loc[review_index, 'rule_id']
                            candidate_mask = candidate_mask.copy()
                            candidate_mask.loc[validation_records.index] = False
                            if rule['stop_on_match']:
                                resolved.loc[validation_records.index] = True
                if not candidate_mask.any():
                    continue
                updated_groups.loc[candidate_mask] = rule['outcome']
                statuses.loc[candidate_mask] = 'organ'
            else:
                updated_groups.loc[candidate_mask] = rule['outcome']
                statuses.loc[candidate_mask] = rule['outcome']
            sources.loc[candidate_mask] = column
            source_texts.loc[candidate_mask] = text_series.loc[candidate_mask]
            rule_ids.loc[candidate_mask] = rule['rule_id']
            if rule['stop_on_match']:
                resolved = resolved | candidate_mask
    for column in normalization_columns:
        if column not in normalized_columns:
            continue
        text_series = normalized_columns[column]
        fallback_mask = (fallback_texts == '') & (~ignored_columns[column])
        if fallback_mask.any():
            fallback_sources.loc[fallback_mask] = column
            fallback_texts.loc[fallback_mask] = text_series.loc[fallback_mask]
    unknown_mask = (statuses == 'unchanged') & (fallback_texts != '')
    if unknown_mask.any():
        statuses.loc[unknown_mask] = 'unknown'
        sources.loc[unknown_mask] = fallback_sources.loc[unknown_mask]
        source_texts.loc[unknown_mask] = fallback_texts.loc[unknown_mask]
    out_df['sample_group_original'] = original_groups
    out_df['sample_group'] = updated_groups
    out_df['sample_group_normalization_status'] = statuses
    out_df['sample_group_normalization_source'] = sources
    out_df['sample_group_normalization_text'] = source_texts
    out_df['sample_group_normalization_rule_id'] = rule_ids
    out_df = review_select_cross_field_tissue_conflicts(out_df, normalize_rules)
    return out_df


def review_select_cross_field_tissue_conflicts(df, normalize_rules):
    safe_rule_ids = {
        rule['rule_id']
        for rule in normalize_rules
        if rule.get('skip_validation', False)
    }
    strong_columns = [
        column
        for column in SELECT_CROSS_FIELD_STRONG_TISSUE_COLUMNS
        if column in df.columns
    ]
    if len(strong_columns) == 0:
        return df
    out_df = df.copy(deep=True)
    assigned_groups = out_df['sample_group'].fillna('').astype(str).str.strip()
    eligible_mask = assigned_groups.isin(SELECT_NORMALIZE_ORGANS)
    if not eligible_mask.any():
        return out_df
    normalization_sources = out_df['sample_group_normalization_source'].fillna('').astype(str)
    normalization_rule_ids = out_df['sample_group_normalization_rule_id'].fillna('').astype(str)
    conflict_mask = pandas.Series(False, index=out_df.index)
    conflict_sources = pandas.Series('', index=out_df.index, dtype=str)
    conflict_texts = pandas.Series('', index=out_df.index, dtype=str)
    for column in strong_columns:
        text_series = normalize_select_series(out_df[column])
        nonempty_mask = ~text_series.str.lower().isin(SELECT_IGNORE_VALUES)
        if not nonempty_mask.any():
            continue
        safe_source_mask = (normalization_sources == column) & normalization_rule_ids.isin(safe_rule_ids)
        for organ, pattern in SELECT_CROSS_FIELD_ORGAN_PATTERNS.items():
            organ_mask = text_series.str.contains(pattern, regex=True).fillna(False)
            current_mask = eligible_mask & nonempty_mask & organ_mask & (assigned_groups != organ)
            current_mask = current_mask & (~safe_source_mask)
            new_mask = current_mask & (~conflict_mask)
            if not new_mask.any():
                continue
            conflict_mask.loc[new_mask] = True
            conflict_sources.loc[new_mask] = column
            previous_rule_ids = normalization_rule_ids.loc[new_mask]
            conflict_texts.loc[new_mask] = (
                'assigned=' + assigned_groups.loc[new_mask]
                + '; conflict=' + column + '=' + text_series.loc[new_mask]
                + '; previous_rule_id=' + previous_rule_ids
            )
    if not conflict_mask.any():
        return out_df
    out_df.loc[conflict_mask, 'sample_group'] = 'review'
    out_df.loc[conflict_mask, 'sample_group_normalization_status'] = 'review'
    out_df.loc[conflict_mask, 'sample_group_normalization_source'] = conflict_sources.loc[conflict_mask]
    out_df.loc[conflict_mask, 'sample_group_normalization_text'] = conflict_texts.loc[conflict_mask]
    out_df.loc[conflict_mask, 'sample_group_normalization_rule_id'] = SELECT_CROSS_FIELD_TISSUE_CONFLICT_RULE_ID
    return out_df


def prepare_select_metadata(metadata, select_rules):
    metadata.remove_specialchars()
    metadata.df = apply_select_aggregate_rules(metadata.df, select_rules)
    metadata.df = normalize_select_metadata_frame(metadata.df, select_rules)
    return metadata


def build_select_rule_match_mask(df, rule):
    return build_select_rule_match_mask_with_cache(df=df, rule=rule, text_cache=None)


def get_select_text_series(df, column, text_cache):
    if text_cache is None:
        return df[column].fillna('').astype(str)
    if column not in text_cache:
        text_cache[column] = df[column].fillna('').astype(str)
    return text_cache[column]


def build_select_rule_match_mask_with_cache(df, rule, text_cache=None):
    matched = pandas.Series(False, index=df.index)
    compiled_pattern = rule.get('regex', rule['pattern'])
    for column in rule['columns']:
        if column not in df.columns:
            continue
        text_col = get_select_text_series(df, column, text_cache=text_cache)
        matched = matched | text_col.str.contains(compiled_pattern, regex=True).fillna(False)
    return matched


def get_select_scope_key_series(df, scope_columns, scope_key_cache):
    cache_key = tuple(scope_columns)
    if scope_key_cache is None:
        return build_select_scope_key_series(df, list(cache_key))
    if cache_key not in scope_key_cache:
        scope_key_cache[cache_key] = build_select_scope_key_series(df, list(cache_key))
    return scope_key_cache[cache_key]


def build_select_control_rule_groups(control_rules):
    grouped = {}
    for rule in control_rules:
        scope_columns = tuple(rule.get('scope_columns', parse_select_rule_columns(rule.get('scope_column', SELECT_DEFAULT_SCOPE_COLUMN))))
        target_column = rule['target_column'] if rule['target_column'] != '' else SELECT_DEFAULT_TARGET_COLUMN
        outcome = rule['outcome'] if rule['outcome'] != '' else 'non_control'
        aggregate_key = (scope_columns, target_column, outcome, rule['stop_on_match'])
        group = grouped.setdefault(
            aggregate_key,
            {
                'scope_columns': scope_columns,
                'target_column': target_column,
                'outcome': outcome,
                'stop_on_match': rule['stop_on_match'],
                'column_patterns': {},
            },
        )
        for column in rule['columns']:
            group['column_patterns'].setdefault(column, []).append(rule['pattern'])
    for aggregate_key, group in grouped.items():
        compiled_by_column = {}
        for column, patterns in group['column_patterns'].items():
            combined_pattern = '|'.join(['(?:{})'.format(pattern) for pattern in patterns])
            compiled_by_column[column] = compile_select_rule_pattern(
                combined_pattern,
                'control_group:{}:{}'.format(','.join(aggregate_key[0]), column),
            )
        group['compiled_by_column'] = compiled_by_column
    return list(grouped.values())


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
    text_cache = {}
    total_matched = 0
    total_wrote = 0
    rules_with_writes = 0
    for rule in exclude_rules:
        target_column = rule['target_column'] if rule['target_column'] != '' else SELECT_DEFAULT_TARGET_COLUMN
        metadata.df = ensure_select_target_column(metadata.df, target_column)
        matched = build_select_rule_match_mask_with_cache(metadata.df, rule, text_cache=text_cache)
        write_mask = resolve_select_rule_write_mask(metadata.df, target_column, matched, rule['stop_on_match'])
        matched_count = int(matched.sum())
        wrote_count = int(write_mask.sum())
        metadata.df.loc[write_mask, target_column] = rule['outcome']
        total_matched += matched_count
        total_wrote += wrote_count
        if wrote_count > 0:
            rules_with_writes += 1
    print(
        '{}: Exclude rules complete: rules={}, matched_total={:,}, wrote_total={:,}, wrote_rules={}'.format(
            datetime.datetime.now(),
            len(exclude_rules),
            total_matched,
            total_wrote,
            rules_with_writes,
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
    text_cache = {}
    scope_key_cache = {}
    control_groups = build_select_control_rule_groups(control_rules)
    total_control = 0
    total_matched = 0
    aggregated_matches = dict()
    for group in control_groups:
        scope_columns = list(group['scope_columns'])
        missing_scope_columns = [column for column in scope_columns if column not in metadata.df.columns]
        if len(missing_scope_columns) > 0:
            print(
                '{}: Skipping control rule group because scope column(s) "{}" are missing'.format(
                    datetime.datetime.now(),
                    ','.join(missing_scope_columns),
                ),
                flush=True,
            )
            continue
        target_column = group['target_column']
        metadata.df = ensure_select_target_column(metadata.df, target_column)
        matched = pandas.Series(False, index=metadata.df.index)
        for column, compiled_pattern in group['compiled_by_column'].items():
            if column not in metadata.df.columns:
                continue
            text_col = get_select_text_series(metadata.df, column, text_cache=text_cache)
            matched = matched | text_col.str.contains(compiled_pattern, regex=True).fillna(False)
        scope_key_series = get_select_scope_key_series(metadata.df, scope_columns, scope_key_cache=scope_key_cache)
        num_control = 0
        scope_values = scope_key_series.loc[matched]
        for scope_value in scope_values.unique().tolist():
            if not is_valid_select_scope_key(scope_value):
                continue
            is_scope = scope_key_series == scope_value
            num_control += int((is_scope & matched).sum())
        total_control += num_control
        total_matched += int(matched.sum())
        aggregate_key = (tuple(scope_columns), target_column, group['outcome'], group['stop_on_match'])
        aggregated_matches[aggregate_key] = matched
    total_protected = 0
    total_marked = 0
    for (scope_columns, target_column, outcome, stop_on_match), matched_any in aggregated_matches.items():
        scope_key_series = get_select_scope_key_series(metadata.df, list(scope_columns), scope_key_cache=scope_key_cache)
        scope_values = scope_key_series.loc[matched_any]
        num_control = 0
        num_treatment = 0
        for scope_value in scope_values.unique().tolist():
            if not is_valid_select_scope_key(scope_value):
                continue
            is_scope = scope_key_series == scope_value
            candidate_mask = is_scope & (~matched_any)
            write_mask = resolve_select_rule_write_mask(metadata.df, target_column, candidate_mask, stop_on_match)
            metadata.df.loc[write_mask, target_column] = outcome
            num_control += int((is_scope & matched_any).sum())
            num_treatment += int(write_mask.sum())
        total_protected += num_control
        total_marked += num_treatment
    print(
        '{}: Control rules complete: rules={}, groups={}, matched_total={:,}, control_total={:,}, protected_total={:,}, marked_total={:,}, scope_groups={}'.format(
            datetime.datetime.now(),
            len(control_rules),
            len(control_groups),
            total_matched,
            total_control,
            total_protected,
            total_marked,
            len(aggregated_matches),
        ),
        flush=True,
    )
    return metadata


def resolve_select_rule_runtime_column(column_template, parameter_value):
    if '{parameter}' not in column_template:
        return column_template
    return column_template.format(parameter=parameter_value)


def apply_select_filter_rules(metadata, args, select_rules):
    filter_rules = [rule for rule in select_rules if rule['stage'] == 'filter']
    if len(filter_rules) == 0:
        return metadata
    print('{}: Applying select filter rules'.format(datetime.datetime.now()), flush=True)
    metadata.df = ensure_select_target_column(metadata.df, SELECT_DEFAULT_TARGET_COLUMN)
    total_marked = 0
    for rule in filter_rules:
        target_column = rule['target_column'] if rule['target_column'] != '' else SELECT_DEFAULT_TARGET_COLUMN
        metadata.df = ensure_select_target_column(metadata.df, target_column)
        marked_mask = pandas.Series(False, index=metadata.df.index)
        if rule['action'] == 'exclude_if_lt_parameter':
            threshold = int(resolve_select_runtime_parameter(args, rule['parameter_name'], rule['rule_id']))
            column = rule['columns'][0]
            if column not in metadata.df.columns:
                raise ValueError(
                    'Column "{}" required by filter rule "{}" was not found in metadata.'.format(
                        column,
                        rule['rule_id'],
                    )
                )
            numeric_series = pandas.to_numeric(metadata.df[column], errors='coerce').fillna(0).astype(int)
            if threshold > 0:
                marked_mask = (numeric_series != 0) & (numeric_series < threshold)
        elif rule['action'] == 'exclude_if_missing_selected_rank':
            selected_rank = resolve_select_runtime_parameter(args, rule['parameter_name'], rule['rule_id'])
            if selected_rank != 'none':
                column = resolve_select_rule_runtime_column(rule['columns'][0], selected_rank)
                if column not in metadata.df.columns:
                    raise ValueError(
                        'Column "{}" required by filter rule "{}" was not found in metadata.'.format(
                            column,
                            rule['rule_id'],
                        )
                    )
                marked_mask = metadata.df[column].isna()
        elif rule['action'] == 'exclude_if_empty':
            empty_masks = []
            for column in rule['columns']:
                if column not in metadata.df.columns:
                    raise ValueError(
                        'Column "{}" required by filter rule "{}" was not found in metadata.'.format(
                            column,
                            rule['rule_id'],
                        )
                    )
                empty_masks.append(metadata.df[column].fillna('').astype(str).str.strip() == '')
            if len(empty_masks) > 0:
                marked_mask = empty_masks[0].copy()
                for extra_mask in empty_masks[1:]:
                    marked_mask = marked_mask & extra_mask
        else:
            raise ValueError(
                'Unsupported filter rule action "{}" in "{}".'.format(rule['action'], rule['rule_id'])
            )
        if bool(marked_mask.any()):
            metadata.df.loc[marked_mask, target_column] = rule['outcome']
            total_marked += int(marked_mask.sum())
    print(
        '{}: Filter rules complete: rules={}, marked_total={:,}'.format(
            datetime.datetime.now(),
            len(filter_rules),
            total_marked,
        ),
        flush=True,
    )
    return metadata


def apply_select_redundant_biosample_filter(metadata, rule, enabled):
    if not enabled:
        return metadata
    required_columns = list(rule['columns']) + [rule['target_column']]
    missing_columns = [column for column in required_columns if column not in metadata.df.columns]
    if len(missing_columns) > 0:
        raise ValueError(
            'Column(s) required for dedup rule "{}" were not found in metadata: {}'.format(
                rule['rule_id'],
                ', '.join(missing_columns),
            )
        )
    print('{}: Applying dedup rule "{}"'.format(datetime.datetime.now(), rule['rule_id']), flush=True)
    exclusion_series = metadata.df['exclusion'].fillna('no').astype(str).str.strip().str.lower()
    key_series = []
    valid_key_mask = pandas.Series(True, index=metadata.df.index)
    for column in rule['columns']:
        series = metadata.df[column].fillna('').astype(str).str.strip()
        key_series.append(series)
        valid_key_mask = valid_key_mask & (series != '')
    eligible_mask = (exclusion_series == 'no') & valid_key_mask
    if not bool(eligible_mask.any()):
        print(
            '{}: Dedup rule "{}" complete: candidates=0, redundant_marked=0'.format(
                datetime.datetime.now(),
                rule['rule_id'],
            ),
            flush=True,
        )
        return metadata
    ranking_series = pandas.to_numeric(metadata.df[rule['target_column']], errors='coerce').fillna(0).astype(int)
    if 'run' in metadata.df.columns:
        run_series = metadata.df['run'].fillna('').astype(str).str.strip()
    else:
        run_series = pandas.Series('', index=metadata.df.index, dtype=str)
    ranking_df = pandas.DataFrame({
        'orig_index': metadata.df.index,
        'run': run_series,
        'ranking_value': ranking_series,
    }).loc[eligible_mask, :].copy()
    subset_columns = []
    sort_columns = []
    ascending = []
    for column, series in zip(rule['columns'], key_series):
        ranking_df[column] = series.loc[eligible_mask]
        subset_columns.append(column)
        sort_columns.append(column)
        ascending.append(True)
    ranking_df = ranking_df.sort_values(
        sort_columns + ['ranking_value', 'run', 'orig_index'],
        ascending=ascending + [False, True, True],
        kind='mergesort',
    )
    redundant_mask = ranking_df.duplicated(subset=subset_columns, keep='first')
    redundant_indices = ranking_df.loc[redundant_mask, 'orig_index']
    if len(redundant_indices) > 0:
        metadata.df.loc[redundant_indices, 'exclusion'] = rule['outcome']
    print(
        '{}: Dedup rule "{}" complete: candidates={:,}, redundant_marked={:,}'.format(
            datetime.datetime.now(),
            rule['rule_id'],
            int(eligible_mask.sum()),
            int(len(redundant_indices)),
        ),
        flush=True,
    )
    return metadata


def apply_select_filters(metadata, args, select_rules):
    metadata = apply_select_exclude_rules(metadata, select_rules)
    metadata = apply_select_control_rules(metadata, select_rules)
    metadata = apply_select_filter_rules(metadata, args, select_rules)
    dedup_rules = [rule for rule in select_rules if rule['stage'] == 'dedup']
    for rule in dedup_rules:
        enabled = True
        if rule['parameter_name'] != '':
            enabled = bool(resolve_select_runtime_parameter(args, rule['parameter_name'], rule['rule_id']))
        metadata = apply_select_redundant_biosample_filter(metadata, rule, enabled)
    metadata.label_sampled_data(
        args.max_sample,
        sampling_strategy=getattr(args, 'sampling_strategy', 'maximize_bioproject_diversity'),
    )
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


def normalize_select_series(series):
    return (
        series.fillna('')
        .astype(str)
        .str.replace('_', ' ', regex=False)
        .str.replace(r'\s+', ' ', regex=True)
        .str.strip()
    )


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
    organs = ['flower', 'leaf', 'root']
    for organ in organs:
        organ_mask = (sample_group_series == organ)
        row['{}_rows'.format(organ)] = int(organ_mask.sum())
        row['{}_qualified_yes'.format(organ)] = int((organ_mask & (qualified_series == 'yes')).sum())
        row['{}_sampled_yes'.format(organ)] = int((organ_mask & (sampled_series == 'yes')).sum())
    row['qualified_yes_total'] = int((qualified_series == 'yes').sum())
    row['sampled_yes_total'] = int((sampled_series == 'yes').sum())
    strict_threshold, moderate_threshold = resolve_select_queue_thresholds(max_sample)
    row['any_tissues_ge1'] = any(row['{}_sampled_yes'.format(organ)] >= 1 for organ in organs)
    row['all_tissues_ge1'] = all(row['{}_sampled_yes'.format(organ)] >= 1 for organ in organs)
    row['all_tissues_ge30'] = all(row['{}_sampled_yes'.format(organ)] >= strict_threshold for organ in organs)
    row['all_tissues_ge3'] = all(row['{}_sampled_yes'.format(organ)] >= moderate_threshold for organ in organs)
    row['queue_tier'] = resolve_select_queue_tier(row)
    row['selected_metadata_path'] = selected_metadata_path
    return row


def resolve_select_queue_thresholds(max_sample):
    max_sample_int = int(max_sample)
    return min(max_sample_int, 30), min(max_sample_int, 3)


def resolve_select_queue_tier(row):
    if row['all_tissues_ge30']:
        return 'all_tissues_ge30'
    if row['all_tissues_ge3']:
        return 'all_tissues_ge3'
    if row['all_tissues_ge1']:
        return 'all_tissues_ge1'
    if row['any_tissues_ge1']:
        return 'any_tissues_ge1'
    return 'no_tissues_ge1'


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
        'any_tissues_ge1',
        'all_tissues_ge1',
        'all_tissues_ge30',
        'all_tissues_ge3',
        'queue_tier',
        'selected_metadata_path',
    ]].copy(deep=True)
    queue_df['queue_sort'] = queue_df['queue_tier'].map({
        'all_tissues_ge30': 0,
        'all_tissues_ge3': 1,
        'all_tissues_ge1': 2,
        'any_tissues_ge1': 3,
        'no_tissues_ge1': 4,
    }).fillna(9)
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
            'all_tissues_ge30': 'external_run_now',
            'all_tissues_ge3': 'external_run_if_capacity',
            'all_tissues_ge1': 'external_run_all_tissues_ge1_panel',
            'any_tissues_ge1': 'external_run_any_tissues_ge1_panel',
            'no_tissues_ge1': 'metadata_gap_review',
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
        base + '_all_tissues_ge30' + ext,
        base + '_all_tissues_ge3' + ext,
        base + '_all_tissues_ge1' + ext,
        base + '_any_tissues_ge1' + ext,
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


def _resolve_select_species_jobs(args, task_count):
    species_jobs, _ = resolve_worker_allocation(
        requested_workers=getattr(args, 'internal_jobs', 'auto'),
        requested_threads=getattr(args, 'threads', 'auto'),
        internal_cpu_budget=getattr(args, 'internal_cpu_budget', 'auto'),
        worker_option_name='internal_jobs',
        context='select:',
    )
    if species_jobs > int(task_count):
        species_jobs = int(task_count)
    if species_jobs <= 0:
        species_jobs = 1
    return species_jobs


def _run_select_batch_task(task, args, species_count, select_rules, select_parameters, metadata_specieswise_dir):
    index, species_name, species_token = task
    print('[{}/{}] {}'.format(index, species_count, species_token), flush=True)
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
    species_out_dir = os.path.join(args.out_dir, species_token)
    metadata_dir = os.path.join(species_out_dir, 'metadata')
    normalized_metadata_path = os.path.join(metadata_dir, 'metadata.tsv')
    normalization_row = summarize_select_normalization(
        species_name=species_name,
        species_token=species_token,
        df=normalized_df,
        workspace_dir=species_out_dir,
        metadata_path=normalized_metadata_path,
    )
    runtime_args = apply_select_config_parameters(
        clone_namespace(args, out_dir=species_out_dir),
        select_parameters,
    )
    metadata_original_df = metadata.df.copy(deep=True)
    metadata = filter_metadata_by_sample_group(metadata, getattr(runtime_args, 'sample_group', None))
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
    summary_row = summarize_selected_metadata(
        species_name=species_name,
        species_token=species_token,
        metadata=metadata,
        selected_metadata_path=path_metadata_table,
        max_sample=runtime_args.max_sample,
    )
    return {
        'normalization_row': normalization_row,
        'summary_row': summary_row,
    }


def select_batch_main(args):
    out_dir = os.path.realpath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    select_rules_tsv = resolve_select_rules_tsv(args, out_dir=out_dir)
    select_config = read_select_config(select_rules_tsv)
    select_rules = select_config['rules']
    select_parameters = select_config['parameters']
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
    print('select batch: processing {} species'.format(species_df.shape[0]), flush=True)
    task_items = [
        (index + 1, row['scientific_name'], row['species_token'])
        for index, row in species_df.iterrows()
    ]
    species_jobs = _resolve_select_species_jobs(args=args, task_count=len(task_items))
    task_results = {}
    if (species_jobs == 1) or (len(task_items) <= 1):
        for task in task_items:
            task_results[task] = _run_select_batch_task(
                task=task,
                args=args,
                species_count=species_df.shape[0],
                select_rules=select_rules,
                select_parameters=select_parameters,
                metadata_specieswise_dir=metadata_specieswise_dir,
            )
    else:
        print(
            'select: running {:,} species with {:,} parallel job(s).'.format(
                len(task_items),
                species_jobs,
            ),
            flush=True,
        )
        task_results, failures = run_tasks_with_optional_threads(
            task_items=task_items,
            task_fn=lambda task: _run_select_batch_task(
                task=task,
                args=args,
                species_count=species_df.shape[0],
                select_rules=select_rules,
                select_parameters=select_parameters,
                metadata_specieswise_dir=metadata_specieswise_dir,
            ),
            max_workers=species_jobs,
        )
        if failures:
            details = '; '.join(
                [
                    '{}: {}'.format(task[2] if isinstance(task, tuple) and len(task) >= 3 else task, err)
                    for task, err in failures
                ]
            )
            raise RuntimeError(
                'select batch failed for {}/{} species. {}'.format(
                    len(failures),
                    len(task_items),
                    details,
                )
            )
    normalization_rows = [task_results[task]['normalization_row'] for task in task_items]
    summary_rows = [task_results[task]['summary_row'] for task in task_items]
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
            'any_tissues_ge1',
            'all_tissues_ge1',
            'all_tissues_ge30',
            'all_tissues_ge3',
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
    ge30_manifest_tsv, ge3_manifest_tsv, all_ge1_manifest_tsv, any_ge1_manifest_tsv = manifest_sidecar_paths(manifest_tsv)
    atomic_write_dataframe(manifest_df.loc[manifest_df['queue_tier'] == 'all_tissues_ge30', :], ge30_manifest_tsv, sep='\t', index=False)
    atomic_write_dataframe(manifest_df.loc[manifest_df['queue_tier'] == 'all_tissues_ge3', :], ge3_manifest_tsv, sep='\t', index=False)
    atomic_write_dataframe(manifest_df.loc[manifest_df['queue_tier'] == 'all_tissues_ge1', :], all_ge1_manifest_tsv, sep='\t', index=False)
    atomic_write_dataframe(manifest_df.loc[manifest_df['queue_tier'] == 'any_tissues_ge1', :], any_ge1_manifest_tsv, sep='\t', index=False)
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
    select_config = read_select_config(select_rules_tsv)
    select_rules = select_config['rules']
    runtime_args = apply_select_config_parameters(runtime_args, select_config['parameters'])
    path_metadata_table = os.path.join(metadata_dir, 'metadata.tsv')
    path_metadata_original = os.path.join(metadata_dir, 'metadata_original.tsv')

    metadata = load_metadata(runtime_args)
    metadata = prepare_select_metadata(metadata, select_rules)
    metadata_original_df = metadata.df.copy(deep=True)
    metadata = filter_metadata_by_sample_group(metadata, getattr(runtime_args, 'sample_group', None))
    metadata = apply_select_filters(metadata, runtime_args, select_rules)
    write_select_outputs(
        path_metadata_original=path_metadata_original,
        path_metadata_table=path_metadata_table,
        metadata_dir=metadata_dir,
        metadata=metadata,
        metadata_original_df=metadata_original_df,
    )
