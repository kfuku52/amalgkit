import json

import pandas

from amalgkit.batch_effect_common import BatchEffectResult


_DCF_FIELD_TYPES = {
    'stable': 'bool',
    'sva_stable': 'bool',
    'group_model_used': 'bool',
    'group_fallback_used': 'bool',
    'ruv_fallback_used': 'bool',
    'latent_converged': 'bool',
    'corrected_run_ids': 'str_list',
    'uncorrected_run_ids': 'str_list',
    'trace_method': 'str_list',
    'trace_B': 'int_list',
    'trace_nsv': 'int_list',
    'counts_shape': 'int_list',
    'resolved_sva_nsv': 'int',
    'resolved_sva_B': 'int',
    'resolved_ruv_k': 'int',
    'resolved_ruv_controls': 'int',
    'resolved_latent_k': 'int',
    'latent_iterations': 'int',
    'negative_values_before_clip': 'int',
    'negative_values_after_clip': 'int',
    'ruv_nb_fallback_genes': 'int',
    'ruv_anova_failure_genes': 'int',
    'metadata_rows': 'int',
    'latent_objective': 'float',
    'ruv_baseline_score': 'float',
    'ruv_selected_score': 'float',
    'ruv_selected_penalized_score': 'float',
    'ruv_penalty': 'float',
    'latent_family': 'nullable_str',
    'group_error_message': 'nullable_str',
    'ruv_residual_method': 'nullable_str',
    'ruv_pvalue_method': 'nullable_str',
    'ruv_fallback_reason': 'nullable_str',
}


def read_expression_matrix_tsv(path):
    return pandas.read_csv(path, sep='\t', index_col=0)


def write_expression_matrix_tsv(df, path):
    df.to_csv(path, sep='\t', index=True)


def read_metadata_tsv(path):
    return pandas.read_csv(path, sep='\t', low_memory=False)


def write_backend_summary_json(summary, path):
    if isinstance(summary, BatchEffectResult):
        payload = summary.to_jsonable()
    else:
        payload = dict(summary)
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write('\n')


def read_backend_summary_json(path):
    with open(path, encoding='utf-8') as handle:
        return json.load(handle)


def write_backend_summary_dcf(summary, path):
    if isinstance(summary, BatchEffectResult):
        payload = summary.to_jsonable()
    else:
        payload = dict(summary)
    with open(path, 'w', encoding='utf-8') as handle:
        for key, value in payload.items():
            if value is None:
                text = ''
            elif isinstance(value, bool):
                text = 'TRUE' if value else 'FALSE'
            elif isinstance(value, (list, tuple)):
                text = '|'.join([str(item) for item in value])
            else:
                text = str(value)
            handle.write('{}: {}\n'.format(str(key), text))


def _parse_backend_summary_dcf_value(key, value):
    field_type = _DCF_FIELD_TYPES.get(key)
    if field_type is None:
        return value
    if field_type == 'nullable_str':
        return None if value == '' else value
    if field_type == 'bool':
        if value == '':
            return None
        normalized = value.upper()
        if normalized == 'TRUE':
            return True
        if normalized == 'FALSE':
            return False
        raise ValueError('Invalid boolean value for DCF field {}: {!r}'.format(key, value))
    if field_type == 'str_list':
        return [] if value == '' else value.split('|')
    if field_type == 'int_list':
        if value == '':
            return []
        try:
            return [int(item) for item in value.split('|')]
        except ValueError as exc:
            raise ValueError('Invalid integer-list value for DCF field {}: {!r}'.format(key, value)) from exc
    if value == '':
        return None
    try:
        if field_type == 'int':
            return int(value)
        if field_type == 'float':
            return float(value)
    except ValueError as exc:
        raise ValueError('Invalid {} value for DCF field {}: {!r}'.format(field_type, key, value)) from exc
    raise RuntimeError('Unsupported DCF field type: {}'.format(field_type))


def read_backend_summary_dcf(path):
    payload = {}
    with open(path, encoding='utf-8') as handle:
        for raw_line in handle:
            line = raw_line.rstrip('\n')
            if line == '':
                continue
            key, value = line.split(':', 1)
            value = value.lstrip()
            payload[key] = _parse_backend_summary_dcf_value(key=key, value=value)
    return payload
