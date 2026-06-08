import json

import pandas

from amalgkit.batch_effect_common import BatchEffectResult


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


def read_backend_summary_dcf(path):
    payload = {}
    with open(path, encoding='utf-8') as handle:
        for raw_line in handle:
            line = raw_line.rstrip('\n')
            if line == '':
                continue
            key, value = line.split(':', 1)
            payload[key] = value.lstrip()
    return payload
