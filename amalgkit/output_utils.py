import os
import re
import tempfile
from contextlib import contextmanager

import pandas


_DELIMITED_TEXT_UNSAFE_PATTERN = re.compile(r'[\r\n\t]+')


@contextmanager
def atomic_output_path(outpath, prefix='amalgkit_atomic_', suffix=None):
    real_outpath = os.path.realpath(outpath)
    parent_dir = os.path.dirname(real_outpath)
    if parent_dir == '':
        parent_dir = '.'
    if os.path.exists(parent_dir) and (not os.path.isdir(parent_dir)):
        raise NotADirectoryError('Output parent path exists but is not a directory: {}'.format(parent_dir))
    os.makedirs(parent_dir, exist_ok=True)
    if os.path.exists(real_outpath) and (not os.path.isfile(real_outpath)):
        raise IsADirectoryError('Output path exists but is not a file: {}'.format(real_outpath))
    tmp_suffix = suffix
    if tmp_suffix is None:
        tmp_suffix = os.path.splitext(real_outpath)[1]
    fd, tmp_path = tempfile.mkstemp(prefix=prefix, suffix=tmp_suffix, dir=parent_dir)
    os.close(fd)
    committed = False
    try:
        yield tmp_path
        os.replace(tmp_path, real_outpath)
        committed = True
    finally:
        if (not committed) and os.path.exists(tmp_path):
            os.remove(tmp_path)


def atomic_write_dataframe(df, outpath, **to_csv_kwargs):
    sep = to_csv_kwargs.get('sep', ',')
    df_to_write = df
    if sep == '\t':
        df_to_write = sanitize_dataframe_for_tsv(df)
    with atomic_output_path(outpath=outpath, suffix='.tsv') as tmp_path:
        df_to_write.to_csv(tmp_path, **to_csv_kwargs)


def sanitize_delimited_text_cell(value):
    if not isinstance(value, str):
        return value
    if value == '':
        return value
    sanitized = _DELIMITED_TEXT_UNSAFE_PATTERN.sub(' ', value)
    return sanitized.strip()


def sanitize_dataframe_for_tsv(df):
    sanitized = df.copy()
    text_columns = sanitized.select_dtypes(include=['object', 'string']).columns
    for column in text_columns:
        sanitized[column] = sanitized[column].map(sanitize_delimited_text_cell)
    return sanitized
