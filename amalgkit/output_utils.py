import os
import re
import tempfile
import uuid
from contextlib import contextmanager

_DELIMITED_TEXT_UNSAFE_PATTERN = re.compile(r'[\r\n\t]+')


def get_default_creation_mode(parent_dir, is_directory):
    probe_path = os.path.join(
        parent_dir,
        '.amalgkit_mode_probe_{}'.format(uuid.uuid4().hex),
    )
    if is_directory:
        os.mkdir(probe_path, mode=0o777)
        try:
            return os.stat(probe_path, follow_symlinks=False).st_mode & 0o777
        finally:
            os.rmdir(probe_path)
    descriptor = os.open(probe_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o666)
    os.close(descriptor)
    try:
        return os.stat(probe_path, follow_symlinks=False).st_mode & 0o777
    finally:
        os.remove(probe_path)


@contextmanager
def atomic_output_path(outpath, prefix='amalgkit_atomic_', suffix=None):
    absolute_outpath = os.path.abspath(outpath)
    if os.path.lexists(absolute_outpath) and os.path.islink(absolute_outpath):
        raise ValueError('Refusing to replace symbolic-link output path: {}'.format(absolute_outpath))
    parent_dir = os.path.realpath(os.path.dirname(absolute_outpath))
    if parent_dir == '':
        parent_dir = '.'
    real_outpath = os.path.join(parent_dir, os.path.basename(absolute_outpath))
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
    output_mode = get_default_creation_mode(parent_dir, is_directory=False)
    if os.path.isfile(real_outpath):
        output_mode = os.stat(real_outpath, follow_symlinks=False).st_mode & 0o777
    committed = False
    try:
        yield tmp_path
        os.chmod(tmp_path, output_mode)
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
