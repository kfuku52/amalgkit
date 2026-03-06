import os
import tempfile
from contextlib import contextmanager


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
    with atomic_output_path(outpath=outpath, suffix='.tsv') as tmp_path:
        df.to_csv(tmp_path, **to_csv_kwargs)
