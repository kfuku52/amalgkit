import os
import tempfile
from contextlib import contextmanager


def _serialize_r_config_value(value):
    if value is None:
        return ''
    if isinstance(value, bool):
        return '1' if value else '0'
    if isinstance(value, (list, tuple)):
        return '|'.join([_serialize_r_config_value(item) for item in value])
    return str(value).replace('\n', '\\n')


@contextmanager
def temporary_r_config(config_map, prefix='amalgkit_r_'):
    fd, config_path = tempfile.mkstemp(prefix=prefix, suffix='.dcf')
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            for key, value in config_map.items():
                handle.write('{}: {}\n'.format(str(key), _serialize_r_config_value(value)))
        yield config_path
    finally:
        if os.path.isfile(config_path):
            os.remove(config_path)
