import os
import sys

try:
    import importlib.resources as ir
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as ir

def check_directory(args):
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    path_config = os.path.join(out_dir, 'config_'+args.config)
    if os.path.exists(path_config):
        if not os.path.isdir(path_config):
            raise NotADirectoryError('Output config path exists but is not a directory: {}'.format(path_config))
        print('Output directory already exists: {}'.format(path_config))
        if args.overwrite:
            print('--overwrite is set to "yes". Any config files will be overwritten in: {}'.format(path_config))
        else:
            print('--overwrite is set to "no". Exiting.')
            sys.exit()
    else:
        os.makedirs(path_config, exist_ok=True)
    return path_config


def list_available_config_sets():
    root = ir.files('amalgkit.config_dir')
    config_sets = []
    for entry in root.iterdir():
        if not entry.is_dir():
            continue
        if entry.name.startswith('_'):
            continue
        has_config_file = any(
            child.is_file() and child.name.endswith('.config')
            for child in entry.iterdir()
        )
        if has_config_file:
            config_sets.append(entry.name)
    return sorted(config_sets)


def validate_config_set(config_name):
    available = list_available_config_sets()
    if config_name in available:
        return
    raise ValueError(
        'Unknown config set "{}". Available config sets: {}'.format(
            config_name,
            ', '.join(available),
        )
    )

def create_config_from_package(args, path_config):
    config_base = 'amalgkit.config_dir' + '.' + args.config
    try:
        config_files = ir.files(config_base).rglob('*.config')
    except ModuleNotFoundError as exc:
        available = ', '.join(list_available_config_sets())
        raise ValueError(
            'Unknown config set "{}". Available config sets: {}'.format(args.config, available)
        ) from exc
    for config_file in config_files:
        file_content = config_file.read_bytes()
        print('Copying from {} to {}'.format(config_file, path_config))
        dst_path = os.path.join(path_config, config_file.name)
        if os.path.exists(dst_path) and (not os.path.isfile(dst_path)):
            raise IsADirectoryError('Config output path exists but is not a file: {}'.format(dst_path))
        with open(dst_path, mode='wb') as f:
            f.write(file_content)

def config_main(args):
    validate_config_set(args.config)
    path_config = check_directory(args)
    create_config_from_package(args, path_config)
