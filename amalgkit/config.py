import os
import sys

from amalgkit.exceptions import AmalgkitExit

SELECT_RULES_FILENAME = 'select_rules.tsv'

try:
    import importlib.resources as ir
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as ir

def check_directory(args):
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    os.makedirs(out_dir, exist_ok=True)
    path_config = os.path.join(out_dir, SELECT_RULES_FILENAME)
    if os.path.exists(path_config):
        if not os.path.isfile(path_config):
            raise IsADirectoryError('Config output path exists but is not a file: {}'.format(path_config))
        print('Output select rules already exists: {}'.format(path_config))
        if args.overwrite:
            print('--overwrite is set to "yes". The select rules file will be overwritten: {}'.format(path_config))
        else:
            raise AmalgkitExit('--overwrite is set to "no". Exiting.', exit_code=0, use_stderr=False)
    return path_config


def list_available_config_sets():
    root = ir.files('amalgkit.config_dir')
    config_sets = []
    for entry in root.iterdir():
        if not entry.is_dir():
            continue
        if entry.name.startswith('_'):
            continue
        if entry.joinpath(SELECT_RULES_FILENAME).is_file():
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
        config_root = ir.files(config_base)
    except ModuleNotFoundError as exc:
        available = ', '.join(list_available_config_sets())
        raise ValueError(
            'Unknown config set "{}". Available config sets: {}'.format(args.config, available)
        ) from exc
    rules_file = config_root.joinpath(SELECT_RULES_FILENAME)
    if not rules_file.is_file():
        raise FileNotFoundError(
            'Bundled select_rules.tsv not found for config set "{}".'.format(args.config)
        )
    file_content = rules_file.read_bytes()
    print('Copying from {} to {}'.format(rules_file, path_config))
    with open(path_config, mode='wb') as f:
        f.write(file_content)

def config_main(args):
    validate_config_set(args.config)
    path_config = check_directory(args)
    create_config_from_package(args, path_config)
