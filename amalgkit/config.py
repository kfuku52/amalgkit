import os
import sys

try:
    import importlib.resources as ir
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as ir

def check_directory(args):
    path_config = os.path.join(args.out_dir, 'config_'+args.config)
    if os.path.exists(path_config):
        print('Output directory already exists: {}'.format(path_config))
        if args.overwrite:
            print('--overwrite is set to "yes". Any config files will be overwritten in: {}'.format(path_config))
        else:
            print('--overwrite is set to "no". Exiting.')
            sys.exit()
    else:
        os.makedirs(path_config, exist_ok=True)
    return path_config

def create_config_from_package(args, path_config):
    config_base = 'amalgkit.config_dir' + '.' + args.config
    config_files = ir.files(config_base).rglob('*.config')
    for config_file in config_files:
        file_content = config_file.read_bytes()
        print('Copying from {} to {}'.format(config_file, path_config))
        with open(os.path.join(path_config, config_file.name), mode='wb') as f:
            f.write(file_content)

def config_main(args):
    path_config = check_directory(args)
    create_config_from_package(args, path_config)
