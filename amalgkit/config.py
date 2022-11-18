import os
import sys
from glob import glob


def check_directory(args):
    if args.config_dir == 'inferred':
        config_dir = args.config
    else:
        config_dir = args.config_dir
    work_dir = os.path.join(args.out_dir, 'config')
    path_config_dir = os.path.join(work_dir, config_dir)

    if not os.path.exists(os.path.join(work_dir)):
        os.makedirs(os.path.join(work_dir))

    print('Checking config directory ...')
    if os.path.exists(path_config_dir):
        print(path_config_dir, ' already exists.')
        if args.overwrite:
            print('--overwrite is set! Any config files in ', path_config_dir, ' will be overwritten.')
            print('Continuing.')
        else:
            print('Exiting.')
            sys.exit()
    else:
        os.makedirs(path_config_dir)

    return


def check_files(args, config_file_list, ext='.config'):
    work_dir = os.path.join(args.out_dir, 'config')
    print('Checking if all config files are present in ', os.path.join(work_dir, args.config_dir), '...')
    missing_files = []
    for config_file in config_file_list:
        path_config_file = os.path.join(work_dir, args.config_dir, config_file + ext)
        if os.path.exists(path_config_file):
            print('found ', config_file + ext, '...')
        else:
            print('could not find ', config_file + ext, 'in ', os.path.join(work_dir, args.config_dir))
            print('make sure amalgkit config ran correctly')
            missing_files.append(config_file + ext)
    if len(missing_files) > 0:
        print('some files are missing: ')
        for missing_file in missing_files:
            print(missing_file)
    else:
        print('All config files are present!')
    return


def create_config_from_package(args):

    if args.config_dir == 'inferred':
        config_dir = args.config
    else:
        config_dir = args.config_dir

    work_dir = os.path.join(args.out_dir, 'config')
    path_config = os.path.join(work_dir, config_dir)

    try:
        import importlib.resources as ir
    except ImportError:
        # Try backported to PY<37 `importlib_resources`.
        import importlib_resources as ir

    config_base = 'config_dir' + '.' + args.config
    config_files = ir.files(config_base).rglob('*.config')

    for config_file in config_files:
        try:
            # extracts the contents of a single file from the config_dir within the amalgkit package
            file_content = ir.files(config_base).joinpath(os.path.basename(config_file)).read_bytes()
        except:
            continue
        print('Writing {} for the dataset {}'.format(config_file, args.config))
        with open(os.path.join(path_config, os.path.basename(config_file)), mode='wb') as f:
            f.write(file_content)


def config_main(args):
    check_directory(args)
    create_config_from_package(args)