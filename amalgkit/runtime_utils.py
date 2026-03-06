import os
import shutil
import subprocess
import sys

from amalgkit.subprocess_utils import probe_dependency_command


def check_rscript(runner=subprocess.run):
    try:
        probe_dependency_command(
            command=['Rscript', '--help'],
            label='Rscript',
            runner=runner,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError('R (Rscript) is not installed.') from exc


def check_config_dir(dir_path, mode):
    mode_to_files = {
        'select': [
            'group_attribute.config',
            'exclude_keyword.config',
            'control_term.config',
        ],
    }
    asserted_files = mode_to_files.get(mode)
    if asserted_files is None:
        raise ValueError('Unsupported config check mode: {}'.format(mode))
    if not os.path.exists(dir_path):
        raise FileNotFoundError('Config directory not found: {}'.format(dir_path))
    if not os.path.isdir(dir_path):
        raise NotADirectoryError('Config path exists but is not a directory: {}'.format(dir_path))
    missing_count = 0
    for af in asserted_files:
        af_path = os.path.join(dir_path, af)
        if os.path.isfile(af_path):
            print('Config file found: {}'.format(af))
        elif os.path.exists(af_path):
            sys.stderr.write('Config entry exists but is not a file: {}\n'.format(af))
            missing_count += 1
        else:
            sys.stderr.write('Config file not found: {}\n'.format(af))
            missing_count += 1
    if (missing_count > 0):
        txt = 'Please refer to the AMALGKIT Wiki for more info: https://github.com/kfuku52/amalgkit/wiki/amalgkit-metadata\n'
        sys.stderr.write(txt)


def cleanup_tmp_amalgkit_files(work_dir='.'):
    try:
        with os.scandir(work_dir) as entries:
            for entry in entries:
                if not entry.name.startswith('tmp.amalgkit.'):
                    continue
                path = entry.path
                try:
                    if entry.is_dir(follow_symlinks=False):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
                except FileNotFoundError:
                    continue
    except FileNotFoundError:
        return


def get_getfastq_run_dir(args, sra_id):
    amalgkit_out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(amalgkit_out_dir) and (not os.path.isdir(amalgkit_out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(amalgkit_out_dir))
    getfastq_dir = os.path.join(amalgkit_out_dir, 'getfastq')
    if os.path.exists(getfastq_dir) and (not os.path.isdir(getfastq_dir)):
        raise NotADirectoryError('getfastq path exists but is not a directory: {}'.format(getfastq_dir))
    os.makedirs(getfastq_dir, exist_ok=True)
    path_run = os.path.join(getfastq_dir, sra_id)
    if os.path.exists(path_run) and (not os.path.isdir(path_run)):
        raise NotADirectoryError('Run path exists but is not a directory: {}'.format(path_run))
    os.makedirs(path_run, exist_ok=True)
    return path_run
