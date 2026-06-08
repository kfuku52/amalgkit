import os
import shutil
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
