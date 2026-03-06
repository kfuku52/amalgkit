import errno
import os
import time
from contextlib import contextmanager

import ete4


DOWNLOAD_LOCK_POLL_SECONDS = 5
DOWNLOAD_LOCK_TIMEOUT_SECONDS = 3600


def resolve_download_dir(args):
    out_dir = getattr(args, 'out_dir', './')
    inferred = os.path.join(os.path.realpath(out_dir), 'downloads')
    raw_dir = getattr(args, 'download_dir', 'inferred')
    if raw_dir is None:
        return inferred
    normalized = str(raw_dir).strip()
    if normalized.lower() in ['', 'inferred']:
        return inferred
    return os.path.realpath(normalized)


def resolve_ete_data_dir(args, resolve_download_dir_fn=None):
    if resolve_download_dir_fn is None:
        resolve_download_dir_fn = resolve_download_dir
    return os.path.join(resolve_download_dir_fn(args), 'ete4')


def _assert_lock_path_is_regular_file(lock_path, lock_label='Lock'):
    if not os.path.lexists(lock_path):
        return
    if os.path.islink(lock_path) or (not os.path.isfile(lock_path)):
        raise IsADirectoryError('{} path exists but is not a file: {}'.format(lock_label, lock_path))


def _try_create_lock_file(lock_path):
    _assert_lock_path_is_regular_file(lock_path)
    try:
        fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, 'w') as lock_handle:
        lock_handle.write('{}\n'.format(os.getpid()))
    return True


def _read_lock_owner_pid(lock_path):
    try:
        with open(lock_path) as lock_handle:
            first_line = lock_handle.readline().strip()
    except OSError:
        return None
    if first_line == '':
        return None
    try:
        pid = int(first_line)
    except ValueError:
        return None
    if pid <= 0:
        return None
    return pid


def _is_process_alive(pid):
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError as exc:
        if getattr(exc, 'errno', None) == errno.ESRCH:
            return False
        return True
    return True


def _break_stale_lock_if_needed(lock_path, lock_label='Lock'):
    if not os.path.lexists(lock_path):
        return False
    _assert_lock_path_is_regular_file(lock_path, lock_label=lock_label)
    try:
        stat_before = os.stat(lock_path)
    except FileNotFoundError:
        return False
    owner_pid = _read_lock_owner_pid(lock_path)
    if owner_pid is None:
        stale_reason = 'missing/invalid owner PID'
    elif _is_process_alive(owner_pid):
        return False
    else:
        stale_reason = 'owner PID {} is not running'.format(owner_pid)
    try:
        stat_now = os.stat(lock_path)
    except FileNotFoundError:
        return False
    if (
        (stat_before.st_ino != stat_now.st_ino)
        or (stat_before.st_size != stat_now.st_size)
        or (stat_before.st_mtime_ns != stat_now.st_mtime_ns)
    ):
        return False
    try:
        os.remove(lock_path)
    except FileNotFoundError:
        return False
    print(
        'Removed stale {} lock: {} ({})'.format(
            lock_label,
            lock_path,
            stale_reason,
        ),
        flush=True,
    )
    return True


@contextmanager
def acquire_exclusive_lock(
    lock_path,
    lock_label='Lock',
    poll_seconds=DOWNLOAD_LOCK_POLL_SECONDS,
    timeout_seconds=DOWNLOAD_LOCK_TIMEOUT_SECONDS,
):
    poll_seconds = int(poll_seconds)
    timeout_seconds = int(timeout_seconds)
    if poll_seconds <= 0:
        raise ValueError('poll_seconds must be > 0.')
    if timeout_seconds <= 0:
        raise ValueError('timeout_seconds must be > 0.')
    lock_path = os.path.realpath(lock_path)
    lock_dir = os.path.dirname(lock_path)
    if lock_dir != '':
        if os.path.exists(lock_dir) and (not os.path.isdir(lock_dir)):
            raise NotADirectoryError('Lock parent path exists but is not a directory: {}'.format(lock_dir))
        os.makedirs(lock_dir, exist_ok=True)
    wait_start = time.time()
    has_reported_wait = False
    while True:
        if _try_create_lock_file(lock_path):
            try:
                yield
            finally:
                if os.path.lexists(lock_path):
                    _assert_lock_path_is_regular_file(lock_path, lock_label=lock_label)
                    os.remove(lock_path)
            return
        _assert_lock_path_is_regular_file(lock_path, lock_label=lock_label)
        if _break_stale_lock_if_needed(lock_path=lock_path, lock_label=lock_label):
            continue
        elapsed = time.time() - wait_start
        if not has_reported_wait:
            print(
                'Another process holds {}. Waiting for lock release: {}'.format(
                    lock_label,
                    lock_path,
                ),
                flush=True,
            )
            has_reported_wait = True
        if elapsed > timeout_seconds:
            raise TimeoutError(
                'Timed out after {:,} sec waiting for {} lock: {}'.format(
                    timeout_seconds,
                    lock_label,
                    lock_path,
                )
            )
        time.sleep(poll_seconds)


def get_ete_ncbitaxa(
    args=None,
    acquire_exclusive_lock_fn=None,
    ncbitaxa_cls=None,
    resolve_ete_data_dir_fn=None,
):
    if acquire_exclusive_lock_fn is None:
        acquire_exclusive_lock_fn = acquire_exclusive_lock
    if ncbitaxa_cls is None:
        ncbitaxa_cls = ete4.NCBITaxa
    if resolve_ete_data_dir_fn is None:
        resolve_ete_data_dir_fn = resolve_ete_data_dir
    if args is None:
        return ncbitaxa_cls()
    ete_data_dir = resolve_ete_data_dir_fn(args)
    os.makedirs(ete_data_dir, exist_ok=True)
    lock_path = os.path.join(ete_data_dir, '.ete4_taxonomy.lock')
    with acquire_exclusive_lock_fn(lock_path=lock_path, lock_label='ETE4 taxonomy DB'):
        return ncbitaxa_cls(
            dbfile=os.path.join(ete_data_dir, 'taxa.sqlite'),
            taxdump_file=os.path.join(ete_data_dir, 'taxdump.tar.gz'),
        )
