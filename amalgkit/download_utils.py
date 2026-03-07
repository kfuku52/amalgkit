import errno
import os
import time
import urllib.request
from contextlib import contextmanager

import ete4


DOWNLOAD_LOCK_POLL_SECONDS = 5
DOWNLOAD_LOCK_TIMEOUT_SECONDS = 3600
NCBI_TAXDUMP_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'


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


def _assert_regular_file_or_absent(path, label='Path'):
    if not os.path.lexists(path):
        return
    if os.path.islink(path) or (not os.path.isfile(path)):
        raise IsADirectoryError('{} exists but is not a file: {}'.format(label, path))


def _assert_lock_path_is_regular_file(lock_path, lock_label='Lock'):
    _assert_regular_file_or_absent(lock_path, label='{} path'.format(lock_label))


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


def ensure_ete_taxdump_file(taxdump_path, urlretrieve_fn=None):
    if urlretrieve_fn is None:
        urlretrieve_fn = urllib.request.urlretrieve
    taxdump_path = os.path.realpath(taxdump_path)
    taxdump_dir = os.path.dirname(taxdump_path)
    if taxdump_dir != '':
        if os.path.exists(taxdump_dir) and (not os.path.isdir(taxdump_dir)):
            raise NotADirectoryError(
                'ETE4 taxdump parent path exists but is not a directory: {}'.format(taxdump_dir)
            )
        os.makedirs(taxdump_dir, exist_ok=True)
    _assert_regular_file_or_absent(taxdump_path, label='ETE4 taxdump path')
    if os.path.exists(taxdump_path):
        return taxdump_path
    tmp_path = taxdump_path + '.tmp'
    if os.path.lexists(tmp_path):
        _assert_regular_file_or_absent(tmp_path, label='ETE4 taxdump temporary path')
        os.remove(tmp_path)
    print('Downloading ETE4 taxonomy dump: {}'.format(NCBI_TAXDUMP_URL), flush=True)
    try:
        urlretrieve_fn(NCBI_TAXDUMP_URL, tmp_path)
        os.replace(tmp_path, taxdump_path)
    except Exception:
        if os.path.lexists(tmp_path):
            _assert_regular_file_or_absent(tmp_path, label='ETE4 taxdump temporary path')
            os.remove(tmp_path)
        raise
    return taxdump_path


def should_refresh_custom_ete_taxonomy_db(dbfile, is_taxadb_up_to_date_fn=None):
    dbfile = os.path.realpath(dbfile)
    _assert_regular_file_or_absent(dbfile, label='ETE4 taxonomy DB path')
    if not os.path.exists(dbfile):
        return True
    if is_taxadb_up_to_date_fn is None:
        is_taxadb_up_to_date_fn = getattr(ete4, 'is_taxadb_up_to_date', None)
    if is_taxadb_up_to_date_fn is None:
        return False
    try:
        return (not is_taxadb_up_to_date_fn(dbfile))
    except Exception:
        return True


def get_ete_ncbitaxa(
    args=None,
    acquire_exclusive_lock_fn=None,
    ncbitaxa_cls=None,
    resolve_ete_data_dir_fn=None,
    urlretrieve_fn=None,
    is_taxadb_up_to_date_fn=None,
):
    if acquire_exclusive_lock_fn is None:
        acquire_exclusive_lock_fn = acquire_exclusive_lock
    if ncbitaxa_cls is None:
        ncbitaxa_cls = ete4.NCBITaxa
    if resolve_ete_data_dir_fn is None:
        resolve_ete_data_dir_fn = resolve_ete_data_dir
    if urlretrieve_fn is None:
        urlretrieve_fn = urllib.request.urlretrieve
    if is_taxadb_up_to_date_fn is None:
        is_taxadb_up_to_date_fn = getattr(ete4, 'is_taxadb_up_to_date', None)
    if args is None:
        return ncbitaxa_cls()
    ete_data_dir = resolve_ete_data_dir_fn(args)
    os.makedirs(ete_data_dir, exist_ok=True)
    lock_path = os.path.join(ete_data_dir, '.ete4_taxonomy.lock')
    dbfile = os.path.join(ete_data_dir, 'taxa.sqlite')
    taxdump_file = os.path.join(ete_data_dir, 'taxdump.tar.gz')
    with acquire_exclusive_lock_fn(lock_path=lock_path, lock_label='ETE4 taxonomy DB'):
        kwargs = {'dbfile': dbfile}
        if should_refresh_custom_ete_taxonomy_db(
            dbfile=dbfile,
            is_taxadb_up_to_date_fn=is_taxadb_up_to_date_fn,
        ):
            kwargs['taxdump_file'] = ensure_ete_taxdump_file(
                taxdump_path=taxdump_file,
                urlretrieve_fn=urlretrieve_fn,
            )
        else:
            kwargs['update'] = False
        return ncbitaxa_cls(**kwargs)
