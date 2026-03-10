import errno
import json
import os
import socket
import threading
import time
import urllib.request
from contextlib import contextmanager

import ete4


DOWNLOAD_LOCK_POLL_SECONDS = 5
DOWNLOAD_LOCK_TIMEOUT_SECONDS = 3600
DOWNLOAD_LOCK_HEARTBEAT_SECONDS = 60
DOWNLOAD_LOCK_STALE_SECONDS = 900
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


def _resolve_local_boot_id():
    boot_id_path = '/proc/sys/kernel/random/boot_id'
    if not os.path.isfile(boot_id_path):
        return None
    try:
        with open(boot_id_path) as handle:
            boot_id = handle.readline().strip()
    except OSError:
        return None
    if boot_id == '':
        return None
    return boot_id


def _build_lock_metadata():
    return {
        'format': 'amalgkit-lock-v2',
        'pid': os.getpid(),
        'hostname': socket.gethostname(),
        'boot_id': _resolve_local_boot_id(),
        'created_at': time.time(),
    }


def _serialize_lock_metadata(metadata):
    return json.dumps(metadata, sort_keys=True) + '\n'


def _read_lock_metadata(lock_path):
    try:
        with open(lock_path) as lock_handle:
            raw = lock_handle.read().strip()
    except OSError:
        return None
    if raw == '':
        return None
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        parsed = None
    if isinstance(parsed, dict):
        return parsed
    first_line = raw.splitlines()[0].strip()
    try:
        pid = int(first_line)
    except ValueError:
        return None
    if pid <= 0:
        return None
    return {
        'format': 'legacy-pid',
        'pid': pid,
    }


def _describe_lock_owner(metadata):
    if not isinstance(metadata, dict):
        return 'owner=unknown'
    parts = []
    hostname = metadata.get('hostname')
    pid = metadata.get('pid')
    created_at = metadata.get('created_at')
    if hostname:
        parts.append('host={}'.format(hostname))
    if pid:
        parts.append('pid={}'.format(pid))
    if isinstance(created_at, (float, int)):
        parts.append('created_at={:.0f}'.format(created_at))
    if len(parts) == 0:
        return 'owner=unknown'
    return ', '.join(parts)


def _try_create_lock_file(lock_path):
    _assert_lock_path_is_regular_file(lock_path)
    try:
        fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, 'w') as lock_handle:
        lock_handle.write(_serialize_lock_metadata(_build_lock_metadata()))
    return True


def _read_lock_owner_pid(lock_path):
    metadata = _read_lock_metadata(lock_path)
    if not isinstance(metadata, dict):
        return None
    pid = metadata.get('pid')
    if pid is None:
        return None
    try:
        pid = int(pid)
    except (TypeError, ValueError):
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


def _is_same_host_lock_owner(metadata):
    if not isinstance(metadata, dict):
        return False
    owner_host = metadata.get('hostname')
    local_host = socket.gethostname()
    if (owner_host is None) or (owner_host != local_host):
        return False
    owner_boot_id = metadata.get('boot_id')
    local_boot_id = _resolve_local_boot_id()
    if owner_boot_id and local_boot_id:
        return owner_boot_id == local_boot_id
    return False


def _stale_lock_heartbeat_expired(stat_result, stale_seconds):
    if stale_seconds <= 0:
        return False
    heartbeat_age = time.time() - stat_result.st_mtime
    return heartbeat_age > stale_seconds


def _start_lock_heartbeat(lock_path, interval_seconds):
    interval_seconds = float(interval_seconds)
    if interval_seconds <= 0:
        raise ValueError('interval_seconds must be > 0.')
    stop_event = threading.Event()

    def heartbeat():
        while True:
            if stop_event.wait(interval_seconds):
                return
            try:
                os.utime(lock_path, None)
            except FileNotFoundError:
                return
            except OSError:
                return

    thread = threading.Thread(
        target=heartbeat,
        name='amalgkit-lock-heartbeat',
        daemon=True,
    )
    thread.start()
    return stop_event, thread


def _break_stale_lock_if_needed(lock_path, lock_label='Lock', stale_seconds=DOWNLOAD_LOCK_STALE_SECONDS):
    if not os.path.lexists(lock_path):
        return False
    _assert_lock_path_is_regular_file(lock_path, lock_label=lock_label)
    try:
        stat_before = os.stat(lock_path)
    except FileNotFoundError:
        return False
    metadata = _read_lock_metadata(lock_path)
    owner_pid = _read_lock_owner_pid(lock_path)
    stale_reason = None
    if _is_same_host_lock_owner(metadata) and (owner_pid is not None):
        if not _is_process_alive(owner_pid):
            stale_reason = 'same-host owner PID {} is not running'.format(owner_pid)
    elif _stale_lock_heartbeat_expired(stat_before, stale_seconds):
        stale_reason = 'heartbeat expired after {:.0f} sec ({})'.format(
            time.time() - stat_before.st_mtime,
            _describe_lock_owner(metadata),
        )
    if stale_reason is None:
        return False
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
            heartbeat_stop, heartbeat_thread = _start_lock_heartbeat(
                lock_path=lock_path,
                interval_seconds=DOWNLOAD_LOCK_HEARTBEAT_SECONDS,
            )
            try:
                yield
            finally:
                heartbeat_stop.set()
                heartbeat_thread.join(timeout=max(1.0, float(DOWNLOAD_LOCK_HEARTBEAT_SECONDS)))
                if os.path.lexists(lock_path):
                    _assert_lock_path_is_regular_file(lock_path, lock_label=lock_label)
                    os.remove(lock_path)
            return
        _assert_lock_path_is_regular_file(lock_path, lock_label=lock_label)
        if _break_stale_lock_if_needed(
            lock_path=lock_path,
            lock_label=lock_label,
            stale_seconds=DOWNLOAD_LOCK_STALE_SECONDS,
        ):
            continue
        elapsed = time.time() - wait_start
        if not has_reported_wait:
            owner_metadata = _read_lock_metadata(lock_path)
            print(
                'Another process holds {}. Waiting for lock release: {} ({})'.format(
                    lock_label,
                    lock_path,
                    _describe_lock_owner(owner_metadata),
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
