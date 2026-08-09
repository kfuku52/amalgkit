import json
import os
import socket
import time
from io import BytesIO

import pytest
from types import SimpleNamespace

from amalgkit import download_utils
from amalgkit.util import (
    acquire_exclusive_lock,
)


# ---------------------------------------------------------------------------
# strtobool
# ---------------------------------------------------------------------------


class TestBoundedDownload:
    def test_writes_nonempty_regular_file_with_bounded_open_timeout(self, tmp_path):
        output_path = tmp_path / 'download.bin'
        observed = []

        def fake_urlopen(url, timeout):
            observed.append((url, timeout))
            return BytesIO(b'payload')

        result = download_utils.download_url_to_regular_file(
            url='https://example.test/download.bin',
            output_path=output_path,
            timeout_seconds=7,
            urlopen_fn=fake_urlopen,
            chunk_size=2,
        )

        assert result == str(output_path)
        assert output_path.read_bytes() == b'payload'
        assert observed == [('https://example.test/download.bin', 7.0)]

    def test_rejects_empty_download_and_removes_partial_output(self, tmp_path):
        output_path = tmp_path / 'download.bin'

        with pytest.raises(ValueError, match='Downloaded file is empty'):
            download_utils.download_url_to_regular_file(
                url='https://example.test/download.bin',
                output_path=output_path,
                timeout_seconds=7,
                urlopen_fn=lambda _url, timeout: BytesIO(b''),
            )

        assert not output_path.exists()

    def test_total_deadline_removes_partial_output(self, tmp_path, monkeypatch):
        output_path = tmp_path / 'download.bin'
        monotonic_values = iter([0.0, 2.0])
        monkeypatch.setattr(
            'amalgkit.download_utils.time.monotonic',
            lambda: next(monotonic_values),
        )

        with pytest.raises(TimeoutError, match='Download exceeded'):
            download_utils.download_url_to_regular_file(
                url='https://example.test/download.bin',
                output_path=output_path,
                timeout_seconds=1,
                urlopen_fn=lambda _url, timeout: BytesIO(b'payload'),
            )

        assert not output_path.exists()


class TestDownloadLockRecovery:
    def test_download_lock_path_honors_custom_lock_directory(self, tmp_path):
        args = SimpleNamespace(
            out_dir=str(tmp_path / 'out'),
            download_dir=str(tmp_path / 'downloads'),
            download_lock_dir=str(tmp_path / 'custom-locks'),
        )

        lock_path = download_utils.resolve_download_lock_path(args, 'silva_refs')

        assert lock_path == os.path.join(
            os.path.realpath(args.download_lock_dir),
            'silva_refs.lock',
        )

    def test_replaced_lock_is_not_touched_or_removed_by_old_owner(self, tmp_path):
        lock_path = tmp_path / 'download.lock'
        owner_state = download_utils._try_create_lock_file(str(lock_path))
        heartbeat_stop, heartbeat_thread = download_utils._start_lock_heartbeat(
            lock_path=str(lock_path),
            interval_seconds=0.01,
            owner_state=owner_state,
        )
        lock_path.unlink()
        replacement_metadata = {
            'format': 'amalgkit-lock-v3',
            'hostname': socket.gethostname(),
            'pid': os.getpid(),
            'created_at': time.time(),
            'owner_token': 'replacement-owner',
        }
        lock_path.write_text(json.dumps(replacement_metadata) + '\n')
        fixed_mtime = time.time() - 5
        os.utime(lock_path, (fixed_mtime, fixed_mtime))
        replacement_mtime = os.stat(lock_path).st_mtime_ns
        time.sleep(0.05)

        download_utils._release_heartbeat_lock(
            lock_path=str(lock_path),
            heartbeat_stop=heartbeat_stop,
            heartbeat_thread=heartbeat_thread,
            owner_state=owner_state,
            lock_label='test lock',
        )

        assert lock_path.exists()
        assert os.stat(lock_path).st_mtime_ns == replacement_mtime
        assert json.loads(lock_path.read_text())['owner_token'] == 'replacement-owner'

    def test_acquire_exclusive_lock_rejects_symlink_lock_path(self, tmp_path):
        lock_path = tmp_path / 'download.lock'
        os.symlink(tmp_path / 'dangling.lock.target', lock_path)

        with pytest.raises(IsADirectoryError, match='test lock path exists but is not a file'):
            with acquire_exclusive_lock(str(lock_path), lock_label='test lock', poll_seconds=1, timeout_seconds=2):
                pass

    def test_acquire_exclusive_lock_reclaims_stale_same_host_lock(self, tmp_path, monkeypatch):
        lock_path = tmp_path / 'download.lock'
        lock_path.write_text(json.dumps({
            'format': 'amalgkit-lock-v2',
            'hostname': socket.gethostname(),
            'boot_id': 'boot-1',
            'pid': 999999,
            'created_at': time.time(),
        }) + '\n')

        def fake_kill(_pid, _sig):
            raise ProcessLookupError()

        monkeypatch.setattr('amalgkit.download_utils._resolve_local_boot_id', lambda: 'boot-1')
        monkeypatch.setattr('amalgkit.util.os.kill', fake_kill)

        with acquire_exclusive_lock(str(lock_path), lock_label='test lock', poll_seconds=1, timeout_seconds=2):
            assert lock_path.exists()
            metadata = json.loads(lock_path.read_text())
            assert metadata['pid'] == os.getpid()
            assert metadata['hostname'] == socket.gethostname()

        assert not lock_path.exists()

    def test_reclaims_heartbeat_expired_same_host_lock_even_when_pid_is_alive(
        self,
        tmp_path,
        monkeypatch,
    ):
        lock_path = tmp_path / 'download.lock'
        lock_path.write_text(json.dumps({
            'format': 'amalgkit-lock-v3',
            'hostname': socket.gethostname(),
            'boot_id': 'boot-alive',
            'pid': os.getpid(),
            'created_at': time.time() - 30,
            'owner_token': 'stalled-owner',
        }) + '\n')
        stale_at = time.time() - 10
        os.utime(lock_path, (stale_at, stale_at))
        monkeypatch.setattr(
            'amalgkit.download_utils._resolve_local_boot_id',
            lambda: 'boot-alive',
        )
        monkeypatch.setattr(
            'amalgkit.download_utils._is_process_alive',
            lambda _pid: True,
        )

        removed = download_utils._break_stale_lock_if_needed(
            str(lock_path),
            lock_label='test lock',
            stale_seconds=1,
        )

        assert removed is True
        assert not lock_path.exists()

    def test_stale_breaker_does_not_remove_lock_when_owner_token_changes(
        self,
        tmp_path,
        monkeypatch,
    ):
        lock_path = tmp_path / 'download.lock'
        old_metadata = {
            'format': 'amalgkit-lock-v3',
            'hostname': socket.gethostname(),
            'boot_id': 'boot-race',
            'pid': os.getpid(),
            'created_at': time.time() - 30,
            'owner_token': 'old-owner',
        }
        replacement_metadata = dict(old_metadata, owner_token='replacement-owner')
        lock_path.write_text(json.dumps(old_metadata) + '\n')
        stale_at = time.time() - 10
        os.utime(lock_path, (stale_at, stale_at))
        metadata_reads = iter([
            old_metadata,
            old_metadata,
            replacement_metadata,
        ])
        monkeypatch.setattr(
            'amalgkit.download_utils._read_lock_metadata',
            lambda _path: next(metadata_reads),
        )
        monkeypatch.setattr(
            'amalgkit.download_utils._resolve_local_boot_id',
            lambda: 'boot-race',
        )
        monkeypatch.setattr(
            'amalgkit.download_utils._is_process_alive',
            lambda _pid: True,
        )

        removed = download_utils._break_stale_lock_if_needed(
            str(lock_path),
            lock_label='test lock',
            stale_seconds=1,
        )

        assert removed is False
        assert lock_path.exists()

    def test_acquire_exclusive_lock_times_out_when_same_host_owner_pid_is_alive(self, tmp_path, monkeypatch):
        lock_path = tmp_path / 'download.lock'
        lock_path.write_text(json.dumps({
            'format': 'amalgkit-lock-v2',
            'hostname': socket.gethostname(),
            'boot_id': 'boot-2',
            'pid': 12345,
            'created_at': time.time(),
        }) + '\n')

        monkeypatch.setattr('amalgkit.download_utils._resolve_local_boot_id', lambda: 'boot-2')
        monkeypatch.setattr('amalgkit.util.os.kill', lambda _pid, _sig: None)
        fake_now = [time.time()]
        monkeypatch.setattr('amalgkit.download_utils.time.time', lambda: fake_now[0])
        monkeypatch.setattr(
            'amalgkit.download_utils.time.sleep',
            lambda seconds: fake_now.__setitem__(0, fake_now[0] + seconds),
        )

        with pytest.raises(TimeoutError, match='Timed out'):
            with acquire_exclusive_lock(str(lock_path), lock_label='test lock', poll_seconds=1, timeout_seconds=1):
                pass

        assert lock_path.exists()

    def test_acquire_exclusive_lock_reclaims_stale_cross_host_lock_after_heartbeat_timeout(self, tmp_path, monkeypatch):
        lock_path = tmp_path / 'download.lock'
        lock_path.write_text(json.dumps({
            'format': 'amalgkit-lock-v2',
            'hostname': 'remote-node-01',
            'boot_id': 'remote-boot',
            'pid': 321,
            'created_at': time.time() - 30,
        }) + '\n')
        stale_at = time.time() - 10
        os.utime(lock_path, (stale_at, stale_at))

        monkeypatch.setattr('amalgkit.download_utils.DOWNLOAD_LOCK_STALE_SECONDS', 1)

        with acquire_exclusive_lock(str(lock_path), lock_label='test lock', poll_seconds=1, timeout_seconds=2):
            metadata = json.loads(lock_path.read_text())
            assert metadata['pid'] == os.getpid()

        assert not lock_path.exists()

    def test_acquire_exclusive_lock_waits_for_fresh_cross_host_lock(self, tmp_path, monkeypatch):
        lock_path = tmp_path / 'download.lock'
        lock_path.write_text(json.dumps({
            'format': 'amalgkit-lock-v2',
            'hostname': 'remote-node-02',
            'boot_id': 'remote-boot',
            'pid': 654,
            'created_at': time.time(),
        }) + '\n')

        monkeypatch.setattr('amalgkit.download_utils.DOWNLOAD_LOCK_STALE_SECONDS', 60)
        fake_now = [time.time()]
        monkeypatch.setattr('amalgkit.download_utils.time.time', lambda: fake_now[0])
        monkeypatch.setattr(
            'amalgkit.download_utils.time.sleep',
            lambda seconds: fake_now.__setitem__(0, fake_now[0] + seconds),
        )

        with pytest.raises(TimeoutError, match='Timed out'):
            with acquire_exclusive_lock(str(lock_path), lock_label='test lock', poll_seconds=1, timeout_seconds=1):
                pass

        assert lock_path.exists()

    @pytest.mark.slow
    def test_acquire_exclusive_lock_updates_lock_heartbeat(self, tmp_path, monkeypatch):
        lock_path = tmp_path / 'download.lock'

        monkeypatch.setattr('amalgkit.download_utils.DOWNLOAD_LOCK_HEARTBEAT_SECONDS', 0.05)

        with acquire_exclusive_lock(str(lock_path), lock_label='test lock', poll_seconds=1, timeout_seconds=2):
            first_mtime = os.stat(lock_path).st_mtime_ns
            time.sleep(0.2)
            second_mtime = os.stat(lock_path).st_mtime_ns
            assert second_mtime > first_mtime


class TestDownloadSemaphoreRecovery:
    def test_default_wait_does_not_expire_after_one_hour(self):
        assert download_utils.DOWNLOAD_LOCK_TIMEOUT_SECONDS is None

    def test_acquire_counting_semaphore_wait_false_returns_none_when_all_slots_busy(self, tmp_path):
        from amalgkit.download_utils import acquire_counting_semaphore

        semaphore_dir = tmp_path / 'semaphore'
        semaphore_dir.mkdir()
        slot_path = semaphore_dir / 'slot-0001.lock'
        slot_path.write_text(json.dumps({
            'format': 'amalgkit-lock-v2',
            'hostname': socket.gethostname(),
            'pid': os.getpid(),
            'created_at': time.time(),
        }) + '\n')

        with acquire_counting_semaphore(
            semaphore_dir=str(semaphore_dir),
            max_concurrency=1,
            lock_label='test semaphore',
            poll_seconds=1,
            timeout_seconds=2,
            wait=False,
        ) as acquired_slot_path:
            assert acquired_slot_path is None

        assert slot_path.exists()

    def test_acquire_counting_semaphore_reclaims_stale_same_host_slot(self, tmp_path, monkeypatch):
        from amalgkit.download_utils import acquire_counting_semaphore

        semaphore_dir = tmp_path / 'semaphore'
        semaphore_dir.mkdir()
        slot_path = semaphore_dir / 'slot-0001.lock'
        slot_path.write_text(json.dumps({
            'format': 'amalgkit-lock-v2',
            'hostname': socket.gethostname(),
            'boot_id': 'boot-semaphore',
            'pid': 999999,
            'created_at': time.time(),
        }) + '\n')

        def fake_kill(_pid, _sig):
            raise ProcessLookupError()

        monkeypatch.setattr('amalgkit.download_utils._resolve_local_boot_id', lambda: 'boot-semaphore')
        monkeypatch.setattr('amalgkit.download_utils.os.kill', fake_kill)

        with acquire_counting_semaphore(
            semaphore_dir=str(semaphore_dir),
            max_concurrency=1,
            lock_label='test semaphore',
            poll_seconds=1,
            timeout_seconds=2,
        ) as acquired_slot_path:
            assert acquired_slot_path == str(slot_path)
            metadata = json.loads(slot_path.read_text())
            assert metadata['pid'] == os.getpid()
            assert metadata['hostname'] == socket.gethostname()

        assert not slot_path.exists()

    @pytest.mark.slow
    def test_acquire_counting_semaphore_updates_slot_heartbeat(self, tmp_path, monkeypatch):
        from amalgkit.download_utils import acquire_counting_semaphore

        semaphore_dir = tmp_path / 'semaphore'
        monkeypatch.setattr('amalgkit.download_utils.DOWNLOAD_LOCK_HEARTBEAT_SECONDS', 0.05)

        with acquire_counting_semaphore(
            semaphore_dir=str(semaphore_dir),
            max_concurrency=1,
            lock_label='test semaphore',
            poll_seconds=1,
            timeout_seconds=2,
        ) as acquired_slot_path:
            first_mtime = os.stat(acquired_slot_path).st_mtime_ns
            time.sleep(0.2)
            second_mtime = os.stat(acquired_slot_path).st_mtime_ns
            assert second_mtime > first_mtime
