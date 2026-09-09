"""Transport regression coverage from the public GSA smoke test."""
import http.server
import shutil
import socket
import ssl
import subprocess
import threading
from types import SimpleNamespace

import pytest

import amalgkit.getfastq as getfastq


def args():
    return SimpleNamespace(sra_download_transfer_timeout_seconds=30, sra_download_method='curl')


def install_failed_transfer(monkeypatch, output, *, code=18, status=200, detail='connection closed'):
    calls = []
    monkeypatch.setattr(getfastq.shutil, 'which', lambda _: '/usr/bin/curl')
    monkeypatch.setattr(getfastq.time, 'sleep', lambda _: None)

    def run(command, **kwargs):
        start = output.stat().st_size if output.exists() else 0
        calls.append(start)
        header = command[command.index('--dump-header') + 1]
        http_status = 206 if start and status == 200 else status
        headers = 'HTTP/1.1 {} Response\r\n'.format(http_status)
        if http_status == 206:
            headers += 'Content-Range: bytes {}-9/10\r\n'.format(start)
        if http_status in (200, 206):
            with output.open('ab') as f:
                f.write(b'x')
        with open(header, 'w') as f:
            f.write(headers + '\r\n')
        return subprocess.CompletedProcess(command, code,
            stdout='AMALGKIT_CURL_STATUS:{}\nAMALGKIT_CURL_REDIRECT:\n'.format(http_status).encode(),
            stderr=detail.encode())

    monkeypatch.setattr(getfastq.subprocess, 'run', run)
    return calls


def test_transient_retries_are_bounded_and_keep_partial_without_urllib(tmp_path, monkeypatch, capsys):
    output = tmp_path / 'input.fastq.gz.part'
    calls = install_failed_transfer(monkeypatch, output)
    monkeypatch.setattr(getfastq, 'download_with_urllib', lambda **_: pytest.fail('must keep the resumable input'))
    assert not getfastq._download_file_from_source_without_semaphore(
        sra_id='CRR0001', source_url_original='https://download.cncb.ac.cn/gsa/run.fastq.gz', output_path=str(output),
        args=args(), sra_source_name='GSA', artifact_label='GSA original FASTQ', resume_existing=True,
    )
    assert calls == [0, 1, 2, 3]
    assert output.read_bytes() == b'xxxx'
    assert 'exit=18, HTTP=206' in capsys.readouterr().err


def test_permanent_http_failure_is_not_retried_and_diagnostics_redact_urls(tmp_path, monkeypatch, capsys):
    output = tmp_path / 'input.fastq.gz.part'
    calls = install_failed_transfer(monkeypatch, output, code=22, status=404,
        detail='not found: https://user:secret@download.cncb.ac.cn/run?token=private#fragment')
    assert not getfastq.download_with_curl(
        'https://download.cncb.ac.cn/gsa/run.fastq.gz', str(output), args(), 'GSA', resume_existing=True,
    )
    assert calls == [0]
    error = capsys.readouterr().err
    assert 'exit=22, HTTP=404' in error and 'not found' in error
    assert all(secret not in error for secret in ('secret', 'private', 'fragment', 'user:'))


def test_automatic_retry_obeys_original_transfer_deadline(tmp_path, monkeypatch):
    output = tmp_path / 'input.fastq.gz.part'
    calls = install_failed_transfer(monkeypatch, output)
    ticks = iter([0, 0, 30])
    monkeypatch.setattr(getfastq.time, 'monotonic', lambda: next(ticks))
    assert not getfastq.download_with_curl(
        'https://download.cncb.ac.cn/gsa/run.fastq.gz', str(output), args(), 'GSA', resume_existing=True,
    )
    assert calls == [0]


@pytest.mark.integration
def test_real_curl_recovers_truncated_https_response(tmp_path, monkeypatch):
    """An actual TLS socket closes mid-body; the next request must resume."""
    if shutil.which('curl') is None or shutil.which('openssl') is None:
        pytest.skip('curl and openssl are required for the local TLS transport test')
    cert, key = tmp_path/'cert.pem', tmp_path/'key.pem'
    subprocess.run([
        'openssl', 'req', '-x509', '-newkey', 'rsa:2048', '-nodes', '-days', '1',
        '-subj', '/CN=127.0.0.1', '-addext', 'subjectAltName=IP:127.0.0.1',
        '-keyout', str(key), '-out', str(cert),
    ], check=True, capture_output=True, timeout=15)
    body = b'ACGT' * 8192
    offsets = []

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *_):
            pass

        def do_GET(self):
            offset = int(self.headers.get('Range', 'bytes=0-').split('=')[1].split('-')[0])
            offsets.append(offset)
            self.send_response(206 if offset else 200)
            self.send_header('Content-Length', str(len(body)-offset))
            if offset:
                self.send_header('Content-Range', 'bytes {}-{}/{}'.format(offset, len(body)-1, len(body)))
            self.end_headers()
            if len(offsets) == 1:
                self.wfile.write(body[:1024])
                self.wfile.flush()
                self.connection.shutdown(socket.SHUT_RDWR)
                self.close_connection = True
            else:
                self.wfile.write(body[offset:])

    server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
    context.load_cert_chain(cert, key)
    server.socket = context.wrap_socket(server.socket, server_side=True)
    worker = threading.Thread(target=server.serve_forever, daemon=True)
    worker.start()
    url = 'https://127.0.0.1:{}/test.fastq.gz'.format(server.server_port)
    monkeypatch.setattr(getfastq, 'assert_allowed_download_url', lambda value: value == url or pytest.fail('unexpected URL'))
    monkeypatch.setenv('CURL_CA_BUNDLE', str(cert))
    monkeypatch.setenv('NO_PROXY', '127.0.0.1')
    monkeypatch.setattr(getfastq.time, 'sleep', lambda _: None)
    output = tmp_path/'input.fastq.gz.part'
    try:
        assert getfastq.download_with_curl(url, str(output), args(), 'GSA', resume_existing=True)
        assert offsets == [0, 1024]
        assert output.read_bytes() == body
    finally:
        server.shutdown()
        server.server_close()
        worker.join(timeout=5)


def test_interrupted_mismatched_range_cannot_poison_retained_prefix(tmp_path, monkeypatch):
    output = tmp_path/'input.fastq.gz.part'
    output.write_bytes(b'original')
    monkeypatch.setattr(getfastq.shutil, 'which', lambda _: '/usr/bin/curl')
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        with open(command[command.index('--dump-header') + 1], 'w') as f:
            f.write('HTTP/1.1 206 Partial Content\r\nContent-Range: bytes 0-99/100\r\n\r\n')
        with output.open('ab') as f:
            f.write(b'wrong offset')
        return subprocess.CompletedProcess(command, 18,
            stdout=b'AMALGKIT_CURL_STATUS:206\n', stderr=b'transfer closed')

    monkeypatch.setattr(getfastq.subprocess, 'run', run)
    assert not getfastq.download_with_curl(
        'https://download.cncb.ac.cn/gsa/run.fastq.gz', str(output), args(), 'GSA', resume_existing=True,
    )
    assert len(calls) == 1 and output.read_bytes() == b'original'


def test_urllib_failure_logs_exception_type_and_redacted_reason(tmp_path, monkeypatch, capsys):
    import urllib.error
    def fail(**_):
        raise urllib.error.URLError('connection reset: https://name:secret@download.cncb.ac.cn/run?token=private')
    monkeypatch.setattr(getfastq, 'download_with_urllib', fail)
    options = args()
    options.sra_download_method = 'urllib'
    assert not getfastq._download_file_from_source_without_semaphore(
        sra_id='CRR0001', source_url_original='https://download.cncb.ac.cn/gsa/run.fastq.gz',
        output_path=str(tmp_path/'input'), args=options, sra_source_name='GSA', artifact_label='FASTQ',
    )
    error = capsys.readouterr().err
    assert 'URLError' in error and 'connection reset' in error
    assert 'secret' not in error and 'private' not in error


def test_missing_curl_still_allows_urllib_to_replace_partial(tmp_path, monkeypatch):
    output = tmp_path/'input.part'
    output.write_bytes(b'old partial')
    monkeypatch.setattr(getfastq.shutil, 'which', lambda _: None)
    def download(**kwargs):
        with open(kwargs['output_path'], 'wb') as f:
            f.write(b'complete input')
    monkeypatch.setattr(getfastq, 'download_with_urllib', download)
    assert getfastq._download_file_from_source_without_semaphore(
        sra_id='CRR0001', source_url_original='https://download.cncb.ac.cn/gsa/run.fastq.gz',
        output_path=str(output), args=args(), sra_source_name='GSA', artifact_label='FASTQ', resume_existing=True,
    )
    assert output.read_bytes() == b'complete input'
