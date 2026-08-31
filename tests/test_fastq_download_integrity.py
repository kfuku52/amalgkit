import gzip
import hashlib
from pathlib import Path
from types import SimpleNamespace

import pytest

from amalgkit import getfastq
from amalgkit.fastq_download_integrity import FastqDownloadIntegrityError, validate_download


def entry(payload, filename='SRR1.fastq.gz'):
    return {
        'filename': filename,
        'sources': [{'source_name': 'ENA', 'url': 'https://ftp.sra.ebi.ac.uk/' + filename}],
        'expected_md5': hashlib.md5(payload, usedforsecurity=False).hexdigest(),
        'expected_bytes': len(payload),
    }


def test_crc_and_late_truncation_are_rejected_without_published_checksum(tmp_path):
    good = gzip.compress(b'@read\nACGT\n+\n!!!!\n')
    corrupt = bytearray(good)
    corrupt[-8] ^= 1
    target = tmp_path / 'download'
    for payload in (bytes(corrupt), good[:-3], good + b'not-gzip', gzip.compress(b'')):
        target.write_bytes(payload)
        with pytest.raises(FastqDownloadIntegrityError):
            validate_download(target, {})
    target.write_bytes(good + good)
    validate_download(target, {})  # All concatenated gzip members are valid.


def test_md5_and_size_are_checked_even_for_valid_gzip(tmp_path):
    good = gzip.compress(b'@read\nACGT\n+\n!!!!\n', mtime=0)
    target = tmp_path / 'download'
    target.write_bytes(good)
    validate_download(target, entry(good))
    with pytest.raises(FastqDownloadIntegrityError, match='MD5'):
        validate_download(target, {**entry(good), 'expected_md5': '0' * 32})
    with pytest.raises(FastqDownloadIntegrityError, match='byte count'):
        validate_download(target, {**entry(good), 'expected_bytes': len(good) + 1})


@pytest.mark.parametrize('recover', [True, False])
def test_corrupt_download_is_retried_once_before_publication(tmp_path, monkeypatch, recover):
    good = gzip.compress(b'@read\nACGT\n+\n!!!!\n', mtime=0)
    calls = []

    def download(**kwargs):
        path = Path(kwargs['output_path'])
        assert path.parent != tmp_path
        assert not (tmp_path / 'SRR1.fastq.gz').exists()
        calls.append(path)
        path.write_bytes(good if recover and len(calls) == 2 else b'broken transfer')
        return True

    monkeypatch.setattr(getfastq, 'download_file_from_candidate_sources', download)
    result = getfastq.download_public_original_fastq_files(
        {'sra_id': 'SRR1', 'getfastq_sra_dir': str(tmp_path), 'total_spot': 1},
        SimpleNamespace(), 1, 1, public_original_fastqs=[entry(good)], source_description='ENA',
    )
    assert len(calls) == 2
    if recover:
        assert result is not None
        assert (tmp_path / 'SRR1.fastq.gz').read_bytes() == good
    else:
        assert result is None
        assert not (tmp_path / 'SRR1.fastq.gz').exists()
    assert not any(path.exists() for path in calls)


def test_metadata_reaches_transfer_without_guessing_file_identity():
    good = gzip.compress(b'@r\nA\n+\n!\n')
    url = 'https://ftp.sra.ebi.ac.uk/SRR1.fastq.gz'
    identity = {'expected_md5': entry(good)['expected_md5'], 'expected_bytes': len(good)}
    assert getfastq._build_ena_original_fastq_entries([url], {url: identity}) == [entry(good)]
