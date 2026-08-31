"""Validate downloaded public FASTQ bytes before atomic publication."""

import gzip
import os
import re
import zlib

from amalgkit.download_utils import calculate_file_md5


class FastqDownloadIntegrityError(ValueError):
    """Downloaded bytes do not match the provider or complete gzip stream."""


def validate_integrity_metadata(entry):
    md5 = entry.get('expected_md5')
    size = entry.get('expected_bytes')
    if md5 is not None and (not isinstance(md5, str) or re.fullmatch(r'[0-9a-f]{32}', md5) is None):
        raise ValueError('Original FASTQ expected MD5 is invalid.')
    if size is not None and (isinstance(size, bool) or not isinstance(size, int) or size <= 0):
        raise ValueError('Original FASTQ expected byte count is invalid.')


def validate_download(path, entry):
    """Stream through all gzip members/trailers; memory use is bounded."""
    validate_integrity_metadata(entry)
    if entry.get('expected_bytes') is not None and os.path.getsize(path) != entry['expected_bytes']:
        raise FastqDownloadIntegrityError('Original FASTQ byte count differs from ENA metadata.')
    if entry.get('expected_md5') is not None and calculate_file_md5(path) != entry['expected_md5']:
        raise FastqDownloadIntegrityError('Original FASTQ MD5 differs from ENA metadata.')
    try:
        uncompressed_bytes = 0
        with gzip.open(path, 'rb') as handle:
            while block := handle.read(1024 * 1024):
                uncompressed_bytes += len(block)
        if not uncompressed_bytes:
            raise FastqDownloadIntegrityError('Original FASTQ gzip stream is empty.')
    except (OSError, EOFError, zlib.error) as exc:
        raise FastqDownloadIntegrityError('Original FASTQ gzip integrity validation failed.') from exc

