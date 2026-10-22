"""Streaming hash identities retain mutation detection and cache ownership."""

import hashlib
import os

import pytest

from amalgkit import quant_provenance


@pytest.mark.parametrize('size', [0, 1, 262143, 262144, 262145, 1048575, 1048576, 1048577])
def test_file_identity_matches_independent_digest_at_buffer_boundaries(tmp_path, size):
    path = tmp_path / '参照.bin'
    content = (bytes(range(256)) * (size // 256 + 1))[:size]
    path.write_bytes(content)
    expected = {'size': size, 'sha256': hashlib.sha256(content).hexdigest()}
    assert quant_provenance.file_identity(str(path), {}) == expected
    assert path.read_bytes() == content


def test_identity_cache_shares_aliases_and_returns_independent_dicts(tmp_path, monkeypatch):
    path, alias = tmp_path / 'input.bin', tmp_path / 'alias.bin'
    path.write_bytes(b'original')
    alias.symlink_to(path)
    cache = {}
    original = quant_provenance.file_identity(str(path), cache)
    original['sha256'] = 'caller mutation'

    def do_not_hash(*args, **kwargs):
        raise AssertionError('unchanged source must use cached digest')

    with monkeypatch.context() as patch:
        patch.setattr(quant_provenance.hashlib, 'file_digest', do_not_hash)
        assert quant_provenance.file_identity(str(alias), cache)['sha256'] == hashlib.sha256(b'original').hexdigest()
    assert len(cache) == 1
    before = path.stat()
    path.write_bytes(b'changed!')
    os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert quant_provenance.file_identity(str(path), cache)['sha256'] == hashlib.sha256(b'changed!').hexdigest()


def test_mutation_during_file_digest_does_not_publish_cached_identity(tmp_path, monkeypatch):
    path = tmp_path / 'input.bin'
    path.write_bytes(b'original')
    real_digest = hashlib.file_digest

    def mutate(handle, algorithm):
        result = real_digest(handle, algorithm)
        path.write_bytes(b'mutated!')
        return result

    monkeypatch.setattr(quant_provenance.hashlib, 'file_digest', mutate)
    cache = {}
    with pytest.raises(ValueError, match='Quant input changed while fingerprinting'):
        quant_provenance.file_identity(str(path), cache)
    assert cache == {}
