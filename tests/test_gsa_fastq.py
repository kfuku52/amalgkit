import bz2
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pandas
import pytest

from amalgkit import gsa_fastq
from amalgkit.metadata_utils import Metadata
from tests.support.gsa import manifest_row, payloads_for_row, fastq_bytes


def downloader(payloads, calls):
    def download(**kwargs):
        name = kwargs['source_candidates'][0]['url'].rsplit('/', 1)[-1]
        calls.append(name)
        Path(kwargs['output_path']).write_bytes(payloads[name])
        return True
    return download


@pytest.mark.parametrize('paired', [True, False])
@pytest.mark.parametrize('groups', [1, 2])
def test_full_counts_and_nonoverlapping_ranges_reuse_native_cache(tmp_path, paired, groups):
    row = manifest_row(paired, groups)
    payloads = payloads_for_row(row)
    calls = []
    args = SimpleNamespace(out_dir=str(tmp_path))
    stats = gsa_fastq.prepare_run(args, row, downloader(payloads, calls))
    assert stats['total_spots'] == 4 * groups
    assert stats['total_bases'] == 4 * groups * (20 if paired else 10)
    assert stats['read_count_status'] == 'measured'
    assert gsa_fastq.prepare_run(args, row, downloader(payloads, calls)) == stats
    assert len(calls) == groups * (2 if paired else 1)
    row.update(stats)
    output = tmp_path / 'output'
    output.mkdir()
    assert gsa_fastq.extract_run(args, row, str(output), 2, 3) == (2, 40 if paired else 20)
    first = gzip.decompress(next(output.glob('*1.fastq.gz' if paired else '*.fastq.gz')).read_bytes())
    assert b'@read1' in first and b'@read2' in first and b'@read0' not in first
    assert gsa_fastq.extract_run(args, row, str(output), 4, 4) == (1, 20 if paired else 10)
    second = gzip.decompress(next(output.glob('*1.fastq.gz' if paired else '*.fastq.gz')).read_bytes())
    assert b'@read3' in second and b'@read2' not in second
    assert len(calls) == groups * (2 if paired else 1)


def test_range_crossing_lane_boundary_keeps_pair_order(tmp_path):
    row = manifest_row(groups=2)
    args = SimpleNamespace(out_dir=str(tmp_path))
    row.update(gsa_fastq.prepare_run(args, row, downloader(payloads_for_row(row), [])))
    output = tmp_path / 'output'
    output.mkdir()
    assert gsa_fastq.extract_run(args, row, str(output), 4, 5) == (2, 40)
    for mate in (1, 2):
        data = gzip.decompress((output / f'CRR0001_{mate}.fastq.gz').read_bytes())
        assert [line for line in data.splitlines() if line.startswith(b'@')] == [f'@read3/{mate}'.encode(), f'@read4/{mate}'.encode()]


@pytest.mark.parametrize('compression', ['gz', 'bz2', 'plain'])
def test_gzip_bzip2_and_plain_fastq_inputs(tmp_path, compression):
    row = manifest_row(False)
    files = json.loads(row['gsa_fastq_files'])
    entry = files[0]
    name = 'CRR0001.fq' + ('.' + compression if compression != 'plain' else '')
    entry['filename'] = name
    entry['sources'][0]['url'] = entry['sources'][0]['url'].rsplit('/', 1)[0] + '/' + name
    row['gsa_fastq_files'] = json.dumps(files)
    raw = fastq_bytes()
    payload = gzip.compress(raw) if compression == 'gz' else bz2.compress(raw) if compression == 'bz2' else raw
    stats = gsa_fastq.prepare_run(SimpleNamespace(out_dir=str(tmp_path)), row, downloader({name: payload}, []))
    assert stats['total_spots'] == 4
    assert stats['total_bases'] == 40


def test_provider_checksum_and_corruption_retry(tmp_path):
    row = manifest_row(False)
    files = json.loads(row['gsa_fastq_files'])
    payloads = payloads_for_row(row)
    name = files[0]['filename']
    files[0]['expected_md5'] = hashlib.md5(payloads[name]).hexdigest()
    files[0]['expected_bytes'] = len(payloads[name])
    row['gsa_fastq_files'] = json.dumps(files)
    calls = []
    def transfer(**kwargs):
        calls.append(kwargs)
        Path(kwargs['output_path']).write_bytes(b'bad bytes' if len(calls) == 1 else payloads[name])
        return True
    stats = gsa_fastq.prepare_run(SimpleNamespace(out_dir=str(tmp_path)), row, transfer)
    assert stats['total_spots'] == 4
    assert len(calls) == 2


def test_truncated_compression_is_not_published_as_valid_cache(tmp_path):
    row = manifest_row(False)
    payloads = payloads_for_row(row)
    payloads = {name: data[:-8] for name, data in payloads.items()}
    args = SimpleNamespace(out_dir=str(tmp_path))
    calls = []
    with pytest.raises(gsa_fastq.GsaCorruptInputError, match='compression'):
        gsa_fastq.prepare_run(args, row, downloader(payloads, calls))
    assert len(calls) == 2
    assert not (Path(gsa_fastq.cache_directory(args, row)) / 'validated.json').exists()


@pytest.mark.parametrize('problem', ['ids', 'count'])
def test_mismatched_pairs_fail_before_any_processed_fastq(tmp_path, problem):
    row = manifest_row()
    payloads = payloads_for_row(row)
    second = json.loads(row['gsa_fastq_files'])[1]['filename']
    payloads[second] = gzip.compress(fastq_bytes(3 if problem == 'count' else 4, mate=2, offset=1 if problem == 'ids' else 0))
    args = SimpleNamespace(out_dir=str(tmp_path))
    with pytest.raises(ValueError, match='counts differ|IDs/order differ'):
        gsa_fastq.prepare_run(args, row, downloader(payloads, []))
    assert not (Path(gsa_fastq.cache_directory(args, row)) / 'validated.json').exists()
    assert not (tmp_path / 'getfastq').exists()


def test_interrupted_download_preserves_part_and_retries_on_next_invocation(tmp_path):
    row = manifest_row(False)
    args = SimpleNamespace(out_dir=str(tmp_path))
    payload = next(iter(payloads_for_row(row).values()))
    def interrupt(**kwargs):
        Path(kwargs['output_path']).write_bytes(payload[:15])
        return False
    with pytest.raises(OSError, match='download failed'):
        gsa_fastq.prepare_run(args, row, interrupt)
    def resume(**kwargs):
        assert kwargs['resume_existing'] is True
        assert Path(kwargs['output_path']).read_bytes() == payload[:15]
        with open(kwargs['output_path'], 'ab') as handle:
            handle.write(payload[15:])
        return True
    assert gsa_fastq.prepare_run(args, row, resume)['total_spots'] == 4


def test_changed_manifest_is_a_separate_input_cache(tmp_path):
    row = manifest_row(False)
    args = SimpleNamespace(out_dir=str(tmp_path))
    first = gsa_fastq.cache_directory(args, row)
    files = json.loads(row['gsa_fastq_files'])
    files[0]['expected_bytes'] = 12
    row['gsa_fastq_files'] = json.dumps(files)
    assert gsa_fastq.cache_directory(args, row) != first


def test_changed_cached_bytes_trigger_validation_and_recovery(tmp_path):
    row = manifest_row(False)
    payloads = payloads_for_row(row)
    args = SimpleNamespace(out_dir=str(tmp_path))
    calls = []
    gsa_fastq.prepare_run(args, row, downloader(payloads, calls))
    path = Path(gsa_fastq.cache_directory(args, row)) / next(iter(payloads))
    path.write_bytes(b'bad cache')
    assert gsa_fastq.prepare_run(args, row, downloader(payloads, calls))['total_spots'] == 4
    assert len(calls) == 2


@pytest.mark.parametrize('mutate', [
    lambda files: files[0].update(filename='../escape.fastq.gz'),
    lambda files: files[0]['sources'][0].update(url='file:///etc/passwd'),
    lambda files: files[0]['sources'][0].update(url='https://other.example/CRR0001/input.fastq.gz'),
    lambda files: files[0].update(mate=2),
    lambda files: files[0].update(group=10),
    lambda files: files[0].update(expected_md5='bad'),
])
def test_invalid_manifest_rejected_before_download(mutate):
    row = manifest_row()
    files = json.loads(row['gsa_fastq_files'])
    mutate(files)
    row['gsa_fastq_files'] = json.dumps(files)
    with pytest.raises((ValueError, TypeError)):
        gsa_fastq.read_manifest(row)


def test_symlink_payload_never_overwrites_target(tmp_path):
    row = manifest_row(False)
    args = SimpleNamespace(out_dir=str(tmp_path))
    directory = Path(gsa_fastq.cache_directory(args, row))
    target = tmp_path / 'user-file'
    target.write_text('preserve')
    (directory / json.loads(row['gsa_fastq_files'])[0]['filename']).symlink_to(target)
    with pytest.raises(ValueError, match='regular file'):
        gsa_fastq.prepare_run(args, row, downloader(payloads_for_row(row), []))
    assert target.read_text() == 'preserve'


def test_measured_counts_replace_blanks_without_private_flag(tmp_path):
    row = manifest_row()
    metadata = Metadata.from_DataFrame(pandas.DataFrame([row]))
    result = gsa_fastq.prepare_gsa_metadata(SimpleNamespace(out_dir=str(tmp_path)), metadata, downloader(payloads_for_row(row), []))
    assert result.df.loc[0, 'total_spots'] == 4
    assert result.df.loc[0, 'total_bases'] == 80
    assert result.df.loc[0, 'spot_length'] == 20
    assert result.df.loc[0, 'private_file'] == 'no'


def test_explicit_header_mates_must_match_manifest(tmp_path):
    row = manifest_row()
    payloads = payloads_for_row(row)
    second = json.loads(row['gsa_fastq_files'])[1]['filename']
    payloads[second] = gzip.compress(fastq_bytes(mate=1))
    with pytest.raises(ValueError, match='contradicts its mate'):
        gsa_fastq.prepare_run(SimpleNamespace(out_dir=str(tmp_path)), row, downloader(payloads, []))


def test_same_size_read_counts_do_not_hide_changed_input_contents(tmp_path):
    row = manifest_row(False)
    payloads = payloads_for_row(row)
    args = SimpleNamespace(out_dir=str(tmp_path))
    initial = gsa_fastq.prepare_run(args, row, downloader(payloads, []))
    path = Path(gsa_fastq.cache_directory(args, row)) / next(iter(payloads))
    path.write_bytes(gzip.compress(fastq_bytes().replace(b'AAAAAAAAAA', b'CCCCCCCCCC'), mtime=0))
    updated = gsa_fastq.prepare_run(args, row, downloader(payloads, []))
    assert initial['total_spots'] == updated['total_spots']
    assert initial['total_bases'] == updated['total_bases']
    assert initial['gsa_input_fingerprint'] != updated['gsa_input_fingerprint']
