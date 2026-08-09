import os
import json
import pathlib
import builtins
from contextlib import contextmanager

import pytest
import pandas

from amalgkit.output_utils import atomic_output_path
from amalgkit.quant import (
    _is_ready_index_file,
    _safely_remove_quant_fastq_files,
    _write_index_ready_marker,
    quant_output_exists,
    purge_existing_quant_outputs,
    parse_quant_option_args,
    adapt_oarfish_outputs,
)
from amalgkit.util import Metadata


def write_valid_quant_outputs(output_dir, sra_id='SRR001', target_id='tx1'):
    os.makedirs(output_dir, exist_ok=True)
    pandas.DataFrame([{
        'target_id': target_id,
        'length': 100,
        'eff_length': 90,
        'est_counts': 5,
        'tpm': 1.0,
    }]).to_csv(os.path.join(output_dir, sra_id + '_abundance.tsv'), sep='\t', index=False)
    with open(os.path.join(output_dir, sra_id + '_run_info.json'), 'w', encoding='utf-8') as handle:
        json.dump({'p_pseudoaligned': 50.0}, handle)


def write_safely_removed_fastq_marker(out_dir, sra_id='SRR001'):
    run_dir = pathlib.Path(out_dir) / 'getfastq' / sra_id
    run_dir.mkdir(parents=True, exist_ok=True)
    marker = run_dir / (sra_id + '.fastq.gz.safely_removed')
    marker.write_text(
        'This fastq file was safely removed after `amalgkit quant`.',
        encoding='utf-8',
    )
    return marker


def make_single_run_quant_metadata():
    return Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['SRR001'],
        'scientific_name': ['Species A'],
        'lib_layout': ['single'],
        'total_spots': [10],
        'spot_length': [100],
        'total_bases': [1000],
        'nominal_length': [200],
        'exclusion': ['no'],
    }))


# ---------------------------------------------------------------------------
# quant_output_exists (checks for abundance.tsv in output directory)
# ---------------------------------------------------------------------------


class TestQuantOutputExists:
    def test_output_exists(self, tmp_path):
        """Returns True when both abundance.tsv and run_info.json exist."""
        write_valid_quant_outputs(tmp_path)
        assert quant_output_exists('SRR001', str(tmp_path)) is True

    def test_output_missing(self, tmp_path):
        """Returns False when required output files are missing."""
        assert quant_output_exists('SRR001', str(tmp_path)) is False

    def test_output_missing_run_info(self, tmp_path):
        (tmp_path / 'SRR001_abundance.tsv').write_text('target_id\tdata\n')
        assert quant_output_exists('SRR001', str(tmp_path)) is False

    def test_ignores_directory_named_as_abundance_file(self, tmp_path):
        (tmp_path / 'SRR001_abundance.tsv').mkdir()
        assert quant_output_exists('SRR001', str(tmp_path)) is False

    def test_rejects_header_only_abundance_table(self, tmp_path):
        (tmp_path / 'SRR001_abundance.tsv').write_text(
            'target_id\tlength\teff_length\test_counts\ttpm\n',
            encoding='utf-8',
        )
        (tmp_path / 'SRR001_run_info.json').write_text(
            '{"p_pseudoaligned": 50}',
            encoding='utf-8',
        )

        assert quant_output_exists('SRR001', str(tmp_path)) is False

    def test_rejects_duplicate_target_ids(self, tmp_path):
        write_valid_quant_outputs(tmp_path)
        abundance_path = tmp_path / 'SRR001_abundance.tsv'
        abundance = pandas.read_csv(abundance_path, sep='\t')
        pandas.concat([abundance, abundance], ignore_index=True).to_csv(
            abundance_path,
            sep='\t',
            index=False,
        )

        assert quant_output_exists('SRR001', str(tmp_path)) is False

    @pytest.mark.parametrize('pseudoaligned', [-1, 101, 'NaN'])
    def test_rejects_invalid_pseudoalignment_percentage(self, tmp_path, pseudoaligned):
        write_valid_quant_outputs(tmp_path)
        (tmp_path / 'SRR001_run_info.json').write_text(
            json.dumps({'p_pseudoaligned': pseudoaligned}),
            encoding='utf-8',
        )

        assert quant_output_exists('SRR001', str(tmp_path)) is False


@pytest.mark.parametrize(
    ('option_name', 'option_string'),
    [
        ('--kallisto_options', '--index=other.idx'),
        ('--kallisto_options', '-o other'),
        ('--kallisto_options', '--threads 99'),
        ('--kallisto_options', '-iother.idx'),
        ('--kallisto_options', '-oelsewhere'),
        ('--kallisto_options', '-t8'),
        ('--kallisto_options', '-l200'),
        ('--kallisto_options', '-s20'),
        ('--oarfish_options', '--reads=other.fastq'),
        ('--oarfish_options', '--seq-tech pacbio'),
        ('--oarfish_options', '--output elsewhere'),
        ('--oarfish_options', '-j8'),
        ('--oarfish_options', '-oelsewhere'),
    ],
)
def test_parse_quant_options_rejects_managed_options(option_name, option_string):
    with pytest.raises(ValueError, match='must not override'):
        parse_quant_option_args(option_string, option_name)


class TestPurgeExistingQuantOutputs:
    def test_removes_known_stale_outputs(self, tmp_path):
        for name in [
            'SRR001_abundance.tsv',
            'SRR001_run_info.json',
            'SRR001_abundance.h5',
            'abundance.tsv',
            'run_info.json',
            'abundance.h5',
        ]:
            (tmp_path / name).write_text('x')

        purge_existing_quant_outputs('SRR001', str(tmp_path))

        for name in [
            'SRR001_abundance.tsv',
            'SRR001_run_info.json',
            'SRR001_abundance.h5',
            'abundance.tsv',
            'run_info.json',
            'abundance.h5',
        ]:
            assert not (tmp_path / name).exists()

    def test_raises_when_stale_output_is_directory(self, tmp_path):
        (tmp_path / 'SRR001_abundance.tsv').mkdir()
        with pytest.raises(IsADirectoryError, match='not a file'):
            purge_existing_quant_outputs('SRR001', str(tmp_path))


def test_safe_fastq_cleanup_rolls_back_when_marker_write_fails(tmp_path, monkeypatch):
    fastq_paths = [tmp_path / 'R1.fastq.gz', tmp_path / 'R2.fastq.gz']
    for path in fastq_paths:
        path.write_bytes(b'reads')
    real_atomic_output_path = atomic_output_path
    calls = {'count': 0}

    @contextmanager
    def fail_second_marker(*args, **kwargs):
        calls['count'] += 1
        if calls['count'] == 2:
            raise OSError('marker failure')
        with real_atomic_output_path(*args, **kwargs) as temporary_path:
            yield temporary_path

    monkeypatch.setattr('amalgkit.quant.atomic_output_path', fail_second_marker)

    with pytest.raises(OSError, match='marker failure'):
        _safely_remove_quant_fastq_files([str(path) for path in fastq_paths])

    assert all(path.read_bytes() == b'reads' for path in fastq_paths)
    assert all(not (tmp_path / (path.name + '.safely_removed')).exists() for path in fastq_paths)


def test_index_ready_marker_invalidates_when_reference_changes(tmp_path):
    index_path = tmp_path / 'Species_A.idx'
    fasta_path = tmp_path / 'Species_A.fa'
    index_path.write_text('index', encoding='utf-8')
    fasta_path.write_text('>a\nAAAA\n', encoding='utf-8')
    _write_index_ready_marker(
        str(index_path),
        backend='kallisto',
        fasta_path=str(fasta_path),
    )

    assert _is_ready_index_file(
        str(index_path),
        backend='kallisto',
        require_ready_marker=True,
        fasta_path=str(fasta_path),
    )

    fasta_path.write_text('>a\nCCCC\n', encoding='utf-8')

    assert not _is_ready_index_file(
        str(index_path),
        backend='kallisto',
        require_ready_marker=True,
        fasta_path=str(fasta_path),
    )


@pytest.mark.parametrize(
    ('quant_rows', 'total_spot', 'error_match'),
    [
        ('tx1\t100\tNaN\n', 10, 'non-finite'),
        ('tx1\t100\t-1\n', 10, 'negative'),
        ('tx1\t100\t6\ntx1\t100\t1\n', 10, 'duplicate target_id'),
        ('tx1\t100\t11\n', 10, 'exceeds total_spot'),
    ],
)
def test_adapt_oarfish_outputs_rejects_invalid_raw_values(
    tmp_path,
    quant_rows,
    total_spot,
    error_match,
):
    output_prefix = tmp_path / 'SRR001'
    output_prefix.with_suffix('.quant').write_text(
        'tname\tlen\tnum_reads\n' + quant_rows,
        encoding='utf-8',
    )
    output_prefix.with_suffix('.meta_info.json').write_text('{}', encoding='utf-8')

    with pytest.raises(ValueError, match=error_match):
        adapt_oarfish_outputs(
            output_dir=str(tmp_path),
            sra_id='SRR001',
            sra_stat={'total_spot': total_spot},
            output_prefix=str(output_prefix),
            seq_tech='ont-cdna',
        )


def test_adapt_oarfish_outputs_uses_utf8_for_json_io(tmp_path, monkeypatch):
    output_prefix = tmp_path / 'SRR001'
    output_prefix.with_suffix('.quant').write_text(
        'tname\tlen\tnum_reads\ntranscript葉\t100\t1\n',
        encoding='utf-8',
    )
    meta_info_path = output_prefix.with_suffix('.meta_info.json')
    meta_info_path.write_text(
        json.dumps({'note': '葉'}, ensure_ascii=False),
        encoding='utf-8',
    )
    run_info_path = tmp_path / 'SRR001_run_info.json'
    original_open = builtins.open
    observed = []

    def non_utf8_locale_open(file, *args, **kwargs):
        resolved = os.path.realpath(os.fspath(file))
        if resolved in {os.path.realpath(meta_info_path), os.path.realpath(run_info_path)}:
            observed.append((resolved, kwargs.get('encoding')))
            if 'encoding' not in kwargs:
                kwargs['encoding'] = 'ascii'
        return original_open(file, *args, **kwargs)

    monkeypatch.setattr(builtins, 'open', non_utf8_locale_open)

    adapt_oarfish_outputs(
        output_dir=str(tmp_path),
        sra_id='SRR001',
        sra_stat={'total_spot': 2},
        output_prefix=str(output_prefix),
        seq_tech='ont-cdna',
    )

    assert json.loads(run_info_path.read_text(encoding='utf-8'))['note'] == '葉'
    assert [encoding for _path, encoding in observed] == ['utf-8', 'utf-8']
