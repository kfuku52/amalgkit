"""Real-tool comparisons; skipped when kallisto is not installed.

This controlled fixture tests wrapper equivalence and length sensitivity. It
does not establish a generally correct fragment prior for biological libraries.
"""

import json
import shutil
import subprocess
from types import SimpleNamespace

import numpy
import pandas
import pytest

from amalgkit.fragment_length import PROVENANCE_KEY
from amalgkit.metadata_utils import Metadata
from amalgkit.main import build_main_parser
from amalgkit.quant import call_kallisto, quant_main

pytestmark = [pytest.mark.integration, pytest.mark.skipif(shutil.which('kallisto') is None, reason='kallisto not installed')]


def make_fragment_fixture(directory, seed=17, mean=150, sd=7, nreads=20000):
    """Uniform molecule abundances, normal lengths conditioned to fit each transcript.

    Select transcripts in proportion to their expected available start positions;
    within each transcript, sample length with its conditioned normal probability
    and then a uniform start. Read length is always 75 and never used as a prior.
    """
    rng = numpy.random.default_rng(seed)
    lengths = numpy.array([160, 180, 250, 500, 1000, 2000])
    sequences = [''.join(rng.choice(list('ACGT'), size=int(length))) for length in lengths]
    # Two isoforms share a prefix, introducing nontrivial multi-mapping.
    sequences[2] = sequences[1][:120] + sequences[2][120:]
    fasta = directory / 'transcripts.fa'
    fasta.write_text(''.join(f'>t{i}\n{sequence}\n' for i, sequence in enumerate(sequences)))
    distributions = []
    effective = []
    for length in lengths:
        fragment_lengths = numpy.arange(75, length + 1)
        weights = numpy.exp(-0.5 * ((fragment_lengths - mean) / sd) ** 2)
        weights /= weights.sum()
        distributions.append((fragment_lengths, weights))
        effective.append(float(numpy.sum((length - fragment_lengths + 1) * weights)))
    frequencies = numpy.asarray(effective) / numpy.sum(effective)
    transcripts = rng.choice(len(lengths), size=nreads, p=frequencies)
    counts = numpy.bincount(transcripts, minlength=len(lengths))
    fastq = directory / 'reads.fastq'
    with fastq.open('w') as handle:
        for read_id, transcript in enumerate(transcripts):
            fragment_lengths, weights = distributions[transcript]
            fragment_length = int(rng.choice(fragment_lengths, p=weights))
            start = int(rng.integers(0, lengths[transcript] - fragment_length + 1))
            read = sequences[transcript][start:start + 75]
            handle.write(f'@r{read_id}\n{read}\n+\n{"I" * 75}\n')
    return fasta, fastq, {'lengths': lengths.tolist(), 'molecule_tpm': [1e6 / len(lengths)] * len(lengths),
                          'fragment_counts': counts.tolist(), 'seed': seed, 'mean': mean, 'sd': sd}


def quantify_fixture(directory, seed=17, mean=150, sd=7):
    fasta, fastq, truth = make_fragment_fixture(directory, seed=seed, mean=mean, sd=sd)
    index = directory / 'reference.idx'
    subprocess.run(['kallisto', 'index', '-i', str(index), str(fasta)], check=True, capture_output=True, timeout=60)
    results = {}
    for label, used_mean, used_sd in [('specified', mean, sd), ('old', max(200, mean), max(200, mean) / 10),
                                      ('mean_only', mean, 20), ('sd_only', 200, sd)]:
        output = directory / label
        subprocess.run(['kallisto', 'quant', '-i', str(index), '-o', str(output), '--threads', '1',
                        '--single', '-l', str(used_mean), '-s', str(used_sd), '--plaintext', str(fastq)],
                       check=True, capture_output=True, timeout=60)
        results[label] = pandas.read_csv(output / 'abundance.tsv', sep='\t').set_index('target_id')
    wrapped = directory / 'amalgkit'
    wrapped.mkdir()
    data = Metadata.from_DataFrame(pandas.DataFrame([{'run': 'R1', 'lib_layout': 'single',
        'fragment_length_mean': mean, 'fragment_length_sd': sd,
        'fragment_length_source': 'simulated', 'fragment_length_source_detail': f'normal input distribution; seed {seed}'}]))
    call_kallisto(SimpleNamespace(threads=1, kallisto_options='--plaintext'), [str(fastq)], data,
                  {'sra_id': 'R1', 'layout': 'single'}, str(wrapped), str(index))
    results['amalgkit'] = pandas.read_csv(wrapped / 'R1_abundance.tsv', sep='\t').set_index('target_id')
    model = json.loads((wrapped / 'R1_run_info.json').read_text())[PROVENANCE_KEY]
    return results, truth, model


def test_real_kallisto_preserves_mean_sd_and_matches_direct_invocation(tmp_path):
    results, truth, model = quantify_fixture(tmp_path)
    pandas.testing.assert_frame_equal(results['amalgkit'], results['specified'], check_exact=True)
    assert model['mean']['value'] == 150 and model['sd']['value'] == 7
    assert model['kallisto_version'] != 'unknown'
    # These checks establish sensitivity, not universal accuracy improvements.
    assert not numpy.allclose(results['specified']['eff_length'], results['old']['eff_length'])
    assert not numpy.allclose(results['specified']['eff_length'], results['mean_only']['eff_length'])
    assert not numpy.allclose(results['specified']['tpm'], results['old']['tpm'])
    for frame in results.values():
        assert numpy.isfinite(frame.to_numpy()).all()
        assert (frame['eff_length'] > 0).all()
        assert (frame['est_counts'] >= 0).all()
        assert frame['tpm'].sum() == pytest.approx(1e6, rel=1e-5)
    assert sum(truth['fragment_counts']) == 20000


def test_real_quant_batch_uses_shared_fragment_file(tmp_path):
    results, _, _ = quantify_fixture(tmp_path)
    rows = [{'run': run, 'scientific_name': 'Test species', 'is_sampled': 'yes',
             'lib_layout': 'single', 'total_spots': 20000, 'spot_length': 75, 'total_bases': 1500000}
            for run in ['R1', 'R2']]
    metadata_path = tmp_path / 'metadata.tsv'
    pandas.DataFrame(rows).to_csv(metadata_path, sep='\t', index=False)
    fragment_path = tmp_path / 'fragments.tsv'
    fragment_path.write_text('run\tfragment_length_mean\tfragment_length_sd\tsource\tsource_detail\n'
                             'R1\t180\t20\tuser\tunused batch\nR2\t150\t7\tuser\tsimulated R2\n')
    index_dir = tmp_path / 'index'
    index_dir.mkdir()
    shutil.copyfile(tmp_path / 'reference.idx', index_dir / 'Test_species.idx')
    read_dir = tmp_path / 'getfastq' / 'R2'
    read_dir.mkdir(parents=True)
    shutil.copyfile(tmp_path / 'reads.fastq', read_dir / 'R2.fastq')
    runtime = build_main_parser().parse_args(['quant', '--out_dir', str(tmp_path), '--metadata', str(metadata_path),
        '--fragment_length_file', str(fragment_path), '--batch', '2', '--threads', '1', '--quant_backend', 'kallisto',
        '--clean_fastq', 'no', '--kallisto_options=--plaintext'])
    quant_main(runtime)
    actual = pandas.read_csv(tmp_path / 'quant' / 'R2' / 'R2_abundance.tsv', sep='\t').set_index('target_id')
    pandas.testing.assert_frame_equal(actual, results['specified'], check_exact=True)
    assert not (tmp_path / 'quant' / 'R1').exists()
