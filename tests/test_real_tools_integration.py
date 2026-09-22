"""Small, network-free real-tool contracts also required by nightly CI."""

import gzip
import hashlib
import json
import shutil

import numpy
import pandas
import pytest

from amalgkit.fastq_utils import validate_fastq_structure
from amalgkit.getfastq import getfastq_main
from amalgkit.main import build_main_parser
from amalgkit.quant import quant_main

pytestmark = pytest.mark.integration


def write_reads(path, sequences):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt') as handle:
        for index, sequence in enumerate(sequences):
            handle.write(f'@r{index}\n{sequence}\n+\n' + 'I' * len(sequence) + '\n')


@pytest.mark.skipif(any(shutil.which(tool) is None for tool in ('fastp', 'seqkit', 'fasterq-dump')),
                    reason='fastp, seqkit and fasterq-dump are required')
@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('compressed', [False, True])
def test_real_default_fastp_private_sources_and_resume(tmp_path, paired, compressed):
    rng = numpy.random.default_rng(73)
    sequences = [''.join(rng.choice(list('ACGT'), size=100)) for _ in range(24)]
    sources = []
    for mate in range(2 if paired else 1):
        path = tmp_path / f'original{mate}.fastq{ ".gz" if compressed else ""}'
        reads = sequences if mate == 0 else [s.translate(str.maketrans('ACGT', 'TGCA'))[::-1] for s in sequences]
        write_reads(path, reads)
        sources.append(path)
    hashes = [hashlib.sha256(path.read_bytes()).hexdigest() for path in sources]
    metadata = tmp_path / 'metadata.tsv'
    pandas.DataFrame([dict(run='R1', scientific_name='Test species', taxid=9606,
                           lib_layout='paired' if paired else 'single', private_file='yes',
                           read1_path=str(sources[0]), read2_path=str(sources[-1]) if paired else '',
                           total_spots=24, total_bases=2400 * len(sources), spot_length=100 * len(sources),
                           size=sum(p.stat().st_size for p in sources), exclusion='no', is_sampled='yes')]).to_csv(
        metadata, sep='\t', index=False)
    args = build_main_parser().parse_args(['getfastq', '--out_dir', str(tmp_path), '--metadata', str(metadata),
                                          '--threads', '1'])
    assert args.fastp
    getfastq_main(args)
    outputs = sorted((tmp_path / 'getfastq' / 'R1').glob('*.amalgkit.fastq.gz'))
    assert len(outputs) == len(sources)
    assert all(validate_fastq_structure(str(path)) == 24 for path in outputs)
    mtimes = [path.stat().st_mtime_ns for path in outputs]
    getfastq_main(args)
    assert [path.stat().st_mtime_ns for path in outputs] == mtimes
    assert [hashlib.sha256(path.read_bytes()).hexdigest() for path in sources] == hashes


@pytest.mark.skipif(shutil.which('oarfish') is None, reason='oarfish is required')
@pytest.mark.parametrize('compressed', [False, True])
def test_real_oarfish_filtered_read_denominator_and_reuse(tmp_path, compressed):
    rng = numpy.random.default_rng(94)
    sequence = ''.join(rng.choice(list('ACGT'), size=2400))
    fasta_dir = tmp_path / 'fasta'
    fasta_dir.mkdir()
    (fasta_dir / 'Test_species.fa').write_text('>t1\n' + sequence + '\n')
    read_dir = tmp_path / 'getfastq' / 'R1'
    read_dir.mkdir(parents=True)
    reads = read_dir / f'R1.fastq{ ".gz" if compressed else ""}'
    write_reads(reads, [sequence] * 12)
    metadata = tmp_path / 'metadata.tsv'
    pandas.DataFrame([dict(run='R1', scientific_name='Test species', lib_layout='single',
                           platform='OXFORD_NANOPORE', total_spots=1200, total_bases=2880000,
                           spot_length=2400, exclusion='no', is_sampled='yes')]).to_csv(metadata, sep='\t', index=False)
    args = build_main_parser().parse_args(['quant', '--out_dir', str(tmp_path), '--metadata', str(metadata),
                                          '--threads', '1', '--build_index', 'yes', '--quant_backend', 'oarfish',
                                          '--oarfish_seq_tech', 'ont-cdna'])
    quant_main(args)
    info_path = tmp_path / 'quant' / 'R1' / 'R1_run_info.json'
    info = json.loads(info_path.read_text())
    assert info['num_processed'] == 12
    assert info['num_pseudoaligned'] == pytest.approx(12)
    assert info['p_pseudoaligned'] == pytest.approx(100)
    assert not reads.exists()
    before = info_path.read_bytes(), info_path.stat().st_mtime_ns
    quant_main(args)
    assert (info_path.read_bytes(), info_path.stat().st_mtime_ns) == before
    args.oarfish_options = '--best-n 25'
    with pytest.raises(ValueError, match='oarfish options differ'):
        quant_main(args)
