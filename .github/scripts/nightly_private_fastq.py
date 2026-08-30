"""Prepare/verify plain and gzip, single and paired private inputs for real E2E."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy
import pandas


def prepare(root):
    parts = []
    with gzip.open(root / 'fasta' / 'Saccharomyces_cerevisiae.fa.gz', 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                if parts:
                    break
            else:
                parts.append(line.strip())
    fragment = ''.join(parts)[:200]
    if len(fragment) < 150:
        raise ValueError('E2E reference needs a transcript of at least 150 bases')
    reads = [fragment[:75], fragment[-75:].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]]
    inputs = root / 'private-fastq' / 'Saccharomyces_cerevisiae'
    inputs.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for compressed in (False, True):
        for paired in (False, True):
            run = ('gzip' if compressed else 'plain') + ('_paired' if paired else '_single')
            for mate in range(2 if paired else 1):
                name = run + (f'_{mate + 1}' if paired else '') + ('.fastq.gz' if compressed else '.fastq')
                path = inputs / name
                opener = gzip.open if compressed else open
                with opener(path, 'wt', encoding='utf-8') as handle:
                    for index in range(20):
                        handle.write(f'@nightly{index}/{mate + 1}\n{reads[mate]}\n+\n' + 'I' * 75 + '\n')
                manifest[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    (root / 'private-input-sha256.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')


def verify(root):
    manifest = json.loads((root / 'private-input-sha256.json').read_text())
    for relative, expected in manifest.items():
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'Private source was modified: {relative}')
    markers = list((root / 'getfastq').glob('*/*.amalgkit.fastq.gz.safely_removed'))
    if len(markers) != 6 or list((root / 'getfastq').glob('*/*.amalgkit.fastq.gz')):
        raise RuntimeError('Default cleanup did not retire all six managed FASTQ entries')
    counts = pandas.read_csv(root / 'merge' / 'Saccharomyces_cerevisiae' / 'Saccharomyces_cerevisiae_est_counts.tsv', sep='\t')
    expected_runs = {'plain_single', 'plain_paired', 'gzip_single', 'gzip_paired'}
    if set(counts.columns) != {'target_id', *expected_runs}:
        raise RuntimeError(f'Missing E2E runs in merge: {counts.columns.tolist()}')
    if not (counts[list(expected_runs)].sum(axis=0) > 0).all():
        raise RuntimeError('E2E reads did not quantify')
    metadata = pandas.read_csv(root / 'cstmm' / 'metadata.tsv', sep='\t')
    factors = metadata['tmm_normalization_factor'].to_numpy(dtype=float)
    if not numpy.isfinite(factors).all() or not (factors > 0).all():
        raise RuntimeError('E2E normalization factors are invalid')
    print('Verified four layouts/formats, six unchanged private sources, cleanup, merge and TMM')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'verify'])
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    (prepare if args.stage == 'prepare' else verify)(args.root)
