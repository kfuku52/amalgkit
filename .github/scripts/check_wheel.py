"""Run with the wheel environment's Python -I, outside the source checkout."""

import gzip
import csv
import subprocess
import sys
import tempfile
from importlib import resources
from pathlib import Path

import amalgkit

installed_module = Path(amalgkit.__file__).resolve()
if not installed_module.is_relative_to(Path(sys.prefix).resolve()):
    raise RuntimeError(f'Wheel smoke test imported outside its environment: {installed_module}')
data = resources.files('amalgkit')
if not data.joinpath('select_rule_sets/base/select_rules.tsv').is_file():
    raise RuntimeError('Wheel is missing its selection rules')
with tempfile.TemporaryDirectory(prefix='amalgkit-wheel-data-') as scratch:
    for args in (['--version'], ['dataset', '--list'], ['dataset', '--name', 'yeast', '--out_dir', scratch]):
        subprocess.run(  # noqa: S603 - fixed CLI commands in an isolated environment
            [sys.executable, '-I', '-m', 'amalgkit', *args], cwd=scratch, check=True,
        )
    with gzip.open(Path(scratch) / 'fasta' / 'Saccharomyces_cerevisiae.fa.gz', 'rt') as handle:
        if not handle.readline().startswith('>'):
            raise RuntimeError('Extracted wheel dataset is not a FASTA')
    metadata_path = Path(scratch) / 'metadata' / 'metadata.tsv'
    metadata_path.write_text(
        'run\tscientific_name\ttaxid\ttaxid_species\tsample_group\ttotal_spots\tbioproject\tbiosample\texclusion\n'
        '0001\tSaccharomyces cerevisiae\t4932\t4932\twild_type\t2000000\tP1\tB1\tno\n'
        'NA\tSaccharomyces cerevisiae\t4932\t4932\theat_stress\t2000000\tP1\tB2\tno\n',
        encoding='utf-8',
    )
    subprocess.run(  # noqa: S603 - fixed CLI command on synthetic, network-free metadata
        [sys.executable, '-I', '-m', 'amalgkit', 'select', '--out_dir', scratch], cwd=scratch, check=True,
    )
    with metadata_path.open() as handle:
        selected = {row['run']: row['sample_group'] for row in csv.DictReader(handle, delimiter='\t')
                    if row['exclusion'] == 'no' and row['is_sampled'] == 'yes'}
    if selected != {'0001': 'wild type', 'NA': 'heat stress'}:
        raise RuntimeError(f'Packaged yeast rules lost reviewed groups or lexical run IDs: {selected}')
print(f'Wheel import, bundled dataset and yeast selection verified: {installed_module}')
