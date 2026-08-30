"""Run with the wheel environment's Python -I, outside the source checkout."""

import gzip
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
print(f'Wheel import and bundled dataset verified: {installed_module}')
