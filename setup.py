import os
import re
import ast

from setuptools import setup, find_packages

with open(os.path.join('amalgkit', '__init__.py')) as f:
    match = re.search(r'__version__\s+=\s+(.*)', f.read())
version = str(ast.literal_eval(match.group(1)))

setup(
    name                    = 'amalgkit',
    version                 = version,
    description             = 'Tools for transcriptome amalgamation',
    license                 = "MIT License",
    author                  = "Kenji Fukushima, Matthias Freund",
    author_email            = 'kfuku52@gmail.com, matthias_freund@outlook.com',
    url                     = 'https://github.com/kfuku52/amalgkit.git',
    keywords                = 'transcriptome amalgamation',
    packages                = find_packages(exclude=('*.tests', 'tests', 'tests.*', '*.test', '*.test.*', '*.__pycache__', '__pycache__')),
    package_data            = {
        'amalgkit': [
            'config_dir/base/*.config',
            'config_dir/plantae/*.config',
            'config_dir/test/*.config',
            'config_dir/vertebrate/*.config',
            'datasets/yeast/*.fa.gz',
            'datasets/yeast/*.tsv',
            'datasets/yeast/*.config',
        ],
    },
    install_requires        = [
        'numpy',
        'pandas',
        'scipy',
        'matplotlib',
        'statsmodels',
        'scikit-learn',
        'inmoose',
        'biopython',
        'ete4',
    ],
    scripts                 = ['amalgkit/amalgkit',],
    include_package_data    = False,
    python_requires         = '>=3.9',
)
