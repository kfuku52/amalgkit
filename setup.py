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
            'select_rule_sets/base/select_rules.tsv',
            'select_rule_sets/plantae/select_rules.tsv',
            'select_rule_sets/test/select_rules.tsv',
            'select_rule_sets/vertebrate/select_rules.tsv',
            'datasets/yeast/*.fa.gz',
            'datasets/yeast/*.tsv',
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
    entry_points            = {
        'console_scripts': [
            'amalgkit=amalgkit.cli_entry:main',
        ],
    },
    include_package_data    = False,
    python_requires         = '>=3.9',
)
