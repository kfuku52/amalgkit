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
    license                 = "BSD 3-clause License",
    author                  = "Kenji Fukushima, Matthias Freund",
    author_email            = 'kfuku52@gmail.com, matthias_freund@outlook.com',
    url                     = 'https://github.com/kfuku52/amalgkit.git',
    keywords                = 'transcriptome amalgamation',
    packages                = find_packages(),
    install_requires        = ['numpy','pandas','biopython','lxml','nltk','obonet'],
    scripts                 = ['amalgkit/amalgkit',],
    include_package_data    = True
)