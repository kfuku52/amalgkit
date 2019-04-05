from setuptools import setup, find_packages

setup(
    name             = 'amalgkit',
    version          = "0.1",
    description      = 'Tools for transcriptome amalgamation',
    license          = "BSD 3-clause License",
    author           = "Kenji Fukushima",
    author_email     = 'kfuku52@gmail.com',
    url              = 'https://github.com/kfuku52/amalgkit.git',
    keywords         = '',
    packages         = find_packages(),
    install_requires = ['numpy','pandas','biopython',],
    scripts          = ['amalgkit/amalgkit',],
)