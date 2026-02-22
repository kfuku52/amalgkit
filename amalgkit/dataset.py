"""
amalgkit dataset: Extract bundled test datasets.

This module provides functionality to extract pre-bundled datasets
for testing and demonstration purposes.
"""

import os
import shutil
import sys

# Available datasets and their descriptions
DATASETS = {
    'yeast': {
        'description': 'Two yeast species (S. cerevisiae + S. pombe) with BUSCO genes only (~3 MB)',
        'species': ['Saccharomyces_cerevisiae', 'Schizosaccharomyces_pombe'],
        'files': {
            'fasta': [
                'Saccharomyces_cerevisiae.fa.gz',
                'Schizosaccharomyces_pombe.fa.gz',
            ],
            'busco': [
                'Saccharomyces_cerevisiae_busco.tsv',
                'Schizosaccharomyces_pombe_busco.tsv',
            ],
            'config': [
                'control_term.config',
                'exclude_keyword.config',
                'group_attribute.config',
            ],
        },
    },
}


def get_dataset_dir():
    """Return the path to the bundled datasets directory."""
    return os.path.join(os.path.dirname(__file__), 'datasets')


def list_datasets():
    """List available datasets with descriptions."""
    print('Available datasets:')
    for name, info in DATASETS.items():
        print(f'  {name}: {info["description"]}')
        print(f'    Species: {", ".join(info["species"])}')


def validate_dataset_name(name):
    if name in DATASETS:
        return
    available = ', '.join(DATASETS.keys())
    sys.stderr.write(f'Error: Unknown dataset "{name}". Available: {available}\n')
    sys.exit(1)


def resolve_dataset_source_dir(name):
    dataset_src = os.path.join(get_dataset_dir(), name)
    if os.path.isdir(dataset_src):
        return dataset_src
    sys.stderr.write(f'Error: Dataset directory not found: {dataset_src}\n')
    sys.exit(1)


def build_extracted_dirs(out_dir):
    out_dir = os.path.realpath(out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    extracted_dirs = {
        'fasta': os.path.join(out_dir, 'fasta'),
        'busco': os.path.join(out_dir, 'busco'),
        'config': os.path.join(out_dir, 'config'),
    }
    for path_dir in extracted_dirs.values():
        os.makedirs(path_dir, exist_ok=True)
    return extracted_dirs


def validate_dataset_source_files(dataset, dataset_src):
    missing_sources = []
    for file_type, files in dataset['files'].items():
        for filename in files:
            src = os.path.join(dataset_src, filename)
            if not os.path.isfile(src):
                missing_sources.append(src)
    if missing_sources:
        raise FileNotFoundError(
            'Dataset source file(s) not found: {}'.format(', '.join(missing_sources))
        )


def copy_dataset_files(dataset, dataset_src, extracted_dirs, overwrite=False):
    for file_type, files in dataset['files'].items():
        dest_dir = extracted_dirs.get(file_type)
        if dest_dir is None:
            continue
        for filename in files:
            src = os.path.join(dataset_src, filename)
            dst = os.path.join(dest_dir, filename)
            if os.path.exists(dst) and (not os.path.isfile(dst)):
                raise NotADirectoryError(
                    'Destination path exists but is not a file: {}'.format(dst)
                )
            if os.path.exists(dst) and not overwrite:
                print(f'  Skipping (exists): {dst}')
                continue
            shutil.copy2(src, dst)
            print(f'  Copied: {filename} -> {dest_dir}/')


def extract_dataset(name, out_dir, overwrite=False):
    """
    Extract a bundled dataset to the specified output directory.

    Parameters
    ----------
    name : str
        Name of the dataset (e.g., 'yeast')
    out_dir : str
        Output directory path
    overwrite : bool
        Whether to overwrite existing files

    Returns
    -------
    dict
        Paths to extracted directories: {'fasta': ..., 'busco': ..., 'config': ...}
    """
    validate_dataset_name(name)
    dataset = DATASETS[name]
    dataset_src = resolve_dataset_source_dir(name)
    validate_dataset_source_files(dataset, dataset_src)
    extracted_dirs = build_extracted_dirs(out_dir)
    copy_dataset_files(dataset, dataset_src, extracted_dirs, overwrite=overwrite)
    return extracted_dirs


def dataset_main(args):
    """Main entry point for the dataset command."""

    if args.list:
        list_datasets()
        return

    if args.name is None:
        sys.stderr.write('Error: Please specify --name or use --list to see available datasets.\n')
        sys.exit(1)

    print(f'Extracting dataset "{args.name}" to: {args.out_dir}')

    paths = extract_dataset(
        name=args.name,
        out_dir=args.out_dir,
        overwrite=args.overwrite,
    )

    print('')
    print('Dataset extracted successfully!')
    print(f'  FASTA files:  {paths["fasta"]}')
    print(f'  BUSCO files:  {paths["busco"]}')
    print(f'  Config files: {paths["config"]}')
    print('')
    print('Example usage with this dataset:')
    print(f'  amalgkit config --out_dir {args.out_dir} --config base --overwrite yes')
    print(f'  cp {paths["config"]}/*.config {args.out_dir}/config_base/')
